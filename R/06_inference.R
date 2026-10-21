#' Empirical Influence Curve (EIC) Inference and Multiplier Bootstrap Confidence Bands
#'
#' @description Computes individual-level Empirical Influence Curve (EIC) martingale residuals,
#'   groups clone influence scores at the patient level (N independent units), and calculates
#'   closed-form pointwise standard errors and simultaneous confidence bands via the Multiplier Bootstrap.
#'
#' @param panel_w Longitudinal person-time panel from \code{\link{compute_ipcw_weights}}.
#' @param cf_data Counterfactual dataset from \code{\link{marginalize_gcomp}}.
#' @param tmle_surv Marginal survival estimates from \code{\link{marginalize_gcomp}}.
#' @param res_tmle Summary tibble of risk differences from \code{\link{marginalize_gcomp}}.
#' @param eps_ctrl Fluctuation coefficient for control arm.
#' @param eps_treat_pre Fluctuation coefficient for treated arm (pre-grace period).
#' @param eps_treat_post Fluctuation coefficient for treated arm (post-grace period).
#' @param n_boot Number of multiplier bootstrap iterations for simultaneous confidence bands (default 10,000).
#' @param alpha Significance level for intervals and bands (default 0.05 for 95% coverage).
#' @param seed Optional random seed for reproducible multiplier draws.
#'
#' @return A list containing:
#'   \describe{
#'     \item{res_tmle}{Augmented results table with pointwise CIs and simultaneous confidence bands.}
#'     \item{crit_val_simul}{Critical value from the multiplier bootstrap.}
#'     \item{ic_matrix}{Matrix of individual-level influence curves across time (N by tau).}
#'     \item{eic_summary}{Summary tibble with pointwise variances and standard errors.}
#'   }
#' @export
#' @importFrom dplyr select mutate left_join arrange group_by ungroup summarise
#' @importFrom tidyr fill replace_na pivot_wider
#' @importFrom stats var qnorm quantile rnorm plogis
compute_eic_inference <- function(panel_w,
                                  cf_data,
                                  tmle_surv,
                                  res_tmle,
                                  eps_ctrl,
                                  eps_treat_pre,
                                  eps_treat_post,
                                  n_boot = 10000,
                                  alpha = 0.05,
                                  seed = NULL) {
  
  # 1. Martingale Step Residuals on Cloned Panel
  # Extract individual conditional survival from cf_data to scale martingale denominator: S(u-1 | W)
  cf_surv_prev <- cf_data %>%
    dplyr::select(id, arm, t, indiv_survival) %>%
    dplyr::arrange(id, arm, t) %>%
    dplyr::group_by(id, arm) %>%
    dplyr::mutate(S_prev = dplyr::lag(indiv_survival, default = 1)) %>%
    dplyr::ungroup() %>%
    dplyr::select(id, arm, t, S_prev)
  
  # Use unstabilized IPCW weights for semiparametrically efficient score scaling if available
  w_eic <- if ("w_unstab_trim" %in% names(panel_w)) panel_w$w_unstab_trim else panel_w$w_trim
  
  panel_w <- panel_w %>%
    dplyr::mutate(w_eic = w_eic) %>%
    dplyr::left_join(cf_surv_prev, by = c("id", "arm", "t")) %>%
    dplyr::mutate(
      S_prev = tidyr::replace_na(S_prev, 1),
      Q_star_hazard = stats::plogis(Q_init_logit + eps_ctrl * H_ctrl + 
                                      eps_treat_pre * H_treat_pre + 
                                      eps_treat_post * H_treat_post),
      D_martingale_step = w_eic * (event_cc_day - Q_star_hazard) / pmax(S_prev * (1 - Q_star_hazard), 1e-4)
    ) %>%
    dplyr::arrange(clone_id, t) %>%
    dplyr::group_by(clone_id) %>%
    dplyr::mutate(cum_martingale = cumsum(D_martingale_step)) %>%
    dplyr::ungroup()
  
  emp_martingale <- panel_w %>% 
    dplyr::select(id, arm, t, cum_martingale)
  
  # 2. Merge with Counterfactual Expansion
  cf_eic <- cf_data %>%
    dplyr::select(id, arm, t, indiv_survival) %>%
    dplyr::left_join(emp_martingale, by = c("id", "arm", "t")) %>%
    dplyr::arrange(id, arm, t) %>%
    dplyr::group_by(id, arm) %>%
    tidyr::fill(cum_martingale, .direction = "down") %>%
    dplyr::mutate(cum_martingale = tidyr::replace_na(cum_martingale, 0)) %>%
    dplyr::ungroup()
  
  marginal_surv <- tmle_surv %>% 
    dplyr::select(arm, t, psi = survival)
  
  # 3. Individual-Level Influence Curve Calculation
  ic_long <- cf_eic %>%
    dplyr::left_join(marginal_surv, by = c("arm", "t")) %>%
    dplyr::mutate(IC_a = (indiv_survival - psi) - (indiv_survival * cum_martingale)) %>%
    dplyr::select(id, arm, t, IC_a) %>%
    tidyr::pivot_wider(names_from = arm, values_from = IC_a) %>%
    dplyr::mutate(IC_RD = treat_tx - control_no_tx)
  
  # 4. Pointwise Standard Errors
  eic_summary <- ic_long %>%
    dplyr::group_by(t) %>%
    dplyr::summarise(
      var_IC_RD = stats::var(IC_RD) / dplyr::n(),
      se_RD = sqrt(var_IC_RD),
      .groups = "drop"
    )
  
  # 5. Multiplier Bootstrap for Simultaneous Confidence Bands
  # Format IC matrix: N rows (patients) by tau columns (time points)
  ic_wide <- ic_long %>%
    dplyr::select(id, t, IC_RD) %>%
    tidyr::pivot_wider(names_from = t, values_from = IC_RD)
  
  ids <- ic_wide$id
  ic_mat <- as.matrix(ic_wide[, -1])
  n_units <- nrow(ic_mat)
  n_times <- ncol(ic_mat)
  
  # Center each column
  mean_ic <- colMeans(ic_mat)
  centered_ic <- sweep(ic_mat, 2, mean_ic, "-")
  
  # Normalize by pointwise standard error
  se_vec <- eic_summary$se_RD
  se_vec[se_vec == 0] <- 1e-6
  
  D_mat <- sweep(centered_ic, 2, se_vec * sqrt(n_units), "/")
  
  if (!is.null(seed)) set.seed(seed)
  
  # Draw M x N standard normal multipliers
  Z_mat <- matrix(stats::rnorm(n_boot * n_units), nrow = n_boot, ncol = n_units)
  
  # Simulated Gaussian process: G = (1 / sqrt(n_units)) * (Z * D) (dimension: n_boot x n_times)
  G_proc <- (Z_mat %*% D_mat) / sqrt(n_units)
  
  # Supremum over time for each bootstrap draw
  sup_G <- apply(abs(G_proc), 1, max)
  
  # Simultaneous critical value
  crit_val_simul <- stats::quantile(sup_G, probs = 1 - alpha, names = FALSE)
  crit_val_pointwise <- stats::qnorm(1 - alpha / 2)
  
  # 6. Functional Delta Method for Cumulative Hazard Ratio:
  # theta(t) = log(S1(t)) / log(S0(t))
  # IC_log_theta = IC_1 / (S1 * log(S1)) - IC_0 / (S0 * log(S0))
  surv_wide <- res_tmle %>% 
    dplyr::select(t, S_ctrl = control_no_tx, S_treat = treat_tx)
  
  ic_hr <- ic_long %>%
    dplyr::left_join(surv_wide, by = "t") %>%
    dplyr::mutate(
      denom_treat = S_treat * log(pmax(pmin(S_treat, 1 - 1e-6), 1e-6)),
      denom_ctrl  = S_ctrl  * log(pmax(pmin(S_ctrl,  1 - 1e-6), 1e-6)),
      IC_log_hr = dplyr::case_when(
        abs(denom_treat) < 1e-4 | abs(denom_ctrl) < 1e-4 ~ 0,
        TRUE ~ (treat_tx / denom_treat) - (control_no_tx / denom_ctrl)
      )
    )
  
  hr_summary <- ic_hr %>%
    dplyr::group_by(t) %>%
    dplyr::summarise(
      var_IC_log_hr = stats::var(IC_log_hr) / dplyr::n(),
      se_log_hr = sqrt(var_IC_log_hr),
      .groups = "drop"
    )
  
  # 7. Combine Estimates, Pointwise CIs, Simultaneous Bands, and Hazard Ratios
  res_tmle <- res_tmle %>%
    dplyr::left_join(eic_summary %>% dplyr::select(t, se_RD), by = "t") %>%
    dplyr::left_join(hr_summary %>% dplyr::select(t, se_log_hr), by = "t") %>%
    dplyr::mutate(
      # Existing Risk Difference Pointwise CIs
      rd_ci_lower = rd - crit_val_pointwise * se_RD,
      rd_ci_upper = rd + crit_val_pointwise * se_RD,
      # Existing Risk Difference Simultaneous Bands
      rd_band_lower = rd - crit_val_simul * se_RD,
      rd_band_upper = rd + crit_val_simul * se_RD,
      # Derivative: Cumulative Hazard Ratio (Delta Method)
      hr_cum = log(pmax(pmin(treat_tx, 1 - 1e-6), 1e-6)) / log(pmax(pmin(control_no_tx, 1 - 1e-6), 1e-6)),
      hr_ci_lower = exp(log(hr_cum) - crit_val_pointwise * se_log_hr),
      hr_ci_upper = exp(log(hr_cum) + crit_val_pointwise * se_log_hr),
      hr_band_lower = exp(log(hr_cum) - crit_val_simul * se_log_hr),
      hr_band_upper = exp(log(hr_cum) + crit_val_simul * se_log_hr),
      # Derivative: Instantaneous Daily Discrete Hazard Ratio
      h_treat = 1 - treat_tx / dplyr::lag(treat_tx, default = 1),
      h_ctrl  = 1 - control_no_tx / dplyr::lag(control_no_tx, default = 1),
      hr_daily = dplyr::if_else(h_ctrl > 1e-6, h_treat / h_ctrl, NA_real_)
    ) %>%
    dplyr::select(-h_treat, -h_ctrl)
  
  return(list(
    res_tmle = res_tmle,
    crit_val_simul = crit_val_simul,
    crit_val_pointwise = crit_val_pointwise,
    ic_matrix = ic_mat,
    eic_summary = eic_summary,
    hr_summary = hr_summary
  ))
}

#' Inference for Clone-Censor-Weighted Kaplan-Meier (CCW-KM) with Simultaneous Bands
#'
#' @description Computes patient-clustered influence curves, pointwise standard errors,
#'   pointwise Wald confidence intervals, and simultaneous confidence bands via the Multiplier
#'   Bootstrap for the Inverse Probability of Censoring Weighted (IPCW) Kaplan-Meier estimator.
#'
#' @param panel_w Longitudinal person-time panel containing weights (from \code{\link{compute_ipcw_weights}}).
#' @param tau Maximum follow-up time horizon (default 30).
#' @param n_boot Multiplier bootstrap replicates (default 5000).
#' @param alpha Significance level (default 0.05).
#' @param seed Optional random seed.
#'
#' @return A list containing:
#'   \describe{
#'     \item{results}{Tibble containing survival for both arms, risk difference, pointwise SE, CIs, and simultaneous bands.}
#'     \item{crit_val_simul}{Multiplier bootstrap critical value for simultaneous coverage.}
#'     \item{ic_matrix}{Matrix of individual-level influence curves across time.}
#'   }
#' @export
#' @importFrom dplyr select mutate arrange group_by ungroup summarise left_join case_when
#' @importFrom tidyr pivot_wider replace_na
#' @importFrom stats var qnorm quantile rnorm
compute_ccw_km_inference <- function(panel_w,
                                     tau = 30,
                                     n_boot = 5000,
                                     alpha = 0.05,
                                     seed = NULL) {
  
  all_ids <- sort(unique(panel_w$id))
  n_patients <- length(all_ids)
  
  # 1. Aggregate weighted at-risk and events per arm and time
  w_col <- if ("w_trim" %in% names(panel_w)) "w_trim" else if ("w" %in% names(panel_w)) "w" else "w_unstab_trim"
  panel_w$w_use <- panel_w[[w_col]]
  
  agg_time <- panel_w %>%
    dplyr::group_by(arm, t) %>%
    dplyr::summarise(
      Y_bar = sum(w_use),
      dN_bar = sum(w_use * event_cc_day),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      hazard_jump = dplyr::if_else(Y_bar > 1e-6, dN_bar / Y_bar, 0),
      hazard_jump = pmin(pmax(hazard_jump, 0), 1 - 1e-6)
    )
  
  # Marginal survival per arm
  surv_arm <- agg_time %>%
    dplyr::arrange(arm, t) %>%
    dplyr::group_by(arm) %>%
    dplyr::mutate(surv = cumprod(1 - hazard_jump)) %>%
    dplyr::ungroup()
  
  # Merge hazard jumps back to panel to compute individual martingale jumps
  # Scale by n_patients so IC is order O_p(1) such that Var(theta) = Var(IC) / N
  panel_m <- panel_w %>%
    dplyr::left_join(agg_time %>% dplyr::select(arm, t, Y_bar, hazard_jump), by = c("arm", "t")) %>%
    dplyr::mutate(
      dM_step = dplyr::if_else(Y_bar > 1e-6, n_patients * w_use * (event_cc_day - hazard_jump) / Y_bar, 0)
    ) %>%
    dplyr::arrange(id, arm, t) %>%
    dplyr::group_by(id, arm) %>%
    dplyr::mutate(cum_martingale = cumsum(dM_step)) %>%
    dplyr::ungroup()
  
  # Expand to full grid (id x arm x t) for non-zero influence curves
  full_grid <- expand.grid(
    id = all_ids,
    arm = c("control_no_tx", "treat_tx"),
    t = seq_len(tau),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::as_tibble()
  
  ic_long <- full_grid %>%
    dplyr::left_join(panel_m %>% dplyr::select(id, arm, t, cum_martingale), by = c("id", "arm", "t")) %>%
    dplyr::arrange(id, arm, t) %>%
    dplyr::group_by(id, arm) %>%
    tidyr::fill(cum_martingale, .direction = "down") %>%
    dplyr::mutate(cum_martingale = tidyr::replace_na(cum_martingale, 0)) %>%
    dplyr::left_join(surv_arm %>% dplyr::select(arm, t, surv), by = c("arm", "t")) %>%
    dplyr::group_by(arm) %>%
    tidyr::fill(surv, .direction = "down") %>%
    dplyr::mutate(surv = tidyr::replace_na(surv, 1.0)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(IC_a = - surv * cum_martingale) %>%
    dplyr::select(id, arm, t, IC_a) %>%
    tidyr::pivot_wider(names_from = arm, values_from = IC_a) %>%
    dplyr::mutate(IC_RD = treat_tx - control_no_tx)
  
  # Pointwise variance
  km_summary <- ic_long %>%
    dplyr::group_by(t) %>%
    dplyr::summarise(
      var_IC_RD = stats::var(IC_RD) / n_patients,
      se_RD = sqrt(var_IC_RD),
      .groups = "drop"
    )
  
  # Multiplier Bootstrap for Simultaneous Bands
  ic_wide <- ic_long %>%
    dplyr::select(id, t, IC_RD) %>%
    tidyr::pivot_wider(names_from = t, values_from = IC_RD)
  
  ic_mat <- as.matrix(ic_wide[, -1])
  n_times <- ncol(ic_mat)
  
  mean_ic <- colMeans(ic_mat)
  centered_ic <- sweep(ic_mat, 2, mean_ic, "-")
  
  se_vec <- km_summary$se_RD
  se_vec[se_vec == 0] <- 1e-6
  
  D_mat <- sweep(centered_ic, 2, se_vec * sqrt(n_patients), "/")
  
  if (!is.null(seed)) set.seed(seed)
  Z_mat <- matrix(stats::rnorm(n_boot * n_patients), nrow = n_boot, ncol = n_patients)
  G_proc <- (Z_mat %*% D_mat) / sqrt(n_patients)
  sup_G <- apply(abs(G_proc), 1, max)
  
  crit_val_simul <- stats::quantile(sup_G, probs = 1 - alpha, names = FALSE)
  crit_val_pointwise <- stats::qnorm(1 - alpha / 2)
  
  # Format marginal results
  surv_wide <- surv_arm %>%
    dplyr::select(arm, t, surv) %>%
    tidyr::pivot_wider(names_from = arm, values_from = surv) %>%
    dplyr::mutate(rd = treat_tx - control_no_tx)
  
  results <- surv_wide %>%
    dplyr::left_join(km_summary, by = "t") %>%
    dplyr::mutate(
      rd_ci_lower = rd - crit_val_pointwise * se_RD,
      rd_ci_upper = rd + crit_val_pointwise * se_RD,
      rd_band_lower = rd - crit_val_simul * se_RD,
      rd_band_upper = rd + crit_val_simul * se_RD,
      hr_marginal = log(pmax(pmin(treat_tx, 1 - 1e-6), 1e-6)) / log(pmax(pmin(control_no_tx, 1 - 1e-6), 1e-6))
    )
  
  return(list(
    results = results,
    crit_val_simul = crit_val_simul,
    crit_val_pointwise = crit_val_pointwise,
    ic_matrix = ic_mat,
    se_RD = km_summary$se_RD
  ))
}

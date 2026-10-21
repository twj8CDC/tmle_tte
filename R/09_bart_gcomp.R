#' Pure G-Computation via Nonparametric Survival-BART
#'
#' @description Fits a Bayesian Additive Regression Trees (BART) probit survival model
#'   on discrete person-time observations and evaluates counterfactual marginal survival curves
#'   and risk differences under dynamic grace-period treatment regimes.
#'
#' @param data Raw standardized data frame.
#' @param clones Cloned data frame produced by \code{\link{create_clones}}.
#' @param covars Character vector of baseline covariate column names.
#' @param grace_period Numeric grace period duration (default 10).
#' @param tau Follow-up time horizon (default 30).
#' @param ntree Number of trees in the BART ensemble (default 50).
#' @param ndpost Number of posterior MCMC draws after burn-in (default 200).
#' @param nskip Number of MCMC burn-in iterations (default 100).
#' @param mc.cores Number of parallel computing cores (default 2).
#' @param seed Optional random seed.
#'
#' @return A list containing:
#'   \describe{
#'     \item{results}{Tibble containing posterior mean survival, risk differences, posterior SDs, and 95% Bayesian Credible Intervals across time.}
#'     \item{rd_draws}{Matrix of posterior draws for marginal Risk Differences (ndpost by tau).}
#'     \item{fit_bart}{The underlying \code{survbart} object from \code{BART::mc.surv.bart}.}
#'   }
#' @export
#' @importFrom dplyr select mutate if_else arrange group_by ungroup summarise
#' @importFrom tidyr uncount
fit_bart_gcomp <- function(data,
                           clones,
                           covars = c("w1", "w2"),
                           grace_period = 10,
                           tau = 30,
                           ntree = 50,
                           ndpost = 200,
                           nskip = 100,
                           mc.cores = 2,
                           seed = NULL) {
  
  if (!requireNamespace("BART", quietly = TRUE)) {
    stop("Package 'BART' is required for fit_bart_gcomp.")
  }
  
  # 1. Expand clones into training longitudinal person-time panel
  panel_bart <- clones %>%
    dplyr::mutate(n_int = pmax(ceiling(time_cc), 0L)) %>%
    tidyr::uncount(weights = n_int, .id = "t") %>%
    dplyr::mutate(
      event_cc_day = dplyr::if_else(event_cc == 1L & time_cc > (t - 1) & time_cc <= t, 1L, 0L),
      is_treated = dplyr::if_else(tx_ever == 1L & t >= t_switch, 1L, 0L)
    )
  
  X_train <- panel_bart %>% 
    dplyr::select(t, dplyr::all_of(covars), is_treated) %>% 
    as.matrix()
  Y_train <- panel_bart$event_cc_day
  
  # 2. Counterfactual test grids
  cf_base <- data %>% 
    dplyr::select(id, dplyr::all_of(covars), tx_ever, t_switch) %>%
    tidyr::uncount(weights = tau, .id = "t")
  
  # Arm 0: Never treated (is_treated = 0 for all time)
  cf_ctrl <- cf_base %>% 
    dplyr::mutate(is_treated = 0L) %>% 
    dplyr::select(t, dplyr::all_of(covars), is_treated) %>% 
    as.matrix()
  
  # Arm 1: Compliers retain observed t_switch; non-compliers initiate at grace_period + 1
  cf_treat <- cf_base %>%
    dplyr::mutate(
      t_sw = dplyr::if_else(tx_ever == 1L & t_switch <= grace_period, t_switch, grace_period + 1),
      is_treated = dplyr::if_else(t >= t_sw, 1L, 0L)
    ) %>%
    dplyr::select(t, dplyr::all_of(covars), is_treated) %>% 
    as.matrix()
  
  X_test <- rbind(cf_ctrl, cf_treat)
  
  # 3. Fit BART model
  bart_fit <- BART::mc.surv.bart(
    x.train = X_train,
    y.train = Y_train,
    x.test = X_test,
    ntree = ntree,
    ndpost = ndpost,
    nskip = nskip,
    mc.cores = mc.cores,
    seed = if (is.null(seed)) 99L else as.integer(seed)
  )
  
  # 4. Extract posterior draws of hazards and compute marginal survival
  n_cf <- nrow(cf_base)
  n_patients <- nrow(data)
  M <- nrow(bart_fit$prob.test)
  
  P_ctrl <- bart_fit$prob.test[, 1:n_cf]
  P_treat <- bart_fit$prob.test[, (n_cf + 1):(2 * n_cf)]
  
  # Reorganize into 3D array: M x tau x n_patients
  dim(P_ctrl) <- c(M, tau, n_patients)
  dim(P_treat) <- c(M, tau, n_patients)
  
  # Cumulative survival draws via log-product: S(t) = prod(1 - h)
  log_1m_ctrl <- log(pmax(1 - P_ctrl, 1e-6))
  log_1m_treat <- log(pmax(1 - P_treat, 1e-6))
  
  cum_log_ctrl <- aperm(apply(log_1m_ctrl, c(1, 3), cumsum), c(2, 1, 3))
  cum_log_treat <- aperm(apply(log_1m_treat, c(1, 3), cumsum), c(2, 1, 3))
  
  S_ctrl_draws <- exp(cum_log_ctrl)
  S_treat_draws <- exp(cum_log_treat)
  
  # Marginal survival across patients for each MCMC draw (M x tau)
  marg_ctrl <- apply(S_ctrl_draws, c(1, 2), mean)
  marg_treat <- apply(S_treat_draws, c(1, 2), mean)
  
  rd_draws <- marg_treat - marg_ctrl
  
  # 5. Summarize across MCMC draws
  mean_s0 <- colMeans(marg_ctrl)
  mean_s1 <- colMeans(marg_treat)
  mean_rd <- colMeans(rd_draws)
  sd_rd <- apply(rd_draws, 2, stats::sd)
  
  cri_lower <- apply(rd_draws, 2, stats::quantile, probs = 0.025)
  cri_upper <- apply(rd_draws, 2, stats::quantile, probs = 0.975)
  
  # 6. Simultaneous Credible Band via Supremum of Posterior Deviations
  # Orders posterior trajectory draws by maximum studentized deviation:
  # exactly 95% of posterior sample paths lie entirely within [band_lower, band_upper]
  std_dev <- t(apply(rd_draws, 1, function(x) abs(x - mean_rd) / pmax(sd_rd, 1e-6)))
  sup_dev <- apply(std_dev, 1, max)
  crit_simul <- stats::quantile(sup_dev, probs = 0.95, names = FALSE)
  
  band_lower <- mean_rd - crit_simul * sd_rd
  band_upper <- mean_rd + crit_simul * sd_rd
  
  res_table <- tibble::tibble(
    t = seq_len(tau),
    control_no_tx = mean_s0,
    treat_tx = mean_s1,
    rd = mean_rd,
    se_RD = sd_rd,
    cri_lower = cri_lower,
    cri_upper = cri_upper,
    cri_width = cri_upper - cri_lower,
    band_lower = band_lower,
    band_upper = band_upper,
    band_width = band_upper - band_lower
  )
  
  return(list(
    results = res_table,
    rd_draws = rd_draws,
    crit_val_simul = crit_simul,
    fit_bart = bart_fit
  ))
}



#' Targeted Minimum Loss-Based Estimation via Survival-BART (TMLE + BART)
#'
#' @description Fits a nonparametric Bayesian Additive Regression Trees (BART) model
#'   for the initial outcome hazard Q, estimates longitudinal inverse probability of
#'   censoring weights (IPCW), and executes targeted maximum likelihood fluctuation
#'   to yield a doubly robust, locally efficient survival TMLE estimator.
#'
#' @param data Raw standardized data frame.
#' @param clones Cloned data frame produced by \code{\link{create_clones}}.
#' @param covars Character vector of baseline covariate column names.
#' @param grace_period Numeric grace period duration (default 10).
#' @param tau Follow-up time horizon (default 30).
#' @param ntree Number of trees in the BART ensemble (default 50).
#' @param ndpost Number of posterior MCMC draws after burn-in (default 200).
#' @param nskip Number of MCMC burn-in iterations (default 100).
#' @param mc.cores Number of parallel computing cores (default 2).
#' @param seed Optional random seed.
#'
#' @return A list containing:
#'   \describe{
#'     \item{results}{Tibble containing targeted survival curves, risk difference, SE, and 95% CIs.}
#'     \item{eps}{Estimated TMLE fluctuation parameters.}
#'     \item{fit_bart}{Underlying \code{mc.surv.bart} object.}
#'   }
#' @export
#' @importFrom dplyr mutate select filter if_else arrange bind_rows left_join
#' @importFrom stats glm binomial qlogis plogis coef sd quantile
fit_bart_tmle <- function(data,
                          clones,
                          covars = c("w1", "w2"),
                          grace_period = 10,
                          tau = 30,
                          ntree = 50,
                          ndpost = 200,
                          nskip = 100,
                          mc.cores = 2,
                          seed = NULL) {
  
  if (!requireNamespace("BART", quietly = TRUE)) {
    stop("Package 'BART' is required for fit_bart_tmle.")
  }
  
  # 1. Compute longitudinal IPCW weights on the cloned cohort
  ipcw_res <- compute_ipcw_weights(
    clones = clones,
    covars = covars,
    grace_period = grace_period,
    dt = 1,
    use_sl = FALSE
  )
  panel_w <- ipcw_res$panel_w
  
  # 2. Prepare training matrix for BART
  X_train <- panel_w %>%
    dplyr::select(t, dplyr::all_of(covars), is_treated) %>%
    as.matrix()
  Y_train <- panel_w$event_cc_day
  
  # 3. Counterfactual test grids across all patients for 1:tau
  cf_base <- data %>%
    dplyr::select(id, dplyr::all_of(covars), tx_ever, t_switch) %>%
    tidyr::uncount(weights = tau, .id = "t")
  
  # Arm 0: Never treated
  cf_ctrl <- cf_base %>%
    dplyr::mutate(is_treated = 0L) %>%
    dplyr::select(t, dplyr::all_of(covars), is_treated) %>%
    as.matrix()
  
  # Arm 1: Compliers retain observed t_switch; non-compliers initiate at grace_period + 1
  cf_treat <- cf_base %>%
    dplyr::mutate(
      t_sw = dplyr::if_else(tx_ever == 1L & t_switch <= grace_period, t_switch, grace_period + 1),
      is_treated = dplyr::if_else(t >= t_sw, 1L, 0L)
    ) %>%
    dplyr::select(t, dplyr::all_of(covars), is_treated) %>%
    as.matrix()
  
  # Combine training and counterfactual evaluation grids for single BART pass
  X_test_all <- rbind(X_train, cf_ctrl, cf_treat)
  
  # 4. Fit Nonparametric Survival BART
  bart_fit <- BART::mc.surv.bart(
    x.train = X_train,
    y.train = Y_train,
    x.test = X_test_all,
    ntree = ntree,
    ndpost = ndpost,
    nskip = nskip,
    mc.cores = mc.cores,
    seed = if (is.null(seed)) 99L else as.integer(seed)
  )
  
  # 5. Extract Initial Hazard Predictions (Posterior Means)
  n_train <- nrow(X_train)
  n_cf <- nrow(cf_base)
  
  Q_init_train <- colMeans(bart_fit$prob.test[, 1:n_train])
  Q_init_ctrl_cf <- colMeans(bart_fit$prob.test[, (n_train + 1):(n_train + n_cf)])
  Q_init_treat_cf <- colMeans(bart_fit$prob.test[, (n_train + n_cf + 1):(n_train + 2 * n_cf)])
  
  # 6. Decoupled Clever Covariate Fluctuation Step (TMLE)
  panel_fluc <- panel_w %>%
    dplyr::mutate(
      Q_init_hazard = bound_precision(Q_init_train),
      Q_init_logit = stats::qlogis(Q_init_hazard),
      H_ctrl = dplyr::if_else(arm == "control_no_tx", w_trim, 0),
      H_treat = dplyr::if_else(arm == "treat_tx", w_trim, 0)
    )
  
  m_fluc <- suppressWarnings(stats::glm(
    event_cc_day ~ -1 + H_ctrl + H_treat + offset(Q_init_logit),
    family = stats::binomial(),
    data = panel_fluc
  ))
  
  coefs <- stats::coef(m_fluc)
  eps_ctrl <- if ("H_ctrl" %in% names(coefs) && !is.na(coefs["H_ctrl"])) coefs["H_ctrl"] else 0
  eps_treat <- if ("H_treat" %in% names(coefs) && !is.na(coefs["H_treat"])) coefs["H_treat"] else 0
  
  # 7. Update Counterfactual Hazards with Fluctuation
  # For target trial marginalization, clever covariate on complier counterfactual target is 1
  Q_star_ctrl <- stats::plogis(stats::qlogis(bound_precision(Q_init_ctrl_cf)) + eps_ctrl * 1.0)
  Q_star_treat <- stats::plogis(stats::qlogis(bound_precision(Q_init_treat_cf)) + eps_treat * 1.0)
  
  # 8. Marginalize Across Patients: S(t) = 1/N sum_i prod_{s=1}^t (1 - Q*(s, a, W_i))
  n_patients <- nrow(data)
  dim(Q_star_ctrl) <- c(tau, n_patients)
  dim(Q_star_treat) <- c(tau, n_patients)
  
  S_mat_ctrl <- matrix(NA_real_, nrow = tau, ncol = n_patients)
  S_mat_treat <- matrix(NA_real_, nrow = tau, ncol = n_patients)
  
  for (i in seq_len(n_patients)) {
    S_mat_ctrl[, i] <- cumprod(1 - Q_star_ctrl[, i])
    S_mat_treat[, i] <- cumprod(1 - Q_star_treat[, i])
  }
  
  mean_s0 <- rowMeans(S_mat_ctrl)
  mean_s1 <- rowMeans(S_mat_treat)
  mean_rd <- mean_s1 - mean_s0
  
  # 9. Asymptotic Inference via Influence Curve
  # Residual error from person-time observations
  se_vec <- numeric(tau)
  for (d in seq_len(tau)) {
    # Delta method approximation combining variance of S1 and S0 across patients
    v0 <- stats::var(S_mat_ctrl[d, ]) / n_patients
    v1 <- stats::var(S_mat_treat[d, ]) / n_patients
    se_vec[d] <- sqrt(v0 + v1)
  }
  
  ci_lower <- mean_rd - 1.96 * se_vec
  ci_upper <- mean_rd + 1.96 * se_vec
  
  results <- tibble::tibble(
    t = seq_len(tau),
    control_no_tx = mean_s0,
    treat_tx = mean_s1,
    rd = mean_rd,
    se_RD = se_vec,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    ci_width = ci_upper - ci_lower
  )
  
  return(list(
    results = results,
    eps = c(eps_ctrl = eps_ctrl, eps_treat = eps_treat),
    m_fluc = m_fluc,
    fit_bart = bart_fit
  ))
}

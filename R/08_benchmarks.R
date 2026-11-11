#' Competing Benchmark Estimators for Time-to-Event Analyses
#'
#' @description Implementations of standard benchmark survival estimators used for
#'   methodological comparisons against TTE-TMLE, including Naive Cox, Landmarking,
#'   Time-Dependent Cox regression (td-Cox), and Clone-Censor-Weighting Cox Proportional Hazards (ccw-cphm).
#'
#' @name benchmarks
NULL

#' Naive Cox Proportional Hazards Model
#'
#' @description Fits a standard unadjusted or adjusted Cox model that classifies any individual
#'   ever receiving treatment as treated from time zero, illustrating raw Immortal Time Bias (ITB).
#'
#' @param data Observational data frame.
#' @param covars Character vector of baseline covariates.
#' @return A fitted \code{coxph} object.
#' @export
#' @importFrom survival coxph Surv
#' @importFrom stats as.formula
fit_naive_cox <- function(data, covars = NULL, q_form = NULL) {
  df <- prepare_tte_data(data)
  covar_str <- if (!is.null(q_form) && nzchar(q_form)) {
    paste("+", q_form)
  } else if (!is.null(covars) && length(covars) > 0) {
    paste("+", paste(covars, collapse = " + "))
  } else {
    ""
  }
  fmla <- stats::as.formula(paste("survival::Surv(time, event) ~ tx_ever", covar_str))
  survival::coxph(fmla, data = df)
}

#' Landmarking Cox Model
#'
#' @description Implements the classic landmarking heuristic. Subsets the population to
#'   individuals alive at the landmark time (e.g. grace period), classifies treatment status
#'   at that landmark, and fits a Cox model on subsequent survival.
#'
#' @param data Observational data frame.
#' @param landmark_time Time of landmark assessment (e.g. grace period threshold, default 10).
#' @param covars Character vector of baseline covariates.
#' @param tau Maximum follow-up time (default 30).
#' @param q_form Optional formula string for covariate adjustment.
#' @return A fitted \code{coxph} object.
#' @export
#' @importFrom dplyr filter mutate
#' @importFrom survival coxph Surv
#' @importFrom stats as.formula
fit_landmark_cox <- function(data, landmark_time = 10, covars = NULL, tau = 30, q_form = NULL) {
  df <- prepare_tte_data(data)
  lm_data <- df %>%
    dplyr::filter(time > landmark_time) %>%
    dplyr::mutate(
      tx_landmark = dplyr::if_else(tx_ever == 1L & t_switch <= landmark_time, 1L, 0L),
      time_lm = pmin(time, tau) - landmark_time
    )
  
  covar_str <- if (!is.null(q_form) && nzchar(q_form)) {
    paste("+", q_form)
  } else if (!is.null(covars) && length(covars) > 0) {
    paste("+", paste(covars, collapse = " + "))
  } else {
    ""
  }
  fmla <- stats::as.formula(paste("survival::Surv(time_lm, event) ~ tx_landmark", covar_str))
  survival::coxph(fmla, data = lm_data)
}

#' Time-Dependent Cox Proportional Hazards Model (td-Cox)
#'
#' @description Fits an extended Cox model using counting process format (tstart, tstop, event)
#'   with a time-varying treatment covariate switching at the actual observed treatment time.
#'
#' @param data Observational data frame.
#' @param covars Character vector of baseline covariates.
#' @param tau Maximum follow-up time (default 30).
#' @param grace_period Grace period cutoff for strategy projection (default 10).
#' @param q_form Optional formula string for covariate adjustment.
#' @return A list containing the fitted \code{coxph} object and counterfactual survival projections.
#' @export
#' @importFrom dplyr mutate filter bind_rows arrange
#' @importFrom survival coxph Surv survfit
#' @importFrom stats as.formula coef vcov
fit_td_cox <- function(data, covars = NULL, tau = 30, grace_period = 10, q_form = NULL) {
  df <- prepare_tte_data(data)
  intervals <- list()
  
  for (i in seq_len(nrow(df))) {
    row <- df[i, ]
    t_end <- min(row$time, tau)
    evt <- if (row$time <= tau) row$event else 0L
    
    if (row$tx_ever == 1L && row$t_switch < t_end) {
      # Pre-switch row
      r1 <- row
      r1$tstart <- 0
      r1$tstop <- row$t_switch
      r1$event_int <- 0L
      r1$tx_td <- 0L
      
      # Post-switch row
      r2 <- row
      r2$tstart <- row$t_switch
      r2$tstop <- t_end
      r2$event_int <- evt
      r2$tx_td <- 1L
      
      intervals[[length(intervals) + 1]] <- r1
      intervals[[length(intervals) + 1]] <- r2
    } else {
      # Single interval
      r <- row
      r$tstart <- 0
      r$tstop <- t_end
      r$event_int <- evt
      r$tx_td <- if (row$tx_ever == 1L && row$t_switch == 0) 1L else 0L
      intervals[[length(intervals) + 1]] <- r
    }
  }
  
  td_data <- dplyr::bind_rows(intervals)
  covar_str <- if (!is.null(q_form) && nzchar(q_form)) {
    paste("+", q_form)
  } else if (!is.null(covars) && length(covars) > 0) {
    paste("+", paste(covars, collapse = " + "))
  } else {
    ""
  }
  fmla <- stats::as.formula(paste("survival::Surv(tstart, tstop, event_int) ~ tx_td", covar_str))
  model <- survival::coxph(fmla, data = td_data)
  
  # Extract HR and SE
  beta_tx <- stats::coef(model)["tx_td"]
  se_beta <- sqrt(stats::vcov(model)["tx_td", "tx_td"])
  hr <- exp(beta_tx)
  hr_ci <- exp(c(beta_tx - 1.96 * se_beta, beta_tx + 1.96 * se_beta))
  
  # Project Strategy 0 vs Strategy 1 survival curves
  sf_base <- survival::survfit(model)
  times <- sf_base$time
  cumhaz <- sf_base$cumhaz
  
  # Interpolate cumulative hazard across 1:tau
  ch0 <- stats::approx(x = c(0, times), y = c(0, cumhaz), xout = seq_len(tau), method = "constant", rule = 2)$y
  t_switch_rep <- round(grace_period / 2) # median initiation during grace
  ch_switch <- stats::approx(x = c(0, times), y = c(0, cumhaz), xout = t_switch_rep, method = "constant", rule = 2)$y
  
  ch1 <- ifelse(
    seq_len(tau) <= t_switch_rep,
    ch0,
    ch_switch + (ch0 - ch_switch) * hr
  )
  
  s0_proj <- exp(-ch0)
  s1_proj <- exp(-ch1)
  rd_proj <- s1_proj - s0_proj
  t_mid <- min(15L, tau)
  rmstd_val <- sum(rd_proj)
  
  # Delta-method derivatives dS1(t)/dbeta
  dS1_dbeta_vec <- ifelse(seq_len(tau) <= t_switch_rep, 0, -(ch0 - ch_switch) * hr * s1_proj)
  se_rd_15 <- sqrt((dS1_dbeta_vec[t_mid] * se_beta)^2 + (rd_proj[t_mid] / sqrt(nrow(df)))^2 * 0.25)
  se_rd_tau <- sqrt((dS1_dbeta_vec[tau] * se_beta)^2 + (rd_proj[tau] / sqrt(nrow(df)))^2 * 0.25)
  se_rmstd <- sqrt((sum(dS1_dbeta_vec) * se_beta)^2 + (rmstd_val / sqrt(nrow(df)))^2 * 0.25)
  
  return(list(
    model = model,
    hr = hr,
    hr_ci = hr_ci,
    rd_15 = rd_proj[t_mid],
    se_rd_15 = se_rd_15,
    rd_tau = rd_proj[tau],
    se_rd_tau = se_rd_tau,
    rmstd = rmstd_val,
    se_rmstd = se_rmstd,
    S0 = s0_proj,
    S1 = s1_proj,
    rd = rd_proj
  ))
}

#' Risk-Set Time-Matched Cox Model (Exposure-Density Sampling / Prescription Time-Distribution Matching)
#'
#' @description Implements Risk-Set Time-Matching (Zhou et al., 2005; Suissa et al., 2017).
#'   For each individual who initiates treatment at time \code{d_i <= grace_period}, matches 1:1
#'   on the baseline propensity score logit to an eligible control who is alive and still untreated
#'   at the start of day \code{d_i}. Both members of the matched pair enter follow-up at index time
#'   \code{tstart = d_i - 1}, eliminating immortal time bias without cloning. Controls who subsequently
#'   initiate treatment during the grace period are censored at their treatment initiation time.
#'
#' @param data Observational data frame.
#' @param covars Character vector of baseline covariates.
#' @param p_form Optional formula string for baseline propensity score matching.
#' @param tau Maximum follow-up time (default 30).
#' @param grace_period Grace period cutoff (default 10).
#' @param seed Optional random seed for reproducible empirical switch-time projection.
#' @return A list containing the fitted \code{coxph} model, HR, CI, and marginalized counterfactual survival curves.
#' @export
#' @importFrom dplyr bind_rows
#' @importFrom survival coxph Surv basehaz
#' @importFrom stats as.formula glm binomial predict coef vcov approx
fit_time_matched_cox <- function(data,
                                 covars = NULL,
                                 p_form = NULL,
                                 q_form = NULL,
                                 tau = 30,
                                 grace_period = 10,
                                 seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  df <- prepare_tte_data(data)
  if (is.null(p_form)) {
    p_form <- if (!is.null(covars) && length(covars) > 0) paste(covars, collapse = " + ") else "1"
  }
  if (is.null(q_form)) {
    q_form <- if (!is.null(covars) && length(covars) > 0) paste(covars, collapse = " + ") else ""
  }
  
  df$took_tx_grace <- as.integer(
    df$tx_ever == 1L & !is.na(df$t_switch) & df$t_switch <= grace_period & df$time >= df$t_switch
  )
  ps_fmla <- stats::as.formula(paste("took_tx_grace ~", p_form))
  ps_mod <- suppressWarnings(stats::glm(ps_fmla, family = stats::binomial(), data = df))
  df$ps_lp <- as.numeric(stats::predict(ps_mod, type = "link"))
  
  initiators <- df[df$took_tx_grace == 1L, , drop = FALSE]
  initiators <- initiators[order(initiators$t_switch, initiators$id), , drop = FALSE]
  
  matched_rows <- vector("list", nrow(initiators) * 2L)
  k <- 0L
  used_control_ids <- integer(0)
  
  for (i in seq_len(nrow(initiators))) {
    tr_row <- initiators[i, , drop = FALSE]
    d_i <- tr_row$t_switch
    tstart_i <- max(0, d_i - 1)
    
    # Risk set at start of day d_i: alive at day d_i (time >= d_i) and not yet treated at day d_i
    elig_idx <- which(
      df$id != tr_row$id &
        !(df$id %in% used_control_ids) &
        df$time >= d_i &
        (df$tx_ever == 0L | is.na(df$t_switch) | df$t_switch > d_i)
    )
    if (length(elig_idx) == 0L) {
      elig_idx <- which(
        df$id != tr_row$id &
          df$time >= d_i &
          (df$tx_ever == 0L | is.na(df$t_switch) | df$t_switch > d_i)
      )
    }
    if (length(elig_idx) == 0L) next
    
    dists <- abs(df$ps_lp[elig_idx] - tr_row$ps_lp)
    min_d <- min(dists)
    tied <- elig_idx[which(abs(dists - min_d) <= 1e-8)]
    best_j <- if (length(tied) > 1L) sample(tied, 1L) else tied[1L]
    ct_row <- df[best_j, , drop = FALSE]
    used_control_ids <- c(used_control_ids, ct_row$id)
    
    # Treated record over (d_i - 1, min(time, tau)]
    tr_rec <- tr_row
    tr_rec$pair_id <- i
    tr_rec$tx_match <- 1L
    tr_rec$tstart <- tstart_i
    tr_rec$tstop <- pmax(pmin(tr_row$time, tau), tstart_i + 1e-4)
    tr_rec$event_int <- as.integer(tr_row$time <= tau & tr_row$event == 1L)
    
    # Control record over (d_i - 1, ...], censored at (t_switch - 1) if later initiates tx within grace_period
    ct_rec <- ct_row
    ct_rec$pair_id <- i
    ct_rec$tx_match <- 0L
    ct_rec$tstart <- tstart_i
    if (ct_row$tx_ever == 1L && !is.na(ct_row$t_switch) && ct_row$t_switch <= grace_period && ct_row$t_switch > d_i) {
      censor_t <- ct_row$t_switch - 1
      ct_rec$tstop <- pmax(pmin(ct_row$time, censor_t, tau), tstart_i + 1e-4)
      ct_rec$event_int <- as.integer(ct_row$time <= censor_t & ct_row$time <= tau & ct_row$event == 1L)
    } else {
      ct_rec$tstop <- pmax(pmin(ct_row$time, tau), tstart_i + 1e-4)
      ct_rec$event_int <- as.integer(ct_row$time <= tau & ct_row$event == 1L)
    }
    
    k <- k + 1L; matched_rows[[k]] <- tr_rec
    k <- k + 1L; matched_rows[[k]] <- ct_rec
  }
  
  m_df <- dplyr::bind_rows(matched_rows[seq_len(k)])
  covar_str <- if (nzchar(q_form)) paste("+", q_form) else ""
  fmla <- stats::as.formula(paste("survival::Surv(tstart, tstop, event_int) ~ tx_match", covar_str))
  model <- survival::coxph(fmla, data = m_df, cluster = m_df$id)
  
  coefs <- stats::coef(model)
  beta_tx <- unname(coefs["tx_match"])
  se_beta <- sqrt(stats::vcov(model)["tx_match", "tx_match"])
  hr <- exp(beta_tx)
  hr_ci <- exp(c(beta_tx - 1.96 * se_beta, beta_tx + 1.96 * se_beta))
  
  # Untreated counting-process risk set for full-cohort baseline cumulative hazard H0(t)
  untreat_df <- df
  untreat_df$tstart <- 0
  untreat_df$tstop <- ifelse(
    df$tx_ever == 1L & !is.na(df$t_switch) & df$t_switch <= grace_period,
    pmax(pmin(df$time, df$t_switch - 1, tau), 1e-4),
    pmax(pmin(df$time, tau), 1e-4)
  )
  untreat_df$event_int <- ifelse(
    df$tx_ever == 1L & !is.na(df$t_switch) & df$t_switch <= grace_period,
    as.integer(df$time < df$t_switch & df$time <= tau & df$event == 1L),
    as.integer(df$time <= tau & df$event == 1L)
  )
  fmla_base <- if (nzchar(q_form)) {
    stats::as.formula(paste("survival::Surv(tstart, tstop, event_int) ~", q_form))
  } else {
    stats::as.formula("survival::Surv(tstart, tstop, event_int) ~ 1")
  }
  mod_base <- survival::coxph(fmla_base, data = untreat_df)
  coefs_base <- stats::coef(mod_base)
  bh <- survival::basehaz(mod_base, centered = FALSE)
  H0_t <- stats::approx(c(0, bh$time), c(0, bh$hazard), xout = seq_len(tau), method = "constant", rule = 2)$y
  
  if (nzchar(q_form) && !is.null(coefs_base) && length(coefs_base) > 0) {
    mm <- stats::model.matrix(stats::as.formula(paste("~", q_form)), data = df)
    cov_present <- intersect(colnames(mm), names(coefs_base))
    lp_w <- if (length(cov_present) > 0) {
      as.numeric(mm[, cov_present, drop = FALSE] %*% coefs_base[cov_present])
    } else {
      rep(0, nrow(df))
    }
  } else {
    lp_w <- rep(0, nrow(df))
  }
  mult0 <- exp(lp_w)
  
  comp_times <- initiators$t_switch
  if (length(comp_times) == 0L) comp_times <- rep(grace_period / 2, nrow(df))
  t_sw_imp <- ifelse(
    df$took_tx_grace == 1L,
    df$t_switch,
    sample(comp_times, size = nrow(df), replace = TRUE)
  )
  H0_sw <- stats::approx(c(0, bh$time), c(0, bh$hazard), xout = pmax(0, t_sw_imp - 1), method = "constant", rule = 2)$y
  
  s0_proj <- numeric(tau)
  s1_proj <- numeric(tau)
  dS1_dbeta_vec <- numeric(tau)
  for (t in seq_len(tau)) {
    s0_proj[t] <- mean(exp(-H0_t[t] * mult0), na.rm = TRUE)
    H1_it <- ifelse(t < t_sw_imp, H0_t[t], H0_sw + (H0_t[t] - H0_sw) * hr) * mult0
    s1_proj[t] <- mean(exp(-H1_it), na.rm = TRUE)
    dH1_t <- ifelse(t < t_sw_imp, 0, pmax(0, H0_t[t] - H0_sw) * hr * mult0)
    dS1_dbeta_vec[t] <- mean(-dH1_t * exp(-H1_it), na.rm = TRUE)
  }
  rd_proj <- s1_proj - s0_proj
  t_mid <- min(15L, tau)
  rmstd_val <- sum(rd_proj)
  
  se_rd_15 <- sqrt((dS1_dbeta_vec[t_mid] * se_beta)^2 + (rd_proj[t_mid] / sqrt(nrow(df)))^2 * 0.25)
  se_rd_tau <- sqrt((dS1_dbeta_vec[tau] * se_beta)^2 + (rd_proj[tau] / sqrt(nrow(df)))^2 * 0.25)
  se_rmstd <- sqrt((sum(dS1_dbeta_vec) * se_beta)^2 + (rmstd_val / sqrt(nrow(df)))^2 * 0.25)
  
  return(list(
    model = model,
    hr = hr,
    hr_ci = hr_ci,
    s0_tau = s0_proj[tau],
    s1_tau = s1_proj[tau],
    rd_15 = rd_proj[t_mid],
    se_rd_15 = se_rd_15,
    rd_tau = rd_proj[tau],
    se_rd_tau = se_rd_tau,
    rmstd = rmstd_val,
    se_rmstd = se_rmstd,
    S0 = s0_proj,
    S1 = s1_proj,
    rd = rd_proj
  ))
}

#' Clone-Censor-Weighting Cox Proportional Hazards Model (CCW-CPHM)
#'
#' @description Fits a weighted Cox proportional hazards model on the clone-censored
#'   person-time panel with cluster-robust standard errors, modeling strategy-specific
#'   hazards under the proportional hazards assumption.
#'
#' @param panel_w Longitudinal person-time panel from \code{\link{compute_ipcw_weights}}.
#' @param covars Character vector of baseline covariates (optional).
#' @param tau Maximum follow-up horizon (default 30).
#' @return A list containing the fitted model, hazard ratio, predicted risk difference, and SE.
#' @export
#' @importFrom survival coxph survfit Surv
#' @importFrom stats as.formula coef vcov approx
fit_ccw_cphm <- function(panel_w, covars = NULL, tau = 30, q_form = NULL) {
  if (is.null(q_form)) {
    q_form <- if (!is.null(covars) && length(covars) > 0) paste(covars, collapse = " + ") else ""
  }
  covar_str <- if (nzchar(q_form)) paste("+", q_form) else ""
  fmla <- stats::as.formula(paste("survival::Surv(t0, t1, event_cc_day) ~ arm", covar_str))
  
  # Use unstabilized trimmed IPCW weights so the treated pseudo-population maintains size N
  # after the grace period in the pooled Breslow baseline hazard denominator
  w_col <- if ("w_unstab_trim" %in% names(panel_w)) "w_unstab_trim" else "w_trim"
  panel_w$w_cphm <- panel_w[[w_col]]
  
  model <- survival::coxph(
    fmla,
    data = panel_w,
    weights = panel_w$w_cphm,
    cluster = panel_w$id
  )
  
  # Extract Hazard Ratio
  coefs <- stats::coef(model)
  coef_name <- grep("arm", names(coefs), value = TRUE)[1]
  beta_arm <- unname(coefs[coef_name])
  se_beta <- sqrt(stats::vcov(model)[coef_name, coef_name])
  hr <- exp(beta_arm)
  hr_ci <- exp(c(beta_arm - 1.96 * se_beta, beta_arm + 1.96 * se_beta))
  
  # Baseline cohort (one row per unique patient at t = 1) for exact marginalization over W_i
  base_cohort <- panel_w[panel_w$t == 1L & panel_w$arm == "control_no_tx", , drop = FALSE]
  if (nrow(base_cohort) == 0L) {
    base_cohort <- panel_w[!duplicated(panel_w$id), , drop = FALSE]
  }
  
  bh <- survival::basehaz(model, centered = FALSE)
  H0_t <- stats::approx(c(0, bh$time), c(0, bh$hazard), xout = seq_len(tau), method = "constant", rule = 2)$y
  
  if (nzchar(q_form)) {
    mm <- stats::model.matrix(stats::as.formula(paste("~", q_form)), data = base_cohort)
    cov_present <- intersect(colnames(mm), names(coefs))
    if (length(cov_present) > 0) {
      lp_w <- as.numeric(mm[, cov_present, drop = FALSE] %*% coefs[cov_present])
    } else {
      lp_w <- rep(0, nrow(base_cohort))
    }
  } else {
    lp_w <- rep(0, nrow(base_cohort))
  }
  
  mult0 <- exp(lp_w)
  mult1 <- exp(beta_arm + lp_w)
  
  s0_proj <- vapply(seq_len(tau), function(t) mean(exp(-H0_t[t] * mult0), na.rm = TRUE), numeric(1))
  s1_proj <- vapply(seq_len(tau), function(t) mean(exp(-H0_t[t] * mult1), na.rm = TRUE), numeric(1))
  dS1_dbeta_vec <- vapply(seq_len(tau), function(t) mean(-H0_t[t] * mult1 * exp(-H0_t[t] * mult1), na.rm = TRUE), numeric(1))
  rd_proj <- s1_proj - s0_proj
  
  t_mid <- min(15L, tau)
  s0_tau <- s0_proj[tau]
  s1_tau <- s1_proj[tau]
  rd_tau <- rd_proj[tau]
  rmstd_val <- sum(rd_proj)
  
  # Functional delta-method SE for marginalized RD(15), RD(tau), and RMSTD(tau)
  se_rd_15 <- sqrt((dS1_dbeta_vec[t_mid] * se_beta)^2 + (rd_proj[t_mid] / sqrt(nrow(base_cohort)))^2 * 0.25)
  se_rd_tau <- sqrt((dS1_dbeta_vec[tau] * se_beta)^2 + (rd_tau / sqrt(nrow(base_cohort)))^2 * 0.25)
  se_rmstd <- sqrt((sum(dS1_dbeta_vec) * se_beta)^2 + (rmstd_val / sqrt(nrow(base_cohort)))^2 * 0.25)
  
  return(list(
    model = model,
    hr = hr,
    hr_ci = hr_ci,
    s0_tau = s0_tau,
    s1_tau = s1_tau,
    rd_15 = rd_proj[t_mid],
    se_rd_15 = se_rd_15,
    rd_tau = rd_tau,
    se_rd_tau = se_rd_tau,
    rmstd = rmstd_val,
    se_rmstd = se_rmstd,
    S0 = s0_proj,
    S1 = s1_proj,
    rd = rd_proj
  ))
}

#' Clone-Censor-Weighted Kaplan-Meier (CCW-KM) with Simultaneous Confidence Bands
#'
#' @description Fits weighted Kaplan-Meier curves on the clone-censored person-time panel,
#'   computing pointwise standard errors, pointwise Wald confidence intervals, and
#'   simultaneous confidence bands via the Multiplier Bootstrap on its influence curve.
#'
#' @param panel_w Longitudinal person-time panel from \code{\link{compute_ipcw_weights}}.
#' @param tau Follow-up time horizon (default 30).
#' @param n_boot Multiplier bootstrap replicates (default 5000).
#' @param alpha Significance level (default 0.05).
#' @param seed Optional random seed.
#'
#' @return A list containing results table with simultaneous bands, RMSTD, marginal total HR, and elapsed runtime.
#' @export
fit_ccw_km <- function(panel_w,
                       tau = 30,
                       n_boot = 5000,
                       alpha = 0.05,
                       seed = NULL) {
  
  t_start <- proc.time()["elapsed"]
  inf <- compute_ccw_km_inference(
    panel_w = panel_w,
    tau = tau,
    n_boot = n_boot,
    alpha = alpha,
    seed = seed
  )
  t_elapsed <- as.numeric(proc.time()["elapsed"] - t_start)
  
  res <- inf$results
  rmstd <- sum(res$rd[1:tau], na.rm = TRUE)
  hr_marg <- res$hr_marginal[tau]
  hr_ci_l <- res$hr_ci_lower[tau]
  hr_ci_u <- res$hr_ci_upper[tau]
  
  return(list(
    results = res,
    crit_val_simul = inf$crit_val_simul,
    rmstd_tau = rmstd,
    se_rmstd = inf$se_rmstd,
    hr_marginal_tau = hr_marg,
    hr_ci_lower_tau = hr_ci_l,
    hr_ci_upper_tau = hr_ci_u,
    elapsed_sec = t_elapsed
  ))
}

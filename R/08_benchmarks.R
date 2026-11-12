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
#' @description Fits a naive baseline Cox model that classifies any individual ever
#'   receiving treatment as treated from time zero, illustrating raw Immortal Time Bias (ITB).
#'   Provides only the naive conditional Cox hazard ratio; does not define valid baseline
#'   marginal survival differences or RMSTD.
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
#' @description Implements classic landmarking by subsetting to individuals alive at
#'   \code{landmark_time} (e.g. end of the grace period) and modeling subsequent survival.
#'   Because individuals with events on or before \code{landmark_time} are excluded,
#'   this estimates a post-landmark conditional hazard ratio rather than baseline
#'   marginal survival differences or 30-day RMSTD from time zero.
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
#'   with a time-varying treatment covariate switching at the observed treatment time.
#'   Returns the conditional instantaneous hazard ratio and Wald CI; static grace-period
#'   survival differences and RMSTD are returned as \code{NA_real_} because counting-process
#'   Cox regression does not natively identify a time-zero sustained treatment strategy.
#'
#' @param data Observational data frame.
#' @param covars Character vector of baseline covariates.
#' @param tau Maximum follow-up time (default 30).
#' @param grace_period Unused for point estimation (kept for signature compatibility).
#' @param q_form Optional formula string for covariate adjustment.
#' @return A list containing the fitted \code{coxph} object, conditional HR/CI, and \code{NA_real_}
#'   placeholders for non-identifiable baseline marginal survival metrics.
#' @export
#' @importFrom dplyr bind_rows
#' @importFrom survival coxph Surv
#' @importFrom stats as.formula coef vcov
fit_td_cox <- function(data, covars = NULL, tau = 30, grace_period = 10, q_form = NULL) {
  df <- prepare_tte_data(data)
  intervals <- list()

  for (i in seq_len(nrow(df))) {
    row <- df[i, ]
    t_end <- min(row$time, tau)
    evt <- if (row$time <= tau) row$event else 0L

    if (row$tx_ever == 1L && row$t_switch < t_end) {
      r1 <- row
      r1$tstart <- 0
      r1$tstop <- row$t_switch
      r1$event_int <- 0L
      r1$tx_td <- 0L

      r2 <- row
      r2$tstart <- row$t_switch
      r2$tstop <- t_end
      r2$event_int <- evt
      r2$tx_td <- 1L

      intervals[[length(intervals) + 1]] <- r1
      intervals[[length(intervals) + 1]] <- r2
    } else {
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

  beta_tx <- unname(stats::coef(model)["tx_td"])
  se_beta <- sqrt(stats::vcov(model)["tx_td", "tx_td"])
  hr_cond <- exp(beta_tx)
  hr_cond_ci <- exp(c(beta_tx - 1.96 * se_beta, beta_tx + 1.96 * se_beta))

  list(
    model = model,
    hr = hr_cond,
    hr_ci = hr_cond_ci,
    hr_cond = hr_cond,
    hr_cond_ci = hr_cond_ci,
    rd_15 = NA_real_,
    se_rd_15 = NA_real_,
    rd_tau = NA_real_,
    se_rd_tau = NA_real_,
    rmstd = NA_real_,
    se_rmstd = NA_real_,
    S0 = rep(NA_real_, tau),
    S1 = rep(NA_real_, tau),
    rd = rep(NA_real_, tau)
  )
}

# Internal single-pass fitter for Risk-Set Time-Matched Cox + G-standardized marginal curves
.fit_time_matched_point <- function(df,
                                    p_form,
                                    q_form,
                                    tau = 30,
                                    grace_period = 10,
                                    robust_vcov = TRUE) {
  n <- nrow(df)
  df$took_tx_grace <- as.integer(
    df$tx_ever == 1L & !is.na(df$t_switch) & df$t_switch <= grace_period & df$time >= df$t_switch
  )

  ps_fmla <- stats::as.formula(paste("took_tx_grace ~", p_form))
  ps_mod <- suppressWarnings(stats::glm(ps_fmla, family = stats::binomial(), data = df))
  ps_lp <- as.numeric(stats::predict(ps_mod, type = "link"))

  init_pos <- which(df$took_tx_grace == 1L)
  if (length(init_pos) == 0L) return(NULL)
  init_pos <- init_pos[order(df$t_switch[init_pos], df$id[init_pos])]
  n_init <- length(init_pos)

  avail <- rep(TRUE, n)
  src_idx <- integer(2L * n_init)
  pair_id_vec <- integer(2L * n_init)
  tx_match_vec <- integer(2L * n_init)
  tstart_vec <- numeric(2L * n_init)
  tstop_vec <- numeric(2L * n_init)
  event_int_vec <- integer(2L * n_init)
  k <- 0L

  time_v <- df$time
  event_v <- df$event
  tx_ever_v <- df$tx_ever
  tsw_v <- df$t_switch

  elig_by_day <- lapply(seq_len(grace_period), function(d) {
    time_v >= d & (tx_ever_v == 0L | is.na(tsw_v) | tsw_v > d)
  })

  for (i in seq_len(n_init)) {
    tr_pos <- init_pos[i]
    d_i <- max(1L, min(as.integer(grace_period), as.integer(tsw_v[tr_pos])))
    tstart_i <- max(0, d_i - 1)

    elig_day <- elig_by_day[[d_i]]
    elig_day[tr_pos] <- FALSE

    elig_idx <- which(avail & elig_day)
    if (length(elig_idx) == 0L) {
      elig_idx <- which(elig_day)
    }
    if (length(elig_idx) == 0L) next

    dists <- abs(ps_lp[elig_idx] - ps_lp[tr_pos])
    min_d <- min(dists)
    tied <- elig_idx[abs(dists - min_d) <= 1e-8]
    best_j <- if (length(tied) > 1L) tied[sample.int(length(tied), 1L)] else tied[1L]
    avail[best_j] <- FALSE

    # Treated record over (d_i - 1, min(time, tau)]
    k <- k + 1L
    src_idx[k] <- tr_pos
    pair_id_vec[k] <- i
    tx_match_vec[k] <- 1L
    tstart_vec[k] <- tstart_i
    tstop_vec[k] <- max(min(time_v[tr_pos], tau), tstart_i + 1e-4)
    event_int_vec[k] <- as.integer(time_v[tr_pos] <= tau && event_v[tr_pos] == 1L)

    # Control record over (d_i - 1, ...], censored at (t_switch - 1) if later initiates tx within grace_period
    k <- k + 1L
    src_idx[k] <- best_j
    pair_id_vec[k] <- i
    tx_match_vec[k] <- 0L
    tstart_vec[k] <- tstart_i
    if (tx_ever_v[best_j] == 1L && !is.na(tsw_v[best_j]) && tsw_v[best_j] <= grace_period && tsw_v[best_j] > d_i) {
      censor_t <- tsw_v[best_j] - 1
      tstop_vec[k] <- max(min(time_v[best_j], censor_t, tau), tstart_i + 1e-4)
      event_int_vec[k] <- as.integer(time_v[best_j] <= censor_t && time_v[best_j] <= tau && event_v[best_j] == 1L)
    } else {
      tstop_vec[k] <- max(min(time_v[best_j], tau), tstart_i + 1e-4)
      event_int_vec[k] <- as.integer(time_v[best_j] <= tau && event_v[best_j] == 1L)
    }
  }

  if (k < 2L) return(NULL)
  seq_k <- seq_len(k)
  m_df <- df[src_idx[seq_k], , drop = FALSE]
  m_df$pair_id <- pair_id_vec[seq_k]
  m_df$tx_match <- tx_match_vec[seq_k]
  m_df$tstart <- tstart_vec[seq_k]
  m_df$tstop <- tstop_vec[seq_k]
  m_df$event_int <- event_int_vec[seq_k]

  if (length(unique(m_df$tx_match)) < 2L || sum(m_df$event_int) == 0L) return(NULL)

  covar_str <- if (nzchar(q_form)) paste("+", q_form) else ""
  fmla <- stats::as.formula(paste("survival::Surv(tstart, tstop, event_int) ~ tx_match", covar_str))
  model <- if (robust_vcov) {
    suppressWarnings(survival::coxph(fmla, data = m_df, cluster = m_df$id))
  } else {
    suppressWarnings(survival::coxph(fmla, data = m_df))
  }

  coefs <- stats::coef(model)
  beta_tx <- unname(coefs["tx_match"])
  if (is.na(beta_tx) || !is.finite(beta_tx)) return(NULL)
  hr_cond <- exp(beta_tx)

  # Untreated counting-process risk set for full-cohort baseline cumulative hazard H0(t)
  untreat_df <- df
  untreat_df$tstart <- 0
  sw_flag <- tx_ever_v == 1L & !is.na(tsw_v) & tsw_v <= grace_period
  untreat_df$tstop <- ifelse(sw_flag, pmax(pmin(time_v, tsw_v - 1, tau), 1e-4), pmax(pmin(time_v, tau), 1e-4))
  untreat_df$event_int <- ifelse(
    sw_flag,
    as.integer(time_v < tsw_v & time_v <= tau & event_v == 1L),
    as.integer(time_v <= tau & event_v == 1L)
  )
  fmla_base <- if (nzchar(q_form)) {
    stats::as.formula(paste("survival::Surv(tstart, tstop, event_int) ~", q_form))
  } else {
    stats::as.formula("survival::Surv(tstart, tstop, event_int) ~ 1")
  }
  mod_base <- suppressWarnings(survival::coxph(fmla_base, data = untreat_df))
  coefs_base <- stats::coef(mod_base)
  bh <- survival::basehaz(mod_base, centered = FALSE)
  H0_t <- stats::approx(c(0, bh$time), c(0, bh$hazard), xout = seq_len(tau), method = "constant", rule = 2)$y

  if (nzchar(q_form) && !is.null(coefs_base) && length(coefs_base) > 0) {
    coefs_base[is.na(coefs_base)] <- 0
    mm <- stats::model.matrix(stats::as.formula(paste("~", q_form)), data = df)
    cov_present <- intersect(colnames(mm), names(coefs_base))
    lp_w <- if (length(cov_present) > 0) {
      as.numeric(mm[, cov_present, drop = FALSE] %*% coefs_base[cov_present])
    } else {
      rep(0, n)
    }
  } else {
    lp_w <- rep(0, n)
  }
  mult0 <- exp(lp_w)

  # Analytical integration across empirical grace-period switch-day distribution
  comp_days <- tsw_v[init_pos]
  tab_d <- table(comp_days)
  d_vals <- as.integer(names(tab_d))
  p_vals <- as.numeric(tab_d) / sum(tab_d)
  H0_d <- stats::approx(c(0, bh$time), c(0, bh$hazard), xout = pmax(0, d_vals - 1), method = "constant", rule = 2)$y

  t_seq <- seq_len(tau)
  s0_proj <- colMeans(exp(-outer(mult0, H0_t, "*")), na.rm = TRUE)

  init_mask <- df$took_tx_grace == 1L
  has_noninit <- any(!init_mask)
  s1_mat <- matrix(0, nrow = n, ncol = tau)

  for (j in seq_along(d_vals)) {
    d <- d_vals[j]
    H_d <- H0_d[j]
    h_vec_d <- ifelse(t_seq < d, H0_t, H_d + (H0_t - H_d) * hr_cond)
    s_d <- exp(-outer(mult0, h_vec_d, "*"))

    match_j <- init_mask & (tsw_v == d)
    if (any(match_j)) {
      s1_mat[match_j, ] <- s_d[match_j, , drop = FALSE]
    }
    if (has_noninit) {
      s1_mat[!init_mask, ] <- s1_mat[!init_mask, , drop = FALSE] + p_vals[j] * s_d[!init_mask, , drop = FALSE]
    }
  }

  s1_proj <- colMeans(s1_mat, na.rm = TRUE)
  rd_proj <- s1_proj - s0_proj

  s0_tau_c <- min(max(s0_proj[tau], 1e-6), 1 - 1e-6)
  s1_tau_c <- min(max(s1_proj[tau], 1e-6), 1 - 1e-6)
  hr_marg <- log(s1_tau_c) / log(s0_tau_c)
  t_mid <- min(15L, tau)

  out <- list(
    rd_15 = rd_proj[t_mid],
    rd_tau = rd_proj[tau],
    rmstd = sum(rd_proj),
    hr_marg = hr_marg,
    hr_cond = hr_cond,
    s0_tau = s0_proj[tau],
    s1_tau = s1_proj[tau],
    S0 = s0_proj,
    S1 = s1_proj,
    rd = rd_proj
  )
  if (robust_vcov) {
    se_beta <- sqrt(stats::vcov(model)["tx_match", "tx_match"])
    out$model <- model
    out$hr_cond_ci <- exp(c(beta_tx - 1.96 * se_beta, beta_tx + 1.96 * se_beta))
  }
  out
}

#' Risk-Set Time-Matched Cox Model (Exposure-Density Sampling / Prescription Time-Distribution Matching)
#'
#' @description Implements Risk-Set Time-Matching (Zhou et al., 2005; Suissa et al., 2017).
#'   For each individual who initiates treatment at time \code{d_i <= grace_period}, matches 1:1
#'   on the baseline propensity score logit to an eligible control who is alive and still untreated
#'   at the start of day \code{d_i}. Both members of the matched pair enter follow-up at index time
#'   \code{tstart = d_i - 1}, eliminating immortal time bias without cloning. Controls who subsequently
#'   initiate treatment during the grace period are censored at their treatment initiation time.
#'   Counterfactual marginal survival curves (\code{S0}, \code{S1}), \code{RD(t)}, \code{RMSTD},
#'   and cumulative \code{Marginal HR(tau)} are obtained via baseline untreated-hazard standardization
#'   integrated deterministically over the empirical initiation-time distribution, with patient-level
#'   nonparametric bootstrap standard errors.
#'
#' @param data Observational data frame.
#' @param covars Character vector of baseline covariates.
#' @param p_form Optional formula string for baseline propensity score matching.
#' @param q_form Optional formula string for outcome Cox regression adjustment.
#' @param tau Maximum follow-up time (default 30).
#' @param grace_period Grace period cutoff (default 10).
#' @param n_boot Number of patient-level nonparametric bootstrap replicates for inference (default 200).
#' @param seed Optional random seed for reproducible tie-breaking and bootstrap inference.
#' @return A list containing the fitted \code{coxph} model, marginal and conditional HRs, CIs,
#'   G-standardized counterfactual survival curves, and bootstrap standard errors.
#' @export
#' @importFrom survival coxph Surv basehaz
#' @importFrom stats as.formula glm binomial predict coef vcov approx sd
fit_time_matched_cox <- function(data,
                                 covars = NULL,
                                 p_form = NULL,
                                 q_form = NULL,
                                 tau = 30,
                                 grace_period = 10,
                                 n_boot = 200L,
                                 seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  df <- prepare_tte_data(data)
  if (is.null(p_form)) {
    p_form <- if (!is.null(covars) && length(covars) > 0) paste(covars, collapse = " + ") else "1"
  }
  if (is.null(q_form)) {
    q_form <- if (!is.null(covars) && length(covars) > 0) paste(covars, collapse = " + ") else ""
  }

  pt <- .fit_time_matched_point(
    df = df,
    p_form = p_form,
    q_form = q_form,
    tau = tau,
    grace_period = grace_period,
    robust_vcov = TRUE
  )
  if (is.null(pt)) {
    stop("fit_time_matched_cox failed to construct matched risk sets.")
  }

  n <- nrow(df)
  n_boot <- as.integer(n_boot)
  if (n_boot > 0L) {
    boot_mat <- matrix(NA_real_, nrow = n_boot, ncol = 4L)
    for (b in seq_len(n_boot)) {
      b_idx <- sample.int(n, size = n, replace = TRUE)
      b_df <- df[b_idx, , drop = FALSE]
      b_df$id <- seq_len(n)
      res_b <- tryCatch(
        .fit_time_matched_point(
          df = b_df,
          p_form = p_form,
          q_form = q_form,
          tau = tau,
          grace_period = grace_period,
          robust_vcov = FALSE
        ),
        error = function(e) NULL
      )
      if (!is.null(res_b) && is.finite(res_b$hr_marg) && res_b$hr_marg > 0) {
        boot_mat[b, ] <- c(res_b$rd_15, res_b$rd_tau, res_b$rmstd, log(res_b$hr_marg))
      }
    }
    valid_n <- sum(stats::complete.cases(boot_mat))
    se_rd_15 <- if (valid_n >= 2L) stats::sd(boot_mat[, 1L], na.rm = TRUE) else NA_real_
    se_rd_tau <- if (valid_n >= 2L) stats::sd(boot_mat[, 2L], na.rm = TRUE) else NA_real_
    se_rmstd <- if (valid_n >= 2L) stats::sd(boot_mat[, 3L], na.rm = TRUE) else NA_real_
    se_log_hr <- if (valid_n >= 2L) stats::sd(boot_mat[, 4L], na.rm = TRUE) else NA_real_
    hr_ci <- if (is.finite(se_log_hr) && is.finite(pt$hr_marg) && pt$hr_marg > 0) {
      exp(log(pt$hr_marg) + c(-1.96, 1.96) * se_log_hr)
    } else {
      c(NA_real_, NA_real_)
    }
  } else {
    valid_n <- 0L
    se_rd_15 <- NA_real_
    se_rd_tau <- NA_real_
    se_rmstd <- NA_real_
    se_log_hr <- NA_real_
    hr_ci <- c(NA_real_, NA_real_)
  }

  list(
    model = pt$model,
    hr = pt$hr_marg,
    hr_ci = hr_ci,
    se_log_hr = se_log_hr,
    hr_cond = pt$hr_cond,
    hr_cond_ci = pt$hr_cond_ci,
    s0_tau = pt$s0_tau,
    s1_tau = pt$s1_tau,
    rd_15 = pt$rd_15,
    se_rd_15 = se_rd_15,
    rd_tau = pt$rd_tau,
    se_rd_tau = se_rd_tau,
    rmstd = pt$rmstd,
    se_rmstd = se_rmstd,
    S0 = pt$S0,
    S1 = pt$S1,
    rd = pt$rd,
    n_boot_valid = valid_n
  )
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

  # Baseline cohort (one row per unique patient at t = 1) extracted before zero-weight filtering
  base_cohort <- panel_w[panel_w$t == 1L & panel_w$arm == "control_no_tx", , drop = FALSE]
  if (nrow(base_cohort) == 0L) {
    base_cohort <- panel_w[!duplicated(panel_w$id), , drop = FALSE]
  }

  # coxph requires strictly positive weights; drop structural zero-weight rows (e.g. end-of-grace non-compliers)
  panel_fit <- panel_w[panel_w$w_cphm > 0, , drop = FALSE]

  model <- survival::coxph(
    fmla,
    data = panel_fit,
    weights = panel_fit$w_cphm,
    cluster = panel_fit$id
  )
  
  # Extract Hazard Ratio
  coefs <- stats::coef(model)
  coef_name <- grep("arm", names(coefs), value = TRUE)[1]
  beta_arm <- unname(coefs[coef_name])
  se_beta <- sqrt(stats::vcov(model)[coef_name, coef_name])
  hr <- exp(beta_arm)
  hr_ci <- exp(c(beta_arm - 1.96 * se_beta, beta_arm + 1.96 * se_beta))
  
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

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
fit_naive_cox <- function(data, covars = NULL) {
  df <- prepare_tte_data(data)
  covar_str <- if (!is.null(covars) && length(covars) > 0) paste("+", paste(covars, collapse = " + ")) else ""
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
#' @return A fitted \code{coxph} object.
#' @export
#' @importFrom dplyr filter mutate
#' @importFrom survival coxph Surv
#' @importFrom stats as.formula
fit_landmark_cox <- function(data, landmark_time = 10, covars = NULL, tau = 30) {
  df <- prepare_tte_data(data)
  lm_data <- df %>%
    dplyr::filter(time > landmark_time) %>%
    dplyr::mutate(
      tx_landmark = dplyr::if_else(tx_ever == 1L & t_switch <= landmark_time, 1L, 0L),
      time_lm = pmin(time, tau) - landmark_time
    )
  
  covar_str <- if (!is.null(covars) && length(covars) > 0) paste("+", paste(covars, collapse = " + ")) else ""
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
#' @return A list containing the fitted \code{coxph} object and counterfactual survival projections.
#' @export
#' @importFrom dplyr mutate filter bind_rows arrange
#' @importFrom survival coxph Surv survfit
#' @importFrom stats as.formula coef vcov
fit_td_cox <- function(data, covars = NULL, tau = 30, grace_period = 10) {
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
  covar_str <- if (!is.null(covars) && length(covars) > 0) paste("+", paste(covars, collapse = " + ")) else ""
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
  
  # Approximate SE via delta method
  se_rd_tau <- abs(rd_proj[tau]) * se_beta
  
  return(list(
    model = model,
    hr = hr,
    hr_ci = hr_ci,
    rd_tau = rd_proj[tau],
    se_rd_tau = se_rd_tau,
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
fit_ccw_cphm <- function(panel_w, covars = NULL, tau = 30) {
  covar_str <- if (!is.null(covars) && length(covars) > 0) paste("+", paste(covars, collapse = " + ")) else ""
  fmla <- stats::as.formula(paste("survival::Surv(t0, t1, event_cc_day) ~ arm", covar_str))
  
  model <- survival::coxph(
    fmla,
    data = panel_w,
    weights = panel_w$w_trim,
    cluster = panel_w$id
  )
  
  # Extract Hazard Ratio
  coef_name <- grep("arm", names(stats::coef(model)), value = TRUE)[1]
  beta_arm <- stats::coef(model)[coef_name]
  se_beta <- sqrt(stats::vcov(model)[coef_name, coef_name])
  hr <- exp(beta_arm)
  hr_ci <- exp(c(beta_arm - 1.96 * se_beta, beta_arm + 1.96 * se_beta))
  
  # Predict counterfactual survival for both arms
  dummy <- data.frame(arm = c("control_no_tx", "treat_tx"), id = c("1", "2"))
  if (!is.null(covars) && length(covars) > 0) {
    for (cv in covars) {
      dummy[[cv]] <- mean(panel_w[[cv]], na.rm = TRUE)
    }
  }
  
  sf <- survival::survfit(model, newdata = dummy)
  s_tau <- summary(sf, times = tau)
  
  s0_tau <- s_tau$surv[1]
  s1_tau <- s_tau$surv[2]
  se0_tau <- s_tau$std.err[1]
  se1_tau <- s_tau$std.err[2]
  
  rd_tau <- s1_tau - s0_tau
  se_rd_tau <- sqrt(se0_tau^2 + se1_tau^2)
  
  return(list(
    model = model,
    hr = hr,
    hr_ci = hr_ci,
    s0_tau = s0_tau,
    s1_tau = s1_tau,
    rd_tau = rd_tau,
    se_rd_tau = se_rd_tau
  ))
}

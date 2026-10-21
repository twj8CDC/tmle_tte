#' Competing Benchmark Estimators for Time-to-Event Analyses
#'
#' @description Implementations of standard benchmark survival estimators used for
#'   methodological comparisons against TTE-TMLE, including Naive Cox, Landmarking,
#'   and Time-Dependent Cox regression.
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

#' Time-Dependent Cox Proportional Hazards Model
#'
#' @description Fits an extended Cox model using counting process format (tstart, tstop, event)
#'   with a time-varying treatment covariate.
#'
#' @param data Observational data frame.
#' @param covars Character vector of baseline covariates.
#' @param tau Maximum follow-up time (default 30).
#' @return A fitted \code{coxph} object.
#' @export
#' @importFrom dplyr mutate filter bind_rows arrange
#' @importFrom survival coxph Surv
#' @importFrom stats as.formula
fit_td_cox <- function(data, covars = NULL, tau = 30) {
  df <- prepare_tte_data(data)
  # Split observations into pre-switch and post-switch intervals
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
  survival::coxph(fmla, data = td_data)
}

#' Data-Generating Process for Scenario 1: Structural Baseline (Pure ITB + Grace Period)
#'
#' @description Generates observational time-to-event data with baseline confounding,
#'   treatment initiation delay (grace period), and potential immortal time bias.
#'
#' @param n Sample size (e.g., 300 or 1000).
#' @param tau Maximum follow-up time (default 30).
#' @param grace_period Grace period duration in days (default 10).
#' @param shape Weibull shape parameter (default 1.5).
#' @param base_scale Baseline scale parameter (default 40).
#' @param coef_w Confounder effect on baseline hazard (default 0.5).
#' @param hr_target True hazard ratio benefit of treatment (default 0.35).
#' @param p_tx_intercept Intercept for treatment propensity logistic model (default 0.0).
#' @param p_tx_coef Coefficient for confounder in treatment propensity model (default 0.6).
#' @param seed Optional random seed.
#'
#' @return A tibble with standardized observational columns:
#'   \describe{
#'     \item{id}{Patient identifier.}
#'     \item{x1}{Baseline confounder (W).}
#'     \item{tx_obs}{Observed treatment initiation indicator (1 if initiated, 0 otherwise).}
#'     \item{tx_time_obs}{Observed initiation time (or Inf if never initiated).}
#'     \item{event_time}{Observed follow-up time (min(event_time, tau)).}
#'     \item{event}{Binary event indicator (1 if event observed by tau, 0 otherwise).}
#'   }

library(dplyr)
library(tibble)

generate_scenario1_data <- function(n = 500,
                                    tau = 30,
                                    grace_period = 10,
                                    shape = 1.5,
                                    base_scale = 40,
                                    coef_w = 0.5,
                                    hr_target = 0.35,
                                    p_tx_intercept = 0.0,
                                    p_tx_coef = 0.6,
                                    seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # 1. Baseline Confounder
  w <- rnorm(n, mean = 0, sd = 1)
  
  # 2. Treatment Initiation Propensity
  prob_tx <- plogis(p_tx_intercept + p_tx_coef * w)
  tx_obs <- rbinom(n, 1, prob_tx)
  
  # 3. Treatment Initiation Time (Complier distribution among initiators)
  tx_time_obs <- ifelse(
    tx_obs == 1L,
    pmax(1, pmin(grace_period, round(rbeta(n, 4, 2) * grace_period))),
    Inf
  )
  
  # 4. Weibull Scales
  scale_0 <- base_scale * exp(- (coef_w * w) / shape)
  scale_t <- scale_0 * (hr_target)^(-1 / shape)
  
  # 5. Event Times (Piecewise Weibull Hazard Inversion)
  u <- runif(n)
  h_target <- -log(u)
  h_at_switch <- ifelse(is.finite(tx_time_obs), (tx_time_obs / scale_0)^shape, Inf)
  
  raw_event_time <- ifelse(
    h_target <= h_at_switch,
    scale_0 * (h_target)^(1 / shape),
    ((h_target - h_at_switch) * (scale_t^shape) + (tx_time_obs^shape))^(1 / shape)
  )
  
  # 6. Discrete follow-up and administrative censoring at tau
  event_time_discrete <- pmax(1, pmin(tau, ceiling(raw_event_time)))
  event_indicator <- ifelse(raw_event_time <= tau, 1L, 0L)
  
  df <- tibble(
    id = seq_len(n),
    x1 = w,
    tx_obs = tx_obs,
    tx_time_obs = tx_time_obs,
    event_time = event_time_discrete,
    event = event_indicator
  )
  
  return(df)
}


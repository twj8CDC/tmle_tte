#' Data-Generating Process for Scenario 5: Positivity Strain & Weight Shock
#'
#' @description Generates observational time-to-event data featuring severe overlap strain
#'   (near-positivity violations) where a frail patient subgroup has an extremely low propensity
#'   of initiating treatment during the grace period (p < 0.03).
#'   Evaluates how pure weighting (CCW-KM) suffers from extreme weight inflation and instability,
#'   whereas TTE-TMLE and Pure G-computation maintain stability and bounded inference.
#'
#' @param n Sample size (default 1000).
#' @param tau Maximum follow-up time (default 30).
#' @param grace_period Grace period duration in days (default 10).
#' @param shape Weibull shape parameter (default 1.5).
#' @param base_scale Baseline scale parameter (default 45).
#' @param coef_w1 Confounder 1 effect on baseline hazard (default 0.5).
#' @param coef_w2 Confounder 2 effect on baseline hazard (default 0.4).
#' @param hr_target True hazard ratio benefit of treatment (default 0.40).
#' @param p_tx_intercept Intercept for treatment propensity (default -1.6, inducing positivity strain).
#' @param p_tx_coef1 Coefficient for x1 in propensity model (default 1.2).
#' @param p_tx_coef2 Coefficient for x2 in propensity model (default 1.0).
#' @param seed Optional random seed.
#'
#' @return A tibble with standardized observational columns:
#'   \describe{
#'     \item{id}{Patient identifier.}
#'     \item{x1}{Baseline clinical confounder (e.g. ECOG / functional score).}
#'     \item{x2}{Baseline frailty index (e.g. comorbidity / organ score).}
#'     \item{tx_obs}{Observed treatment initiation indicator (1 if initiated, 0 otherwise).}
#'     \item{tx_time_obs}{Observed initiation time (or Inf if never initiated).}
#'     \item{event_time}{Observed follow-up time (min(event_time, tau)).}
#'     \item{event}{Binary event indicator (1 if event observed by tau, 0 otherwise).}
#'   }
#' @export
#' @importFrom dplyr tibble mutate
#' @importFrom stats rnorm rbinom runif plogis rbeta

generate_scenario5_data <- function(n = 1000,
                                    tau = 30,
                                    grace_period = 10,
                                    shape = 1.5,
                                    base_scale = 45,
                                    coef_w1 = 0.5,
                                    coef_w2 = 0.4,
                                    hr_target = 0.40,
                                    p_tx_intercept = -1.6,
                                    p_tx_coef1 = 1.2,
                                    p_tx_coef2 = 1.0,
                                    seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # 1. Baseline Covariates
  x1 <- stats::rnorm(n, mean = 0, sd = 1)
  x2 <- stats::rnorm(n, mean = 0, sd = 1)
  
  # 2. Treatment Initiation Propensity with severe positivity strain in lower tails
  # Logit = -1.6 + 1.2*x1 + 1.0*x2
  # For x1 = -1.5, x2 = -1.5 -> logit = -4.9 -> prob ~ 0.007
  prob_tx <- stats::plogis(p_tx_intercept + p_tx_coef1 * x1 + p_tx_coef2 * x2)
  tx_obs <- stats::rbinom(n, 1, prob_tx)
  
  # 3. Treatment Initiation Time (Complier distribution among initiators)
  tx_time_obs <- ifelse(
    tx_obs == 1L,
    pmax(1, pmin(grace_period, round(stats::rbeta(n, 4, 2) * grace_period))),
    Inf
  )
  
  # 4. Weibull Scales
  lin_pred_0 <- coef_w1 * x1 + coef_w2 * x2
  scale_0 <- base_scale * exp(- lin_pred_0 / shape)
  scale_t <- scale_0 * (hr_target)^(-1 / shape)
  
  # 5. Event Times (Piecewise Weibull Hazard Inversion)
  u <- stats::runif(n)
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
  
  df <- dplyr::tibble(
    id = seq_len(n),
    x1 = x1,
    x2 = x2,
    tx_obs = tx_obs,
    tx_time_obs = tx_time_obs,
    event_time = event_time_discrete,
    event = event_indicator
  )
  
  return(df)
}


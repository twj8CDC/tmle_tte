#' Data-Generating Process for Scenario 8: Depletion of Susceptibles & Survivor Collider Bias
#'
#' @description Generates observational time-to-event data with unmeasured frailty (U)
#'   influencing survival time. Even when the true individual hazard ratio is constant (HR = 0.35),
#'   the differential survival of frail individuals induces an artificial time-waning in the
#'   observable marginal hazard ratio (collider conditioning on T >= t).
#'   Demonstrates how hazard-ratio estimators (td-Cox, ccw-cphm) report biased, uninterpretable summaries,
#'   whereas TMLE targeting prospective causal Risk Differences RD(t) remains completely valid.
#'
#' @param n Sample size (default 1000).
#' @param tau Maximum follow-up time (default 30).
#' @param grace_period Grace period duration in days (default 10).
#' @param shape Weibull shape parameter (default 1.5).
#' @param base_scale Baseline scale parameter (default 45).
#' @param coef_w Confounder effect on baseline hazard (default 0.4).
#' @param frailty_var Variance of unmeasured frailty U (default 0.6).
#' @param hr_target True biological individual hazard ratio benefit (default 0.35).
#' @param p_tx_intercept Intercept for treatment propensity (default 0.0).
#' @param p_tx_coef Coefficient for confounder in treatment propensity (default 0.5).
#' @param seed Optional random seed.
#'
#' @return A tibble with standardized observational columns:
#'   \describe{
#'     \item{id}{Patient identifier.}
#'     \item{x1}{Measured baseline confounder.}
#'     \item{tx_obs}{Observed treatment initiation indicator.}
#'     \item{tx_time_obs}{Observed initiation time (or Inf if never initiated).}
#'     \item{event_time}{Observed follow-up time.}
#'     \item{event}{Binary event indicator.}
#'   }
#' @export
#' @importFrom dplyr tibble mutate
#' @importFrom stats rnorm rbinom runif plogis rbeta rgamma

generate_scenario8_data <- function(n = 1000,
                                    tau = 30,
                                    grace_period = 10,
                                    shape = 1.5,
                                    base_scale = 45,
                                    coef_w = 0.4,
                                    frailty_var = 0.6,
                                    hr_target = 0.35,
                                    p_tx_intercept = 0.0,
                                    p_tx_coef = 0.5,
                                    seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # 1. Measured Confounder and Unmeasured Frailty
  x1 <- stats::rnorm(n, mean = 0, sd = 1)
  # Unmeasured frailty U ~ Gamma(shape = 1/frailty_var, scale = frailty_var), mean = 1
  k_shape <- 1 / frailty_var
  u_frailty <- stats::rgamma(n, shape = k_shape, scale = frailty_var)
  
  # 2. Treatment Propensity (Depends ONLY on measured x1, NOT unmeasured U)
  prob_tx <- stats::plogis(p_tx_intercept + p_tx_coef * x1)
  tx_obs <- stats::rbinom(n, 1, prob_tx)
  
  # 3. Treatment Initiation Time (Complier distribution among initiators)
  tx_time_obs <- ifelse(
    tx_obs == 1L,
    pmax(1, pmin(grace_period, round(stats::rbeta(n, 4, 2) * grace_period))),
    Inf
  )
  
  # 4. Weibull Scales with multiplicative frailty
  # scale_0 is scaled by (u_frailty)^(-1/shape)
  scale_0 <- (base_scale * exp(- (coef_w * x1) / shape)) * (u_frailty)^(-1 / shape)
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
    tx_obs = tx_obs,
    tx_time_obs = tx_time_obs,
    event_time = event_time_discrete,
    event = event_indicator
  )
  
  return(df)
}


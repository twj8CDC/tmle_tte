#' Data-Generating Process for Scenario 6: Confounder-Treatment Feedback (G-Null / Time-Varying Confounding)
#'
#' @description Generates observational time-to-event data featuring a time-varying intermediate
#'   biomarker L(t) (e.g., acute organ injury at Day 3) affected by early treatment (A_early -> L),
#'   which subsequently confounds late treatment decisions (L -> A_late) and directly affects survival (L -> Y).
#'   Demonstrates how conditioning on L in standard Cox models induces mediator-collider stratification bias,
#'   while ignoring L induces unmeasured confounding. Only Target Trial Emulation (CCW / TTE-TMLE) succeeds.
#'
#' @param n Sample size (default 1000).
#' @param tau Maximum follow-up time (default 30).
#' @param grace_period Grace period duration in days (default 10).
#' @param t_med Time of intermediate biomarker assessment (default 3).
#' @param shape Weibull shape parameter (default 1.5).
#' @param base_scale Baseline scale parameter (default 45).
#' @param coef_w Confounder effect on baseline hazard (default 0.4).
#' @param coef_l Intermediate biomarker effect on baseline hazard (default 0.8, HR ~ 2.2).
#' @param hr_target True hazard ratio benefit of treatment (default 0.35).
#' @param seed Optional random seed.
#'
#' @return A tibble with standardized observational columns:
#'   \describe{
#'     \item{id}{Patient identifier.}
#'     \item{x1}{Baseline clinical confounder.}
#'     \item{l_med}{Intermediate biomarker status assessed at Day 3 (1 if organ injury, 0 otherwise).}
#'     \item{tx_obs}{Observed treatment initiation indicator (1 if initiated, 0 otherwise).}
#'     \item{tx_time_obs}{Observed initiation time (or Inf if never initiated).}
#'     \item{event_time}{Observed follow-up time (min(event_time, tau)).}
#'     \item{event}{Binary event indicator (1 if event observed by tau, 0 otherwise).}
#'   }
#' @export
#' @importFrom dplyr tibble mutate
#' @importFrom stats rnorm rbinom runif plogis rbeta

generate_scenario6_data <- function(n = 1000,
                                    tau = 30,
                                    grace_period = 10,
                                    t_med = 3,
                                    shape = 1.5,
                                    base_scale = 45,
                                    coef_w = 0.4,
                                    coef_l = 0.8,
                                    hr_target = 0.35,
                                    seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # 1. Baseline Confounder
  x1 <- stats::rnorm(n, mean = 0, sd = 1)
  
  # 2. Early Treatment Initiation (t <= t_med, e.g. Days 1-2)
  prob_tx_early <- stats::plogis(-0.8 + 0.6 * x1)
  tx_early <- stats::rbinom(n, 1, prob_tx_early)
  t_switch_early <- ifelse(tx_early == 1L, sample(1:(t_med - 1), size = n, replace = TRUE), Inf)
  
  # 3. Intermediate Time-Varying Biomarker L at Day 3
  # Treatment early reduces probability of organ injury (A -> L)
  # Baseline sickness increases probability (x1 -> L)
  prob_l <- stats::plogis(-0.2 + 0.7 * x1 - 1.2 * tx_early)
  l_med <- stats::rbinom(n, 1, prob_l)
  
  # 4. Late Treatment Initiation (t_med < t <= grace_period, e.g. Days 4-10)
  # Clinicians escalate therapy if patient develops organ injury (L -> A)
  prob_tx_late <- stats::plogis(-1.2 + 0.5 * x1 + 1.4 * l_med)
  tx_late <- ifelse(tx_early == 1L, 0L, stats::rbinom(n, 1, prob_tx_late))
  t_switch_late <- ifelse(tx_late == 1L, sample((t_med + 1):grace_period, size = n, replace = TRUE), Inf)
  
  # Combined Treatment Initiation
  tx_obs <- ifelse(tx_early == 1L | tx_late == 1L, 1L, 0L)
  tx_time_obs <- pmin(t_switch_early, t_switch_late)
  
  # 5. Event Times (Piecewise Weibull Hazard Inversion)
  # Baseline hazard scale depends on x1 and L(t)
  lin_pred_0 <- coef_w * x1 + coef_l * l_med
  scale_0 <- base_scale * exp(- lin_pred_0 / shape)
  scale_t <- scale_0 * (hr_target)^(-1 / shape)
  
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
    l_med = l_med,
    tx_obs = tx_obs,
    tx_time_obs = tx_time_obs,
    event_time = event_time_discrete,
    event = event_indicator
  )
  
  return(df)
}


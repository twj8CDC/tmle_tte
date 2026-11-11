#' Data-Generating Process for Scenario 6: Confounder-Treatment Feedback (Time-Varying Confounding)
#'
#' @description Builds directly off Scenario 2's complex non-linear (w1, w2) confounding DGP,
#'   adding a time-varying intermediate biomarker L(t) at Day 3 affected by early treatment
#'   (A_early -> L), which subsequently confounds late treatment decisions (L -> A_late)
#'   and directly affects survival (L -> Y).
#'
#' @param n Sample size (default 1000).
#' @param tau Maximum follow-up time (default 30).
#' @param grace_period Grace period duration in days (default 10).
#' @param t_med Time of intermediate biomarker assessment (default 3).
#' @param shape Weibull shape parameter (default 1.5).
#' @param base_scale Baseline scale parameter (default 45).
#' @param coef_l Intermediate biomarker effect on baseline hazard (default 0.8).
#' @param hr_target True hazard ratio benefit of treatment (default 0.35).
#' @param seed Optional random seed.
#' @return A tibble with standardized observational columns.

library(dplyr)
library(tibble)

generate_scenario6_data <- function(n = 1000,
                                    tau = 30,
                                    grace_period = 10,
                                    t_med = 3,
                                    shape = 1.5,
                                    base_scale = 45,
                                    coef_l = 0.8,
                                    hr_target = 0.35,
                                    seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # 1. Baseline Confounders (Scenario 2 Complex Non-Linear Base)
  w1 <- stats::rnorm(n, mean = 3, sd = 1.5)
  w2 <- stats::rbinom(n, 1, prob = 0.4)
  w1_std <- (w1 - 3) / 1.5
  gate <- ifelse(w1 > 3.8 & w2 == 1L, 1.0, 0.0)
  nl_score <- -0.35 + 0.35 * exp(0.70 * w1_std) + 0.25 * (w1_std^2) + 0.50 * w2 + 0.45 * gate

  # 2. Early Treatment Initiation (Days 1..t_med)
  prob_tx_early <- pmin(pmax(stats::plogis(-1.1 + 0.85 * nl_score), 0.03), 0.90)
  tx_early <- stats::rbinom(n, 1, prob_tx_early)
  t_switch_early <- ifelse(tx_early == 1L, sample(seq_len(t_med), size = n, replace = TRUE), Inf)

  # 3. Intermediate Time-Varying Biomarker L at Day 3 (A_early -> L, W -> L)
  prob_l <- stats::plogis(-0.3 + 0.80 * nl_score - 1.25 * tx_early)
  l_med <- stats::rbinom(n, 1, prob_l)

  # 4. Late Treatment Initiation (Days t_med+1..grace_period) (L -> A_late, W -> A_late)
  prob_tx_late <- pmin(pmax(stats::plogis(-1.3 + 0.70 * nl_score + 1.40 * l_med), 0.03), 0.95)
  tx_late <- ifelse(tx_early == 1L, 0L, stats::rbinom(n, 1, prob_tx_late))
  t_switch_late <- ifelse(tx_late == 1L, sample((t_med + 1L):grace_period, size = n, replace = TRUE), Inf)

  tx_obs <- ifelse(tx_early == 1L | tx_late == 1L, 1L, 0L)
  tx_time_obs <- pmin(t_switch_early, t_switch_late)

  # 5. Event Times (Piecewise Weibull Hazard Inversion; continuous switch time = d - 1)
  lin_pred_0 <- nl_score + coef_l * l_med
  scale_0 <- base_scale * exp(-lin_pred_0 / shape)
  indiv_hr <- hr_target * exp(0.15 * w1_std)
  scale_t <- scale_0 * (indiv_hr)^(-1 / shape)

  u <- stats::runif(n)
  h_target <- -log(u)
  t_sw_cont <- ifelse(is.finite(tx_time_obs), pmax(0, tx_time_obs - 1), Inf)
  h_at_switch <- ifelse(is.finite(t_sw_cont), (t_sw_cont / scale_0)^shape, Inf)

  raw_event_time <- ifelse(
    h_target <= h_at_switch,
    scale_0 * (h_target)^(1 / shape),
    ((h_target - h_at_switch) * (scale_t^shape) + (t_sw_cont^shape))^(1 / shape)
  )

  event_time_discrete <- pmax(1, pmin(tau, ceiling(raw_event_time)))
  event_indicator <- ifelse(raw_event_time <= tau, 1L, 0L)

  dplyr::tibble(
    id = seq_len(n),
    w1 = w1,
    w2 = w2,
    l_med = l_med,
    tx_obs = tx_obs,
    tx_time_obs = tx_time_obs,
    time = event_time_discrete,
    event_time = event_time_discrete,
    tx_ever = tx_obs,
    event = event_indicator
  )
}

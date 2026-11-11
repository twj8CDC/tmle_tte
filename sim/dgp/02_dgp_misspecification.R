#' Data-Generating Process for Scenario 2: Complex Non-Linear DGP (Exponential + Quadratic + Threshold)
#'
#' @description Generates observational survival data with strong exponential and quadratic
#'   non-linear confounding in w1 and a threshold biomarker interaction with w2 in both
#'   treatment propensity and outcome hazard. Serves as the foundational non-linear DGP
#'   for Scenarios 2 through 8.
#'
#' @param n Sample size (default 1000).
#' @param tau Maximum follow-up time (default 30).
#' @param grace_period Grace period duration in days (default 10).
#' @param shape Weibull shape parameter (default 1.5).
#' @param base_scale Baseline scale parameter (default 45).
#' @param hr_target Base treatment effect hazard ratio (default 0.35).
#' @param seed Optional random seed.
#' @return A tibble with columns id, w1, w2, tx_obs, tx_time_obs, time, event_time, tx_ever, event.

library(dplyr)
library(tibble)

generate_scenario2_data <- function(n = 1000,
                                    tau = 30,
                                    grace_period = 10,
                                    shape = 1.5,
                                    base_scale = 45,
                                    hr_target = 0.35,
                                    seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # 1. Baseline Confounders: Continuous Frailty (w1) and Biomarker Gate (w2)
  w1 <- rnorm(n, mean = 3, sd = 1.5)
  w2 <- rbinom(n, 1, prob = 0.4)
  w1_std <- (w1 - 3) / 1.5
  gate <- ifelse(w1 > 3.8 & w2 == 1L, 1.0, 0.0)

  # 2. Non-Linear Treatment Propensity (Exponential + Quadratic + Threshold Gate)
  logit_p_tx <- -0.75 + 0.45 * exp(0.70 * w1_std) + 0.30 * (w1_std^2) + 0.40 * w2 + 1.25 * gate
  prob_tx <- pmin(pmax(plogis(logit_p_tx), 0.04), 0.96)
  tx_obs <- rbinom(n, 1, prob_tx)

  # 3. Treatment Initiation Day d in {1, ..., grace_period}
  tx_time_obs <- ifelse(
    tx_obs == 1L,
    pmax(1, pmin(grace_period, round(rbeta(n, 4, 2) * grace_period))),
    Inf
  )

  # 4. Non-Linear Outcome Scale (Exponential + Quadratic + Threshold Gate)
  lin_pred_0 <- -0.35 + 0.35 * exp(0.70 * w1_std) + 0.25 * (w1_std^2) + 0.50 * w2 + 0.45 * gate
  scale_0 <- base_scale * exp(-lin_pred_0 / shape)

  # Treatment effect with frailty interaction
  indiv_hr <- hr_target * exp(0.15 * w1_std)
  scale_t <- scale_0 * (indiv_hr)^(-1 / shape)

  # 5. Event Times via Piecewise Inversion (continuous switch time = d - 1)
  u <- runif(n)
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

  tibble(
    id = seq_len(n),
    w1 = w1,
    w2 = w2,
    tx_obs = tx_obs,
    tx_time_obs = tx_time_obs,
    time = event_time_discrete,
    event_time = event_time_discrete,
    tx_ever = tx_obs,
    event = event_indicator
  )
}

#' Data-Generating Process for Scenario 3: Non-Proportional Hazards & Waning Efficacy
#'
#' @description Builds directly off Scenario 2's complex non-linear (w1, w2) confounding DGP,
#'   adding structural non-proportional hazards: strong early protection (HR = 0.20) for t <= 15 days,
#'   followed by efficacy waning (HR = 1.10) for t > 15 days.
#'
#' @param n Sample size (default 1000).
#' @param tau Follow-up time horizon (default 30).
#' @param grace_period Grace period duration in days (default 10).
#' @param t_wane Time at which treatment efficacy wanes (default 15).
#' @param shape Weibull shape parameter (default 1.5).
#' @param base_scale Baseline scale parameter (default 45).
#' @param hr_early Early treatment hazard ratio (default 0.20).
#' @param hr_late Late/waned treatment hazard ratio (default 1.10).
#' @param seed Optional random seed.
#' @return A tibble with standardized observational columns.

library(dplyr)
library(tibble)

generate_scenario3_data <- function(n = 1000,
                                    tau = 30,
                                    grace_period = 10,
                                    t_wane = 15,
                                    shape = 1.5,
                                    base_scale = 45,
                                    hr_early = 0.20,
                                    hr_late = 1.10,
                                    seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # 1. Baseline Confounders (Scenario 2 Complex Non-Linear Base)
  w1 <- stats::rnorm(n, mean = 3, sd = 1.5)
  w2 <- stats::rbinom(n, 1, prob = 0.4)
  w1_std <- (w1 - 3) / 1.5
  gate <- ifelse(w1 > 3.8 & w2 == 1L, 1.0, 0.0)

  # 2. Non-Linear Treatment Propensity
  logit_p_tx <- -0.75 + 0.45 * exp(0.70 * w1_std) + 0.30 * (w1_std^2) + 0.40 * w2 + 1.25 * gate
  prob_tx <- pmin(pmax(stats::plogis(logit_p_tx), 0.04), 0.96)
  tx_obs <- stats::rbinom(n, 1, prob = prob_tx)

  # 3. Treatment Initiation Day d in {1, ..., grace_period}
  tx_time_obs <- ifelse(
    tx_obs == 1L,
    pmax(1, pmin(grace_period, round(stats::rbeta(n, 4, 2) * grace_period))),
    Inf
  )

  # 4. Non-Linear Baseline Scale
  lin_pred_0 <- -0.35 + 0.35 * exp(0.70 * w1_std) + 0.25 * (w1_std^2) + 0.50 * w2 + 0.45 * gate
  scale_0 <- base_scale * exp(-lin_pred_0 / shape)

  # 5. Piecewise Hazard Inversion (continuous switch time = d - 1)
  u <- stats::runif(n)
  h_target <- -log(u)
  t_sw_cont <- ifelse(is.finite(tx_time_obs), pmax(0, tx_time_obs - 1), Inf)

  h_at_switch <- ifelse(is.finite(t_sw_cont), (t_sw_cont / scale_0)^shape, Inf)
  h_at_wane_base <- (t_wane / scale_0)^shape
  h_at_wane_actual <- ifelse(
    t_sw_cont <= t_wane,
    h_at_switch + hr_early * (h_at_wane_base - h_at_switch),
    h_at_wane_base
  )

  raw_event_time <- numeric(n)
  for (i in seq_len(n)) {
    ht <- h_target[i]
    tsw <- t_sw_cont[i]
    sc0 <- scale_0[i]
    h_sw <- h_at_switch[i]
    h_wn <- h_at_wane_actual[i]

    if (ht <= h_sw) {
      raw_event_time[i] <- sc0 * (ht)^(1 / shape)
    } else if (ht <= h_wn) {
      term <- (tsw / sc0)^shape + (ht - h_sw) / hr_early
      raw_event_time[i] <- sc0 * (term)^(1 / shape)
    } else {
      term <- (t_wane / sc0)^shape + (ht - h_wn) / hr_late
      raw_event_time[i] <- sc0 * (term)^(1 / shape)
    }
  }

  event_time_discrete <- pmax(1, pmin(tau, ceiling(raw_event_time)))
  event_indicator <- ifelse(raw_event_time <= tau, 1L, 0L)

  tibble::tibble(
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

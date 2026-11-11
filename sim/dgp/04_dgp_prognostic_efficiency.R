#' Data-Generating Process for Scenario 4: Prognostic Efficiency & Local Efficiency
#'
#' @description Builds directly off Scenario 2's complex non-linear (w1, w2) confounding DGP,
#'   adding a conditionally independent, strongly prognostic non-linear biomarker (w_prog).
#'   Evaluates how outcome-model marginalization absorbs prognostic variance to shrink SEs
#'   and CI widths relative to pure IPW (CCW-KM), while misspecifying the functional form
#'   (without omitting any variable) degrades parametric precision and accuracy.
#'
#' @param n Sample size (default 1000).
#' @param tau Maximum follow-up time (default 30).
#' @param grace_period Grace period duration in days (default 10).
#' @param shape Weibull shape parameter (default 1.5).
#' @param base_scale Baseline scale parameter (default 45).
#' @param hr_target True hazard ratio benefit of treatment (default 0.35).
#' @param seed Optional random seed.
#' @return A tibble with standardized observational columns.

library(dplyr)
library(tibble)

generate_scenario4_data <- function(n = 1000,
                                    tau = 30,
                                    grace_period = 10,
                                    shape = 1.5,
                                    base_scale = 45,
                                    hr_target = 0.35,
                                    seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # 1. Baseline Confounders (w1, w2) + Independent Prognostic Biomarker (w_prog)
  w1 <- stats::rnorm(n, mean = 3, sd = 1.5)
  w2 <- stats::rbinom(n, 1, prob = 0.4)
  w_prog <- stats::rnorm(n, mean = 0, sd = 1)

  w1_std <- (w1 - 3) / 1.5
  gate <- ifelse(w1 > 3.8 & w2 == 1L, 1.0, 0.0)

  # 2. Treatment Initiation Propensity (Depends on w1, w2, independent of w_prog)
  logit_p_tx <- -0.75 + 0.45 * exp(0.70 * w1_std) + 0.30 * (w1_std^2) + 0.40 * w2 + 1.25 * gate
  prob_tx <- pmin(pmax(stats::plogis(logit_p_tx), 0.04), 0.96)
  tx_obs <- stats::rbinom(n, 1, prob_tx)

  # 3. Treatment Initiation Day d in {1, ..., grace_period}
  tx_time_obs <- ifelse(
    tx_obs == 1L,
    pmax(1, pmin(grace_period, round(stats::rbeta(n, 4, 2) * grace_period))),
    Inf
  )

  # 4. Weibull Scales (Driven by non-linear w1, w2 AND non-linear prognostic w_prog)
  lin_pred_0 <- -0.55 + 0.35 * exp(0.70 * w1_std) + 0.25 * (w1_std^2) + 0.50 * w2 + 0.45 * gate +
    0.55 * exp(0.55 * w_prog) + 0.25 * (w_prog^2)
  scale_0 <- base_scale * exp(-lin_pred_0 / shape)
  indiv_hr <- hr_target * exp(0.15 * w1_std)
  scale_t <- scale_0 * (indiv_hr)^(-1 / shape)

  # 5. Event Times (Piecewise Weibull Hazard Inversion; continuous switch time = d - 1)
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
    w_prog = w_prog,
    tx_obs = tx_obs,
    tx_time_obs = tx_time_obs,
    time = event_time_discrete,
    event_time = event_time_discrete,
    tx_ever = tx_obs,
    event = event_indicator
  )
}

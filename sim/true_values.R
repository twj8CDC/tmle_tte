#' Ground Truth Counterfactual Value Generator (Scenarios 1 to 8)
#'
#' @description Computes exact counterfactual survival curves and causal estimands
#'   by generating a large Monte Carlo oracle population (N = 200,000) under intervention
#'   Strategy 1 (initiate treatment within the grace period following the empirical complier
#'   initiation distribution) and Strategy 0 (never treated), using Common Random Numbers (u1 = u0)
#'   and exact continuous-to-discrete switch boundary alignment (t_sw_cont = d - 1).

library(dplyr)
library(tibble)

#' Ground Truth Counterfactual Value Generator for Scenario 1 (Simple Linear DGP)
#' @export
compute_scenario1_truth <- function(n_oracle = 200000,
                                    tau = 30,
                                    grace_period = 10,
                                    shape = 1.5,
                                    base_scale = 40,
                                    coef_w = 0.5,
                                    hr_target = 0.35,
                                    seed = 99999) {
  if (!is.null(seed)) set.seed(seed)

  w <- stats::rnorm(n_oracle, mean = 0, sd = 1)
  scale_0 <- base_scale * exp(-(coef_w * w) / shape)
  scale_t <- scale_0 * (hr_target)^(-1 / shape)

  u0 <- stats::runif(n_oracle)
  h_target0 <- -log(u0)
  t_event_0 <- scale_0 * (h_target0)^(1 / shape)

  t_switch_1 <- pmax(1, pmin(grace_period, round(stats::rbeta(n_oracle, 4, 2) * grace_period)))
  t_sw_cont <- pmax(0, t_switch_1 - 1)
  h_target1 <- h_target0
  h_at_switch <- (t_sw_cont / scale_0)^shape

  t_event_1 <- ifelse(
    h_target1 <= h_at_switch,
    scale_0 * (h_target1)^(1 / shape),
    ((h_target1 - h_at_switch) * (scale_t^shape) + (t_sw_cont^shape))^(1 / shape)
  )

  days <- seq_len(tau)
  s0_vec <- vapply(days, function(d) mean(t_event_0 > d), numeric(1))
  s1_vec <- vapply(days, function(d) mean(t_event_1 > d), numeric(1))

  tibble(
    t = days,
    S0_true = s0_vec,
    S1_true = s1_vec,
    rd_true = s1_vec - s0_vec,
    hr_true = log(pmax(s1_vec, 1e-6)) / log(pmax(s0_vec, 1e-6))
  )
}

#' Ground Truth Counterfactual Value Generator for Scenario 2 (Complex Non-Linear DGP)
#' @export
compute_scenario2_truth <- function(n_oracle = 200000,
                                    tau = 30,
                                    grace_period = 10,
                                    shape = 1.5,
                                    base_scale = 45,
                                    hr_target = 0.35,
                                    seed = 99999) {
  if (!is.null(seed)) set.seed(seed)

  w1 <- stats::rnorm(n_oracle, mean = 3, sd = 1.5)
  w2 <- stats::rbinom(n_oracle, 1, prob = 0.4)
  w1_std <- (w1 - 3) / 1.5
  gate <- ifelse(w1 > 3.8 & w2 == 1L, 1.0, 0.0)

  lin_pred_0 <- -0.35 + 0.35 * exp(0.70 * w1_std) + 0.25 * (w1_std^2) + 0.50 * w2 + 0.45 * gate
  scale_0 <- base_scale * exp(-lin_pred_0 / shape)
  indiv_hr <- hr_target * exp(0.15 * w1_std)
  scale_t <- scale_0 * (indiv_hr)^(-1 / shape)

  u0 <- stats::runif(n_oracle)
  h_target0 <- -log(u0)
  t_event_0 <- scale_0 * (h_target0)^(1 / shape)

  t_switch_1 <- pmax(1, pmin(grace_period, round(stats::rbeta(n_oracle, 4, 2) * grace_period)))
  t_sw_cont <- pmax(0, t_switch_1 - 1)
  h_target1 <- h_target0
  h_at_switch <- (t_sw_cont / scale_0)^shape

  t_event_1 <- ifelse(
    h_target1 <= h_at_switch,
    scale_0 * (h_target1)^(1 / shape),
    ((h_target1 - h_at_switch) * (scale_t^shape) + (t_sw_cont^shape))^(1 / shape)
  )

  days <- seq_len(tau)
  s0_vec <- vapply(days, function(d) mean(t_event_0 > d), numeric(1))
  s1_vec <- vapply(days, function(d) mean(t_event_1 > d), numeric(1))

  tibble(
    t = days,
    S0_true = s0_vec,
    S1_true = s1_vec,
    rd_true = s1_vec - s0_vec,
    hr_true = log(pmax(s1_vec, 1e-6)) / log(pmax(s0_vec, 1e-6))
  )
}

#' Ground Truth Counterfactual Value Generator for Scenario 3 (Non-Proportional Hazards & Waning)
#' @export
compute_scenario3_truth <- function(n_oracle = 200000,
                                    tau = 30,
                                    grace_period = 10,
                                    t_wane = 15,
                                    shape = 1.5,
                                    base_scale = 45,
                                    hr_early = 0.20,
                                    hr_late = 1.10,
                                    seed = 99999) {
  if (!is.null(seed)) set.seed(seed)

  w1 <- stats::rnorm(n_oracle, mean = 3, sd = 1.5)
  w2 <- stats::rbinom(n_oracle, 1, prob = 0.4)
  w1_std <- (w1 - 3) / 1.5
  gate <- ifelse(w1 > 3.8 & w2 == 1L, 1.0, 0.0)

  lin_pred_0 <- -0.35 + 0.35 * exp(0.70 * w1_std) + 0.25 * (w1_std^2) + 0.50 * w2 + 0.45 * gate
  scale_0 <- base_scale * exp(-lin_pred_0 / shape)

  u0 <- stats::runif(n_oracle)
  h_target0 <- -log(u0)
  t_event_0 <- scale_0 * (h_target0)^(1 / shape)

  t_switch_1 <- pmax(1, pmin(grace_period, round(stats::rbeta(n_oracle, 4, 2) * grace_period)))
  t_sw_cont <- pmax(0, t_switch_1 - 1)
  h_target1 <- h_target0

  h_at_switch <- (t_sw_cont / scale_0)^shape
  h_at_wane_base <- (t_wane / scale_0)^shape
  h_at_wane_actual <- h_at_switch + hr_early * (h_at_wane_base - h_at_switch)

  t_event_1 <- numeric(n_oracle)
  idx_pre <- which(h_target1 <= h_at_switch)
  idx_early <- which(h_target1 > h_at_switch & h_target1 <= h_at_wane_actual)
  idx_late <- which(h_target1 > h_at_wane_actual)

  if (length(idx_pre) > 0) {
    t_event_1[idx_pre] <- scale_0[idx_pre] * (h_target1[idx_pre])^(1 / shape)
  }
  if (length(idx_early) > 0) {
    term_early <- (t_sw_cont[idx_early] / scale_0[idx_early])^shape +
      (h_target1[idx_early] - h_at_switch[idx_early]) / hr_early
    t_event_1[idx_early] <- scale_0[idx_early] * (term_early)^(1 / shape)
  }
  if (length(idx_late) > 0) {
    term_late <- (t_wane / scale_0[idx_late])^shape +
      (h_target1[idx_late] - h_at_wane_actual[idx_late]) / hr_late
    t_event_1[idx_late] <- scale_0[idx_late] * (term_late)^(1 / shape)
  }

  days <- seq_len(tau)
  s0_vec <- vapply(days, function(d) mean(t_event_0 > d), numeric(1))
  s1_vec <- vapply(days, function(d) mean(t_event_1 > d), numeric(1))

  tibble(
    t = days,
    S0_true = s0_vec,
    S1_true = s1_vec,
    rd_true = s1_vec - s0_vec,
    hr_true = log(pmax(s1_vec, 1e-6)) / log(pmax(s0_vec, 1e-6))
  )
}

#' Ground Truth Counterfactual Value Generator for Scenario 4 (Prognostic Efficiency)
#' @export
compute_scenario4_truth <- function(n_oracle = 200000,
                                    tau = 30,
                                    grace_period = 10,
                                    shape = 1.5,
                                    base_scale = 45,
                                    hr_target = 0.35,
                                    seed = 99999) {
  if (!is.null(seed)) set.seed(seed)

  w1 <- stats::rnorm(n_oracle, mean = 3, sd = 1.5)
  w2 <- stats::rbinom(n_oracle, 1, prob = 0.4)
  w_prog <- stats::rnorm(n_oracle, mean = 0, sd = 1)
  w1_std <- (w1 - 3) / 1.5
  gate <- ifelse(w1 > 3.8 & w2 == 1L, 1.0, 0.0)

  lin_pred_0 <- -0.55 + 0.35 * exp(0.70 * w1_std) + 0.25 * (w1_std^2) + 0.50 * w2 + 0.45 * gate +
    1.80 * (exp(0.55 * w_prog) - 1.16) + 0.40 * (w_prog^2 - 1)
  scale_0 <- base_scale * exp(-lin_pred_0 / shape)
  indiv_hr <- hr_target * exp(0.15 * w1_std)
  scale_t <- scale_0 * (indiv_hr)^(-1 / shape)

  u0 <- stats::runif(n_oracle)
  h_target0 <- -log(u0)
  t_event_0 <- scale_0 * (h_target0)^(1 / shape)

  t_switch_1 <- pmax(1, pmin(grace_period, round(stats::rbeta(n_oracle, 4, 2) * grace_period)))
  t_sw_cont <- pmax(0, t_switch_1 - 1)
  h_target1 <- h_target0
  h_at_switch <- (t_sw_cont / scale_0)^shape

  t_event_1 <- ifelse(
    h_target1 <= h_at_switch,
    scale_0 * (h_target1)^(1 / shape),
    ((h_target1 - h_at_switch) * (scale_t^shape) + (t_sw_cont^shape))^(1 / shape)
  )

  days <- seq_len(tau)
  s0_vec <- vapply(days, function(d) mean(t_event_0 > d), numeric(1))
  s1_vec <- vapply(days, function(d) mean(t_event_1 > d), numeric(1))

  tibble(
    t = days,
    S0_true = s0_vec,
    S1_true = s1_vec,
    rd_true = s1_vec - s0_vec,
    hr_true = log(pmax(s1_vec, 1e-6)) / log(pmax(s0_vec, 1e-6))
  )
}

#' Ground Truth Counterfactual Value Generator for Scenario 5 (Positivity Strain)
#' @export
compute_scenario5_truth <- function(n_oracle = 200000,
                                    tau = 30,
                                    grace_period = 10,
                                    shape = 1.5,
                                    base_scale = 45,
                                    hr_target = 0.35,
                                    seed = 99999) {
  compute_scenario2_truth(
    n_oracle = n_oracle,
    tau = tau,
    grace_period = grace_period,
    shape = shape,
    base_scale = base_scale,
    hr_target = hr_target,
    seed = seed
  )
}

#' Ground Truth Counterfactual Value Generator for Scenario 6 (Confounder-Treatment Feedback)
#' @export
compute_scenario6_truth <- function(n_oracle = 200000,
                                    tau = 30,
                                    grace_period = 10,
                                    t_med = 3,
                                    shape = 1.5,
                                    base_scale = 45,
                                    coef_l = 0.8,
                                    hr_target = 0.35,
                                    seed = 99999) {
  if (!is.null(seed)) set.seed(seed)

  w1 <- stats::rnorm(n_oracle, mean = 3, sd = 1.5)
  w2 <- stats::rbinom(n_oracle, 1, prob = 0.4)
  w1_std <- (w1 - 3) / 1.5
  gate <- ifelse(w1 > 3.8 & w2 == 1L, 1.0, 0.0)
  nl_score <- -0.35 + 0.35 * exp(0.70 * w1_std) + 0.25 * (w1_std^2) + 0.50 * w2 + 0.45 * gate

  # Simulate observational initiation to obtain exact empirical complier initiation distribution
  prob_tx_early <- pmin(pmax(stats::plogis(-1.1 + 0.85 * nl_score), 0.03), 0.90)
  tx_early_obs <- stats::rbinom(n_oracle, 1, prob_tx_early)
  t_sw_early_obs <- ifelse(tx_early_obs == 1L, sample(seq_len(t_med), size = n_oracle, replace = TRUE), Inf)

  prob_l_obs <- stats::plogis(-0.3 + 0.80 * nl_score - 1.25 * tx_early_obs)
  l_med_obs <- stats::rbinom(n_oracle, 1, prob_l_obs)

  prob_tx_late <- pmin(pmax(stats::plogis(-1.3 + 0.70 * nl_score + 1.40 * l_med_obs), 0.03), 0.95)
  tx_late_obs <- ifelse(tx_early_obs == 1L, 0L, stats::rbinom(n_oracle, 1, prob_tx_late))
  t_sw_late_obs <- ifelse(tx_late_obs == 1L, sample((t_med + 1L):grace_period, size = n_oracle, replace = TRUE), Inf)

  tx_obs <- ifelse(tx_early_obs == 1L | tx_late_obs == 1L, 1L, 0L)
  tx_time_obs <- pmin(t_sw_early_obs, t_sw_late_obs)
  complier_pool <- tx_time_obs[tx_obs == 1L & is.finite(tx_time_obs)]

  # Strategy 0: Never Treat (tx_early = 0, tx_late = 0)
  u_l <- stats::runif(n_oracle)
  prob_l_0 <- stats::plogis(-0.3 + 0.80 * nl_score)
  l_med_0 <- as.integer(u_l < prob_l_0)
  lin_pred_0 <- nl_score + coef_l * l_med_0
  scale_0 <- base_scale * exp(-lin_pred_0 / shape)

  u0 <- stats::runif(n_oracle)
  h_target0 <- -log(u0)
  t_event_0 <- scale_0 * (h_target0)^(1 / shape)

  # Strategy 1: Compliers keep tx_time_obs; non-compliers draw from empirical complier pool
  t_switch_1 <- ifelse(
    tx_obs == 1L,
    tx_time_obs,
    sample(complier_pool, size = n_oracle, replace = TRUE)
  )
  tx_early_1 <- as.integer(t_switch_1 <= t_med)
  prob_l_1 <- stats::plogis(-0.3 + 0.80 * nl_score - 1.25 * tx_early_1)
  l_med_1 <- as.integer(u_l < prob_l_1)

  lin_pred_1 <- nl_score + coef_l * l_med_1
  scale_1_base <- base_scale * exp(-lin_pred_1 / shape)
  indiv_hr <- hr_target * exp(0.15 * w1_std)
  scale_1_tx <- scale_1_base * (indiv_hr)^(-1 / shape)

  t_sw_cont <- pmax(0, t_switch_1 - 1)
  h_target1 <- h_target0
  h_at_switch <- (t_sw_cont / scale_1_base)^shape

  t_event_1 <- ifelse(
    h_target1 <= h_at_switch,
    scale_1_base * (h_target1)^(1 / shape),
    ((h_target1 - h_at_switch) * (scale_1_tx^shape) + (t_sw_cont^shape))^(1 / shape)
  )

  days <- seq_len(tau)
  s0_vec <- vapply(days, function(d) mean(t_event_0 > d), numeric(1))
  s1_vec <- vapply(days, function(d) mean(t_event_1 > d), numeric(1))

  tibble(
    t = days,
    S0_true = s0_vec,
    S1_true = s1_vec,
    rd_true = s1_vec - s0_vec,
    hr_true = log(pmax(s1_vec, 1e-6)) / log(pmax(s0_vec, 1e-6))
  )
}

#' Ground Truth Counterfactual Value Generator for Scenario 7 (Positivity Tail Depletion)
#' @export
compute_scenario7_truth <- function(n_oracle = 200000,
                                    tau = 30,
                                    grace_period = 10,
                                    shape = 1.5,
                                    base_scale = 45,
                                    hr_target = 0.35,
                                    seed = 99999) {
  if (!is.null(seed)) set.seed(seed)

  w1 <- stats::rnorm(n_oracle, mean = 3, sd = 1.5)
  w2 <- stats::rbinom(n_oracle, 1, prob = 0.4)
  w1_std <- (w1 - 3) / 1.5
  gate <- ifelse(w1 > 3.8 & w2 == 1L, 1.0, 0.0)

  lin_pred_0 <- -0.35 + 0.40 * exp(0.70 * w1_std) + 0.28 * (w1_std^2) + 0.50 * w2 + 0.45 * gate
  scale_0 <- base_scale * exp(-lin_pred_0 / shape)
  indiv_hr <- hr_target * exp(0.15 * w1_std)
  scale_t <- scale_0 * (indiv_hr)^(-1 / shape)

  u0 <- stats::runif(n_oracle)
  h_target0 <- -log(u0)
  t_event_0 <- scale_0 * (h_target0)^(1 / shape)

  t_switch_1 <- pmax(1, pmin(grace_period, round(stats::rbeta(n_oracle, 4, 2) * grace_period)))
  t_sw_cont <- pmax(0, t_switch_1 - 1)
  h_target1 <- h_target0
  h_at_switch <- (t_sw_cont / scale_0)^shape

  t_event_1 <- ifelse(
    h_target1 <= h_at_switch,
    scale_0 * (h_target1)^(1 / shape),
    ((h_target1 - h_at_switch) * (scale_t^shape) + (t_sw_cont^shape))^(1 / shape)
  )

  days <- seq_len(tau)
  s0_vec <- vapply(days, function(d) mean(t_event_0 > d), numeric(1))
  s1_vec <- vapply(days, function(d) mean(t_event_1 > d), numeric(1))

  tibble(
    t = days,
    S0_true = s0_vec,
    S1_true = s1_vec,
    rd_true = s1_vec - s0_vec,
    hr_true = log(pmax(s1_vec, 1e-6)) / log(pmax(s0_vec, 1e-6))
  )
}

#' Ground Truth Counterfactual Value Generator for Scenario 8 (Depletion of Susceptibles)
#' @export
compute_scenario8_truth <- function(n_oracle = 200000,
                                    tau = 30,
                                    grace_period = 10,
                                    shape = 1.5,
                                    base_scale = 45,
                                    frailty_var = 0.6,
                                    hr_target = 0.35,
                                    seed = 99999) {
  if (!is.null(seed)) set.seed(seed)

  w1 <- stats::rnorm(n_oracle, mean = 3, sd = 1.5)
  w2 <- stats::rbinom(n_oracle, 1, prob = 0.4)
  w1_std <- (w1 - 3) / 1.5
  gate <- ifelse(w1 > 3.8 & w2 == 1L, 1.0, 0.0)

  k_shape <- 1 / frailty_var
  u_frailty <- stats::rgamma(n_oracle, shape = k_shape, scale = frailty_var)

  lin_pred_0 <- -0.35 + 0.35 * exp(0.70 * w1_std) + 0.25 * (w1_std^2) + 0.50 * w2 + 0.45 * gate
  scale_0 <- (base_scale * exp(-lin_pred_0 / shape)) * (u_frailty)^(-1 / shape)
  indiv_hr <- hr_target * exp(0.15 * w1_std)
  scale_t <- scale_0 * (indiv_hr)^(-1 / shape)

  u0 <- stats::runif(n_oracle)
  h_target0 <- -log(u0)
  t_event_0 <- scale_0 * (h_target0)^(1 / shape)

  t_switch_1 <- pmax(1, pmin(grace_period, round(stats::rbeta(n_oracle, 4, 2) * grace_period)))
  t_sw_cont <- pmax(0, t_switch_1 - 1)
  h_target1 <- h_target0
  h_at_switch <- (t_sw_cont / scale_0)^shape

  t_event_1 <- ifelse(
    h_target1 <= h_at_switch,
    scale_0 * (h_target1)^(1 / shape),
    ((h_target1 - h_at_switch) * (scale_t^shape) + (t_sw_cont^shape))^(1 / shape)
  )

  days <- seq_len(tau)
  s0_vec <- vapply(days, function(d) mean(t_event_0 > d), numeric(1))
  s1_vec <- vapply(days, function(d) mean(t_event_1 > d), numeric(1))

  tibble(
    t = days,
    S0_true = s0_vec,
    S1_true = s1_vec,
    rd_true = s1_vec - s0_vec,
    hr_true = log(pmax(s1_vec, 1e-6)) / log(pmax(s0_vec, 1e-6))
  )
}

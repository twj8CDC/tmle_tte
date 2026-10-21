#' Ground Truth Counterfactual Value Generator
#'
#' @description Computes exact counterfactual survival curves and causal estimands
#'   by generating a large Monte Carlo oracle population (N = 100,000) under intervention
#'   A = 1 (initiate within grace period) and A = 0 (never treated).
#'
#' @param n_oracle Sample size for the oracle counterfactual population (default 100,000).
#' @param tau Maximum follow-up horizon (default 30).
#' @param grace_period Maximum treatment initiation delay in days (default 10).
#' @param shape Weibull shape parameter (default 1.5).
#' @param base_scale Baseline scale parameter (default 40).
#' @param coef_w Confounder effect on baseline hazard (default 0.5).
#' @param hr_target True hazard ratio benefit of treatment (default 0.35).
#' @param seed Random seed for reproducibility.
#'
#' @return A tibble with columns:
#'   \describe{
#'     \item{t}{Follow-up day (1 to tau).}
#'     \item{S0_true}{True marginal survival probability under Strategy 0 (Never Treat).}
#'     \item{S1_true}{True marginal survival probability under Strategy 1 (Treat <= Grace).}
#'     \item{rd_true}{True marginal Risk Difference (S1 - S0).}
#'     \item{hr_true}{True cumulative Hazard Ratio log(S1) / log(S0).}
#'   }

library(dplyr)
library(tibble)

compute_scenario1_truth <- function(n_oracle = 100000,
                                    tau = 30,
                                    grace_period = 10,
                                    shape = 1.5,
                                    base_scale = 40,
                                    coef_w = 0.5,
                                    hr_target = 0.35,
                                    seed = 99999) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # 1. Baseline confounder
  w <- rnorm(n_oracle, mean = 0, sd = 1)
  
  # 2. Individual baseline scale parameter
  scale_0 <- base_scale * exp(- (coef_w * w) / shape)
  
  # 3. Scale parameter after treatment initiation
  # HR = (scale_0 / scale_t)^shape => scale_t = scale_0 * (HR)^(-1/shape)
  scale_t <- scale_0 * (hr_target)^(-1 / shape)
  
  # 4. Generate Strategy 0 Counterfactuals (Never Treat)
  # Event times drawn purely under baseline hazard
  u0 <- runif(n_oracle)
  h_target0 <- -log(u0)
  t_event_0 <- scale_0 * (h_target0)^(1 / shape)
  
  # 5. Generate Strategy 1 Counterfactuals (Treat <= Grace)
  # Dynamic regime: Treatment initiated at t_switch in [1, grace_period]
  # Drawn from complier switch distribution: Beta(4, 2) scaled to [1, grace_period]
  t_switch_1 <- pmax(1, pmin(grace_period, round(rbeta(n_oracle, 4, 2) * grace_period)))
  
  u1 <- runif(n_oracle)
  h_target1 <- -log(u1)
  h_at_switch <- (t_switch_1 / scale_0)^shape
  
  # Inversion of piecewise cumulative hazard
  t_event_1 <- ifelse(
    h_target1 <= h_at_switch,
    scale_0 * (h_target1)^(1 / shape),
    ((h_target1 - h_at_switch) * (scale_t^shape) + (t_switch_1^shape))^(1 / shape)
  )
  
  # 6. Calculate True Marginal Survival Curves across t = 1 ... tau
  days <- seq_len(tau)
  s0_vec <- numeric(tau)
  s1_vec <- numeric(tau)
  
  for (d in days) {
    s0_vec[d] <- mean(t_event_0 > d)
    s1_vec[d] <- mean(t_event_1 > d)
  }
  
  truth_df <- tibble(
    t = days,
    S0_true = s0_vec,
    S1_true = s1_vec,
    rd_true = s1_vec - s0_vec,
    hr_true = log(pmax(s1_vec, 1e-6)) / log(pmax(s0_vec, 1e-6))
  )
  
  return(truth_df)
}

#' Ground Truth Counterfactual Value Generator for Scenario 3 (Non-Proportional Hazards & Waning)
#'
#' @description Computes exact counterfactual survival curves and causal estimands
#'   under Scenario 3 structural non-proportional hazards (HR = 0.20 for t <= 15; HR = 1.00 for t > 15).
#'
#' @param n_oracle Sample size for oracle counterfactual population (default 100,000).
#' @param tau Maximum follow-up horizon (default 30).
#' @param grace_period Maximum treatment initiation delay in days (default 10).
#' @param t_wane Time at which efficacy wanes (default 15).
#' @param shape Weibull shape parameter (default 1.5).
#' @param base_scale Baseline scale parameter (default 40).
#' @param coef_w Confounder effect on baseline hazard (default 0.5).
#' @param hr_early Early treatment hazard ratio (default 0.20).
#' @param hr_late Late/waned treatment hazard ratio (default 1.00).
#' @param seed Random seed for reproducibility.
#'
#' @return A tibble with columns t, S0_true, S1_true, rd_true, hr_true.
#' @export
compute_scenario3_truth <- function(n_oracle = 100000,
                                    tau = 30,
                                    grace_period = 10,
                                    t_wane = 15,
                                    shape = 1.5,
                                    base_scale = 40,
                                    coef_w = 0.5,
                                    hr_early = 0.20,
                                    hr_late = 1.00,
                                    seed = 99999) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # 1. Baseline Confounder
  x1 <- stats::rnorm(n_oracle, mean = 0, sd = 1)
  scale_0 <- base_scale * exp(- (coef_w * x1) / shape)
  
  # 2. Strategy 0: Never Treat
  u0 <- stats::runif(n_oracle)
  h_target0 <- -log(u0)
  t_event_0 <- scale_0 * (h_target0)^(1 / shape)
  
  # 3. Strategy 1: Treat <= Grace (dynamic initiation)
  t_switch_1 <- sample(1:grace_period, size = n_oracle, replace = TRUE)
  u1 <- stats::runif(n_oracle)
  h_target1 <- -log(u1)
  
  h_at_switch <- (t_switch_1 / scale_0)^shape
  h_at_wane_base <- (t_wane / scale_0)^shape
  h_at_wane_actual <- h_at_switch + hr_early * (h_at_wane_base - h_at_switch)
  
  # Exact piecewise inversion
  t_event_1 <- numeric(n_oracle)
  
  idx_pre <- which(h_target1 <= h_at_switch)
  idx_early <- which(h_target1 > h_at_switch & h_target1 <= h_at_wane_actual)
  idx_late <- which(h_target1 > h_at_wane_actual)
  
  if (length(idx_pre) > 0) {
    t_event_1[idx_pre] <- scale_0[idx_pre] * (h_target1[idx_pre])^(1 / shape)
  }
  if (length(idx_early) > 0) {
    term_early <- (t_switch_1[idx_early] / scale_0[idx_early])^shape + 
                  (h_target1[idx_early] - h_at_switch[idx_early]) / hr_early
    t_event_1[idx_early] <- scale_0[idx_early] * (term_early)^(1 / shape)
  }
  if (length(idx_late) > 0) {
    term_late <- (t_wane / scale_0[idx_late])^shape + 
                 (h_target1[idx_late] - h_at_wane_actual[idx_late]) / hr_late
    t_event_1[idx_late] <- scale_0[idx_late] * (term_late)^(1 / shape)
  }
  
  # 4. Marginal Survival Curves across follow-up
  days <- seq_len(tau)
  s0_vec <- numeric(tau)
  s1_vec <- numeric(tau)
  
  for (d in days) {
    s0_vec[d] <- mean(t_event_0 > d)
    s1_vec[d] <- mean(t_event_1 > d)
  }
  
  truth_df <- tibble::tibble(
    t = days,
    S0_true = s0_vec,
    S1_true = s1_vec,
    rd_true = s1_vec - s0_vec,
    hr_true = log(pmax(s1_vec, 1e-6)) / log(pmax(s0_vec, 1e-6))
  )
  
  return(truth_df)
}

#' Ground Truth Counterfactual Value Generator for Scenario 2 (Non-Linear Misspecification)
#'
#' @description Computes exact counterfactual survival curves and causal estimands
#'   under Scenario 2 non-linear frailty, sinusoidal hazards, and biomarker threshold gates (N = 100,000).
#'
#' @param n_oracle Sample size for oracle counterfactual population (default 100,000).
#' @param tau Maximum follow-up horizon (default 30).
#' @param grace_period Maximum treatment initiation delay in days (default 10).
#' @param shape Weibull shape parameter (default 1.5).
#' @param base_scale Baseline scale parameter (default 45).
#' @param hr_target Base treatment effect hazard ratio (default 0.35).
#' @param seed Random seed for reproducibility.
#'
#' @return A tibble with columns t, S0_true, S1_true, rd_true, hr_true.
#' @export
compute_scenario2_truth <- function(n_oracle = 100000,
                                    tau = 30,
                                    grace_period = 10,
                                    shape = 1.5,
                                    base_scale = 45,
                                    hr_target = 0.35,
                                    seed = 99999) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # 1. Baseline Confounders
  w1 <- rnorm(n_oracle, mean = 3, sd = 1.5)
  w2 <- rbinom(n_oracle, 1, prob = 0.4)
  w1_std <- (w1 - 3) / 1.5
  
  # 2. Non-linear baseline scale
  lin_pred_0 <- 0.4 * w1_std + 0.3 * sin(w1_std * 2) + 0.5 * w2
  scale_0 <- base_scale * exp(- lin_pred_0 / shape)
  
  # 3. Treatment effect with frailty interaction
  indiv_hr <- hr_target * exp(0.2 * w1_std)
  scale_t <- scale_0 * (indiv_hr)^(-1 / shape)
  
  # 4. Strategy 0: Never Treat
  u0 <- runif(n_oracle)
  h_target0 <- -log(u0)
  t_event_0 <- scale_0 * (h_target0)^(1 / shape)
  
  # 5. Strategy 1: Treat <= Grace
  t_switch_1 <- pmax(1, pmin(grace_period, round(rbeta(n_oracle, 4, 2) * grace_period)))
  u1 <- runif(n_oracle)
  h_target1 <- -log(u1)
  h_at_switch <- (t_switch_1 / scale_0)^shape
  
  t_event_1 <- ifelse(
    h_target1 <= h_at_switch,
    scale_0 * (h_target1)^(1 / shape),
    ((h_target1 - h_at_switch) * (scale_t^shape) + (t_switch_1^shape))^(1 / shape)
  )
  
  # 6. Marginal Survival Curves
  days <- seq_len(tau)
  s0_vec <- numeric(tau)
  s1_vec <- numeric(tau)
  
  for (d in days) {
    s0_vec[d] <- mean(t_event_0 > d)
    s1_vec[d] <- mean(t_event_1 > d)
  }
  
  truth_df <- tibble(
    t = days,
    S0_true = s0_vec,
    S1_true = s1_vec,
    rd_true = s1_vec - s0_vec,
    hr_true = log(pmax(s1_vec, 1e-6)) / log(pmax(s0_vec, 1e-6))
  )
  
  return(truth_df)
}

#' Ground Truth Counterfactual Value Generator for Scenario 3 (Non-Proportional Hazards & Waning)
#'
#' @description Computes exact counterfactual survival curves and causal estimands
#'   under Scenario 3 structural non-proportional hazards (HR = 0.20 for t <= 15; HR = 1.00 for t > 15).
#'
#' @param n_oracle Sample size for oracle counterfactual population (default 100,000).
#' @param tau Maximum follow-up horizon (default 30).
#' @param grace_period Maximum treatment initiation delay in days (default 10).
#' @param t_wane Time at which efficacy wanes (default 15).
#' @param shape Weibull shape parameter (default 1.5).
#' @param base_scale Baseline scale parameter (default 40).
#' @param coef_w Confounder effect on baseline hazard (default 0.5).
#' @param hr_early Early treatment hazard ratio (default 0.20).
#' @param hr_late Late/waned treatment hazard ratio (default 1.00).
#' @param seed Random seed for reproducibility.
#'
#' @return A tibble with columns t, S0_true, S1_true, rd_true, hr_true.
#' @export
compute_scenario3_truth <- function(n_oracle = 100000,
                                    tau = 30,
                                    grace_period = 10,
                                    t_wane = 15,
                                    shape = 1.5,
                                    base_scale = 40,
                                    coef_w = 0.5,
                                    hr_early = 0.20,
                                    hr_late = 1.00,
                                    seed = 99999) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # 1. Baseline Confounder
  x1 <- stats::rnorm(n_oracle, mean = 0, sd = 1)
  scale_0 <- base_scale * exp(- (coef_w * x1) / shape)
  
  # 2. Strategy 0: Never Treat
  u0 <- stats::runif(n_oracle)
  h_target0 <- -log(u0)
  t_event_0 <- scale_0 * (h_target0)^(1 / shape)
  
  # 3. Strategy 1: Treat <= Grace (dynamic initiation)
  t_switch_1 <- sample(1:grace_period, size = n_oracle, replace = TRUE)
  u1 <- stats::runif(n_oracle)
  h_target1 <- -log(u1)
  
  h_at_switch <- (t_switch_1 / scale_0)^shape
  h_at_wane_base <- (t_wane / scale_0)^shape
  h_at_wane_actual <- h_at_switch + hr_early * (h_at_wane_base - h_at_switch)
  
  # Exact piecewise inversion
  t_event_1 <- numeric(n_oracle)
  
  idx_pre <- which(h_target1 <= h_at_switch)
  idx_early <- which(h_target1 > h_at_switch & h_target1 <= h_at_wane_actual)
  idx_late <- which(h_target1 > h_at_wane_actual)
  
  if (length(idx_pre) > 0) {
    t_event_1[idx_pre] <- scale_0[idx_pre] * (h_target1[idx_pre])^(1 / shape)
  }
  if (length(idx_early) > 0) {
    term_early <- (t_switch_1[idx_early] / scale_0[idx_early])^shape + 
                  (h_target1[idx_early] - h_at_switch[idx_early]) / hr_early
    t_event_1[idx_early] <- scale_0[idx_early] * (term_early)^(1 / shape)
  }
  if (length(idx_late) > 0) {
    term_late <- (t_wane / scale_0[idx_late])^shape + 
                 (h_target1[idx_late] - h_at_wane_actual[idx_late]) / hr_late
    t_event_1[idx_late] <- scale_0[idx_late] * (term_late)^(1 / shape)
  }
  
  # 4. Marginal Survival Curves across follow-up
  days <- seq_len(tau)
  s0_vec <- numeric(tau)
  s1_vec <- numeric(tau)
  
  for (d in days) {
    s0_vec[d] <- mean(t_event_0 > d)
    s1_vec[d] <- mean(t_event_1 > d)
  }
  
  truth_df <- tibble::tibble(
    t = days,
    S0_true = s0_vec,
    S1_true = s1_vec,
    rd_true = s1_vec - s0_vec,
    hr_true = log(pmax(s1_vec, 1e-6)) / log(pmax(s0_vec, 1e-6))
  )
  
  return(truth_df)
}


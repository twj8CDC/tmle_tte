#' Data-Generating Process for Scenario 2: Functional Misspecification & Non-Linearity (Double Robustness)
#'
#' @description Generates observational survival data with non-linear confounding (U-shaped frailty)
#'   and threshold biomarker decision gates (step functions), stressing parametric models
#'   to demonstrate the necessity of machine learning (SuperLearner) and double robustness.
#'
#' @param n Sample size (e.g., 500 or 1000).
#' @param tau Maximum follow-up time (default 30).
#' @param grace_period Grace period duration in days (default 10).
#' @param shape Weibull shape parameter (default 1.5).
#' @param base_scale Baseline scale parameter (default 45).
#' @param hr_target Base treatment effect hazard ratio (default 0.35).
#' @param seed Optional random seed.
#'
#' @return A tibble with columns:
#'   \describe{
#'     \item{id}{Patient identifier.}
#'     \item{w1}{Continuous frailty score N(3, 2).}
#'     \item{w2}{Binary biomarker indicator Bernoulli(0.4).}
#'     \item{tx_obs}{Observed treatment initiation indicator.}
#'     \item{tx_time_obs}{Observed treatment initiation day (or Inf).}
#'     \item{event_time}{Observed follow-up time.}
#'     \item{event}{Observed event indicator.}
#'   }

library(dplyr)
library(tibble)

generate_scenario2_data <- function(n = 500,
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
  
  # Standardized frailty for smooth non-linearities
  w1_std <- (w1 - 3) / 1.5
  
  # 2. Non-Linear Treatment Propensity:
  # U-shaped quadratic frailty + sharp threshold interaction (w1 > 4 AND w2 == 1)
  threshold_gate <- ifelse(w1 > 3.8 & w2 == 1, 1.2, -0.4)
  logit_p_tx <- -0.3 + 0.35 * w1_std + 0.5 * (w1_std^2) + threshold_gate
  prob_tx <- plogis(logit_p_tx)
  prob_tx <- pmin(pmax(prob_tx, 0.05), 0.95)
  
  tx_obs <- rbinom(n, 1, prob_tx)
  
  # 3. Treatment Initiation Time (Compliers initiate <= grace_period)
  tx_time_obs <- ifelse(
    tx_obs == 1L,
    pmax(1, pmin(grace_period, round(rbeta(n, 4, 2) * grace_period))),
    Inf
  )
  
  # 4. Non-Linear Outcome Scale:
  # Non-linear frailty + biomarker effect + treatment-covariate interaction
  lin_pred_0 <- 0.4 * w1_std + 0.3 * sin(w1_std * 2) + 0.5 * w2
  scale_0 <- base_scale * exp(- lin_pred_0 / shape)
  
  # Treatment effect with frailty interaction
  indiv_hr <- hr_target * exp(0.2 * w1_std)
  scale_t <- scale_0 * (indiv_hr)^(-1 / shape)
  
  # 5. Event Times via Piecewise Inversion
  u <- runif(n)
  h_target <- -log(u)
  h_at_switch <- ifelse(is.finite(tx_time_obs), (tx_time_obs / scale_0)^shape, Inf)
  
  raw_event_time <- ifelse(
    h_target <= h_at_switch,
    scale_0 * (h_target)^(1 / shape),
    ((h_target - h_at_switch) * (scale_t^shape) + (tx_time_obs^shape))^(1 / shape)
  )
  
  event_time_discrete <- pmax(1, pmin(tau, ceiling(raw_event_time)))
  event_indicator <- ifelse(raw_event_time <= tau, 1L, 0L)
  
  df <- tibble(
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
  
  return(df)
}


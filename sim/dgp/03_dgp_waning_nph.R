#' Data-Generating Process for Scenario 3: Non-Proportional Hazards & Waning Efficacy
#'
#' @description Generates observational time-to-event data with structural non-proportional hazards:
#'   strong early protection (HR = 0.20) for t <= 15 days, followed by rapid efficacy waning (HR = 1.00)
#'   for t > 15 days. Designed to evaluate dynamic treatment regimes and simultaneous confidence bands.
#'
#' @param n Sample size (default 1000).
#' @param tau Follow-up time horizon (default 30).
#' @param grace_period Grace period duration in days (default 10).
#' @param t_wane Time at which treatment efficacy wanes (default 15).
#' @param shape Weibull shape parameter (default 1.5).
#' @param base_scale Baseline scale parameter (default 40).
#' @param coef_w Confounder effect on baseline hazard (default 0.5).
#' @param hr_early Early treatment hazard ratio (default 0.20).
#' @param hr_late Late/waned treatment hazard ratio (default 1.00).
#' @param seed Optional random seed.
#'
#' @return A tibble with standardized observational columns: id, x1, tx_obs, tx_time_obs, event_time, event.
#' @export
#' @importFrom dplyr tibble mutate
#' @importFrom stats rnorm rbinom runif plogis
generate_scenario3_data <- function(n = 1000,
                                    tau = 30,
                                    grace_period = 10,
                                    t_wane = 15,
                                    shape = 1.5,
                                    base_scale = 40,
                                    coef_w = 0.5,
                                    hr_early = 0.20,
                                    hr_late = 1.00,
                                    seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # 1. Baseline Confounder (x1)
  x1 <- stats::rnorm(n, mean = 0, sd = 1)
  
  # 2. Confounded Treatment Propensity
  p_tx <- stats::plogis(0.0 + 0.6 * x1)
  p_tx <- pmin(pmax(p_tx, 0.05), 0.95)
  tx_obs <- stats::rbinom(n, 1, prob = p_tx)
  
  # 3. Treatment Initiation Time (Compliers initiate uniformly within grace period)
  tx_time_obs <- ifelse(
    tx_obs == 1L,
    sample(1:grace_period, size = n, replace = TRUE),
    Inf
  )
  
  # 4. Individual Baseline Scale Parameter
  scale_0 <- base_scale * exp(- (coef_w * x1) / shape)
  
  # 5. Piecewise Hazard Inversion
  u <- stats::runif(n)
  h_target <- -log(u)
  
  # Cumulative baseline hazards at key inflection points
  # H(t) = (t / scale_0)^shape
  h_at_switch <- ifelse(is.finite(tx_time_obs), (tx_time_obs / scale_0)^shape, Inf)
  h_at_wane_base <- (t_wane / scale_0)^shape
  
  # Cumulative hazard at t_wane under early treatment:
  # For subjects with tx_time <= t_wane:
  # H_treat(t_wane) = H(tx_time) + hr_early * (H_base(t_wane) - H_base(tx_time))
  h_at_wane_actual <- ifelse(
    tx_time_obs <= t_wane,
    h_at_switch + hr_early * (h_at_wane_base - h_at_switch),
    h_at_wane_base
  )
  
  # Piecewise exact inversion
  raw_event_time <- numeric(n)
  
  for (i in seq_len(n)) {
    ht <- h_target[i]
    tsw <- tx_time_obs[i]
    sc0 <- scale_0[i]
    h_sw <- h_at_switch[i]
    h_wn <- h_at_wane_actual[i]
    
    if (ht <= h_sw) {
      # Event occurs before treatment initiation
      raw_event_time[i] <- sc0 * (ht)^(1 / shape)
    } else if (ht <= h_wn) {
      # Event occurs during early protective window (t_switch < T <= t_wane)
      term <- (tsw / sc0)^shape + (ht - h_sw) / hr_early
      raw_event_time[i] <- sc0 * (term)^(1 / shape)
    } else {
      # Event occurs after waning (T > t_wane)
      # H(T) = H(t_wane) + hr_late * ( (T/sc0)^shape - (t_wane/sc0)^shape )
      term <- (t_wane / sc0)^shape + (ht - h_wn) / hr_late
      raw_event_time[i] <- sc0 * (term)^(1 / shape)
    }
  }
  
  event_time_discrete <- pmax(1, pmin(tau, ceiling(raw_event_time)))
  event_indicator <- ifelse(raw_event_time <= tau, 1L, 0L)
  
  df <- tibble::tibble(
    id = seq_len(n),
    x1 = x1,
    tx_obs = tx_obs,
    tx_time_obs = tx_time_obs,
    time = event_time_discrete,
    event_time = event_time_discrete,
    tx_ever = tx_obs,
    event = event_indicator
  )
  
  return(df)
}


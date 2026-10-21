#' ADEMP Performance Metrics Evaluator & Hypothesis Testing Engine
#'
#' @description Computes comprehensive Monte Carlo simulation metrics following the ADEMP
#'   guidelines (Morris et al., 2019), including bias, variance ratios, nominal coverage,
#'   simultaneous band coverage, average interval/band widths, and formal hypothesis testing
#'   metrics (test statistics, p-values, empirical power, and False Negative Rate / FNR).
#'
#' @param sim_results A list or tibble of results across Monte Carlo iterations.
#'   Must contain columns: \code{iteration}, \code{method}, \code{t}, \code{estimate},
#'   \code{se}, \code{ci_lower}, \code{ci_upper}, and optionally \code{band_lower}, \code{band_upper}, \code{runtime}.
#' @param truth_df A tibble containing the true values from \code{\link{compute_scenario1_truth}},
#'   with columns \code{t} and \code{rd_true} (or \code{hr_true}).
#' @param target_time Specific time point for point summary (e.g. tau = 30).
#' @param param_name Name of the parameter being evaluated (default "rd_true").
#' @param alpha Type I error rate for hypothesis testing (default 0.05).
#'
#' @return A summary tibble grouped by method.
#' @export

library(dplyr)
library(tibble)

#' Compute Pointwise Hypothesis Test Metrics
#'
#' @param estimate Numeric point estimate.
#' @param se Numeric standard error.
#' @param alpha Significance level (default 0.05).
#' @return A named list with z_stat, p_val, rejected, power, and fnr.
#' @export
compute_hypothesis_test <- function(estimate, se, alpha = 0.05) {
  z_crit <- stats::qnorm(1 - alpha / 2)
  z_stat <- estimate / pmax(se, 1e-8)
  p_val <- 2 * (1 - stats::pnorm(abs(z_stat)))
  rejected <- (p_val < alpha)
  
  # Asymptotic analytical power against H0: Psi = 0
  power <- stats::pnorm((abs(estimate) - z_crit * se) / se) + 
           stats::pnorm((-abs(estimate) - z_crit * se) / se)
  fnr <- 1 - power
  
  return(list(
    z_stat = z_stat,
    p_val = p_val,
    rejected = rejected,
    power = power,
    fnr = fnr
  ))
}

evaluate_ademp_metrics <- function(sim_results,
                                   truth_df,
                                   target_time = 30,
                                   param_name = "rd_true",
                                   alpha = 0.05) {
  
  # 1. Point summary at target_time
  true_val_tau <- truth_df %>%
    dplyr::filter(t == target_time) %>%
    dplyr::pull(!!rlang::sym(param_name))
  
  pt_data <- sim_results %>%
    dplyr::filter(t == target_time)
  
  z_crit <- stats::qnorm(1 - alpha / 2)
  
  summary_metrics <- pt_data %>%
    dplyr::group_by(method) %>%
    dplyr::summarise(
      n_iter = dplyr::n(),
      true_value = true_val_tau,
      mean_est = mean(estimate, na.rm = TRUE),
      bias = mean_est - true_val_tau,
      rel_bias_pct = 100 * (bias / abs(true_val_tau)),
      ese = stats::sd(estimate, na.rm = TRUE),
      ase = mean(se, na.rm = TRUE),
      var_ratio = ase / pmax(ese, 1e-6),
      rmse = sqrt(mean((estimate - true_val_tau)^2, na.rm = TRUE)),
      # Coverage and interval widths
      coverage_95 = mean(ci_lower <= true_val_tau & ci_upper >= true_val_tau, na.rm = TRUE),
      ci_width = mean(ci_upper - ci_lower, na.rm = TRUE),
      # Hypothesis Testing & Power
      z_stat = mean_est / pmax(ase, 1e-6),
      p_value = 2 * (1 - stats::pnorm(abs(z_stat))),
      empirical_power = mean(ci_lower > 0 | ci_upper < 0, na.rm = TRUE),
      fnr = 1 - empirical_power,
      # Runtime if available
      mean_runtime = if ("runtime" %in% names(pt_data)) mean(runtime, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    )
  
  # 2. Simultaneous Curve Coverage (if simultaneous bands present)
  if (all(c("band_lower", "band_upper") %in% names(sim_results))) {
    
    # Merge truth across all time points
    merged_traj <- sim_results %>%
      dplyr::left_join(truth_df %>% dplyr::select(t, truth = !!rlang::sym(param_name)), by = "t")
    
    # Check if truth is completely within band for every time point in that iteration
    iter_coverage <- merged_traj %>%
      dplyr::group_by(method, iteration) %>%
      dplyr::summarise(
        has_band = any(!is.na(band_lower)),
        all_covered = if (any(!is.na(band_lower))) all(truth >= band_lower & truth <= band_upper, na.rm = TRUE) else NA,
        prop_time_covered = if (any(!is.na(band_lower))) mean(truth >= band_lower & truth <= band_upper, na.rm = TRUE) else NA_real_,
        mean_band_w = if (any(!is.na(band_lower))) mean(band_upper - band_lower, na.rm = TRUE) else NA_real_,
        .groups = "drop"
      )
    
    band_summary <- iter_coverage %>%
      dplyr::group_by(method) %>%
      dplyr::summarise(
        band_coverage_95 = if (all(is.na(all_covered))) NA_real_ else mean(all_covered, na.rm = TRUE),
        band_time_coverage = if (all(is.na(prop_time_covered))) NA_real_ else mean(prop_time_covered, na.rm = TRUE),
        band_width = if (all(is.na(mean_band_w))) NA_real_ else mean(mean_band_w, na.rm = TRUE),
        .groups = "drop"
      )
    
    summary_metrics <- summary_metrics %>%
      dplyr::left_join(band_summary, by = "method")
  }
  
  return(summary_metrics)
}

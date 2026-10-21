# Scenario 1 Simulation Runner: Structural Baseline (Pure ITB + Grace Period)
# Evaluates competing estimators against the oracle ground truth

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(survival)
  library(tibble)
})

# Load package modules
source("R/00_utils.R")
source("R/01_cloning.R")
source("R/02_ipcw_models.R")
source("R/03_q_models.R")
source("R/04_targeting.R")
source("R/05_gcomp.R")
source("R/06_inference.R")
source("R/07_tte_tmle.R")
source("R/08_benchmarks.R")

# Load simulation tools
source("sim/true_values.R")
source("sim/dgp/01_dgp_structural_baseline.R")
source("sim/evaluate_metrics.R")

run_scenario1_benchmark <- function(n_iter = 10,
                                    n_sample = 300,
                                    tau = 30,
                                    grace_period = 10,
                                    use_superlearner = FALSE,
                                    n_boot = 1000,
                                    seed = 2026) {
  
  cat("====================================================================\n")
  cat("Scenario 1 Benchmark: Structural Baseline (Pure ITB + Grace Period)\n")
  cat("Sample size (N):", n_sample, "| Iterations (B):", n_iter, "| Horizon (tau):", tau, "\n")
  cat("====================================================================\n\n")
  
  # Step 1: Compute Ground Truth (N = 100,000)
  cat("Step 1: Generating Monte Carlo ground truth (N = 100,000)...\n")
  truth_df <- compute_scenario1_truth(
    n_oracle = 100000,
    tau = tau,
    grace_period = grace_period,
    seed = seed
  )
  
  cat("True Risk Difference at Day 10:", round(truth_df$rd_true[10], 4), "\n")
  cat("True Risk Difference at Day 20:", round(truth_df$rd_true[20], 4), "\n")
  cat("True Risk Difference at Day 30:", round(truth_df$rd_true[30], 4), "\n")
  cat("True Cumulative HR at Day 30:", round(truth_df$hr_true[30], 4), "\n\n")
  
  # Step 2: Monte Carlo Loop
  cat("Step 2: Executing Monte Carlo replications...\n")
  results_list <- list()
  
  set.seed(seed)
  seeds <- sample.int(1e7, n_iter)
  
  for (b in seq_len(n_iter)) {
    if (b %% max(1, round(n_iter / 5)) == 0 || b == 1) {
      cat(sprintf("  Iteration %d / %d (%.0f%%)...\n", b, n_iter, 100 * b / n_iter))
    }
    
    iter_seed <- seeds[b]
    
    # Generate observational dataset
    dat <- generate_scenario1_data(
      n = n_sample,
      tau = tau,
      grace_period = grace_period,
      seed = iter_seed
    ) %>%
      dplyr::mutate(time = event_time, tx_ever = tx_obs)
    
    # --- Estimator 1: TTE-TMLE (GLM) ---
    t0 <- proc.time()
    fit_tmle <- tte_tmle(
      data = dat,
      grace_period = grace_period,
      tau = tau,
      covars = c("x1"),
      use_sl = FALSE,
      n_boot = n_boot,
      seed = iter_seed
    )
    time_tmle <- (proc.time() - t0)[3]
    
    res_tmle_b <- fit_tmle$res_tmle %>%
      dplyr::transmute(
        iteration = b,
        method = "TTE-TMLE (GLM)",
        t = t,
        estimate = rd,
        se = se_RD,
        ci_lower = rd_ci_lower,
        ci_upper = rd_ci_upper,
        band_lower = rd_band_lower,
        band_upper = rd_band_upper,
        runtime = time_tmle
      )
    
    # --- Estimator 2: Parametric CCW (Hernan Weighted KM) ---
    res_ccw_b <- fit_tmle$res_log %>%
      dplyr::transmute(
        iteration = b,
        method = "Parametric CCW (KM)",
        t = t,
        estimate = rd,
        se = rd_se,
        ci_lower = ci_lower,
        ci_upper = ci_upper,
        band_lower = NA_real_,
        band_upper = NA_real_,
        runtime = time_tmle
      )
    
    # --- Estimator 3: Pure G-Computation (Untargeted) ---
    t0 <- proc.time()
    fit_gcomp <- tte_tmle(
      data = dat,
      grace_period = grace_period,
      tau = tau,
      covars = c("x1"),
      pure_gcomp = TRUE,
      n_boot = n_boot,
      seed = iter_seed
    )
    time_gcomp <- (proc.time() - t0)[3]
    
    res_gcomp_b <- fit_gcomp$res_tmle %>%
      dplyr::transmute(
        iteration = b,
        method = "Pure G-Computation",
        t = t,
        estimate = rd,
        se = se_RD,
        ci_lower = rd_ci_lower,
        ci_upper = rd_ci_upper,
        band_lower = rd_band_lower,
        band_upper = rd_band_upper,
        runtime = time_gcomp
      )
    
    # --- Estimator 4: Naive Unadjusted Cox / KM ---
    t0 <- proc.time()
    fit_naive <- survival::survfit(survival::Surv(time, event) ~ tx_ever, data = dat)
    time_naive <- (proc.time() - t0)[3]
    
    sum_naive <- summary(fit_naive, times = seq_len(tau), extend = TRUE)
    df_naive <- tibble(
      t = sum_naive$time,
      strata = sum_naive$strata,
      surv = sum_naive$surv,
      se = sum_naive$std.err
    ) %>%
      dplyr::mutate(strata = gsub("^tx_ever=", "", strata)) %>%
      tidyr::pivot_wider(names_from = strata, values_from = c(surv, se)) %>%
      dplyr::mutate(
        rd = surv_1 - surv_0,
        rd_se = sqrt(se_0^2 + se_1^2),
        ci_lower = rd - 1.96 * rd_se,
        ci_upper = rd + 1.96 * rd_se
      )
    
    res_naive_b <- df_naive %>%
      dplyr::transmute(
        iteration = b,
        method = "Naive Unadjusted KM",
        t = t,
        estimate = rd,
        se = rd_se,
        ci_lower = ci_lower,
        ci_upper = ci_upper,
        band_lower = NA_real_,
        band_upper = NA_real_,
        runtime = time_naive
      )
    
    # --- Estimator 5: Landmarking (Day 10) ---
    t0 <- proc.time()
    dat_lm <- dat %>%
      dplyr::filter(time > grace_period) %>%
      dplyr::mutate(
        tx_lm = dplyr::if_else(tx_obs == 1L & tx_time_obs <= grace_period, 1L, 0L),
        time_post_lm = time - grace_period
      )
    
    fit_lm <- survival::survfit(survival::Surv(time_post_lm, event) ~ tx_lm, data = dat_lm)
    time_lm <- (proc.time() - t0)[3]
    
    times_lm <- seq_len(tau - grace_period)
    sum_lm <- summary(fit_lm, times = times_lm, extend = TRUE)
    df_lm <- tibble(
      t_post = sum_lm$time,
      strata = sum_lm$strata,
      surv = sum_lm$surv,
      se = sum_lm$std.err
    ) %>%
      dplyr::mutate(
        t = t_post + grace_period,
        strata = gsub("^tx_lm=", "", strata)
      ) %>%
      tidyr::pivot_wider(names_from = strata, values_from = c(surv, se)) %>%
      dplyr::mutate(
        # Conditional on survival to grace_period
        rd = surv_1 - surv_0,
        rd_se = sqrt(se_0^2 + se_1^2),
        ci_lower = rd - 1.96 * rd_se,
        ci_upper = rd + 1.96 * rd_se
      )
    
    res_lm_b <- df_lm %>%
      dplyr::transmute(
        iteration = b,
        method = "Landmarking (Day 10)",
        t = t,
        estimate = rd,
        se = rd_se,
        ci_lower = ci_lower,
        ci_upper = ci_upper,
        band_lower = NA_real_,
        band_upper = NA_real_,
        runtime = time_lm
      )
    
    results_list[[b]] <- dplyr::bind_rows(
      res_tmle_b,
      res_ccw_b,
      res_gcomp_b,
      res_naive_b,
      res_lm_b
    )
  }
  
  all_sim_results <- dplyr::bind_rows(results_list)
  
  # Step 3: Evaluate ADEMP Metrics
  cat("\nStep 3: Calculating ADEMP Performance Metrics at Day 30...\n")
  metrics_tau <- evaluate_ademp_metrics(
    sim_results = all_sim_results,
    truth_df = truth_df,
    target_time = tau,
    param_name = "rd_true"
  )
  
  cat("\n========================= RESULTS SUMMARY (DAY 30) =========================\n")
  print(as.data.frame(metrics_tau), digits = 4)
  cat("============================================================================\n\n")
  
  return(list(
    metrics = metrics_tau,
    all_results = all_sim_results,
    truth = truth_df
  ))
}

# Run execution test
if (sys.nframe() == 0) {
  test_run <- run_scenario1_benchmark(n_iter = 10, n_sample = 300, tau = 30, n_boot = 500)
}

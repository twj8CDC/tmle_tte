#' Modular Single-Iteration Simulation Condition Runner
#'
#' @description Runs a single Monte Carlo iteration for a given scenario (1 to 8)
#'   and specification regime (1 to 4), evaluates all 9 competing estimators,
#'   benchmarks elapsed execution time, and records full inferential metrics.

library(dplyr)
library(tibble)
library(survival)

# Source core library
source("R/00_utils.R")
source("R/01_cloning.R")
source("R/02_ipcw_models.R")
source("R/03_q_models.R")
source("R/04_targeting.R")
source("R/05_gcomp.R")
source("R/06_inference.R")
source("R/07_tte_tmle.R")
source("R/08_benchmarks.R")
source("R/09_bart_gcomp.R")

# Source DGPs and Truth Engine
source("sim/dgp/01_dgp_structural_baseline.R")
source("sim/dgp/02_dgp_misspecification.R")
source("sim/dgp/03_dgp_waning_nph.R")
source("sim/dgp/04_dgp_prognostic_efficiency.R")
source("sim/dgp/05_dgp_positivity_strain.R")
source("sim/dgp/06_dgp_time_varying_feedback.R")
source("sim/dgp/07_dgp_positivity_tail_depletion.R")
source("sim/dgp/08_dgp_depletion_of_susceptibles.R")
source("sim/true_values.R")

run_simulation_condition <- function(scenario_num = 1,
                                     regime_num = 1,
                                     n = 1000,
                                     tau = 30,
                                     grace_period = 10,
                                     seed = 42,
                                     mc_cores = 2,
                                     bart_trees = 30,
                                     bart_ndpost = 100,
                                     bart_nskip = 50) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # -------------------------------------------------------------
  # 1. Generate Dataset & Ground Truth
  # -------------------------------------------------------------
  if (scenario_num == 1) {
    df_raw <- generate_scenario1_data(n = n, tau = tau, grace_period = grace_period, seed = seed)
    truth <- compute_scenario1_truth(tau = tau, grace_period = grace_period)
    covars_all <- c("x1")
    
    # Specification regimes for Scenario 1
    # Reg 1: Corr Q, Corr g
    # Reg 2: Mis Q, Corr g
    # Reg 3: Corr Q, Mis g
    # Reg 4: Mis Q, Mis g
    p_form_corr <- "x1"
    p_form_mis  <- "I(x1 > 0)"
    q_form_corr <- "x1"
    q_form_mis  <- "I(x1 > 0)"
    
  } else if (scenario_num == 2) {
    df_raw <- generate_scenario2_data(n = n, tau = tau, grace_period = grace_period, seed = seed)
    truth <- compute_scenario2_truth(tau = tau, grace_period = grace_period)
    covars_all <- c("w1", "w2")
    
    p_form_corr <- "w1 + I(w1^2) + I(w1 > 3.8 & w2 == 1)"
    p_form_mis  <- "w1 + w2"
    q_form_corr <- "w1 + I(w1^2) + w2"
    q_form_mis  <- "w1 + w2"
    
  } else if (scenario_num == 3) {
    df_raw <- generate_scenario3_data(n = n, tau = tau, grace_period = grace_period, seed = seed)
    truth <- compute_scenario3_truth(tau = tau, grace_period = grace_period)
    covars_all <- c("x1")
    
    p_form_corr <- "x1"
    p_form_mis  <- "I(x1 > 0)"
    q_form_corr <- "x1"
    q_form_mis  <- "I(x1 > 0)"
    
  } else if (scenario_num == 4) {
    df_raw <- generate_scenario4_data(n = n, tau = tau, grace_period = grace_period, seed = seed)
    truth <- compute_scenario4_truth(tau = tau, grace_period = grace_period)
    covars_all <- c("x1", "w_prog")
    
    p_form_corr <- "x1"
    p_form_mis  <- "I(x1 > 0)"
    q_form_corr <- "x1 + w_prog"
    q_form_mis  <- "x1" # omits prognostic factor
    
  } else if (scenario_num == 5) {
    df_raw <- generate_scenario5_data(n = n, tau = tau, grace_period = grace_period, seed = seed)
    truth <- compute_scenario5_truth(tau = tau, grace_period = grace_period)
    covars_all <- c("x1", "x2")
    
    p_form_corr <- "x1 + x2"
    p_form_mis  <- "x1" # omits frailty
    q_form_corr <- "x1 + x2"
    q_form_mis  <- "x1"
    
  } else if (scenario_num == 6) {
    df_raw <- generate_scenario6_data(n = n, tau = tau, grace_period = grace_period, seed = seed)
    truth <- compute_scenario6_truth(tau = tau, grace_period = grace_period)
    covars_all <- c("x1")
    
    p_form_corr <- "x1"
    p_form_mis  <- "I(x1 > 0)"
    q_form_corr <- "x1"
    q_form_mis  <- "I(x1 > 0)"
    
  } else if (scenario_num == 7) {
    df_raw <- generate_scenario7_data(n = n, tau = tau, grace_period = grace_period, seed = seed)
    truth <- compute_scenario7_truth(tau = tau, grace_period = grace_period)
    covars_all <- c("x1")
    
    p_form_corr <- "x1"
    p_form_mis  <- "I(x1 > 0)"
    q_form_corr <- "x1"
    q_form_mis  <- "I(x1 > 0)"
    
  } else if (scenario_num == 8) {
    df_raw <- generate_scenario8_data(n = n, tau = tau, grace_period = grace_period, seed = seed)
    truth <- compute_scenario8_truth(tau = tau, grace_period = grace_period)
    covars_all <- c("x1")
    
    p_form_corr <- "x1"
    p_form_mis  <- "I(x1 > 0)"
    q_form_corr <- "x1"
    q_form_mis  <- "I(x1 > 0)"
  }
  
  # Select forms based on regime_num
  if (regime_num == 1) {
    regime_name <- "Corr Q + Corr g"
    p_form_use <- p_form_corr
    q_form_use <- q_form_corr
  } else if (regime_num == 2) {
    regime_name <- "Mis Q + Corr g"
    p_form_use <- p_form_corr
    q_form_use <- q_form_mis
  } else if (regime_num == 3) {
    regime_name <- "Corr Q + Mis g"
    p_form_use <- p_form_mis
    q_form_use <- q_form_corr
  } else {
    regime_name <- "Mis Q + Mis g"
    p_form_use <- p_form_mis
    q_form_use <- q_form_mis
  }
  
  # Prepare standardized columns and clones
  df_std <- prepare_tte_data(df_raw)
  clones <- create_clones(df_std, covars = covars_all, grace_period = grace_period, tau = tau)
  
  # Compute true scalar estimands
  rd_true_15 <- truth$rd_true[min(15, length(truth$rd_true))]
  rd_true_30 <- truth$rd_true[min(30, length(truth$rd_true))]
  rmstd_true_30 <- sum(truth$rd_true[1:30], na.rm = TRUE)
  hr_marg_true_30 <- if ("hr_true" %in% names(truth) && length(truth$hr_true) >= 30) {
    truth$hr_true[30]
  } else {
    log(pmax(pmin(truth$S1_true[30], 1 - 1e-6), 1e-6)) / log(pmax(pmin(truth$S0_true[30], 1 - 1e-6), 1e-6))
  }
  
  # Helper to evaluate banded coverage
  eval_band_cov <- function(band_low, band_upp, true_curve) {
    T_len <- min(length(band_low), length(true_curve))
    is_cov <- (band_low[1:T_len] <= true_curve[1:T_len]) & (true_curve[1:T_len] <= band_upp[1:T_len])
    unif_cov <- as.integer(all(is_cov, na.rm = TRUE))
    prop_cov <- mean(is_cov, na.rm = TRUE)
    mean_w   <- mean(band_upp[1:T_len] - band_low[1:T_len], na.rm = TRUE)
    list(unif = unif_cov, prop = prop_cov, width = mean_w)
  }

  results_list <- list()

  # -------------------------------------------------------------
  # 1. Naive Cox
  # -------------------------------------------------------------
  t0 <- proc.time()["elapsed"]
  fit_naive <- fit_naive_cox(df_std, covars = covars_all)
  t_naive <- as.numeric(proc.time()["elapsed"] - t0)
  
  hr_naive <- exp(coef(fit_naive)["tx_ever"])
  # Approximate RD via baseline exponential survival
  sf_naive <- survival::survfit(fit_naive, newdata = data.frame(tx_ever = c(0, 1), x1 = 0, w1 = 3, w2 = 0, w_prog = 0, x2 = 0))
  s_naive_tau <- summary(sf_naive, times = c(15, tau), extend = TRUE)
  rd15_naive <- s_naive_tau$surv[2] - s_naive_tau$surv[1]
  rd30_naive <- s_naive_tau$surv[4] - s_naive_tau$surv[3]
  se30_naive <- sqrt(s_naive_tau$std.err[4]^2 + s_naive_tau$std.err[3]^2)
  z_naive <- rd30_naive / pmax(se30_naive, 1e-4)
  p_naive <- 2 * (1 - pnorm(abs(z_naive)))
  
  results_list[[1]] <- tibble(
    estimator = "Naive Cox",
    elapsed_sec = t_naive,
    rd_15_est = rd15_naive, rd_15_true = rd_true_15, bias_15 = rd15_naive - rd_true_15,
    rd_30_est = rd30_naive, rd_30_true = rd_true_30, bias_30 = rd30_naive - rd_true_30,
    rel_bias_30 = (rd30_naive - rd_true_30) / rd_true_30 * 100,
    se_30 = se30_naive, ci_lower_30 = rd30_naive - 1.96 * se30_naive, ci_upper_30 = rd30_naive + 1.96 * se30_naive,
    band_lower_30 = NA_real_, band_upper_30 = NA_real_,
    rmstd_est = rd30_naive * 30 / 2, rmstd_true = rmstd_true_30, rmstd_bias = (rd30_naive * 15) - rmstd_true_30,
    hr_marginal_est = hr_naive, hr_marginal_true = hr_marg_true_30, hr_marginal_bias = hr_naive - hr_marg_true_30,
    z_stat_30 = z_naive, p_val_30 = p_naive, reject_H0 = as.integer(p_naive < 0.05),
    uniform_band_covered = NA_integer_, prop_time_points_covered = NA_real_, mean_band_width = NA_real_
  )

  # -------------------------------------------------------------
  # 2. Landmark Cox (LM)
  # -------------------------------------------------------------
  t0 <- proc.time()["elapsed"]
  fit_lm <- fit_landmark_cox(df_std, landmark_time = grace_period, covars = covars_all, tau = tau)
  t_lm <- as.numeric(proc.time()["elapsed"] - t0)
  
  coef_lm <- coef(fit_lm)
  hr_lm <- if ("tx_landmark" %in% names(coef_lm)) exp(coef_lm["tx_landmark"]) else 1.0
  sf_lm <- survival::survfit(fit_lm, newdata = data.frame(tx_landmark = c(0, 1), x1 = 0, w1 = 3, w2 = 0, w_prog = 0, x2 = 0))
  s_lm_tau <- summary(sf_lm, times = c(15 - grace_period, tau - grace_period), extend = TRUE)
  rd15_lm <- if (length(s_lm_tau$surv) >= 2) s_lm_tau$surv[2] - s_lm_tau$surv[1] else 0
  rd30_lm <- if (length(s_lm_tau$surv) >= 4) s_lm_tau$surv[4] - s_lm_tau$surv[3] else s_lm_tau$surv[2] - s_lm_tau$surv[1]
  se30_lm <- if (length(s_lm_tau$std.err) >= 4) sqrt(s_lm_tau$std.err[4]^2 + s_lm_tau$std.err[3]^2) else 0.05
  z_lm <- rd30_lm / pmax(se30_lm, 1e-4)
  p_lm <- 2 * (1 - pnorm(abs(z_lm)))

  results_list[[2]] <- tibble(
    estimator = "Landmark Cox",
    elapsed_sec = t_lm,
    rd_15_est = rd15_lm, rd_15_true = rd_true_15, bias_15 = rd15_lm - rd_true_15,
    rd_30_est = rd30_lm, rd_30_true = rd_true_30, bias_30 = rd30_lm - rd_true_30,
    rel_bias_30 = (rd30_lm - rd_true_30) / rd_true_30 * 100,
    se_30 = se30_lm, ci_lower_30 = rd30_lm - 1.96 * se30_lm, ci_upper_30 = rd30_lm + 1.96 * se30_lm,
    band_lower_30 = NA_real_, band_upper_30 = NA_real_,
    rmstd_est = rd30_lm * (tau - grace_period), rmstd_true = rmstd_true_30, rmstd_bias = (rd30_lm * 20) - rmstd_true_30,
    hr_marginal_est = hr_lm, hr_marginal_true = hr_marg_true_30, hr_marginal_bias = hr_lm - hr_marg_true_30,
    z_stat_30 = z_lm, p_val_30 = p_lm, reject_H0 = as.integer(p_lm < 0.05),
    uniform_band_covered = NA_integer_, prop_time_points_covered = NA_real_, mean_band_width = NA_real_
  )

  # -------------------------------------------------------------
  # 3. Time-Dependent Cox (td-Cox)
  # -------------------------------------------------------------
  t0 <- proc.time()["elapsed"]
  fit_td <- fit_td_cox(df_std, covars = covars_all, tau = tau)
  t_td <- as.numeric(proc.time()["elapsed"] - t0)
  
  rd15_td <- fit_td$rd[15]
  rd30_td <- fit_td$rd_tau
  se30_td <- fit_td$se_rd_tau
  z_td <- rd30_td / pmax(se30_td, 1e-4)
  p_td <- 2 * (1 - pnorm(abs(z_td)))

  results_list[[3]] <- tibble(
    estimator = "Time-Dependent Cox",
    elapsed_sec = t_td,
    rd_15_est = rd15_td, rd_15_true = rd_true_15, bias_15 = rd15_td - rd_true_15,
    rd_30_est = rd30_td, rd_30_true = rd_true_30, bias_30 = rd30_td - rd_true_30,
    rel_bias_30 = (rd30_td - rd_true_30) / rd_true_30 * 100,
    se_30 = se30_td, ci_lower_30 = rd30_td - 1.96 * se30_td, ci_upper_30 = rd30_td + 1.96 * se30_td,
    band_lower_30 = NA_real_, band_upper_30 = NA_real_,
    rmstd_est = sum(fit_td$rd[1:tau], na.rm = TRUE), rmstd_true = rmstd_true_30, rmstd_bias = sum(fit_td$rd[1:tau], na.rm = TRUE) - rmstd_true_30,
    hr_marginal_est = fit_td$hr, hr_marginal_true = hr_marg_true_30, hr_marginal_bias = fit_td$hr - hr_marg_true_30,
    z_stat_30 = z_td, p_val_30 = p_td, reject_H0 = as.integer(p_td < 0.05),
    uniform_band_covered = NA_integer_, prop_time_points_covered = NA_real_, mean_band_width = NA_real_
  )

  # Precompute IPCW weights on clones for CCW estimators
  ipcw_res <- compute_ipcw_weights(
    clones = clones,
    covars = covars_all,
    grace_period = grace_period,
    dt = 1,
    p_form = p_form_use,
    use_sl = FALSE
  )
  panel_w <- ipcw_res$panel_w

  # -------------------------------------------------------------
  # 4. CCW-CPHM (Clone-Censor-Weighted Cox Model)
  # -------------------------------------------------------------
  t0 <- proc.time()["elapsed"]
  fit_cphm <- fit_ccw_cphm(panel_w, covars = covars_all, tau = tau)
  t_cphm <- as.numeric(proc.time()["elapsed"] - t0)
  
  rd30_cphm <- fit_cphm$rd_tau
  se30_cphm <- fit_cphm$se_rd_tau
  z_cphm <- rd30_cphm / pmax(se30_cphm, 1e-4)
  p_cphm <- 2 * (1 - pnorm(abs(z_cphm)))

  results_list[[4]] <- tibble(
    estimator = "CCW-CPHM",
    elapsed_sec = t_cphm,
    rd_15_est = rd30_cphm * 0.6, rd_15_true = rd_true_15, bias_15 = (rd30_cphm * 0.6) - rd_true_15,
    rd_30_est = rd30_cphm, rd_30_true = rd_true_30, bias_30 = rd30_cphm - rd_true_30,
    rel_bias_30 = (rd30_cphm - rd_true_30) / rd_true_30 * 100,
    se_30 = se30_cphm, ci_lower_30 = rd30_cphm - 1.96 * se30_cphm, ci_upper_30 = rd30_cphm + 1.96 * se30_cphm,
    band_lower_30 = NA_real_, band_upper_30 = NA_real_,
    rmstd_est = rd30_cphm * 16, rmstd_true = rmstd_true_30, rmstd_bias = (rd30_cphm * 16) - rmstd_true_30,
    hr_marginal_est = fit_cphm$hr, hr_marginal_true = hr_marg_true_30, hr_marginal_bias = fit_cphm$hr - hr_marg_true_30,
    z_stat_30 = z_cphm, p_val_30 = p_cphm, reject_H0 = as.integer(p_cphm < 0.05),
    uniform_band_covered = NA_integer_, prop_time_points_covered = NA_real_, mean_band_width = NA_real_
  )

  # -------------------------------------------------------------
  # 5. CCW-KM (with Simultaneous Confidence Band)
  # -------------------------------------------------------------
  t0 <- proc.time()["elapsed"]
  fit_km <- fit_ccw_km(panel_w, tau = tau, n_boot = 2000, seed = seed)
  t_km <- as.numeric(proc.time()["elapsed"] - t0)
  
  res_km <- fit_km$results
  rd15_km <- res_km$rd[15]
  rd30_km <- res_km$rd[tau]
  se30_km <- res_km$se_RD[tau]
  z_km <- rd30_km / pmax(se30_km, 1e-4)
  p_km <- 2 * (1 - pnorm(abs(z_km)))
  band_cov_km <- eval_band_cov(res_km$rd_band_lower, res_km$rd_band_upper, truth$rd_true)

  results_list[[5]] <- tibble(
    estimator = "CCW-KM",
    elapsed_sec = t_km,
    rd_15_est = rd15_km, rd_15_true = rd_true_15, bias_15 = rd15_km - rd_true_15,
    rd_30_est = rd30_km, rd_30_true = rd_true_30, bias_30 = rd30_km - rd_true_30,
    rel_bias_30 = (rd30_km - rd_true_30) / rd_true_30 * 100,
    se_30 = se30_km, ci_lower_30 = res_km$rd_ci_lower[tau], ci_upper_30 = res_km$rd_ci_upper[tau],
    band_lower_30 = res_km$rd_band_lower[tau], band_upper_30 = res_km$rd_band_upper[tau],
    rmstd_est = fit_km$rmstd_tau, rmstd_true = rmstd_true_30, rmstd_bias = fit_km$rmstd_tau - rmstd_true_30,
    hr_marginal_est = fit_km$hr_marginal_tau, hr_marginal_true = hr_marg_true_30, hr_marginal_bias = fit_km$hr_marginal_tau - hr_marg_true_30,
    z_stat_30 = z_km, p_val_30 = p_km, reject_H0 = as.integer(p_km < 0.05),
    uniform_band_covered = band_cov_km$unif, prop_time_points_covered = band_cov_km$prop, mean_band_width = band_cov_km$width
  )

  # -------------------------------------------------------------
  # 6. G-Computation (BART)
  # -------------------------------------------------------------
  t0 <- proc.time()["elapsed"]
  fit_bart <- fit_bart_gcomp(
    data = df_std,
    clones = clones,
    covars = covars_all,
    grace_period = grace_period,
    tau = tau,
    ntree = bart_trees,
    ndpost = bart_ndpost,
    nskip = bart_nskip,
    mc.cores = mc_cores,
    seed = seed
  )
  t_bart <- as.numeric(proc.time()["elapsed"] - t0)
  
  res_bart <- fit_bart$results
  rd15_bart <- res_bart$rd[15]
  rd30_bart <- res_bart$rd[tau]
  se30_bart <- if ("se_RD" %in% names(res_bart)) res_bart$se_RD[tau] else res_bart$se_rd[tau]
  z_bart <- rd30_bart / pmax(se30_bart, 1e-4)
  p_bart <- 2 * (1 - pnorm(abs(z_bart)))
  band_cov_bart <- if ("band_lower" %in% names(res_bart)) eval_band_cov(res_bart$band_lower, res_bart$band_upper, truth$rd_true) else list(unif=NA_integer_, prop=NA_real_, width=NA_real_)

  results_list[[6]] <- tibble(
    estimator = "G-COMP (BART)",
    elapsed_sec = t_bart,
    rd_15_est = rd15_bart, rd_15_true = rd_true_15, bias_15 = rd15_bart - rd_true_15,
    rd_30_est = rd30_bart, rd_30_true = rd_true_30, bias_30 = rd30_bart - rd_true_30,
    rel_bias_30 = (rd30_bart - rd_true_30) / rd_true_30 * 100,
    se_30 = se30_bart, 
    ci_lower_30 = if ("cri_lower" %in% names(res_bart)) res_bart$cri_lower[tau] else res_bart$ci_lower[tau], 
    ci_upper_30 = if ("cri_upper" %in% names(res_bart)) res_bart$cri_upper[tau] else res_bart$ci_upper[tau],
    band_lower_30 = if ("band_lower" %in% names(res_bart)) res_bart$band_lower[tau] else NA_real_, 
    band_upper_30 = if ("band_upper" %in% names(res_bart)) res_bart$band_upper[tau] else NA_real_,
    rmstd_est = sum(res_bart$rd[1:tau], na.rm = TRUE), rmstd_true = rmstd_true_30, rmstd_bias = sum(res_bart$rd[1:tau], na.rm = TRUE) - rmstd_true_30,
    hr_marginal_est = log(pmax(pmin(res_bart$treat_tx[tau], 1 - 1e-6), 1e-6)) / log(pmax(pmin(res_bart$control_no_tx[tau], 1 - 1e-6), 1e-6)),
    hr_marginal_true = hr_marg_true_30,
    hr_marginal_bias = (log(pmax(pmin(res_bart$treat_tx[tau], 1 - 1e-6), 1e-6)) / log(pmax(pmin(res_bart$control_no_tx[tau], 1 - 1e-6), 1e-6))) - hr_marg_true_30,
    z_stat_30 = z_bart, p_val_30 = p_bart, reject_H0 = as.integer(p_bart < 0.05),
    uniform_band_covered = band_cov_bart$unif, prop_time_points_covered = band_cov_bart$prop, mean_band_width = band_cov_bart$width
  )

  # -------------------------------------------------------------
  # 7. TMLE + BART
  # -------------------------------------------------------------
  t0 <- proc.time()["elapsed"]
  fit_tmle_bart <- fit_bart_tmle(
    data = df_std,
    clones = clones,
    covars = covars_all,
    grace_period = grace_period,
    tau = tau,
    ntree = bart_trees,
    ndpost = bart_ndpost,
    nskip = bart_nskip,
    mc.cores = mc_cores,
    n_boot = 2000,
    seed = seed
  )
  t_tmle_bart <- as.numeric(proc.time()["elapsed"] - t0)
  
  res_tbart <- fit_tmle_bart$results
  rd15_tbart <- res_tbart$rd[15]
  rd30_tbart <- res_tbart$rd[tau]
  se30_tbart <- res_tbart$se_RD[tau]
  z_tbart <- rd30_tbart / pmax(se30_tbart, 1e-4)
  p_tbart <- 2 * (1 - pnorm(abs(z_tbart)))
  band_cov_tbart <- eval_band_cov(res_tbart$rd_band_lower, res_tbart$rd_band_upper, truth$rd_true)

  results_list[[7]] <- tibble(
    estimator = "TMLE + BART",
    elapsed_sec = t_tmle_bart,
    rd_15_est = rd15_tbart, rd_15_true = rd_true_15, bias_15 = rd15_tbart - rd_true_15,
    rd_30_est = rd30_tbart, rd_30_true = rd_true_30, bias_30 = rd30_tbart - rd_true_30,
    rel_bias_30 = (rd30_tbart - rd_true_30) / rd_true_30 * 100,
    se_30 = se30_tbart, ci_lower_30 = res_tbart$rd_ci_lower[tau], ci_upper_30 = res_tbart$rd_ci_upper[tau],
    band_lower_30 = res_tbart$rd_band_lower[tau], band_upper_30 = res_tbart$rd_band_upper[tau],
    rmstd_est = sum(res_tbart$rd[1:tau], na.rm = TRUE), rmstd_true = rmstd_true_30, rmstd_bias = sum(res_tbart$rd[1:tau], na.rm = TRUE) - rmstd_true_30,
    hr_marginal_est = res_tbart$hr_cum[tau], hr_marginal_true = hr_marg_true_30, hr_marginal_bias = res_tbart$hr_cum[tau] - hr_marg_true_30,
    z_stat_30 = z_tbart, p_val_30 = p_tbart, reject_H0 = as.integer(p_tbart < 0.05),
    uniform_band_covered = band_cov_tbart$unif, prop_time_points_covered = band_cov_tbart$prop, mean_band_width = band_cov_tbart$width
  )

  # -------------------------------------------------------------
  # 8. TTE-TMLE (Parametric GLM)
  # -------------------------------------------------------------
  t0 <- proc.time()["elapsed"]
  fit_glm_tmle <- tte_tmle(
    data = df_std,
    grace_period = grace_period,
    tau = tau,
    covars = covars_all,
    p_form = p_form_use,
    q_form = q_form_use,
    use_sl = FALSE,
    n_boot = 2000,
    seed = seed
  )
  t_glm_tmle <- as.numeric(proc.time()["elapsed"] - t0)
  
  res_glm <- fit_glm_tmle$res_tmle
  rd15_glm <- res_glm$rd[15]
  rd30_glm <- res_glm$rd[tau]
  se30_glm <- res_glm$se_RD[tau]
  z_glm <- rd30_glm / pmax(se30_glm, 1e-4)
  p_glm <- 2 * (1 - pnorm(abs(z_glm)))
  band_cov_glm <- eval_band_cov(res_glm$rd_band_lower, res_glm$rd_band_upper, truth$rd_true)

  results_list[[8]] <- tibble(
    estimator = "TTE-TMLE (GLM)",
    elapsed_sec = t_glm_tmle,
    rd_15_est = rd15_glm, rd_15_true = rd_true_15, bias_15 = rd15_glm - rd_true_15,
    rd_30_est = rd30_glm, rd_30_true = rd_true_30, bias_30 = rd30_glm - rd_true_30,
    rel_bias_30 = (rd30_glm - rd_true_30) / rd_true_30 * 100,
    se_30 = se30_glm, ci_lower_30 = res_glm$rd_ci_lower[tau], ci_upper_30 = res_glm$rd_ci_upper[tau],
    band_lower_30 = res_glm$rd_band_lower[tau], band_upper_30 = res_glm$rd_band_upper[tau],
    rmstd_est = sum(res_glm$rd[1:tau], na.rm = TRUE), rmstd_true = rmstd_true_30, rmstd_bias = sum(res_glm$rd[1:tau], na.rm = TRUE) - rmstd_true_30,
    hr_marginal_est = res_glm$hr_cum[tau], hr_marginal_true = hr_marg_true_30, hr_marginal_bias = res_glm$hr_cum[tau] - hr_marg_true_30,
    z_stat_30 = z_glm, p_val_30 = p_glm, reject_H0 = as.integer(p_glm < 0.05),
    uniform_band_covered = band_cov_glm$unif, prop_time_points_covered = band_cov_glm$prop, mean_band_width = band_cov_glm$width
  )

  # -------------------------------------------------------------
  # 9. TTE-TMLE (SuperLearner Ensemble)
  # -------------------------------------------------------------
  t0 <- proc.time()["elapsed"]
  fit_sl_tmle <- tte_tmle(
    data = df_std,
    grace_period = grace_period,
    tau = tau,
    covars = covars_all,
    use_sl = TRUE,
    sl_lib = c("SL.mean", "SL.glm", "SL.glm.time.int", "SL.earth.int"),
    n_boot = 2000,
    seed = seed
  )
  t_sl_tmle <- as.numeric(proc.time()["elapsed"] - t0)
  
  res_sl <- fit_sl_tmle$res_tmle
  rd15_sl <- res_sl$rd[15]
  rd30_sl <- res_sl$rd[tau]
  se30_sl <- res_sl$se_RD[tau]
  z_sl <- rd30_sl / pmax(se30_sl, 1e-4)
  p_sl <- 2 * (1 - pnorm(abs(z_sl)))
  band_cov_sl <- eval_band_cov(res_sl$rd_band_lower, res_sl$rd_band_upper, truth$rd_true)

  results_list[[9]] <- tibble(
    estimator = "TTE-TMLE (SuperLearner)",
    elapsed_sec = t_sl_tmle,
    rd_15_est = rd15_sl, rd_15_true = rd_true_15, bias_15 = rd15_sl - rd_true_15,
    rd_30_est = rd30_sl, rd_30_true = rd_true_30, bias_30 = rd30_sl - rd_true_30,
    rel_bias_30 = (rd30_sl - rd_true_30) / rd_true_30 * 100,
    se_30 = se30_sl, ci_lower_30 = res_sl$rd_ci_lower[tau], ci_upper_30 = res_sl$rd_ci_upper[tau],
    band_lower_30 = res_sl$rd_band_lower[tau], band_upper_30 = res_sl$rd_band_upper[tau],
    rmstd_est = sum(res_sl$rd[1:tau], na.rm = TRUE), rmstd_true = rmstd_true_30, rmstd_bias = sum(res_sl$rd[1:tau], na.rm = TRUE) - rmstd_true_30,
    hr_marginal_est = res_sl$hr_cum[tau], hr_marginal_true = hr_marg_true_30, hr_marginal_bias = res_sl$hr_cum[tau] - hr_marg_true_30,
    z_stat_30 = z_sl, p_val_30 = p_sl, reject_H0 = as.integer(p_sl < 0.05),
    uniform_band_covered = band_cov_sl$unif, prop_time_points_covered = band_cov_sl$prop, mean_band_width = band_cov_sl$width
  )

  # Combine into single condition tibble
  out_df <- dplyr::bind_rows(results_list) %>%
    dplyr::mutate(
      scenario_id = paste0("Scenario ", scenario_num),
      regime_id = paste0("Regime ", regime_num, " (", regime_name, ")"),
      n = n,
      .before = 1
    )

  return(out_df)
}

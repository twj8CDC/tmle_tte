suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(parallel)
  library(readr)
})

source("sim/run_simulation_condition.R")

t_start <- proc.time()["elapsed"]

out_rds <- "sim/scratch/check_4estimators_s1_s2_s3_20iter.rds"
out_csv <- "sim/scratch/check_4estimators_s1_s2_s3_20iter.csv"

cat("Precomputing Oracle Ground Truth for Scenarios 1, 2, 3 (N_oracle = 200,000)...\n")
truth_map <- list(
  "1" = compute_scenario1_truth(),
  "2" = compute_scenario2_truth(),
  "3" = compute_scenario3_truth()
)

base_seed <- 2026L
tau <- 30L
grace_period <- 10L
n_obs <- 1000L

run_one_task <- function(task) {
  s <- task$s
  it <- task$it
  job_seed <- base_seed + (s * 100000L) + it * 37L
  set.seed(job_seed)

  truth_s <- truth_map[[as.character(s)]]
  spec <- get_scenario_spec(
    scenario_num = s,
    n = n_obs,
    tau = tau,
    grace_period = grace_period,
    seed = job_seed,
    truth_precomp = truth_s
  )
  df_std <- prepare_tte_data(spec$df_raw)
  covars_all <- spec$covars_all
  clones <- create_clones(df_std, covars = covars_all, grace_period = grace_period, tau = tau)

  rd_true_15 <- truth_s$rd_true[15]
  rd_true_30 <- truth_s$rd_true[tau]
  rmstd_true_30 <- sum(truth_s$rd_true[1:tau], na.rm = TRUE)
  hr_marg_true_30 <- if ("hr_true" %in% names(truth_s) && length(truth_s$hr_true) >= tau) {
    truth_s$hr_true[tau]
  } else {
    log(pmax(pmin(truth_s$S1_true[tau], 1 - 1e-6), 1e-6)) /
      log(pmax(pmin(truth_s$S0_true[tau], 1 - 1e-6), 1e-6))
  }

  ev_trt_days <- sort(ceiling(df_std$time[df_std$tx_ever == 1L & df_std$t_switch <= grace_period & df_std$event == 1L & df_std$time >= df_std$t_switch]))
  t_band_start <- if (length(ev_trt_days) >= 5L) max(1L, as.integer(ev_trt_days[5L])) else as.integer(grace_period)
  idx_band <- seq(from = min(t_band_start, tau), to = tau)

  make_row <- function(est_label, res_tbl, rmstd_est, se_rmstd, elapsed_sec, hr_est = res_tbl$hr_cum[tau],
                       hr_l = res_tbl$hr_ci_lower[tau], hr_u = res_tbl$hr_ci_upper[tau]) {
    is_cov <- (res_tbl$rd_band_lower[idx_band] <= truth_s$rd_true[idx_band]) &
      (truth_s$rd_true[idx_band] <= res_tbl$rd_band_upper[idx_band])
    rmstd_l <- rmstd_est - 1.96 * se_rmstd
    rmstd_u <- rmstd_est + 1.96 * se_rmstd

    tibble(
      scenario_id = paste0("Scenario ", s),
      scenario_num = s,
      iter = it,
      estimator = est_label,
      elapsed_sec = elapsed_sec,
      rd_15_est = res_tbl$rd[15],
      bias_15 = res_tbl$rd[15] - rd_true_15,
      se_15 = res_tbl$se_RD[15],
      covered_15 = as.integer(res_tbl$rd_ci_lower[15] <= rd_true_15 & rd_true_15 <= res_tbl$rd_ci_upper[15]),
      ci_width_15 = res_tbl$rd_ci_upper[15] - res_tbl$rd_ci_lower[15],
      rd_30_est = res_tbl$rd[tau],
      bias_30 = res_tbl$rd[tau] - rd_true_30,
      se_30 = res_tbl$se_RD[tau],
      covered_30 = as.integer(res_tbl$rd_ci_lower[tau] <= rd_true_30 & rd_true_30 <= res_tbl$rd_ci_upper[tau]),
      ci_width_30 = res_tbl$rd_ci_upper[tau] - res_tbl$rd_ci_lower[tau],
      rmstd_est = rmstd_est,
      rmstd_bias = rmstd_est - rmstd_true_30,
      se_rmstd = se_rmstd,
      rmstd_covered = as.integer(rmstd_l <= rmstd_true_30 & rmstd_true_30 <= rmstd_u),
      rmstd_ci_width = rmstd_u - rmstd_l,
      hr_marginal_est = hr_est,
      hr_marginal_bias = hr_est - hr_marg_true_30,
      hr_covered_30 = as.integer(hr_l <= hr_marg_true_30 & hr_marg_true_30 <= hr_u),
      uniform_band_covered = as.integer(all(is_cov, na.rm = TRUE)),
      prop_time_points_covered = mean(is_cov, na.rm = TRUE),
      mean_band_width = mean(res_tbl$rd_band_upper[idx_band] - res_tbl$rd_band_lower[idx_band], na.rm = TRUE)
    )
  }

  # 1. TTE-TMLE (BART) using package defaults (SurvBART.preset_standard(): trees=50, ndpost=250, nskip=250, thin=4, chains=4)
  fit_bart <- fit_bart_tmle(
    data = df_std,
    clones = clones,
    covars = covars_all,
    grace_period = grace_period,
    tau = tau,
    n_boot = 2000L,
    seed = job_seed
  )
  row_bart <- make_row(
    est_label = "TTE-TMLE (BART)",
    res_tbl = fit_bart$results,
    rmstd_est = sum(fit_bart$results$rd[1:tau], na.rm = TRUE),
    se_rmstd = fit_bart$se_rmstd,
    elapsed_sec = fit_bart$elapsed_sec
  )

  # 2. TTE-TMLE (SuperLearner)
  t0_sl <- proc.time()["elapsed"]
  fit_sl <- tte_tmle(
    data = df_std,
    grace_period = grace_period,
    tau = tau,
    covars = covars_all,
    use_sl = TRUE,
    sl_lib = c("SL.mean", "SL.glm.time.int", "SL.hal.time", "SL.earth.int"),
    trim_quantile = c(0.001, 0.999),
    n_boot = 2000L,
    seed = job_seed
  )
  t_sl <- as.numeric(proc.time()["elapsed"] - t0_sl)
  row_sl <- make_row(
    est_label = "TTE-TMLE (SuperLearner)",
    res_tbl = fit_sl$res_tmle,
    rmstd_est = sum(fit_sl$res_tmle$rd[1:tau], na.rm = TRUE),
    se_rmstd = fit_sl$se_rmstd,
    elapsed_sec = t_sl
  )

  # 3. TTE-TMLE (GLM) [Corr Q + Corr g]
  t0_glm <- proc.time()["elapsed"]
  fit_glm <- tte_tmle(
    data = df_std,
    grace_period = grace_period,
    tau = tau,
    covars = covars_all,
    p_form = spec$p_form_corr,
    q_form = spec$q_form_corr,
    use_sl = FALSE,
    n_boot = 2000L,
    seed = job_seed
  )
  t_glm <- as.numeric(proc.time()["elapsed"] - t0_glm)
  row_glm <- make_row(
    est_label = "TTE-TMLE (GLM) [Corr]",
    res_tbl = fit_glm$res_tmle,
    rmstd_est = sum(fit_glm$res_tmle$rd[1:tau], na.rm = TRUE),
    se_rmstd = fit_glm$se_rmstd,
    elapsed_sec = t_glm
  )

  # 4. CCW-KM [Corr]
  t0_km <- proc.time()["elapsed"]
  ipcw_res <- compute_ipcw_weights(
    clones = clones,
    covars = covars_all,
    grace_period = grace_period,
    dt = 1,
    p_form = spec$p_form_corr,
    use_sl = FALSE
  )
  fit_km <- fit_ccw_km(ipcw_res$panel_w, tau = tau, n_boot = 2000L, seed = job_seed)
  t_km <- as.numeric(proc.time()["elapsed"] - t0_km)
  row_km <- make_row(
    est_label = "CCW-KM [Corr]",
    res_tbl = fit_km$results,
    rmstd_est = fit_km$rmstd_tau,
    se_rmstd = fit_km$se_rmstd,
    elapsed_sec = t_km,
    hr_est = fit_km$hr_marginal_tau,
    hr_l = fit_km$hr_ci_lower_tau,
    hr_u = fit_km$hr_ci_upper_tau
  )

  bind_rows(row_bart, row_sl, row_glm, row_km)
}

tasks <- unlist(lapply(1:3, function(s) lapply(1:20, function(it) list(s = s, it = it))), recursive = FALSE)
cat(sprintf(
  "Running %d tasks (Scenarios 1-3 x 20 iters) across 10 cores (BART defaults: ntree=%d, ndpost=%d, nskip=%d, thin=%d, chains=%d)...\n",
  length(tasks),
  formals(fit_bart_pair)$ntree,
  formals(fit_bart_pair)$ndpost,
  formals(fit_bart_pair)$nskip,
  formals(fit_bart_pair)$thin,
  eval(formals(fit_bart_pair)$chains)
))

res_list <- mclapply(tasks, function(tk) {
  tryCatch(run_one_task(tk), error = function(e) {
    warning(sprintf("Scenario %d Iter %d failed: %s", tk$s, tk$it, e$message))
    NULL
  })
}, mc.cores = 10L, mc.preschedule = FALSE)

n_null <- sum(vapply(res_list, is.null, logical(1)))
if (n_null > 0L) stop(sprintf("%d tasks failed!", n_null))

new_df <- bind_rows(res_list)
saveRDS(new_df, out_rds)
write_csv(new_df, out_csv)

summ <- new_df %>%
  group_by(scenario_num, scenario_id, estimator) %>%
  summarise(
    reps = n(),
    bias15 = mean(bias_15),
    rmse15 = sqrt(mean(bias_15^2)),
    cov15 = mean(covered_15) * 100,
    se_r15 = mean(se_15) / sd(rd_15_est),
    bias30 = mean(bias_30),
    rmse30 = sqrt(mean(bias_30^2)),
    cov30 = mean(covered_30) * 100,
    se_r30 = mean(se_30) / sd(rd_30_est),
    band_cov = mean(uniform_band_covered) * 100,
    rmstd_b = mean(rmstd_bias),
    rmstd_c = mean(rmstd_covered) * 100,
    hr_b = mean(hr_marginal_bias),
    hr_c = mean(hr_covered_30) * 100,
    time_s = mean(elapsed_sec),
    .groups = "drop"
  ) %>%
  arrange(scenario_num, estimator)

cat("\n=== Phase 2 + BART Preset Standard: Scenarios 1-3 (20 Iters Each) ===\n")
print(as.data.frame(summ), digits = 4)

if (file.exists("sim/simulation_results_100iter.rds")) {
  orig <- readRDS("sim/simulation_results_100iter.rds") %>%
    filter(
      scenario_id %in% c("Scenario 1", "Scenario 2", "Scenario 3"),
      iter %in% 1:20,
      (estimator %in% c("TTE-TMLE (BART)", "TTE-TMLE (SuperLearner)")) |
        (estimator == "CCW-KM" & spec_group == "Parametric (Correct)") |
        (estimator == "TTE-TMLE (GLM)" & spec_group == "Parametric (Correct)")
    ) %>%
    mutate(
      estimator = case_when(
        estimator == "CCW-KM" ~ "CCW-KM [Corr]",
        estimator == "TTE-TMLE (GLM)" ~ "TTE-TMLE (GLM) [Corr]",
        TRUE ~ estimator
      )
    ) %>%
    group_by(scenario_id, estimator) %>%
    summarise(
      p1_bias15 = mean(bias_15),
      p1_cov15 = mean(covered_15) * 100,
      p1_bias30 = mean(bias_30),
      p1_rmse30 = sqrt(mean(bias_30^2)),
      p1_cov30 = mean(covered_30) * 100,
      p1_band_cov = mean(uniform_band_covered) * 100,
      .groups = "drop"
    )
  comp <- summ %>%
    select(scenario_id, estimator, bias15, cov15, bias30, rmse30, cov30, band_cov) %>%
    left_join(orig, by = c("scenario_id", "estimator"))
  cat("\n=== Paired Delta vs Phase 1 (Exact Same 20 Datasets: Iter 1..20) ===\n")
  print(as.data.frame(comp), digits = 4)
}

cat(sprintf("\nCOMPLETED IN %.2f MINS!\n", (proc.time()["elapsed"] - t_start) / 60))

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(parallel)
  library(readr)
})

source("sim/run_simulation_condition.R")

t_start <- proc.time()["elapsed"]

out_rds <- "sim/scratch/bart_preset_standard_s3_20iter.rds"
out_csv <- "sim/scratch/bart_preset_standard_s3_20iter.csv"

cat("Precomputing Oracle Ground Truth for Scenario 3 (N_oracle = 200,000)...\n")
truth_s3 <- compute_scenario3_truth()

base_seed <- 2026L
s <- 3L
tau <- 30L
grace_period <- 10L
n_obs <- 1000L

# Exact SurvBART.preset_standard() configuration:
# trees = 50, draws = 1000 (ndpost = 250 * thin = 4), tune = 250, chains = 4
std_ntree <- 50L
std_ndpost <- 250L
std_nskip <- 250L
std_thin <- 4L
std_chains <- 4L
std_cores <- 1L

run_s3_iter <- function(it) {
  job_seed <- base_seed + (s * 100000L) + it * 37L
  set.seed(job_seed)

  spec <- get_scenario_spec(
    scenario_num = s,
    n = n_obs,
    tau = tau,
    grace_period = grace_period,
    seed = job_seed,
    truth_precomp = truth_s3
  )
  df_std <- prepare_tte_data(spec$df_raw)
  covars_all <- spec$covars_all
  clones <- create_clones(df_std, covars = covars_all, grace_period = grace_period, tau = tau)

  rd_true_15 <- truth_s3$rd_true[15]
  rd_true_30 <- truth_s3$rd_true[tau]
  rmstd_true_30 <- sum(truth_s3$rd_true[1:tau], na.rm = TRUE)
  hr_marg_true_30 <- if ("hr_true" %in% names(truth_s3) && length(truth_s3$hr_true) >= tau) {
    truth_s3$hr_true[tau]
  } else {
    log(pmax(pmin(truth_s3$S1_true[tau], 1 - 1e-6), 1e-6)) /
      log(pmax(pmin(truth_s3$S0_true[tau], 1 - 1e-6), 1e-6))
  }

  ev_trt_days <- sort(ceiling(df_std$time[df_std$tx_ever == 1L & df_std$t_switch <= grace_period & df_std$event == 1L & df_std$time >= df_std$t_switch]))
  t_band_start <- if (length(ev_trt_days) >= 5L) max(1L, as.integer(ev_trt_days[5L])) else as.integer(grace_period)
  idx_band <- seq(from = min(t_band_start, tau), to = tau)

  bart_pair <- fit_bart_pair(
    data = df_std,
    clones = clones,
    covars = covars_all,
    grace_period = grace_period,
    tau = tau,
    ntree = std_ntree,
    ndpost = std_ndpost,
    nskip = std_nskip,
    thin = std_thin,
    chains = std_chains,
    cores = std_cores,
    n_boot = 2000L,
    seed = job_seed
  )

  # 1. G-COMP (BART-StdPreset)
  fit_bg <- bart_pair$gcomp
  res_bg <- fit_bg$results
  is_cov_bg <- (res_bg$band_lower[idx_band] <= truth_s3$rd_true[idx_band]) &
    (truth_s3$rd_true[idx_band] <= res_bg$band_upper[idx_band])

  row_bg <- tibble(
    scenario_id = "Scenario 3",
    iter = it,
    estimator = "G-COMP (BART-StdPreset)",
    elapsed_sec = fit_bg$elapsed_sec,
    rd_15_est = res_bg$rd[15],
    bias_15 = res_bg$rd[15] - rd_true_15,
    se_15 = res_bg$se_RD[15],
    covered_15 = as.integer(res_bg$cri_lower[15] <= rd_true_15 & rd_true_15 <= res_bg$cri_upper[15]),
    ci_width_15 = res_bg$cri_upper[15] - res_bg$cri_lower[15],
    rd_30_est = res_bg$rd[tau],
    bias_30 = res_bg$rd[tau] - rd_true_30,
    se_30 = res_bg$se_RD[tau],
    covered_30 = as.integer(res_bg$cri_lower[tau] <= rd_true_30 & rd_true_30 <= res_bg$cri_upper[tau]),
    ci_width_30 = res_bg$cri_upper[tau] - res_bg$cri_lower[tau],
    rmstd_est = fit_bg$rmstd,
    rmstd_bias = fit_bg$rmstd - rmstd_true_30,
    se_rmstd = fit_bg$se_rmstd,
    rmstd_covered = as.integer(fit_bg$rmstd_ci_lower <= rmstd_true_30 & rmstd_true_30 <= fit_bg$rmstd_ci_upper),
    rmstd_ci_width = fit_bg$rmstd_ci_upper - fit_bg$rmstd_ci_lower,
    hr_marginal_est = res_bg$hr_cum[tau],
    hr_marginal_bias = res_bg$hr_cum[tau] - hr_marg_true_30,
    hr_covered_30 = as.integer(res_bg$hr_ci_lower[tau] <= hr_marg_true_30 & hr_marg_true_30 <= res_bg$hr_ci_upper[tau]),
    uniform_band_covered = as.integer(all(is_cov_bg, na.rm = TRUE)),
    prop_time_points_covered = mean(is_cov_bg, na.rm = TRUE),
    mean_band_width = mean(res_bg$band_upper[idx_band] - res_bg$band_lower[idx_band], na.rm = TRUE)
  )

  # 2. TTE-TMLE (BART-StdPreset)
  fit_tb <- bart_pair$tmle
  res_tb <- fit_tb$results
  is_cov_tb <- (res_tb$rd_band_lower[idx_band] <= truth_s3$rd_true[idx_band]) &
    (truth_s3$rd_true[idx_band] <= res_tb$rd_band_upper[idx_band])
  rmstd_tb <- sum(res_tb$rd[1:tau], na.rm = TRUE)
  rmstd_l_tb <- rmstd_tb - 1.96 * fit_tb$se_rmstd
  rmstd_u_tb <- rmstd_tb + 1.96 * fit_tb$se_rmstd

  row_tb <- tibble(
    scenario_id = "Scenario 3",
    iter = it,
    estimator = "TTE-TMLE (BART-StdPreset)",
    elapsed_sec = fit_tb$elapsed_sec,
    rd_15_est = res_tb$rd[15],
    bias_15 = res_tb$rd[15] - rd_true_15,
    se_15 = res_tb$se_RD[15],
    covered_15 = as.integer(res_tb$rd_ci_lower[15] <= rd_true_15 & rd_true_15 <= res_tb$rd_ci_upper[15]),
    ci_width_15 = res_tb$rd_ci_upper[15] - res_tb$rd_ci_lower[15],
    rd_30_est = res_tb$rd[tau],
    bias_30 = res_tb$rd[tau] - rd_true_30,
    se_30 = res_tb$se_RD[tau],
    covered_30 = as.integer(res_tb$rd_ci_lower[tau] <= rd_true_30 & rd_true_30 <= res_tb$rd_ci_upper[tau]),
    ci_width_30 = res_tb$rd_ci_upper[tau] - res_tb$rd_ci_lower[tau],
    rmstd_est = rmstd_tb,
    rmstd_bias = rmstd_tb - rmstd_true_30,
    se_rmstd = fit_tb$se_rmstd,
    rmstd_covered = as.integer(rmstd_l_tb <= rmstd_true_30 & rmstd_true_30 <= rmstd_u_tb),
    rmstd_ci_width = rmstd_u_tb - rmstd_l_tb,
    hr_marginal_est = res_tb$hr_cum[tau],
    hr_marginal_bias = res_tb$hr_cum[tau] - hr_marg_true_30,
    hr_covered_30 = as.integer(res_tb$hr_ci_lower[tau] <= hr_marg_true_30 & hr_marg_true_30 <= res_tb$hr_ci_upper[tau]),
    uniform_band_covered = as.integer(all(is_cov_tb, na.rm = TRUE)),
    prop_time_points_covered = mean(is_cov_tb, na.rm = TRUE),
    mean_band_width = mean(res_tb$rd_band_upper[idx_band] - res_tb$rd_band_lower[idx_band], na.rm = TRUE)
  )

  bind_rows(row_bg, row_tb)
}

cat(sprintf(
  "Running Scenario 3 x 20 iterations with SurvBART.preset_standard (trees=%d, draws=%d, tune=%d, thin=%d, chains=%d) across 10 cores...\n",
  std_ntree, std_ndpost * std_thin, std_nskip, std_thin, std_chains
))

res_list <- mclapply(1:20, function(it) {
  tryCatch(run_s3_iter(it), error = function(e) {
    warning(sprintf("Iter %d failed: %s", it, e$message))
    NULL
  })
}, mc.cores = 10L, mc.preschedule = FALSE)

n_null <- sum(vapply(res_list, is.null, logical(1)))
if (n_null > 0L) stop(sprintf("%d jobs failed!", n_null))

new_df <- bind_rows(res_list)
saveRDS(new_df, out_rds)
write_csv(new_df, out_csv)

# Paired comparison on the EXACT same 20 datasets (Scenario 3, iter 1..20)
orig <- readRDS("sim/simulation_results_100iter.rds") %>%
  filter(
    scenario_id == "Scenario 3",
    iter %in% 1:20,
    (estimator %in% c("G-COMP (BART)", "TTE-TMLE (BART)", "TTE-TMLE (SuperLearner)")) |
      (estimator == "CCW-KM" & spec_group == "Parametric (Correct)") |
      (estimator == "TTE-TMLE (GLM)" & spec_group == "Parametric (Correct)")
  ) %>%
  mutate(
    estimator = case_when(
      estimator == "G-COMP (BART)" ~ "G-COMP (BART-Fast)",
      estimator == "TTE-TMLE (BART)" ~ "TTE-TMLE (BART-Fast)",
      estimator == "CCW-KM" ~ "CCW-KM [Corr]",
      estimator == "TTE-TMLE (GLM)" ~ "TTE-TMLE (GLM) [Corr]",
      TRUE ~ estimator
    )
  )

comp <- bind_rows(new_df, orig %>% select(any_of(names(new_df)))) %>%
  group_by(estimator) %>%
  summarise(
    reps = n(),
    bias15 = mean(bias_15),
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
  arrange(bias15)

cat("\n=== Scenario 3 Paired Comparison (Exact Same 20 Datasets: Iter 1..20) ===\n")
print(as.data.frame(comp), digits = 4)
cat(sprintf("\nCOMPLETED IN %.2f MINS!\n", (proc.time()["elapsed"] - t_start) / 60))

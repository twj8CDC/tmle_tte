suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(parallel)
  library(readr)
})

source("sim/run_simulation_condition.R")

t_start <- proc.time()["elapsed"]

out_rds <- "sim/simulation_results_100iter.rds"
out_csv <- "sim/simulation_results_100iter.csv"
bak_rds <- "sim/scratch/simulation_results_100iter_pre_option1_bak.rds"

if (!file.exists(bak_rds)) {
  cat("Backing up pre-Option-1 master RDS to:", bak_rds, "\n")
  existing <- readRDS(out_rds)
  saveRDS(existing, bak_rds)
} else {
  cat("Using existing pre-Option-1 backup:", bak_rds, "\n")
  existing <- readRDS(bak_rds)
}

keep_ests <- c(
  "TTE-TMLE (BART)",
  "G-COMP (BART)",
  "TTE-TMLE (SuperLearner)",
  "TTE-TMLE (GLM)",
  "CCW-KM"
)
keep_existing <- existing %>%
  filter(estimator %in% keep_ests) %>%
  mutate(
    hr_cond_est = NA_real_,
    hr_cond_ci_lower = NA_real_,
    hr_cond_ci_upper = NA_real_,
    hr_cond_ci_width = NA_real_
  )

cat(sprintf(
  "Preserved %d rows (%d BART/SL + %d GLM + %d CCW-KM).\n",
  nrow(keep_existing),
  sum(keep_existing$spec_group == "Non-Parametric"),
  sum(keep_existing$estimator == "TTE-TMLE (GLM)"),
  sum(keep_existing$estimator == "CCW-KM")
))
stopifnot(nrow(keep_existing) == 7200L)

cat("Precomputing oracle ground truth for Scenarios 1-8...\n")
truth_list <- list(
  `1` = compute_scenario1_truth(),
  `2` = compute_scenario2_truth(),
  `3` = compute_scenario3_truth(),
  `4` = compute_scenario4_truth(),
  `5` = compute_scenario5_truth(),
  `6` = compute_scenario6_truth(),
  `7` = compute_scenario7_truth(),
  `8` = compute_scenario8_truth()
)

base_seed <- 2026L
grid <- expand.grid(iter = 1:100, scenario = 1:8, stringsAsFactors = FALSE) %>%
  arrange(scenario, iter)

cat(sprintf("Running %d cox_only dataset jobs across 10 cores...\n", nrow(grid)))

run_job <- function(idx) {
  s <- grid$scenario[idx]
  it <- grid$iter[idx]
  job_seed <- base_seed + (s * 100000L) + it * 37L
  tryCatch({
    df_out <- run_simulation_dataset(
      scenario_num = s,
      n = 1000,
      seed = job_seed,
      cox_only = TRUE,
      truth_precomp = truth_list[[as.character(s)]]
    )
    df_out$iter <- it
    df_out
  }, error = function(e) {
    warning(sprintf("Error Scen %d Iter %d: %s", s, it, e$message))
    NULL
  })
}

res_list <- mclapply(seq_len(nrow(grid)), run_job, mc.cores = 10L, mc.preschedule = FALSE)
n_null <- sum(vapply(res_list, is.null, logical(1)))
if (n_null > 0L) stop(sprintf("Found %d failed jobs!", n_null))

new_cox <- bind_rows(res_list)
cat(sprintf("Completed %d Cox rows in %.2f mins!\n", nrow(new_cox), (proc.time()["elapsed"] - t_start) / 60))
stopifnot(nrow(new_cox) == 4800L)

col_order <- names(new_cox)
keep_aligned <- keep_existing %>% select(all_of(col_order))

merged_df <- bind_rows(keep_aligned, new_cox) %>%
  arrange(scenario_id, iter, spec_group, estimator, regime_id)

stopifnot(nrow(merged_df) == 12000L)
print(table(merged_df$scenario_id, merged_df$estimator))

saveRDS(merged_df, out_rds)
write_csv(merged_df, out_csv)
cat("Saved Option-1 RDS and CSV (12,000 rows).\n")

cat("Regenerating Master Excel: sim/simulation_results_100iter.xlsx...\n")
system2(Sys.which("python3"), args = "sim/export_100iter_results_xlsx.py")

cat("Regenerating 8 publication figures...\n")
source("sim/generate_100iter_plots.R")
generate_100iter_plots(rds_path = out_rds, out_dir = "sim/figures")

# Print verification summary for Time-Matched Cox, Landmark Cox, and Time-Dependent Cox
tm_sum <- merged_df %>%
  filter(estimator == "Time-Matched Cox") %>%
  group_by(scenario_id, spec_group) %>%
  summarise(
    bias15 = mean(bias_15),
    bias30 = mean(bias_30),
    rmse30 = sqrt(mean(bias_30^2)),
    cov30 = mean(covered_30) * 100,
    se_ratio_30 = mean(se_30) / sd(rd_30_est),
    rmstd_bias = mean(rmstd_bias),
    rmstd_cov = mean(rmstd_covered) * 100,
    se_ratio_rmstd = mean(se_rmstd) / sd(rmstd_est),
    hr_bias = mean(hr_marginal_bias),
    hr_cov = mean(hr_covered_30) * 100,
    cond_hr_mean = mean(hr_cond_est),
    .groups = "drop"
  )
cat("\n=== Time-Matched Cox (Bootstrapped Option 1 Summary) ===\n")
print(as.data.frame(tm_sum), digits = 4)

cond_sum <- merged_df %>%
  filter(estimator %in% c("Landmark Cox", "Time-Dependent Cox", "Time-Matched Cox")) %>%
  group_by(scenario_id, estimator, spec_group) %>%
  summarise(
    cond_hr = mean(hr_cond_est, na.rm = TRUE),
    cond_hr_sd = sd(hr_cond_est, na.rm = TRUE),
    cond_hr_w = mean(hr_cond_ci_width, na.rm = TRUE),
    n_non_na_rd30 = sum(!is.na(rd_30_est)),
    .groups = "drop"
  )
cat("\n=== Conditional Cox HRs & Marginal NA Guard Verification ===\n")
print(as.data.frame(cond_sum), digits = 4)

cat(sprintf("\nALL OPTION 1 BENCHMARK UPDATES COMPLETED IN %.2f MINS!\n", (proc.time()["elapsed"] - t_start) / 60))

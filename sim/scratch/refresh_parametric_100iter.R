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
bak_rds <- "sim/scratch/simulation_results_100iter_pre_phase1_bak.rds"

if (!file.exists(bak_rds)) {
  cat("Backing up existing master RDS to:", bak_rds, "\n")
  existing <- readRDS(out_rds)
  saveRDS(existing, bak_rds)
} else {
  cat("Using existing pre-Phase-1 backup:", bak_rds, "\n")
  existing <- readRDS(bak_rds)
}

old_nonpar <- existing %>% filter(spec_group == "Non-Parametric")
cat(sprintf("Preserved %d Non-Parametric rows (G-COMP BART, TTE-TMLE BART, TTE-TMLE SuperLearner).\n", nrow(old_nonpar)))

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

cat(sprintf("Running %d fast_only parametric dataset jobs across 10 cores...\n", nrow(grid)))

run_job <- function(idx) {
  s <- grid$scenario[idx]
  it <- grid$iter[idx]
  job_seed <- base_seed + (s * 100000L) + it * 37L
  res <- tryCatch({
    df_out <- run_simulation_dataset(
      scenario_num = s,
      n = 1000,
      seed = job_seed,
      fast_only = TRUE,
      truth_precomp = truth_list[[as.character(s)]]
    )
    df_out$iter <- it
    df_out
  }, error = function(e) {
    warning(sprintf("Error Scen %d Iter %d: %s", s, it, e$message))
    NULL
  })
  res
}

res_list <- mclapply(seq_len(nrow(grid)), run_job, mc.cores = 10L, mc.preschedule = FALSE)
n_null <- sum(vapply(res_list, is.null, logical(1)))
if (n_null > 0L) stop(sprintf("Found %d failed jobs!", n_null))

new_param <- bind_rows(res_list)
cat(sprintf("Completed %d parametric rows in %.2f mins!\n", nrow(new_param), (proc.time()["elapsed"] - t_start) / 60))

stopifnot(nrow(new_param) == 12800L)

# Match column ordering and merge cleanly
col_order <- names(existing)
new_param <- new_param %>% select(all_of(col_order))

merged_df <- bind_rows(old_nonpar, new_param) %>%
  arrange(scenario_id, iter, spec_group, estimator, regime_id)

stopifnot(nrow(merged_df) == 15200L)
print(table(merged_df$scenario_id, merged_df$estimator))

saveRDS(merged_df, out_rds)
write_csv(merged_df, out_csv)
cat("Saved updated RDS and CSV.\n")

cat("Regenerating Excel workbook...\n")
system2(Sys.which("python3"), args = "sim/export_100iter_results_xlsx.py")

cat("Regenerating 8 publication figures...\n")
source("sim/generate_100iter_plots.R")
generate_100iter_plots(rds_path = out_rds, out_dir = "sim/figures")

cat(sprintf("ALL PHASE 1 BENCHMARK REFRESH COMPLETED IN %.2f MINS!\n", (proc.time()["elapsed"] - t_start) / 60))

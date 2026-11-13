suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(parallel)
})

source("sim/run_simulation_condition.R")

cat("Precomputing oracle truths for Scenarios 4 and 7...\n")
truth4 <- compute_scenario4_truth()
truth7 <- compute_scenario6_truth()

jobs <- expand.grid(scenario = c(4L, 7L), iter = 1:10, KEEP.OUT.ATTRS = FALSE)

cat("Running 10 iterations for Scenarios 4 and 7 across 10 cores...\n")
t0 <- proc.time()["elapsed"]
res_list <- parallel::mclapply(seq_len(nrow(jobs)), function(idx) {
  sc <- jobs$scenario[idx]
  it <- jobs$iter[idx]
  tr <- if (sc == 4L) truth4 else truth7
  seed <- 2026L + sc * 10000L + it
  df <- run_simulation_dataset(
    scenario_num = sc,
    n = 1000,
    tau = 30,
    grace_period = 10,
    seed = seed,
    truth_precomp = tr
  )
  df$iter <- it
  df
}, mc.cores = 10L)

res_df <- dplyr::bind_rows(res_list)
saveRDS(res_df, "scratch/test_10iter_s4_s7.rds")
cat(sprintf("Completed 20 datasets in %.2f mins.\n\n", (proc.time()["elapsed"] - t0) / 60))

summary_tbl <- res_df %>%
  dplyr::filter(estimator %in% c(
    "TMLE + BART", "G-COMP (BART)", "TTE-TMLE (SuperLearner)",
    "TTE-TMLE (GLM)", "CCW-KM", "Time-Matched Cox"
  )) %>%
  dplyr::group_by(scenario_id, estimator, regime_id) %>%
  dplyr::summarise(
    bias_15 = mean(bias_15, na.rm = TRUE),
    bias_30 = mean(bias_30, na.rm = TRUE),
    rmse_30 = sqrt(mean(bias_30^2, na.rm = TRUE)),
    cov_30  = mean(covered_30, na.rm = TRUE),
    ci_w_30 = mean(ci_width_30, na.rm = TRUE),
    band_cov = mean(band_covered, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(scenario_id, estimator, regime_id)

print(summary_tbl, n = 50, width = 140)

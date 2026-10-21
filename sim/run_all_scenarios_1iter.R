#' Master Runner for 1-Iteration Simulation Benchmark Across All 8 Scenarios and 4 Regimes
#'
#' @description Executes a single complete iteration for N = 1000 across all 32 experimental
#'   conditions (8 Scenarios x 4 Specification Regimes), evaluating all 9 competing estimators,
#'   benchmarking computational runtimes, extracting RMSTD, Marginal Total HR, and simultaneous bands.
#'   Exports results to sim/simulation_results_1iter.rds, .csv, and calls sim/export_results_xlsx.py.

library(dplyr)
library(tibble)

source("sim/run_simulation_condition.R")

run_all_1iter_benchmark <- function(n = 1000,
                                    seed = 42,
                                    out_csv = "sim/simulation_results_1iter.csv",
                                    out_rds = "sim/simulation_results_1iter.rds") {
  
  cat("=================================================================================\n")
  cat("STARTING 1-ITERATION BENCHMARK ACROSS 8 SCENARIOS x 4 REGIMES (N = 1000)\n")
  cat("=================================================================================\n")
  
  total_conditions <- 8 * 4
  cond_idx <- 0
  all_results <- list()
  total_start <- proc.time()["elapsed"]
  
  for (s in 1:8) {
    for (r in 1:4) {
      cond_idx <- cond_idx + 1
      regime_desc <- switch(r,
        "1" = "Corr Q + Corr g",
        "2" = "Mis Q + Corr g",
        "3" = "Corr Q + Mis g",
        "4" = "Mis Q + Mis g"
      )
      
      cat(sprintf("[%02d/%02d] Running Scenario %d, Regime %d (%s)... ", 
                  cond_idx, total_conditions, s, r, regime_desc))
      
      t_cond_start <- proc.time()["elapsed"]
      
      # Execute condition
      cond_res <- suppressWarnings(run_simulation_condition(
        scenario_num = s,
        regime_num = r,
        n = n,
        seed = seed + cond_idx * 10,
        bart_trees = 30,
        bart_ndpost = 100,
        bart_nskip = 50,
        mc_cores = 2
      ))
      
      t_cond_elapsed <- as.numeric(proc.time()["elapsed"] - t_cond_start)
      cat(sprintf("DONE in %.2f sec\n", t_cond_elapsed))
      
      all_results[[cond_idx]] <- cond_res
    }
  }
  
  master_df <- dplyr::bind_rows(all_results)
  total_elapsed <- as.numeric(proc.time()["elapsed"] - total_start)
  
  cat("=================================================================================\n")
  cat(sprintf("ALL 32 CONDITIONS COMPLETED IN %.2f SECONDS (%.2f MINUTES)\n", total_elapsed, total_elapsed / 60))
  cat(sprintf("Total Estimator Fits: %d rows\n", nrow(master_df)))
  cat("=================================================================================\n")
  
  # Save RDS and CSV
  saveRDS(master_df, file = out_rds)
  readr::write_csv(master_df, file = out_csv)
  cat(sprintf("Saved RDS to: %s\n", out_rds))
  cat(sprintf("Saved CSV to: %s\n", out_csv))
  
  # Trigger Python Excel Generator
  cat("Generating Master Excel Workbook: sim/simulation_results_1iter.xlsx...\n")
  system("python3 sim/export_results_xlsx.py")
  cat("Excel generation complete!\n")
  
  return(invisible(master_df))
}

if (!interactive()) {
  run_all_1iter_benchmark()
}


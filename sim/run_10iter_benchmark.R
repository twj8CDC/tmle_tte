#' 10-Iteration Monte Carlo Benchmark across 8 Scenarios (Regime 1 & Regime 4, N = 1000)
#'
#' @description Executes B = 10 independent Monte Carlo replicates for N = 1000 across:
#'   - Regime 1: Both Correctly Specified (Corr Q + Corr g)
#'   - Regime 4: Both Misspecified (Mis Q + Mis g)
#'   for all 8 scenarios (16 conditions x 10 iters = 160 runs total).
#'   Leverages parallel::mclapply across available CPU cores, benchmarks runtimes,
#'   evaluates point estimates, RMSTD, Marginal Total HR, CIs, and simultaneous bands.

library(dplyr)
library(tibble)
library(parallel)

source("sim/run_simulation_condition.R")

run_10iter_benchmark <- function(n = 1000,
                                 n_iter = 10,
                                 scenarios = 1:8,
                                 regimes = c(1, 4),
                                 n_cores = 6,
                                 base_seed = 1000,
                                 out_csv = "sim/simulation_results_10iter.csv",
                                 out_rds = "sim/simulation_results_10iter.rds") {
  
  cat("=================================================================================\n")
  cat(sprintf("STARTING 10-ITERATION BENCHMARK (N = %d, B = %d iters)\n", n, n_iter))
  cat(sprintf("Scenarios: 1 to 8 | Regimes: 1 (Corr Q + Corr g) and 4 (Mis Q + Mis g)\n"))
  cat(sprintf("Parallelizing across %d CPU cores\n", n_cores))
  cat("=================================================================================\n")
  
  # Build grid of jobs
  job_grid <- expand.grid(
    iter = seq_len(n_iter),
    regime = regimes,
    scenario = scenarios,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::arrange(scenario, regime, iter)
  
  n_jobs <- nrow(job_grid)
  cat(sprintf("Total jobs to execute: %d (16 conditions x %d iters x 9 estimators = %d fits)\n", 
              n_jobs, n_iter, n_jobs * 9))
  
  t_total_start <- proc.time()["elapsed"]
  
  # Worker function for a single job
  run_single_job <- function(job_idx) {
    row <- job_grid[job_idx, ]
    s <- row$scenario
    r <- row$regime
    it <- row$iter
    
    # Unique reproducible seed per job
    job_seed <- base_seed + (s * 10000) + (r * 1000) + it * 17
    
    res <- tryCatch({
      cond_res <- run_simulation_condition(
        scenario_num = s,
        regime_num = r,
        n = n,
        seed = job_seed,
        bart_trees = 30,
        bart_ndpost = 100,
        bart_nskip = 50,
        mc_cores = 1 # single core per job to prevent oversubscription
      )
      cond_res$iter <- it
      cond_res
    }, error = function(e) {
      warning(sprintf("Error in job %d (Scen %d, Reg %d, Iter %d): %s", job_idx, s, r, it, e$message))
      NULL
    })
    
    return(res)
  }
  
  # Execute parallel jobs
  cat("Executing jobs with parallel::mclapply...\n")
  all_results <- parallel::mclapply(
    seq_len(n_jobs),
    run_single_job,
    mc.cores = n_cores,
    mc.preschedule = FALSE
  )
  
  # Remove failed jobs if any and bind
  valid_results <- all_results[!sapply(all_results, is.null)]
  master_df <- dplyr::bind_rows(valid_results)
  
  t_total_elapsed <- as.numeric(proc.time()["elapsed"] - t_total_start)
  
  cat("=================================================================================\n")
  cat(sprintf("ALL %d JOBS COMPLETED IN %.2f SECONDS (%.2f MINUTES)\n", 
              length(valid_results), t_total_elapsed, t_total_elapsed / 60))
  cat(sprintf("Total model fits collected: %d rows\n", nrow(master_df)))
  cat("=================================================================================\n")
  
  # Save RDS and CSV
  saveRDS(master_df, file = out_rds)
  readr::write_csv(master_df, file = out_csv)
  cat(sprintf("Saved RDS to: %s\n", out_rds))
  cat(sprintf("Saved CSV to: %s\n", out_csv))
  
  # Post-process into master Excel
  cat("Generating Master Excel: sim/simulation_results_10iter.xlsx...\n")
  system("python3 sim/export_10iter_results_xlsx.py")
  cat("10-Iteration Excel workbook generated!\n")
  
  return(invisible(master_df))
}

if (!interactive()) {
  run_10iter_benchmark()
}

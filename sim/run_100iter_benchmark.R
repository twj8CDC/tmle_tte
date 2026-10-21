#' 100-Iteration Monte Carlo Benchmark across 8 Scenarios (Regime 1 & Regime 4, N = 1000)
#'
#' @description Executes B = 100 independent Monte Carlo replicates for N = 1000 across:
#'   - Regime 1: Both Correctly Specified (Corr Q + Corr g)
#'   - Regime 4: Both Misspecified (Mis Q + Mis g)
#'   for all 8 scenarios (16 conditions x 100 iters = 1,600 runs total, 14,400 model fits).
#'   Leverages parallel::mclapply across 10 CPU cores, benchmarks runtimes,
#'   evaluates point estimates, RMSTD, Marginal Total HR, CIs, and simultaneous bands.
#'   Saves checkpoints after every 10-iteration batch to guarantee fault tolerance.

library(dplyr)
library(tibble)
library(parallel)
library(readr)

source("sim/run_simulation_condition.R")

run_100iter_benchmark <- function(n = 1000,
                                  n_iter = 100,
                                  batch_size = 10,
                                  scenarios = 1:8,
                                  regimes = c(1, 4),
                                  n_cores = 10,
                                  base_seed = 2026,
                                  out_csv = "sim/simulation_results_100iter.csv",
                                  out_rds = "sim/simulation_results_100iter.rds",
                                  checkpoint_dir = "sim/checkpoints_100iter") {
  
  if (!dir.exists(checkpoint_dir)) {
    dir.create(checkpoint_dir, recursive = TRUE)
  }
  
  cat("=================================================================================\n")
  cat(sprintf("STARTING 100-ITERATION BENCHMARK (N = %d, B = %d iters)\n", n, n_iter))
  cat(sprintf("Scenarios: 1 to 8 | Regimes: 1 (Corr Q + Corr g) and 4 (Mis Q + Mis g)\n"))
  cat(sprintf("Parallelizing across %d CPU cores with batch size %d iters\n", n_cores, batch_size))
  cat("=================================================================================\n")
  
  n_batches <- ceiling(n_iter / batch_size)
  all_batches_dfs <- list()
  
  t_total_start <- proc.time()["elapsed"]
  
  for (b in seq_len(n_batches)) {
    b_iter_start <- (b - 1) * batch_size + 1
    b_iter_end <- min(b * batch_size, n_iter)
    b_iters <- b_iter_start:b_iter_end
    
    cat(sprintf("\n>>> Starting Batch %d/%d (Iterations %d to %d, %d conditions x %d iters = %d jobs)...\n",
                b, n_batches, b_iter_start, b_iter_end, length(scenarios) * length(regimes), 
                length(b_iters), length(scenarios) * length(regimes) * length(b_iters)))
    
    t_batch_start <- proc.time()["elapsed"]
    
    # Build grid of jobs for this batch
    batch_grid <- expand.grid(
      iter = b_iters,
      regime = regimes,
      scenario = scenarios,
      stringsAsFactors = FALSE
    ) %>%
      dplyr::arrange(scenario, regime, iter)
    
    n_batch_jobs <- nrow(batch_grid)
    
    # Worker function for a single job
    run_batch_job <- function(job_idx) {
      row <- batch_grid[job_idx, ]
      s <- row$scenario
      r <- row$regime
      it <- row$iter
      
      job_seed <- base_seed + (s * 100000) + (r * 10000) + it * 37
      
      res <- tryCatch({
        cond_res <- run_simulation_condition(
          scenario_num = s,
          regime_num = r,
          n = n,
          seed = job_seed,
          bart_trees = 30,
          bart_ndpost = 100,
          bart_nskip = 50,
          mc_cores = 1
        )
        cond_res$iter <- it
        cond_res
      }, error = function(e) {
        warning(sprintf("Error in job %d (Scen %d, Reg %d, Iter %d): %s", job_idx, s, r, e$message))
        NULL
      })
      
      return(res)
    }
    
    # Execute batch jobs in parallel
    batch_results <- parallel::mclapply(
      seq_len(n_batch_jobs),
      run_batch_job,
      mc.cores = n_cores,
      mc.preschedule = FALSE
    )
    
    valid_batch_results <- batch_results[!sapply(batch_results, is.null)]
    batch_df <- dplyr::bind_rows(valid_batch_results)
    all_batches_dfs[[b]] <- batch_df
    
    t_batch_elapsed <- as.numeric(proc.time()["elapsed"] - t_batch_start)
    t_total_so_far <- as.numeric(proc.time()["elapsed"] - t_total_start)
    est_total_time <- (t_total_so_far / b) * n_batches
    est_time_remaining <- est_total_time - t_total_so_far
    
    # Save batch checkpoint
    checkpoint_file <- file.path(checkpoint_dir, sprintf("batch_%02d.rds", b))
    saveRDS(batch_df, file = checkpoint_file)
    
    # Update cumulative master files
    master_df <- dplyr::bind_rows(all_batches_dfs)
    saveRDS(master_df, file = out_rds)
    readr::write_csv(master_df, file = out_csv)
    
    cat(sprintf(">>> Batch %d/%d completed in %.2f mins. Cumulative: %d rows (%d fits). Elapsed: %.2f mins. Est. remaining: %.2f mins.\n",
                b, n_batches, t_batch_elapsed / 60, nrow(master_df), nrow(master_df),
                t_total_so_far / 60, est_time_remaining / 60))
  }
  
  t_total_elapsed <- as.numeric(proc.time()["elapsed"] - t_total_start)
  
  cat("=================================================================================\n")
  cat(sprintf("ALL %d BATCHES COMPLETED IN %.2f MINUTES!\n", n_batches, t_total_elapsed / 60))
  cat(sprintf("Total model fits collected: %d rows\n", nrow(master_df)))
  cat(sprintf("Saved cumulative RDS: %s\n", out_rds))
  cat(sprintf("Saved cumulative CSV: %s\n", out_csv))
  cat("=================================================================================\n")
  
  cat("Generating Master Excel: sim/simulation_results_100iter.xlsx...\n")
  system("python3 sim/export_100iter_results_xlsx.py")
  cat("100-Iteration Excel workbook generated successfully!\n")
  
  return(invisible(master_df))
}

if (!interactive()) {
  run_100iter_benchmark()
}


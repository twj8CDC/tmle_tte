#' 100-Iteration Monte Carlo Benchmark across 8 Scenarios (Single Dataset per Scenario-Iteration)
#'
#' @description Executes B = 100 independent Monte Carlo replicates for N = 1000 across all 8 scenarios
#'   (8 scenarios x 100 datasets = 800 datasets).
#'   On each dataset:
#'   - Non-parametric estimators (TTE-TMLE SuperLearner, G-COMP BART, TMLE + BART) are fit ONCE.
#'   - Parametric benchmarks (Naive Cox, Landmark Cox, TD-Cox, Time-Matched Cox, CCW-CPHM, CCW-KM)
#'     are fit under Correct and Misspecified functional forms (never omitting variables).
#'   - TTE-TMLE (GLM) is fit across all 4 specification regimes (Corr/Mis Q x Corr/Mis g).

library(dplyr)
library(tibble)
library(parallel)
library(readr)

source("sim/run_simulation_condition.R")

run_100iter_benchmark <- function(n = 1000,
                                  n_iter = 100,
                                  batch_size = 10,
                                  scenarios = 1:8,
                                  n_cores = 10,
                                  base_seed = 2026,
                                  out_csv = "sim/simulation_results_100iter.csv",
                                  out_rds = "sim/simulation_results_100iter.rds",
                                  checkpoint_dir = "sim/checkpoints_100iter") {
  if (dir.exists(checkpoint_dir)) {
    unlink(checkpoint_dir, recursive = TRUE, force = TRUE)
  }
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

  cat("=================================================================================\n")
  cat(sprintf("PRECOMPUTING ORACLE GROUND TRUTH FOR SCENARIOS 1 TO 8 (N_oracle = 200,000)...\n"))
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
  cat("Oracle ground truths precomputed successfully!\n")
  cat("=================================================================================\n")
  cat(sprintf("STARTING 100-ITERATION BENCHMARK (N = %d, B = %d datasets per scenario)\n", n, n_iter))
  cat(sprintf("Scenarios: 1 to 8 | 800 total datasets | Parallelizing across %d CPU cores\n", n_cores))
  cat("=================================================================================\n")

  n_batches <- ceiling(n_iter / batch_size)
  all_batches_dfs <- list()
  t_total_start <- proc.time()["elapsed"]

  for (b in seq_len(n_batches)) {
    b_iter_start <- (b - 1) * batch_size + 1
    b_iter_end <- min(b * batch_size, n_iter)
    b_iters <- b_iter_start:b_iter_end

    batch_grid <- expand.grid(
      iter = b_iters,
      scenario = scenarios,
      stringsAsFactors = FALSE
    ) %>%
      dplyr::arrange(scenario, iter)

    n_batch_jobs <- nrow(batch_grid)

    cat(sprintf(
      "\n>>> Starting Batch %d/%d (Iterations %d to %d, %d scenarios x %d iters = %d dataset jobs)...\n",
      b, n_batches, b_iter_start, b_iter_end, length(scenarios), length(b_iters), n_batch_jobs
    ))

    t_batch_start <- proc.time()["elapsed"]

    run_batch_job <- function(job_idx) {
      row <- batch_grid[job_idx, ]
      s <- row$scenario
      it <- row$iter
      job_seed <- base_seed + (s * 100000L) + it * 37L

      res <- tryCatch({
        cond_res <- run_simulation_dataset(
          scenario_num = s,
          n = n,
          seed = job_seed,
          bart_trees = 30,
          bart_ndpost = 100,
          bart_nskip = 100,
          bart_thin = 2,
          truth_precomp = truth_list[[as.character(s)]]
        )
        cond_res$iter <- it
        cond_res
      }, error = function(e) {
        warning(sprintf("Error in job %d (Scen %d, Iter %d): %s", job_idx, s, it, e$message))
        NULL
      })

      return(res)
    }

    batch_results <- parallel::mclapply(
      seq_len(n_batch_jobs),
      run_batch_job,
      mc.cores = n_cores,
      mc.preschedule = FALSE
    )

    valid_batch_results <- batch_results[!vapply(batch_results, is.null, logical(1))]
    batch_df <- dplyr::bind_rows(valid_batch_results)
    all_batches_dfs[[b]] <- batch_df

    t_batch_elapsed <- as.numeric(proc.time()["elapsed"] - t_batch_start)
    t_total_so_far <- as.numeric(proc.time()["elapsed"] - t_total_start)
    est_total_time <- (t_total_so_far / b) * n_batches
    est_time_remaining <- est_total_time - t_total_so_far

    checkpoint_file <- file.path(checkpoint_dir, sprintf("batch_%02d.rds", b))
    saveRDS(batch_df, file = checkpoint_file)

    master_df <- dplyr::bind_rows(all_batches_dfs)
    saveRDS(master_df, file = out_rds)
    readr::write_csv(master_df, file = out_csv)

    cat(sprintf(
      ">>> Batch %d/%d completed in %.2f mins. Cumulative: %d rows. Elapsed: %.2f mins. Est. remaining: %.2f mins.\n",
      b, n_batches, t_batch_elapsed / 60, nrow(master_df),
      t_total_so_far / 60, est_time_remaining / 60
    ))
  }

  t_total_elapsed <- as.numeric(proc.time()["elapsed"] - t_total_start)

  cat("=================================================================================\n")
  cat(sprintf("ALL %d BATCHES COMPLETED IN %.2f MINUTES!\n", n_batches, t_total_elapsed / 60))
  cat(sprintf("Total model fits collected: %d rows\n", nrow(master_df)))
  cat(sprintf("Saved cumulative RDS: %s\n", out_rds))
  cat(sprintf("Saved cumulative CSV: %s\n", out_csv))
  cat("=================================================================================\n")

  cat("Generating Master Excel: sim/simulation_results_100iter.xlsx...\n")
  py_bin <- Sys.which("python3")
  system2(py_bin, args = "sim/export_100iter_results_xlsx.py")
  cat("100-Iteration Excel workbook generated successfully!\n")

  cat("Generating Publication Benchmark Figures: sim/figures/scenario_1..8_benchmark_100iter.png...\n")
  source("sim/generate_100iter_plots.R")
  generate_100iter_plots(rds_path = out_rds, out_dir = "sim/figures")
  cat("All 8 scenario benchmark figures generated and copied to artifact directory!\n")

  return(invisible(master_df))
}

if (!interactive()) {
  run_100iter_benchmark()
}

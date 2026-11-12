suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(parallel)
  library(readr)
})

source("sim/run_simulation_condition.R")

t_start <- proc.time()["elapsed"]

out_rds <- "sim/scratch/gcomp_bart_standard_s1_s2.rds"
out_csv <- "sim/scratch/gcomp_bart_standard_s1_s2.csv"

cat("Precomputing oracle ground truth for Scenarios 1 & 2 (N_oracle = 200,000)...\n")
truth_list <- list(
  `1` = compute_scenario1_truth(),
  `2` = compute_scenario2_truth()
)

base_seed <- 2026L
tau <- 30L
grace_period <- 10L
n_obs <- 1000L

# Full SurvBART standard settings (trees=50, draws=2500/chain, tune=1000, thin=10, chains=4 -> 1000 draws)
std_ntree <- 50L
std_ndpost <- 250L
std_nskip <- 1000L
std_thin <- 10L
std_chains <- 4L
std_cores <- 1L

run_single_gcomp_std <- function(s, it) {
  job_seed <- base_seed + (s * 100000L) + it * 37L
  set.seed(job_seed)
  truth <- truth_list[[as.character(s)]]
  spec <- get_scenario_spec(
    scenario_num = s,
    n = n_obs,
    tau = tau,
    grace_period = grace_period,
    seed = job_seed,
    truth_precomp = truth
  )
  df_std <- prepare_tte_data(spec$df_raw)
  covars_all <- spec$covars_all
  clones <- create_clones(df_std, covars = covars_all, grace_period = grace_period, tau = tau)

  fit_bg <- fit_bart_gcomp(
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
    seed = job_seed
  )
  res_bg <- fit_bg$results

  rd_true_15 <- truth$rd_true[15]
  rd_true_30 <- truth$rd_true[30]
  rmstd_true_30 <- sum(truth$rd_true[1:tau], na.rm = TRUE)
  hr_marg_true_30 <- if ("hr_true" %in% names(truth) && length(truth$hr_true) >= tau) {
    truth$hr_true[tau]
  } else {
    log(pmax(pmin(truth$S1_true[tau], 1 - 1e-6), 1e-6)) /
      log(pmax(pmin(truth$S0_true[tau], 1 - 1e-6), 1e-6))
  }

  ev_trt_days <- sort(ceiling(df_std$time[df_std$tx_ever == 1L & df_std$t_switch <= grace_period & df_std$event == 1L & df_std$time >= df_std$t_switch]))
  t_band_start <- if (length(ev_trt_days) >= 5L) max(1L, as.integer(ev_trt_days[5L])) else as.integer(grace_period)
  idx <- seq(from = min(t_band_start, tau), to = tau)
  is_cov <- (res_bg$band_lower[idx] <= truth$rd_true[idx]) & (truth$rd_true[idx] <= res_bg$band_upper[idx])

  rd15 <- res_bg$rd[15]
  se15 <- res_bg$se_RD[15]
  ci15_l <- res_bg$cri_lower[15]
  ci15_u <- res_bg$cri_upper[15]
  rd30 <- res_bg$rd[tau]
  se30 <- res_bg$se_RD[tau]
  ci30_l <- res_bg$cri_lower[tau]
  ci30_u <- res_bg$cri_upper[tau]
  rmstd <- fit_bg$rmstd
  rmstd_l <- fit_bg$rmstd_ci_lower
  rmstd_u <- fit_bg$rmstd_ci_upper
  hr <- res_bg$hr_cum[tau]
  hr_l <- res_bg$hr_ci_lower[tau]
  hr_u <- res_bg$hr_ci_upper[tau]

  tibble(
    scenario_id = paste0("Scenario ", s),
    iter = it,
    estimator = "G-COMP (BART-Standard)",
    spec_group = "Non-Parametric",
    elapsed_sec = fit_bg$elapsed_sec,
    rd_15_est = rd15,
    rd_15_true = rd_true_15,
    bias_15 = rd15 - rd_true_15,
    se_15 = se15,
    ci_lower_15 = ci15_l,
    ci_upper_15 = ci15_u,
    covered_15 = as.integer(ci15_l <= rd_true_15 & rd_true_15 <= ci15_u),
    ci_width_15 = ci15_u - ci15_l,
    rd_30_est = rd30,
    rd_30_true = rd_true_30,
    bias_30 = rd30 - rd_true_30,
    se_30 = se30,
    ci_lower_30 = ci30_l,
    ci_upper_30 = ci30_u,
    covered_30 = as.integer(ci30_l <= rd_true_30 & rd_true_30 <= ci30_u),
    ci_width_30 = ci30_u - ci30_l,
    rmstd_est = rmstd,
    rmstd_true = rmstd_true_30,
    rmstd_bias = rmstd - rmstd_true_30,
    se_rmstd = fit_bg$se_rmstd,
    rmstd_ci_lower = rmstd_l,
    rmstd_ci_upper = rmstd_u,
    rmstd_covered = as.integer(rmstd_l <= rmstd_true_30 & rmstd_true_30 <= rmstd_u),
    rmstd_ci_width = rmstd_u - rmstd_l,
    hr_marginal_est = hr,
    hr_marginal_true = hr_marg_true_30,
    hr_marginal_bias = hr - hr_marg_true_30,
    hr_ci_lower_30 = hr_l,
    hr_ci_upper_30 = hr_u,
    hr_covered_30 = as.integer(hr_l <= hr_marg_true_30 & hr_marg_true_30 <= hr_u),
    hr_ci_width_30 = hr_u - hr_l,
    uniform_band_covered = as.integer(all(is_cov, na.rm = TRUE)),
    prop_time_points_covered = mean(is_cov, na.rm = TRUE),
    mean_band_width = mean(res_bg$band_upper[idx] - res_bg$band_lower[idx], na.rm = TRUE)
  )
}

grid <- expand.grid(iter = 1:100, scenario = 1:2, stringsAsFactors = FALSE) %>%
  arrange(scenario, iter)

batch_size <- 48L # 12 cores x 4 waves per batch (~10 mins/batch), 5 batches total
n_batches <- ceiling(nrow(grid) / batch_size)
all_batches <- vector("list", n_batches)

cat(sprintf(
  "Starting G-COMP (BART) Standard benchmark: 200 datasets (Scen 1 & 2 x 100 iters)\nConfig: trees=%d, draws=%d (ndpost=%d x thin=%d), tune=%d, chains=%d, 12 cores\n",
  std_ntree, std_ndpost * std_thin, std_ndpost, std_thin, std_nskip, std_chains
))

for (b in seq_len(n_batches)) {
  idx_start <- (b - 1L) * batch_size + 1L
  idx_end <- min(b * batch_size, nrow(grid))
  sub_idx <- seq(idx_start, idx_end)
  t_b <- proc.time()["elapsed"]

  res_b <- mclapply(sub_idx, function(k) {
    tryCatch(
      run_single_gcomp_std(grid$scenario[k], grid$iter[k]),
      error = function(e) {
        warning(sprintf("Error Scen %d Iter %d: %s", grid$scenario[k], grid$iter[k], e$message))
        NULL
      }
    )
  }, mc.cores = 12L, mc.preschedule = FALSE)

  n_null <- sum(vapply(res_b, is.null, logical(1)))
  if (n_null > 0L) stop(sprintf("Batch %d had %d failed jobs!", b, n_null))

  all_batches[[b]] <- bind_rows(res_b)
  cum_df <- bind_rows(all_batches[seq_len(b)])
  saveRDS(cum_df, out_rds)
  write_csv(cum_df, out_csv)

  el_min <- (proc.time()["elapsed"] - t_start) / 60
  b_min <- (proc.time()["elapsed"] - t_b) / 60
  cat(sprintf(
    ">>> Batch %d/%d finished in %.2f mins (%d/200 datasets done, elapsed=%.2f mins).\n",
    b, n_batches, b_min, nrow(cum_df), el_min
  ))
}

# Compare directly against stored fast G-COMP (BART) & TTE-TMLE (BART) rows
orig <- readRDS("sim/simulation_results_100iter.rds") %>%
  filter(scenario_id %in% c("Scenario 1", "Scenario 2"),
         estimator %in% c("G-COMP (BART)", "TTE-TMLE (BART)", "TTE-TMLE (GLM)") &
           spec_group %in% c("Non-Parametric", "Parametric (Correct)")) %>%
  mutate(estimator = ifelse(estimator == "G-COMP (BART)", "G-COMP (BART-Fast)", estimator))

comp <- bind_rows(
  cum_df,
  orig %>% select(any_of(names(cum_df)))
) %>%
  group_by(scenario_id, estimator) %>%
  summarise(
    reps = n(),
    bias15 = mean(bias_15),
    bias30 = mean(bias_30),
    rmse30 = sqrt(mean(bias_30^2)),
    cov30 = mean(covered_30) * 100,
    se_ratio_30 = mean(se_30) / sd(rd_30_est),
    ci_w30 = mean(ci_width_30),
    band_cov = mean(uniform_band_covered) * 100,
    rmstd_b = mean(rmstd_bias),
    rmstd_c = mean(rmstd_covered) * 100,
    hr_b = mean(hr_marginal_bias),
    hr_c = mean(hr_covered_30) * 100,
    time_s = mean(elapsed_sec),
    .groups = "drop"
  ) %>%
  arrange(scenario_id, estimator)

cat("\n=== Scenario 1 & 2: G-COMP (BART-Standard) vs Fast Benchmark ===\n")
print(as.data.frame(comp), digits = 4)
cat(sprintf("\nCOMPLETED IN %.2f MINS!\n", (proc.time()["elapsed"] - t_start) / 60))

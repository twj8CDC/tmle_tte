suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
})

source("R/00_utils.R")
source("R/01_cloning.R")
source("R/02_ipcw_models.R")
source("R/03_q_models.R")
source("R/04_targeting.R")
source("R/05_gcomp.R")
source("R/06_inference.R")
source("R/07_tte_tmle.R")
source("R/09_bart_gcomp.R")
source("sim/dgp/03_dgp_waning_nph.R")
source("sim/true_values.R")

cat("Computing ground truth for Scenario 3 (oracle N = 100,000)...\n")
truth3 <- compute_scenario3_truth(n_oracle = 100000, seed = 99999)
true_rd15 <- truth3$rd_true[15]
true_rd30 <- truth3$rd_true[30]
cat(sprintf("Ground Truth: Day 15 RD = %.4f | Day 30 RD = %.4f\n", true_rd15, true_rd30))

reps <- 5
seeds <- 100 + seq_len(reps)

results_list <- list()

for (r in seq_len(reps)) {
  s <- seeds[r]
  cat(sprintf("\n--- Replicate %d/%d (Seed %d) ---\n", r, reps, s))
  data_r <- generate_scenario3_data(n = 1000, seed = s)
  prep_r <- prepare_tte_data(data_r)
  clones_r <- create_clones(prep_r, covars = c("w1", "w2"), grace_period = 10, tau = 30)

  # 1. GLM Pooled TMLE
  fit_glm_pooled <- tte_tmle(
    data = data_r,
    grace_period = 10,
    tau = 30,
    covars = c("w1", "w2"),
    targeting_method = "pooled_glm",
    tx_strategy = "conditional_complier",
    n_boot = 500,
    seed = s
  )

  # 2. GLM One-Step TMLE
  fit_glm_onestep <- tte_tmle(
    data = data_r,
    grace_period = 10,
    tau = 30,
    covars = c("w1", "w2"),
    targeting_method = "one_step",
    tx_strategy = "conditional_complier",
    n_boot = 500,
    seed = s
  )

  # 3. BART Pooled TMLE
  fit_bart_pooled <- fit_bart_tmle(
    data = prep_r,
    clones = clones_r,
    covars = c("w1", "w2"),
    grace_period = 10,
    tau = 30,
    ntree = 30,
    ndpost = 100,
    nskip = 100,
    thin = 2,
    targeting_method = "pooled_glm",
    tx_strategy = "conditional_complier",
    n_boot = 500,
    seed = s
  )

  # 4. BART One-Step TMLE
  fit_bart_onestep <- fit_bart_tmle(
    data = prep_r,
    clones = clones_r,
    covars = c("w1", "w2"),
    grace_period = 10,
    tau = 30,
    ntree = 30,
    ndpost = 100,
    nskip = 100,
    thin = 2,
    targeting_method = "one_step",
    tx_strategy = "conditional_complier",
    n_boot = 500,
    seed = s
  )

  res_r <- tibble(
    replicate = r,
    seed = s,
    method = c("GLM (Pooled)", "GLM (One-Step)", "BART (Pooled)", "BART (One-Step)"),
    rd15 = c(fit_glm_pooled$res_tmle$rd[15], fit_glm_onestep$res_tmle$rd[15],
             fit_bart_pooled$results$rd[15], fit_bart_onestep$results$rd[15]),
    bias15 = rd15 - true_rd15,
    rd30 = c(fit_glm_pooled$res_tmle$rd[30], fit_glm_onestep$res_tmle$rd[30],
             fit_bart_pooled$results$rd[30], fit_bart_onestep$results$rd[30]),
    bias30 = rd30 - true_rd30
  )
  results_list[[r]] <- res_r
  print(res_r %>% dplyr::select(method, rd15, bias15, rd30, bias30))
}

all_res <- bind_rows(results_list)

summary_table <- all_res %>%
  group_by(method) %>%
  summarise(
    mean_rd15 = mean(rd15),
    bias_d15 = mean(bias15),
    rmse_d15 = sqrt(mean(bias15^2)),
    mean_rd30 = mean(rd30),
    bias_d30 = mean(bias30),
    rmse_d30 = sqrt(mean(bias30^2)),
    .groups = "drop"
  )

cat("\n=======================================================\n")
cat("SUMMARY BENCHMARK: POOLED_GLM VS ONE_STEP (SCENARIO 3)\n")
cat("=======================================================\n")
print(summary_table)
saveRDS(list(all_res = all_res, summary = summary_table), "scratch/summary_onestep_scen3.rds")

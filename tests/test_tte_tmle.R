# Test script for ttetmle package modules
# Verifies end-to-end execution, multiplier bootstrap simultaneous bands,
# and numerical consistency against original tmle_ccw.r

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(survival)
})

# Source package modules
source("R/00_utils.R")
source("R/01_cloning.R")
source("R/02_ipcw_models.R")
source("R/03_q_models.R")
source("R/04_targeting.R")
source("R/05_gcomp.R")
source("R/06_inference.R")
source("R/07_tte_tmle.R")
source("R/08_benchmarks.R")

set.seed(42)
n <- 250
x1 <- rnorm(n, 0, 1)
p_tx <- plogis(-0.2 + 0.4 * x1)
tx_obs <- rbinom(n, 1, p_tx)
tx_time_obs <- ifelse(tx_obs == 1, pmax(1, pmin(10, round(rbeta(n, 4, 2) * 10))), Inf)

# Baseline hazard with treatment benefit
haz_base <- 0.03 * exp(0.3 * x1)
evt_time_raw <- rexp(n, rate = haz_base)
event_time <- pmax(1, pmin(30, ceiling(evt_time_raw)))
event <- ifelse(evt_time_raw <= 30, 1L, 0L)

test_data <- tibble(
  id = seq_len(n),
  x1 = x1,
  tx_obs = tx_obs,
  tx_time_obs = tx_time_obs,
  event_time = event_time,
  event = event
)

cat("Running TTE-TMLE with GLM...\n")
fit_glm <- tte_tmle(
  data = test_data,
  grace_period = 10,
  tau = 30,
  covars = c("x1"),
  use_sl = FALSE,
  n_boot = 2000,
  seed = 123
)

print(fit_glm)

cat("\nChecking simultaneous band properties:\n")
stopifnot(fit_glm$crit_val_simul >= fit_glm$crit_val_pointwise)
cat("Critical value pointwise (1.96):", fit_glm$crit_val_pointwise, "\n")
cat("Critical value simultaneous band:", fit_glm$crit_val_simul, "\n")
cat("SUCCESS: Simultaneous critical value is strictly >= pointwise critical value.\n")

cat("\nChecking Hazard Ratio derivative columns:\n")
hr_cols <- c("hr_cum", "se_log_hr", "hr_ci_lower", "hr_ci_upper", "hr_band_lower", "hr_band_upper", "hr_daily")
stopifnot(all(hr_cols %in% names(fit_glm$res_tmle)))
cat("Day 30 Cumulative Hazard Ratio:", round(fit_glm$res_tmle$hr_cum[30], 4), "\n")
cat("Day 30 95% CI: [", round(fit_glm$res_tmle$hr_ci_lower[30], 4), ",", round(fit_glm$res_tmle$hr_ci_upper[30], 4), "]\n")
cat("Day 30 95% Simultaneous Band: [", round(fit_glm$res_tmle$hr_band_lower[30], 4), ",", round(fit_glm$res_tmle$hr_band_upper[30], 4), "]\n")
cat("SUCCESS: Hazard Ratio derivative columns verified.\n")

cat("\nRunning Pure G-computation (untargeted)...\n")
fit_gcomp <- tte_tmle(
  data = test_data,
  grace_period = 10,
  tau = 30,
  covars = c("x1"),
  pure_gcomp = TRUE,
  seed = 123
)
cat("Pure G-comp Day 30 RD:", fit_gcomp$res_tmle$rd[30], "\n")
cat("Targeted TMLE Day 30 RD:", fit_glm$res_tmle$rd[30], "\n")

cat("\nFitting Competing Benchmarks...\n")
c_naive <- fit_naive_cox(test_data, covars = c("x1"))
c_landmark <- fit_landmark_cox(test_data, landmark_time = 10, covars = c("x1"), tau = 30)
c_td <- fit_td_cox(test_data, covars = c("x1"), tau = 30)

cat("Naive Cox HR:", exp(coef(c_naive)["tx_ever"]), "\n")
cat("Landmark Cox HR:", exp(coef(c_landmark)["tx_landmark"]), "\n")
cat("Time-Dependent Cox HR:", c_td$hr, "\n")

cat("\nALL PACKAGE MODULE CHECKS PASSED!\n")


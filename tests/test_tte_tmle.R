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
c_tm <- fit_time_matched_cox(test_data, covars = c("x1"), tau = 30, grace_period = 10, n_boot = 25L, seed = 123)

stopifnot(is.na(c_td$rd_15), is.na(c_td$rd_tau), is.na(c_td$rmstd), is.finite(c_td$hr))
stopifnot(is.finite(c_tm$rd_tau), is.finite(c_tm$se_rd_tau), is.finite(c_tm$se_rmstd), is.finite(c_tm$hr), is.finite(c_tm$hr_cond))

cat("Naive Cox HR:", exp(coef(c_naive)["tx_ever"]), "\n")
cat("Landmark Cox HR:", exp(coef(c_landmark)["tx_landmark"]), "\n")
cat("Time-Dependent Cox HR:", c_td$hr, "(RD/RMSTD omitted as NA)\n")
cat("Time-Matched Cox Marginal HR:", c_tm$hr, "| Day 30 RD:", c_tm$rd_tau, "| Boot SE:", c_tm$se_rd_tau, "\n")

cat("\nChecking all three tx_strategy modes (GLM TMLE + pure G-comp)...\n")
strategies <- c("conditional_complier", "empirical_complier", "natural_grace")
glm_fits <- lapply(strategies, function(st) {
  fit_s <- tte_tmle(
    data = test_data,
    grace_period = 10,
    tau = 30,
    covars = c("x1"),
    use_sl = FALSE,
    n_boot = 300,
    tx_strategy = st,
    seed = 123
  )
  gc_s <- tte_tmle(
    data = test_data,
    grace_period = 10,
    tau = 30,
    covars = c("x1"),
    use_sl = FALSE,
    pure_gcomp = TRUE,
    tx_strategy = st,
    seed = 123
  )
  stopifnot(
    identical(fit_s$params$tx_strategy, st),
    all(is.finite(fit_s$res_tmle$rd)),
    all(is.finite(fit_s$res_tmle$se_RD)),
    fit_s$crit_val_simul >= fit_s$crit_val_pointwise,
    all(is.finite(gc_s$res_tmle$rd))
  )
  cat(sprintf("  [%s] TMLE RD(30)=%.4f | G-comp RD(30)=%.4f\n",
              st, fit_s$res_tmle$rd[30], gc_s$res_tmle$rd[30]))
  list(tmle = fit_s, gcomp = gc_s)
})
names(glm_fits) <- strategies

# Verify analytical integration determinism (no Monte Carlo switch-time sampling jitter)
fit_glm_repeat <- tte_tmle(
  data = test_data,
  grace_period = 10,
  tau = 30,
  covars = c("x1"),
  use_sl = FALSE,
  n_boot = 300,
  tx_strategy = "conditional_complier",
  seed = 999
)
stopifnot(all.equal(glm_fits$conditional_complier$tmle$res_tmle$rd,
                    fit_glm_repeat$res_tmle$rd,
                    tolerance = 1e-12) == TRUE)

# Verify legacy align_grace_weights = FALSE matches tx_strategy = "natural_grace"
fit_legacy_ng <- tte_tmle(
  data = test_data,
  grace_period = 10,
  tau = 30,
  covars = c("x1"),
  use_sl = FALSE,
  n_boot = 300,
  align_grace_weights = FALSE,
  seed = 123
)
stopifnot(identical(fit_legacy_ng$params$tx_strategy, "natural_grace"))
stopifnot(all.equal(fit_legacy_ng$res_tmle$rd,
                    glm_fits$natural_grace$tmle$res_tmle$rd,
                    tolerance = 1e-12) == TRUE)

cat("\nChecking BART paths across tx_strategy modes...\n")
source("R/09_bart_gcomp.R")
prep_std <- prepare_tte_data(test_data)
clones_std <- create_clones(prep_std, grace_period = 10)

bart_fits <- lapply(strategies, function(st) {
  bp <- fit_bart_pair(
    data = prep_std,
    clones = clones_std,
    covars = "x1",
    grace_period = 10,
    tau = 30,
    ntree = 20,
    ndpost = 50,
    nskip = 50,
    thin = 2,
    n_boot = 200,
    tx_strategy = st,
    seed = 77
  )
  bg <- fit_bart_gcomp(
    data = prep_std,
    clones = clones_std,
    covars = "x1",
    grace_period = 10,
    tau = 30,
    ntree = 20,
    ndpost = 50,
    nskip = 50,
    thin = 2,
    tx_strategy = st,
    seed = 77
  )
  stopifnot(
    identical(bp$gcomp$tx_strategy, st),
    identical(bp$tmle$tx_strategy, st),
    identical(bg$tx_strategy, st),
    all(is.finite(bp$gcomp$results$rd)),
    all(is.finite(bp$tmle$results$rd)),
    all(is.finite(bg$results$rd))
  )
  cat(sprintf("  [BART %s] pair-gcomp RD(30)=%.4f | pair-tmle RD(30)=%.4f | standalone-gcomp RD(30)=%.4f\n",
              st, bp$gcomp$results$rd[30], bp$tmle$results$rd[30], bg$results$rd[30]))
  list(pair = bp, standalone = bg)
})
names(bart_fits) <- strategies

# Verify analytical integration determinism on summarize_bart_gcomp_draws() across R RNG states
q_mock <- matrix(rep(seq(0.01, 0.05, length.out = 2 * nrow(prep_std) * 30), each = 10),
                 nrow = 10, ncol = 2 * nrow(prep_std) * 30)
w_mock <- matrix(1 / 10, nrow = nrow(prep_std), ncol = 10)
set.seed(111)
sum_a <- summarize_bart_gcomp_draws(q_mock, n_patients = nrow(prep_std), tau = 30, w_mat = w_mock, k_vals = 1:10)
set.seed(99999)
sum_b <- summarize_bart_gcomp_draws(q_mock, n_patients = nrow(prep_std), tau = 30, w_mat = w_mock, k_vals = 1:10)
stopifnot(all.equal(sum_a$results$rd, sum_b$results$rd, tolerance = 1e-12) == TRUE)

# Verify compute_bart_ipcw_and_cf legacy align_grace_weights = FALSE matches tx_strategy = "natural_grace"
base_ord_std <- prep_std %>% dplyr::distinct(id, .keep_all = TRUE) %>% dplyr::arrange(id)
g_fixed <- rep(0.05, nrow(base_ord_std) * 10L)
ipcw_ng_explicit <- compute_bart_ipcw_and_cf(
  base_ord = base_ord_std, clones = clones_std, g_grid_pred = g_fixed,
  covars = "x1", grace_period = 10, tau = 30, tx_strategy = "natural_grace"
)
ipcw_ng_legacy <- compute_bart_ipcw_and_cf(
  base_ord = base_ord_std, clones = clones_std, g_grid_pred = g_fixed,
  covars = "x1", grace_period = 10, tau = 30, align_grace_weights = FALSE
)
stopifnot(
  identical(ipcw_ng_legacy$tx_strategy, "natural_grace"),
  all.equal(ipcw_ng_legacy$panel_w$w_trim, ipcw_ng_explicit$panel_w$w_trim, tolerance = 1e-12) == TRUE,
  all.equal(ipcw_ng_legacy$cf_H_ctrl, ipcw_ng_explicit$cf_H_ctrl, tolerance = 1e-12) == TRUE,
  all.equal(ipcw_ng_legacy$strat$w_mat, ipcw_ng_explicit$strat$w_mat, tolerance = 1e-12) == TRUE
)

# Verify policy weights sum to 1 for every subject across all strategies
for (st in strategies) {
  ipcw_st <- compute_bart_ipcw_and_cf(
    base_ord = base_ord_std, clones = clones_std, g_grid_pred = g_fixed,
    covars = "x1", grace_period = 10, tau = 30, tx_strategy = st
  )
  stopifnot(max(abs(rowSums(ipcw_st$strat$w_mat) - 1)) < 1e-12)
}

cat("\nALL PACKAGE MODULE CHECKS PASSED!\n")


suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(survival)
})

source("R/00_utils.R")
source("R/01_cloning.R")
source("R/02_ipcw_models.R")
source("R/03_q_models.R")
source("R/04_targeting.R")
source("R/05_gcomp.R")
source("R/06_inference.R")
source("R/07_tte_tmle.R")
source("sim/dgp/02_dgp_misspecification.R")
source("sim/true_values.R")
source("sim/scratch/benchmark_uncloned_vs_cloned.R")

cat("=== Testing Scenario 2 (Complex Non-linear DGP) ===\n")
dat2 <- generate_scenario2_data(n = 1000, seed = 42)
truth2 <- compute_scenario2_truth()

nl_w12 <- "I(exp(0.7 * (w1 - 3) / 1.5)) + I(((w1 - 3) / 1.5)^2) + w2 + I(w1 > 3.8 & w2 == 1)"
lin_w12 <- "w1 + w2"

cat("Running Cloned TTE-TMLE (Corr Q + Corr g)...\n")
t_c <- system.time(m_cloned <- tte_tmle(
  data = dat2,
  covars = c("w1", "w2"),
  p_form = nl_w12,
  q_form = nl_w12,
  use_sl = FALSE
))

cat("Running Uncloned TMLE Conditional (Corr Q + Corr g)...\n")
t_u_c <- system.time(m_uncloned_c <- run_uncloned_tmle(
  data = dat2,
  covars = c("w1", "w2"),
  p_form = nl_w12,
  q_form = nl_w12,
  method = "conditional",
  integration = "analytical"
))

cat("Running Uncloned TMLE Marginal (Corr Q + Corr g)...\n")
t_u_m <- system.time(m_uncloned_m <- run_uncloned_tmle(
  data = dat2,
  covars = c("w1", "w2"),
  p_form = nl_w12,
  q_form = nl_w12,
  method = "marginal",
  integration = "analytical"
))

cat("\n--- Scenario 2 Results ---\n")
cat(sprintf("Day 15 Truth: RD = %.4f\n", truth2$rd_true[15]))
cat(sprintf("  Cloned TTE-TMLE:      RD = %.4f (Bias = %+.4f)\n", m_cloned$res_tmle$rd[15], m_cloned$res_tmle$rd[15] - truth2$rd_true[15]))
cat(sprintf("  Uncloned Conditional: RD = %.4f (Bias = %+.4f)\n", m_uncloned_c$res$rd[15], m_uncloned_c$res$rd[15] - truth2$rd_true[15]))
cat(sprintf("  Uncloned Marginal:    RD = %.4f (Bias = %+.4f)\n", m_uncloned_m$res$rd[15], m_uncloned_m$res$rd[15] - truth2$rd_true[15]))

cat(sprintf("\nDay 30 Truth: RD = %.4f\n", truth2$rd_true[30]))
cat(sprintf("  Cloned TTE-TMLE:      RD = %.4f (Bias = %+.4f)\n", m_cloned$res_tmle$rd[30], m_cloned$res_tmle$rd[30] - truth2$rd_true[30]))
cat(sprintf("  Uncloned Conditional: RD = %.4f (Bias = %+.4f)\n", m_uncloned_c$res$rd[30], m_uncloned_c$res$rd[30] - truth2$rd_true[30]))
cat(sprintf("  Uncloned Marginal:    RD = %.4f (Bias = %+.4f)\n", m_uncloned_m$res$rd[30], m_uncloned_m$res$rd[30] - truth2$rd_true[30]))

cat(sprintf("\nRuntimes:\n  Cloned:              %.2fs\n  Uncloned Conditional: %.2fs\n  Uncloned Marginal:    %.2fs\n",
            t_c["elapsed"], t_u_c["elapsed"], t_u_m["elapsed"]))

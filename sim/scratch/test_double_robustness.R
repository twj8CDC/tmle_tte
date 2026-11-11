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
source("sim/dgp/01_dgp_structural_baseline.R")
source("sim/true_values.R")

# Load run_uncloned_tmle from our previous script
source("sim/scratch/benchmark_uncloned_vs_cloned.R")

cat("=== Testing Double Robustness across 10 Iterations (Scenario 1) ===\n")

truth <- compute_scenario1_truth()
rd15_true <- truth$rd_true[15]
rd30_true <- truth$rd_true[30]

seeds <- 101:110

results <- list()

for (s in seeds) {
  dat <- generate_scenario1_data(n = 1000, seed = s)
  
  # Regimes:
  # R1: Corr Q ("x1"), Corr g ("x1")
  # R2: Mis Q ("I(x1 > 0)"), Corr g ("x1")
  # R3: Corr Q ("x1"), Mis g ("I(x1 > 0)")
  # R4: Mis Q ("I(x1 > 0)"), Mis g ("I(x1 > 0)")
  
  regimes <- list(
    R1 = list(q = "x1", p = "x1"),
    R2 = list(q = "I(x1 > 0)", p = "x1"),
    R3 = list(q = "x1", p = "I(x1 > 0)"),
    R4 = list(q = "I(x1 > 0)", p = "I(x1 > 0)")
  )
  
  for (rname in names(regimes)) {
    reg <- regimes[[rname]]
    
    # 1. Cloned TTE-TMLE
    fit_cloned <- tte_tmle(
      data = dat,
      covars = "x1",
      p_form = reg$p,
      q_form = reg$q,
      use_sl = FALSE
    )
    res_c <- fit_cloned$res_tmle
    rd15_c <- res_c$rd[15]
    rd30_c <- res_c$rd[30]
    
    # 2. Uncloned TMLE (Marginal Analytical)
    fit_uncloned_marg <- run_uncloned_tmle(
      data = dat,
      covars = "x1",
      p_form = reg$p,
      q_form = reg$q,
      method = "marginal",
      integration = "analytical"
    )
    rd15_u_m <- fit_uncloned_marg$res$rd[15]
    rd30_u_m <- fit_uncloned_marg$res$rd[30]

    # 3. Uncloned TMLE (Conditional Analytical)
    fit_uncloned_cond <- run_uncloned_tmle(
      data = dat,
      covars = "x1",
      p_form = reg$p,
      q_form = reg$q,
      method = "conditional",
      integration = "analytical"
    )
    rd15_u_c <- fit_uncloned_cond$res$rd[15]
    rd30_u_c <- fit_uncloned_cond$res$rd[30]

    results[[length(results) + 1]] <- tibble(
      seed = s,
      regime = rname,
      rd15_cloned = rd15_c,
      rd30_cloned = rd30_c,
      rd15_uncloned_marg = rd15_u_m,
      rd30_uncloned_marg = rd30_u_m,
      rd15_uncloned_cond = rd15_u_c,
      rd30_uncloned_cond = rd30_u_c
    )
  }
}

df_res <- bind_rows(results)

summary_table <- df_res %>%
  group_by(regime) %>%
  summarise(
    # RD 15
    bias15_cloned = mean(rd15_cloned) - rd15_true,
    bias15_uncloned_marg = mean(rd15_uncloned_marg) - rd15_true,
    bias15_uncloned_cond = mean(rd15_uncloned_cond) - rd15_true,
    rmse15_cloned = sqrt(mean((rd15_cloned - rd15_true)^2)),
    rmse15_uncloned_marg = sqrt(mean((rd15_uncloned_marg - rd15_true)^2)),
    rmse15_uncloned_cond = sqrt(mean((rd15_uncloned_cond - rd15_true)^2)),
    # RD 30
    bias30_cloned = mean(rd30_cloned) - rd30_true,
    bias30_uncloned_marg = mean(rd30_uncloned_marg) - rd30_true,
    bias30_uncloned_cond = mean(rd30_uncloned_cond) - rd30_true,
    rmse30_cloned = sqrt(mean((rd30_cloned - rd30_true)^2)),
    rmse30_uncloned_marg = sqrt(mean((rd30_uncloned_marg - rd30_true)^2)),
    rmse30_uncloned_cond = sqrt(mean((rd30_uncloned_cond - rd30_true)^2)),
    .groups = "drop"
  )

print(as.data.frame(summary_table), digits = 4)

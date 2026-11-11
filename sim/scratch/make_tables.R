library(dplyr)
sum_df <- readRDS("sim/comprehensive_summary_metrics.rds")

dir.create("sim/scratch", recursive=TRUE, showWarnings=FALSE)
sink("sim/scratch/scenario_tables.md")

for (s in 1:8) {
  scen_name <- paste0("Scenario ", s)
  sub <- sum_df %>% filter(scenario_id == scen_name)
  
  order_key <- c(
    "TMLE + BART_Non-Parametric",
    "G-COMP (BART)_Non-Parametric",
    "TTE-TMLE (SuperLearner)_Non-Parametric",
    "TTE-TMLE (GLM)_Regime 1 (Corr Q + Corr g)",
    "TTE-TMLE (GLM)_Regime 2 (Mis Q + Corr g)",
    "TTE-TMLE (GLM)_Regime 3 (Corr Q + Mis g)",
    "TTE-TMLE (GLM)_Regime 4 (Mis Q + Mis g)",
    "CCW-KM_Parametric (Correct)",
    "CCW-KM_Parametric (Misspecified)",
    "CCW-CPHM_Parametric (Correct)",
    "CCW-CPHM_Parametric (Misspecified)",
    "Time-Matched Cox_Parametric (Correct)",
    "Time-Matched Cox_Parametric (Misspecified)",
    "Time-Dependent Cox_Parametric (Correct)",
    "Time-Dependent Cox_Parametric (Misspecified)",
    "Landmark Cox_Parametric (Correct)",
    "Landmark Cox_Parametric (Misspecified)",
    "Naive Cox_Parametric (Correct)",
    "Naive Cox_Parametric (Misspecified)"
  )
  
  sub <- sub %>%
    mutate(match_key = paste0(estimator, "_", ifelse(estimator == "TTE-TMLE (GLM)", regime_id, spec_group))) %>%
    arrange(match(match_key, order_key))
  
  cat(sprintf("\n### Scenario %d Table Data\n\n", s))
  cat("| Estimator | Regime / Spec | Bias(15) | Bias(30) | RMSE(30) | 95% Cov(30) | CI Width | Band Cov | RMSTD Bias | RMSTD Cov | Mean Time |\n")
  cat("| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |\n")
  for (i in seq_len(nrow(sub))) {
    r <- sub[i, ]
    band_cov_str <- if (is.nan(r$band_cov) || is.na(r$band_cov)) "—" else sprintf("%.1f%%", r$band_cov * 100)
    rmstd_bias_str <- if (is.nan(r$mean_rmstd_bias) || is.na(r$mean_rmstd_bias)) "—" else sprintf("%+.3f", r$mean_rmstd_bias)
    rmstd_cov_str <- if (is.nan(r$rmstd_cov) || is.na(r$rmstd_cov)) "—" else sprintf("%.1f%%", r$rmstd_cov * 100)
    
    grp_label <- if (r$estimator == "TTE-TMLE (GLM)") {
      gsub("Regime [1-4] \\((.*)\\)", "GLM: \\1", r$regime_id)
    } else {
      gsub("Parametric \\((.*)\\)", "\\1", r$spec_group)
    }
    
    cat(sprintf("| **%s** | %s | %+.4f | %+.4f | %.4f | %.1f%% | %.4f | %s | %s | %s | %.2fs |\n",
                r$estimator, grp_label, r$mean_bias_15, r$mean_bias_30, r$rmse_30,
                r$cov_30 * 100, r$ci_width_30, band_cov_str, rmstd_bias_str, rmstd_cov_str, r$mean_elapsed))
  }
}
sink()
cat("Tables written to sim/scratch/scenario_tables.md successfully!\n")

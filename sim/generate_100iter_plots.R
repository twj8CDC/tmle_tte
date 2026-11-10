#' Generate Publication-Quality Scenario Plots from 100-Iteration Benchmark
#'
#' Creates 8 vertically-stacked paired regime plots (Regime 1 on top, Regime 4 on bottom)
#' displaying all 9 requested performance metrics across all 9 estimators.
#' Groups estimators by methodological family (Cox = Warm Rust/Amber, CCW/TTE = Blues, BART = Greens).

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

generate_100iter_plots <- function(rds_path = "sim/simulation_results_100iter.rds",
                                   out_dir = "sim/figures",
                                   artifact_dir = "/Users/jacobtiegs/.gemini/antigravity/brain/d66d2fba-9f6c-4b16-a430-36f36bd7fc90") {
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  df <- readRDS(rds_path)
  
  # Grouped estimator order:
  # 1. Cox family (3)
  # 2. CCW / TTE family (4)
  # 3. BART family (2, G-COMP (BART) right next to TMLE + BART)
  est_all_levels <- c(
    "Naive Cox",
    "Landmark Cox",
    "Time-Dependent Cox",
    "Time-Matched Cox",
    "CCW-CPHM",
    "CCW-KM",
    "TTE-TMLE (GLM)",
    "TTE-TMLE (SuperLearner)",
    "G-COMP (BART)",
    "TMLE + BART"
  )
  est_levels <- intersect(est_all_levels, unique(df$estimator))
  
  # Family-grouped color palette:
  # - Cox models: Warm Amber / Terracotta / Deep Rust / Burgundy
  # - CCW / TTE models: Sky Blue -> Steel Blue -> Royal Blue -> Deep Navy
  # - BART models: Bright Emerald Green -> Deep Forest Green
  est_colors <- c(
    "Naive Cox"               = "#E9964A",
    "Landmark Cox"            = "#D95F02",
    "Time-Dependent Cox"      = "#9E2A2B",
    "Time-Matched Cox"        = "#6A040F",
    "CCW-CPHM"                = "#6BAED6",
    "CCW-KM"                  = "#3182BD",
    "TTE-TMLE (GLM)"          = "#08519C",
    "TTE-TMLE (SuperLearner)" = "#04254E",
    "G-COMP (BART)"           = "#41AB5D",
    "TMLE + BART"             = "#005A32"
  )[est_levels]
  
  # Family-grouped shapes:
  # - Cox models: Squares (15)
  # - CCW / TTE models: Circles (16)
  # - BART models: Triangles (17)
  est_shapes <- c(
    "Naive Cox"               = 15,
    "Landmark Cox"            = 15,
    "Time-Dependent Cox"      = 15,
    "Time-Matched Cox"        = 15,
    "CCW-CPHM"                = 16,
    "CCW-KM"                  = 16,
    "TTE-TMLE (GLM)"          = 16,
    "TTE-TMLE (SuperLearner)" = 16,
    "G-COMP (BART)"           = 17,
    "TMLE + BART"             = 17
  )[est_levels]
  
  # Aggregate across the 100 iterations per (scenario_id, regime_id, estimator)
  summary_df <- df %>%
    mutate(
      cov_30_flag = as.numeric(ci_lower_30 <= rd_30_true & rd_30_true <= ci_upper_30),
      rd30_ci_w   = ci_upper_30 - ci_lower_30,
      hr30_ci_w   = hr_ci_upper_30 - hr_ci_lower_30
    ) %>%
    group_by(scenario_id, regime_id, estimator) %>%
    summarise(
      `Bias(15)`                    = mean(bias_15, na.rm = TRUE),
      `Bias(30)`                    = mean(bias_30, na.rm = TRUE),
      `RMSTD Bias (Days)`           = mean(rmstd_bias, na.rm = TRUE),
      `Marginal HR Bias`            = mean(hr_marginal_bias, na.rm = TRUE),
      emp_se_30                     = sd(rd_30_est, na.rm = TRUE),
      mod_se_30                     = mean(se_30, na.rm = TRUE),
      `SE Ratio [RD(30)]`           = mean(se_30, na.rm = TRUE) / sd(rd_30_est, na.rm = TRUE),
      `RMSE(15)`                    = sqrt(mean(bias_15^2, na.rm = TRUE)),
      `RMSE(30)`                    = sqrt(mean(bias_30^2, na.rm = TRUE)),
      `RMSTD RMSE (Days)`           = sqrt(mean(rmstd_bias^2, na.rm = TRUE)),
      `Marginal HR RMSE`            = sqrt(mean(hr_marginal_bias^2, na.rm = TRUE)),
      `95% Cov [RD(30)] (%)`        = mean(cov_30_flag, na.rm = TRUE) * 100,
      `RD(30) 95% CI Width`         = mean(rd30_ci_w, na.rm = TRUE),
      `HR 95% Cov [HR(30)] (%)`     = mean(hr_covered_30, na.rm = TRUE) * 100,
      `HR(30) 95% CI Width`         = mean(hr30_ci_w, na.rm = TRUE),
      `Mean Pointwise Band Cov (%)` = if (all(is.na(prop_time_points_covered))) NA_real_ else mean(prop_time_points_covered, na.rm = TRUE) * 100,
      `Mean Band Width`             = if (all(is.na(mean_band_width))) NA_real_ else mean(mean_band_width, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      estimator = factor(estimator, levels = rev(est_levels))
    )
  
  metric_order <- c(
    "Bias(15)",
    "Bias(30)",
    "RMSTD Bias (Days)",
    "Marginal HR Bias",
    "SE Ratio [RD(30)]",
    "RMSE(15)",
    "RMSE(30)",
    "RMSTD RMSE (Days)",
    "Marginal HR RMSE",
    "95% Cov [RD(30)] (%)",
    "RD(30) 95% CI Width",
    "HR 95% Cov [HR(30)] (%)",
    "HR(30) 95% CI Width",
    "Mean Pointwise Band Cov (%)",
    "Mean Band Width"
  )
  
  # Reference target lines per metric
  ref_lines <- tibble(
    metric = factor(metric_order, levels = metric_order),
    ref_val = c(
      0, 0, 0, 0, 1.0,
      0, 0, 0, 0, 95.0,
      NA_real_, 95.0, NA_real_, 95.0, NA_real_
    )
  )
  
  scenario_titles <- c(
    "Scenario 1" = "Scenario 1: Structural Baseline (Well-Behaved Linear Confounding & Moderate Grace Uptake)",
    "Scenario 2" = "Scenario 2: Severe Non-Linear & Interaction Confounding (Polynomial + Threshold Step)",
    "Scenario 3" = "Scenario 3: Waning Treatment Efficacy & Non-Proportional Hazards (Day-15 Hazard Kink)",
    "Scenario 4" = "Scenario 4: Strong Prognostic Covariate (Efficiency & Variance Reduction Benchmark)",
    "Scenario 5" = "Scenario 5: Severe Positivity Strain & Near-Violations (Extreme Propensity Weights)",
    "Scenario 6" = "Scenario 6: Time-Varying Confounder-Treatment Feedback (Day-3 Biomarker Collider Trap)",
    "Scenario 7" = "Scenario 7: Positivity Tail Depletion & Empty-Cell Support Bias (Sparse Frailty Tail)",
    "Scenario 8" = "Scenario 8: Depletion of Susceptibles & Survivor Collider Bias (Unmeasured Gamma Frailty)"
  )
  
  build_regime_panel <- function(sub_long, reg_title, show_legend = FALSE) {
    ggplot(sub_long, aes(x = value, y = estimator, color = estimator, shape = estimator)) +
      geom_hline(yintercept = c(2.5, 6.5), linetype = "dotted", color = "#CCCCCC", linewidth = 0.45) +
      geom_vline(
        data = ref_lines %>% filter(!is.na(ref_val)),
        aes(xintercept = ref_val),
        linetype = "dashed", color = "#444444", linewidth = 0.5, alpha = 0.8
      ) +
      geom_segment(
        data = sub_long %>% filter(!is.na(value) & metric %in% c("Bias(15)", "Bias(30)", "RMSTD Bias (Days)", "Marginal HR Bias", "RMSE(15)", "RMSE(30)", "RMSTD RMSE (Days)", "Marginal HR RMSE")),
        aes(x = 0, xend = value, y = estimator, yend = estimator),
        linewidth = 0.7, alpha = 0.8, na.rm = TRUE
      ) +
      geom_point(size = 2.9, na.rm = TRUE) +
      facet_wrap(~ metric, scales = "free_x", nrow = 3, ncol = 5) +
      scale_color_manual(
        name = "Estimator Family Key",
        values = est_colors,
        breaks = est_levels,
        drop = FALSE
      ) +
      scale_shape_manual(
        name = "Estimator Family Key",
        values = est_shapes,
        breaks = est_levels,
        drop = FALSE
      ) +
      scale_y_discrete(drop = FALSE) +
      labs(
        subtitle = reg_title,
        x = NULL,
        y = NULL
      ) +
      theme_bw(base_size = 10.5) +
      theme(
        plot.subtitle = element_text(size = 12, face = "bold", color = "#1F4E78", margin = margin(b = 6)),
        strip.background = element_rect(fill = "#1F4E78", color = "#1F4E78"),
        strip.text = element_text(color = "white", face = "bold", size = 9.5),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(color = "#F2F2F2", linewidth = 0.35),
        axis.text.y = element_text(size = 8.5, face = "bold", color = "#222222"),
        axis.text.x = element_text(size = 8.5, color = "#333333"),
        legend.position = if (show_legend) "top" else "none",
        legend.title = element_text(face = "bold", size = 9.5),
        legend.text = element_text(size = 8.5),
        legend.background = element_rect(fill = "#F8F9FA", color = "#CCCCCC")
      ) +
      guides(
        color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 3.3)),
        shape = guide_legend(nrow = 2, byrow = TRUE)
      )
  }
  
  for (s_idx in 1:8) {
    s_name <- paste0("Scenario ", s_idx)
    s_df <- summary_df %>%
      filter(scenario_id == s_name) %>%
      select(-emp_se_30, -mod_se_30) %>%
      pivot_longer(
        cols = all_of(metric_order),
        names_to = "metric",
        values_to = "value"
      ) %>%
      mutate(metric = factor(metric, levels = metric_order))
    
    s_reg1 <- s_df %>% filter(grepl("Regime 1", regime_id))
    s_reg4 <- s_df %>% filter(grepl("Regime 4", regime_id))
    
    p1 <- build_regime_panel(s_reg1, "Regime 1: Both Correctly Specified (Corr Q + Corr g)", show_legend = TRUE)
    p4 <- build_regime_panel(s_reg4, "Regime 4: Both Models Misspecified (Mis Q + Mis g)", show_legend = FALSE)
    
    combined <- (p1 / p4) +
      plot_annotation(
        title = scenario_titles[s_name],
        subtitle = "100-Iteration Monte Carlo Benchmark (N = 1,000, tau = 30d, Grace = 10d) | Families: Cox (Amber/Rust/Burgundy), CCW/TTE (Blues), BART (Greens)",
        caption = "NA Disclaimer: Mean Band Width and Mean Pointwise Band Cov (%) are NA (unplotted) for Naive Cox, Landmark Cox, Time-Dependent Cox, Time-Matched Cox, and CCW-CPHM (no simultaneous risk-difference bands).",
        theme = theme(
          plot.title = element_text(size = 14, face = "bold", color = "#0F2C59"),
          plot.subtitle = element_text(size = 10, face = "italic", color = "#444444"),
          plot.caption = element_text(size = 9, face = "italic", color = "#666666", hjust = 0)
        )
      )
    
    fname <- sprintf("scenario_%d_benchmark_100iter.png", s_idx)
    out_file <- file.path(out_dir, fname)
    art_file <- file.path(artifact_dir, fname)
    
    ggsave(out_file, plot = combined, width = 17.5, height = 14.5, dpi = 220, bg = "white")
    file.copy(out_file, art_file, overwrite = TRUE)
    cat(sprintf("Saved %s and copied to artifact directory.\n", out_file))
  }
}

if (!interactive()) {
  generate_100iter_plots()
}

#' Generate Publication-Quality Scenario Plots from 100-Iteration Benchmark
#'
#' Layout per scenario: Unified single figure with a 4x4 Estimand Grid + 1x4 Diagnostic Row
#' - Columns 1-4 (Estimands): RD(15), RD(30), RMSTD(30), Marginal HR(30)
#' - Rows 1-4: Bias, RMSE, 95% Coverage (%), 95% CI Width
#' - Row 5 (Diagnostics): SE Ratio [RD(30)], Uniform Band Cov (%), Mean Pointwise Band Cov (%), Mean Band Width
#' - Y-axis Grouping: Non-Parametric (Top), Parametric [Corr] (Middle, Filled), Parametric [Mis] (Bottom, Open)
#' - Colors & Shapes: Non-Parametric (Green Triangles), TTE-TMLE GLM (Blue Triangles), CCW-KM (Blue Circles), Time-Matched Cox (Red Squares)

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

generate_100iter_plots <- function(rds_path = "sim/simulation_results_100iter.rds",
                                   out_dir = "sim/figures",
                                   artifact_dir = "/Users/jacobtiegs/.gemini/antigravity/brain/d66d2fba-9f6c-4b16-a430-36f36bd7fc90") {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  df <- readRDS(rds_path)

  # Plot only estimators that cleanly target baseline t=0 marginal survival trajectories
  df_plot <- df %>%
    filter(
      spec_group %in% c("Non-Parametric", "Parametric (Correct)", "Parametric (Misspecified)"),
      estimator %in% c(
        "TTE-TMLE (BART)", "G-COMP (BART)", "TTE-TMLE (SuperLearner)",
        "TTE-TMLE (GLM)", "CCW-KM", "Time-Matched Cox"
      )
    )

  # Create unified Y-axis factor
  df_plot <- df_plot %>%
    mutate(
      y_label = case_when(
        spec_group == "Non-Parametric" ~ estimator,
        spec_group == "Parametric (Correct)" ~ paste0(estimator, " [Corr]"),
        spec_group == "Parametric (Misspecified)" ~ paste0(estimator, " [Mis]"),
        TRUE ~ estimator
      ),
      group_cat = factor(
        spec_group,
        levels = c("Parametric (Misspecified)", "Parametric (Correct)", "Non-Parametric")
      )
    )

  y_levels_ordered <- c(
    # Parametric (Misspecified) - Bottom
    "Time-Matched Cox [Mis]",
    "CCW-KM [Mis]",
    "TTE-TMLE (GLM) [Mis]",
    # Parametric (Correct) - Middle
    "Time-Matched Cox [Corr]",
    "CCW-KM [Corr]",
    "TTE-TMLE (GLM) [Corr]",
    # Non-Parametric - Top
    "TTE-TMLE (SuperLearner)",
    "G-COMP (BART)",
    "TTE-TMLE (BART)"
  )

  # Aggregate summary metrics across iterations per scenario and y_label
  summary_df <- df_plot %>%
    group_by(scenario_id, y_label, group_cat) %>%
    summarise(
      # RD(15)
      `RD(15) Bias` = mean(bias_15, na.rm = TRUE),
      `RD(15) RMSE` = sqrt(mean(bias_15^2, na.rm = TRUE)),
      `RD(15) 95% Cov (%)` = mean(covered_15, na.rm = TRUE) * 100,
      `RD(15) 95% CI Width` = mean(ci_width_15, na.rm = TRUE),
      # RD(30)
      `RD(30) Bias` = mean(bias_30, na.rm = TRUE),
      `RD(30) RMSE` = sqrt(mean(bias_30^2, na.rm = TRUE)),
      `RD(30) 95% Cov (%)` = mean(covered_30, na.rm = TRUE) * 100,
      `RD(30) 95% CI Width` = mean(ci_width_30, na.rm = TRUE),
      # RMSTD(30)
      `RMSTD Bias (Days)` = if (all(is.na(rmstd_bias))) NA_real_ else mean(rmstd_bias, na.rm = TRUE),
      `RMSTD RMSE (Days)` = if (all(is.na(rmstd_bias))) NA_real_ else sqrt(mean(rmstd_bias^2, na.rm = TRUE)),
      `RMSTD 95% Cov (%)` = if (all(is.na(rmstd_covered))) NA_real_ else mean(rmstd_covered, na.rm = TRUE) * 100,
      `RMSTD 95% CI Width` = if (all(is.na(rmstd_ci_width))) NA_real_ else mean(rmstd_ci_width, na.rm = TRUE),
      # Marginal HR(30)
      `Marginal HR Bias` = mean(hr_marginal_bias, na.rm = TRUE),
      `Marginal HR RMSE` = sqrt(mean(hr_marginal_bias^2, na.rm = TRUE)),
      `Marginal HR 95% Cov (%)` = mean(hr_covered_30, na.rm = TRUE) * 100,
      `Marginal HR 95% CI Width` = mean(hr_ci_width_30, na.rm = TRUE),
      # Diagnostics
      emp_se_30 = sd(rd_30_est, na.rm = TRUE),
      mod_se_30 = mean(se_30, na.rm = TRUE),
      `SE Ratio [RD(30)]` = mean(se_30, na.rm = TRUE) / pmax(sd(rd_30_est, na.rm = TRUE), 1e-6),
      `Uniform Band Cov (%)` = if (all(is.na(uniform_band_covered))) NA_real_ else mean(uniform_band_covered, na.rm = TRUE) * 100,
      `Mean Pointwise Band Cov (%)` = if (all(is.na(prop_time_points_covered))) NA_real_ else mean(prop_time_points_covered, na.rm = TRUE) * 100,
      `Mean Band Width` = if (all(is.na(mean_band_width))) NA_real_ else mean(mean_band_width, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      y_label = factor(y_label, levels = y_levels_ordered)
    )

  # Scenario descriptive titles
  scenario_titles <- c(
    "Scenario 1" = "Scenario 1: Structural Baseline (Linear Confounding & Grace Period Uptake)",
    "Scenario 2" = "Scenario 2: Severe Non-Linear Confounding (Exponential + Quadratic + Threshold Gate)",
    "Scenario 3" = "Scenario 3: Non-Proportional Hazards & Day-15 Treatment Efficacy Waning",
    "Scenario 4" = "Scenario 4: Strong Prognostic Covariate (Variance Reduction & Efficiency)",
    "Scenario 5" = "Scenario 5: Severe Positivity Strain (Near-Violations in Frail Subgroups)",
    "Scenario 6" = "Scenario 6: Time-Varying Confounder Feedback (Day-3 Biomarker Collider Trap)",
    "Scenario 7" = "Scenario 7: Positivity Tail Depletion & Empty-Cell Support Bias (Sparse Tail)",
    "Scenario 8" = "Scenario 8: Depletion of Susceptibles & Survivor Collider Bias (Unmeasured Frailty)"
  )

  # Color, fill, and shape definitions per y_label
  # Colors: Non-Parametric = Green (#2CA02C), CCW-KM & TTE-TMLE (GLM) = Blue (#1F77B4), Time-Matched Cox = Red (#D62728)
  # Shapes (pch 21=circle, 22=square, 24=triangle):
  #   Triangles (24): 3 Non-Parametrics + TTE-TMLE (GLM)
  #   Circles (21):   CCW-KM
  #   Squares (22):   Time-Matched Cox
  # Fill: Filled with group color for Non-Parametric and [Corr]; "white" (open) for [Mis]
  label_colors <- c(
    "Time-Matched Cox [Mis]"     = "#D62728",
    "CCW-KM [Mis]"               = "#1F77B4",
    "TTE-TMLE (GLM) [Mis]"       = "#1F77B4",
    "Time-Matched Cox [Corr]"    = "#D62728",
    "CCW-KM [Corr]"              = "#1F77B4",
    "TTE-TMLE (GLM) [Corr]"      = "#1F77B4",
    "TTE-TMLE (SuperLearner)"    = "#2CA02C",
    "G-COMP (BART)"              = "#2CA02C",
    "TTE-TMLE (BART)"            = "#2CA02C"
  )

  label_fills <- c(
    "Time-Matched Cox [Mis]"     = "white",
    "CCW-KM [Mis]"               = "white",
    "TTE-TMLE (GLM) [Mis]"       = "white",
    "Time-Matched Cox [Corr]"    = "#D62728",
    "CCW-KM [Corr]"              = "#1F77B4",
    "TTE-TMLE (GLM) [Corr]"      = "#1F77B4",
    "TTE-TMLE (SuperLearner)"    = "#2CA02C",
    "G-COMP (BART)"              = "#2CA02C",
    "TTE-TMLE (BART)"            = "#2CA02C"
  )

  label_shapes <- c(
    "Time-Matched Cox [Mis]"     = 22,
    "CCW-KM [Mis]"               = 21,
    "TTE-TMLE (GLM) [Mis]"       = 24,
    "Time-Matched Cox [Corr]"    = 22,
    "CCW-KM [Corr]"              = 21,
    "TTE-TMLE (GLM) [Corr]"      = 24,
    "TTE-TMLE (SuperLearner)"    = 24,
    "G-COMP (BART)"              = 24,
    "TTE-TMLE (BART)"            = 24
  )

  make_cell_plot <- function(data, metric_col, title_str, ref_val = NA_real_, is_leftmost = FALSE) {
    p <- ggplot(data, aes(x = .data[[metric_col]], y = y_label, color = y_label, fill = y_label, shape = y_label)) +
      geom_hline(yintercept = c(3.5, 6.5), linetype = "dashed", color = "#BBBBBB", linewidth = 0.5)

    if (!is.na(ref_val)) {
      p <- p + geom_vline(xintercept = ref_val, linetype = "dotted", color = "#444444", linewidth = 0.6)
    }

    if (grepl("Bias", metric_col)) {
      p <- p + geom_segment(
        aes(x = 0, xend = .data[[metric_col]], y = y_label, yend = y_label),
        linewidth = 0.6, alpha = 0.7, na.rm = TRUE
      )
    }

    p <- p +
      geom_point(size = 2.4, stroke = 0.85, na.rm = TRUE) +
      scale_color_manual(values = label_colors, drop = FALSE) +
      scale_fill_manual(values = label_fills, drop = FALSE) +
      scale_shape_manual(values = label_shapes, drop = FALSE) +
      scale_y_discrete(drop = FALSE) +
      labs(title = title_str, x = NULL, y = NULL) +
      theme_bw(base_size = 9) +
      theme(
        plot.title = element_text(size = 9.5, face = "bold", color = "#1F4E78", hjust = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(color = "#F0F0F0", linewidth = 0.3),
        axis.text.x = element_text(size = 8, color = "#333333"),
        axis.text.y = if (is_leftmost) element_text(size = 8, face = "bold", color = "#222222") else element_blank(),
        axis.ticks.y = if (is_leftmost) element_line(color = "#CCCCCC") else element_blank(),
        legend.position = "none"
      )
    p
  }

  for (s_idx in 1:8) {
    s_name <- paste0("Scenario ", s_idx)
    s_sub <- summary_df %>% filter(scenario_id == s_name)

    # ROW 1: Bias
    p11 <- make_cell_plot(s_sub, "RD(15) Bias", "Bias: RD(15)", ref_val = 0, is_leftmost = TRUE)
    p12 <- make_cell_plot(s_sub, "RD(30) Bias", "Bias: RD(30)", ref_val = 0, is_leftmost = FALSE)
    p13 <- make_cell_plot(s_sub, "RMSTD Bias (Days)", "Bias: RMSTD(30) (Days)", ref_val = 0, is_leftmost = FALSE)
    p14 <- make_cell_plot(s_sub, "Marginal HR Bias", "Bias: Marginal HR(30)", ref_val = 0, is_leftmost = FALSE)

    # ROW 2: RMSE
    p21 <- make_cell_plot(s_sub, "RD(15) RMSE", "RMSE: RD(15)", ref_val = 0, is_leftmost = TRUE)
    p22 <- make_cell_plot(s_sub, "RD(30) RMSE", "RMSE: RD(30)", ref_val = 0, is_leftmost = FALSE)
    p23 <- make_cell_plot(s_sub, "RMSTD RMSE (Days)", "RMSE: RMSTD(30) (Days)", ref_val = 0, is_leftmost = FALSE)
    p24 <- make_cell_plot(s_sub, "Marginal HR RMSE", "RMSE: Marginal HR(30)", ref_val = 0, is_leftmost = FALSE)

    # ROW 3: 95% Coverage (%)
    p31 <- make_cell_plot(s_sub, "RD(15) 95% Cov (%)", "95% Cov (%): RD(15)", ref_val = 95, is_leftmost = TRUE)
    p32 <- make_cell_plot(s_sub, "RD(30) 95% Cov (%)", "95% Cov (%): RD(30)", ref_val = 95, is_leftmost = FALSE)
    p33 <- make_cell_plot(s_sub, "RMSTD 95% Cov (%)", "95% Cov (%): RMSTD(30)", ref_val = 95, is_leftmost = FALSE)
    p34 <- make_cell_plot(s_sub, "Marginal HR 95% Cov (%)", "95% Cov (%): Marginal HR(30)", ref_val = 95, is_leftmost = FALSE)

    # ROW 4: 95% CI Width
    p41 <- make_cell_plot(s_sub, "RD(15) 95% CI Width", "CI Width: RD(15)", ref_val = NA, is_leftmost = TRUE)
    p42 <- make_cell_plot(s_sub, "RD(30) 95% CI Width", "CI Width: RD(30)", ref_val = NA, is_leftmost = FALSE)
    p43 <- make_cell_plot(s_sub, "RMSTD 95% CI Width", "CI Width: RMSTD(30)", ref_val = NA, is_leftmost = FALSE)
    p44 <- make_cell_plot(s_sub, "Marginal HR 95% CI Width", "CI Width: Marginal HR(30)", ref_val = NA, is_leftmost = FALSE)

    # ROW 5: Diagnostics
    p51 <- make_cell_plot(s_sub, "SE Ratio [RD(30)]", "SE Ratio [RD(30)]", ref_val = 1.0, is_leftmost = TRUE)
    p52 <- make_cell_plot(s_sub, "Uniform Band Cov (%)", "Uniform Band Cov (%)", ref_val = 95, is_leftmost = FALSE)
    p53 <- make_cell_plot(s_sub, "Mean Pointwise Band Cov (%)", "Pointwise Band Cov (%)", ref_val = 95, is_leftmost = FALSE)
    p54 <- make_cell_plot(s_sub, "Mean Band Width", "Mean Band Width", ref_val = NA, is_leftmost = FALSE)

    r1 <- p11 + p12 + p13 + p14 + plot_layout(ncol = 4, widths = c(1.4, 1, 1, 1))
    r2 <- p21 + p22 + p23 + p24 + plot_layout(ncol = 4, widths = c(1.4, 1, 1, 1))
    r3 <- p31 + p32 + p33 + p34 + plot_layout(ncol = 4, widths = c(1.4, 1, 1, 1))
    r4 <- p41 + p42 + p43 + p44 + plot_layout(ncol = 4, widths = c(1.4, 1, 1, 1))
    r5 <- p51 + p52 + p53 + p54 + plot_layout(ncol = 4, widths = c(1.4, 1, 1, 1))

    combined <- (r1 / r2 / r3 / r4 / r5) +
      plot_layout(heights = c(1, 1, 1, 1, 1)) +
      plot_annotation(
        title = scenario_titles[s_name],
        subtitle = "100-Iteration Monte Carlo Benchmark (N = 1,000, tau = 30d, Grace = 10d) | 4x4 Estimand Grid (Rows 1-4) + 1x4 Diagnostics (Row 5)\nColors: Non-Parametric (Green), TTE-TMLE (GLM) & CCW-KM (Blue), Time-Matched Cox (Red) | Shapes: TMLE/Non-Param (Triangle), CCW-KM (Circle), Matched Cox (Square) | Fill: [Corr] Filled, [Mis] Open",
        caption = "Notes: Time-Matched Cox uses untreated-hazard standardization, empirical switch integration, and patient-level bootstrap (200 resamples); Cox models do not construct simultaneous bands (NA).\nLandmark Cox and Time-Dependent Cox estimate conditional HRs only (see visual report Section 5) and are omitted from baseline t=0 marginal survival plots.",
        theme = theme(
          plot.title = element_text(size = 13, face = "bold", color = "#0F2C59"),
          plot.subtitle = element_text(size = 9.5, face = "italic", color = "#444444"),
          plot.caption = element_text(size = 8.5, face = "italic", color = "#666666", hjust = 0)
        )
      )

    fname <- sprintf("scenario_%d_benchmark_100iter.png", s_idx)
    out_file <- file.path(out_dir, fname)
    art_file <- file.path(artifact_dir, fname)

    ggsave(out_file, plot = combined, width = 16, height = 18, dpi = 200, bg = "white")
    if (dir.exists(artifact_dir)) {
      file.copy(out_file, art_file, overwrite = TRUE)
    }
    cat(sprintf("Saved %s and copied to artifact directory.\n", out_file))
  }
}

if (sys.nframe() == 0L && !interactive()) {
  generate_100iter_plots()
}

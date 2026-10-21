#' Target Trial Emulation via Time-to-Event TMLE (TTE-TMLE)
#'
#' @description End-to-end implementation of Target Trial Emulation combined with Survival
#'   Targeted Minimum Loss-Based Estimation (TMLE). Addresses immortal time bias across dynamic
#'   treatment grace periods, incorporates machine learning via \code{SuperLearner}, calculates
#'   patient-clustered Empirical Influence Curve (EIC) standard errors, and constructs
#'   simultaneous confidence bands via the Multiplier Bootstrap.
#'
#' @param data A data.frame containing observational survival data. Must contain variables for
#'   individual ID, observed event time, event indicator, treatment indicator, and treatment initiation time.
#' @param grace_period Numeric grace period duration for treatment initiation strategy (default 10).
#' @param tau Follow-up time horizon (default 30).
#' @param covars Character vector of baseline covariate column names (default \code{c("x1")}).
#' @param p_form Optional formula string for compliance/propensity models.
#' @param q_form Optional formula string for outcome hazard models.
#' @param dt Discrete time interval step size (default 1).
#' @param use_sl Logical indicating whether to use \code{SuperLearner} for nuisance models.
#' @param sl_lib Character vector of SuperLearner learners (default \code{c("SL.mean", "SL.glm", "SL.glm.time.int", "SL.earth.int", "SL.ranger")}).
#' @param trim_quantile Numeric vector of length 2 specifying quantiles for IPCW trimming (default \code{c(0.01, 0.99)}).
#' @param n_boot Number of multiplier bootstrap iterations for simultaneous confidence bands (default 10,000).
#' @param alpha Significance level for intervals and bands (default 0.05).
#' @param pure_gcomp Logical. If \code{TRUE}, evaluates pure counterfactual G-computation without the TMLE targeting update.
#' @param seed Optional random seed.
#' @param id_col Column name for individual ID (default \code{"id"}).
#' @param time_col Column name for follow-up time (default \code{"event_time"}).
#' @param event_col Column name for event indicator (default \code{"event"}).
#' @param tx_col Column name for treatment ever indicator (default \code{"tx_obs"}).
#' @param tx_time_col Column name for treatment initiation time (default \code{"tx_time_obs"}).
#'
#' @return An object of class \code{"tte_tmle"} containing:
#'   \describe{
#'     \item{res_tmle}{Tibble containing marginal survival for both arms, risk difference, pointwise CIs, and simultaneous confidence bands.}
#'     \item{res_log}{Tibble containing standard Inverse Probability Weighted Kaplan-Meier estimates for comparison.}
#'     \item{fit_ccw}{Survfit object for the weighted CCW Kaplan-Meier.}
#'     \item{crit_val_simul}{Multiplier bootstrap critical value for simultaneous bands.}
#'     \item{eps}{Estimated targeting coefficients.}
#'     \item{data}{Processed longitudinal person-time data.}
#'   }
#' @export
#' @importFrom dplyr mutate select left_join
#' @importFrom tibble as_tibble
#' @importFrom stats qnorm
tte_tmle <- function(data,
                     grace_period = 10,
                     tau = 30,
                     covars = c("x1"),
                     p_form = NULL,
                     q_form = NULL,
                     dt = 1,
                     use_sl = FALSE,
                     sl_lib = c("SL.mean", "SL.glm", "SL.glm.time.int", "SL.earth.int", "SL.ranger"),
                     trim_quantile = c(0.01, 0.99),
                     n_boot = 10000,
                     alpha = 0.05,
                     pure_gcomp = FALSE,
                     seed = NULL,
                     id_col = "id",
                     time_col = "event_time",
                     event_col = "event",
                     tx_col = "tx_obs",
                     tx_time_col = "tx_time_obs") {
  
  # Standardize column mapping
  data_clean <- prepare_tte_data(
    data, 
    id_col = id_col, 
    time_col = time_col, 
    event_col = event_col, 
    tx_col = tx_col, 
    tx_time_col = tx_time_col
  )
  
  # Step 1: Clone Creation & Protocol Censoring
  clones <- create_clones(
    data = data_clean, 
    covars = covars, 
    grace_period = grace_period, 
    tau = tau
  )
  
  # Step 2: Fit IPCW Models & Apply Weights
  ipcw_res <- compute_ipcw_weights(
    clones = clones,
    covars = covars,
    grace_period = grace_period,
    dt = dt,
    p_form = p_form,
    use_sl = use_sl,
    sl_lib = sl_lib,
    trim_quantile = trim_quantile
  )
  
  panel_w <- ipcw_res$panel_w
  
  # Step 3: Fit Stratified Initial Hazard Models (Q-models)
  q_res <- fit_q_models(
    panel_w = panel_w,
    covars = covars,
    q_form = q_form,
    use_sl = use_sl,
    sl_lib = sl_lib
  )
  
  # Step 4: TMLE Fluctuation / Targeting Step
  if (pure_gcomp) {
    panel_w <- panel_w %>%
      dplyr::mutate(
        Q_init_hazard = bound_precision(q_res$pred_Q_init),
        Q_init_logit = stats::qlogis(Q_init_hazard),
        H_ctrl = dplyr::if_else(arm == "control_no_tx", w_trim, 0),
        H_treat_pre = dplyr::if_else(arm == "treat_tx" & t <= grace_period, w_trim, 0),
        H_treat_post = dplyr::if_else(arm == "treat_tx" & t > grace_period, w_trim, 0)
      )
    eps_ctrl <- 0
    eps_treat_pre <- 0
    eps_treat_post <- 0
    targ_res <- list(panel_w = panel_w, eps_ctrl = 0, eps_treat_pre = 0, eps_treat_post = 0)
  } else {
    targ_res <- target_hazard(
      panel_w = panel_w,
      pred_Q_init = q_res$pred_Q_init,
      grace_period = grace_period
    )
    eps_ctrl <- targ_res$eps_ctrl
    eps_treat_pre <- targ_res$eps_treat_pre
    eps_treat_post <- targ_res$eps_treat_post
  }
  
  # Step 5: Counterfactual G-Computation Marginalization
  gcomp_res <- marginalize_gcomp(
    data = data_clean,
    clones = clones,
    m_Q_ctrl = q_res$m_Q_ctrl,
    m_Q_treat = q_res$m_Q_treat,
    m_denom_ctrl = ipcw_res$m_denom_ctrl,
    m_denom_treat = ipcw_res$m_denom_treat,
    m_num_ctrl = ipcw_res$m_num_ctrl,
    m_num_treat = ipcw_res$m_num_treat,
    eps_ctrl = eps_ctrl,
    eps_treat_pre = eps_treat_pre,
    eps_treat_post = eps_treat_post,
    grace_period = grace_period,
    tau = tau,
    covars = covars,
    q_lim = ipcw_res$q_lim,
    use_sl = use_sl,
    do_targeting = !pure_gcomp
  )
  
  # Step 6: Empirical Influence Curve & Multiplier Bootstrap Inference
  inf_res <- compute_eic_inference(
    panel_w = targ_res$panel_w,
    cf_data = gcomp_res$cf_data,
    tmle_surv = gcomp_res$tmle_surv,
    res_tmle = gcomp_res$res_tmle,
    eps_ctrl = eps_ctrl,
    eps_treat_pre = eps_treat_pre,
    eps_treat_post = eps_treat_post,
    n_boot = n_boot,
    alpha = alpha,
    seed = seed
  )
  
  # Step 7: Standard CCW Weighted KM Comparison
  sum_log <- summary(ipcw_res$fit_ccw, times = 1:tau, extend = TRUE)
  res_log <- tibble::tibble(
    t = sum_log$time, 
    tx = sum_log$strata, 
    sv = sum_log$surv,
    se = sum_log$std.err
  ) %>%
    dplyr::mutate(tx = gsub("^arm=", "", tx)) %>% 
    tidyr::pivot_wider(names_from = tx, values_from = c(sv, se)) %>% 
    dplyr::mutate(
      rd = sv_treat_tx - sv_control_no_tx,
      rd_se = sqrt(se_control_no_tx^2 + se_treat_tx^2),
      ci_lower = rd - stats::qnorm(1 - alpha / 2) * rd_se,
      ci_upper = rd + stats::qnorm(1 - alpha / 2) * rd_se
    )
  
  out <- list(
    res_tmle = inf_res$res_tmle,
    res_log = res_log,
    fit_ccw = ipcw_res$fit_ccw,
    crit_val_simul = inf_res$crit_val_simul,
    crit_val_pointwise = inf_res$crit_val_pointwise,
    eps = c(eps_ctrl = eps_ctrl, eps_treat_pre = eps_treat_pre, eps_treat_post = eps_treat_post),
    eic_summary = inf_res$eic_summary,
    ic_matrix = inf_res$ic_matrix,
    hr_summary = inf_res$hr_summary,
    data = panel_w,
    panel_w = panel_w,
    params = list(
      grace_period = grace_period,
      tau = tau,
      alpha = alpha,
      use_sl = use_sl,
      pure_gcomp = pure_gcomp
    )
  )
  
  class(out) <- "tte_tmle"
  return(out)
}

#' Backwards compatibility alias for run_ccw_tmle_analysis
#' @export
run_ccw_tmle_analysis <- tte_tmle

#' Print method for tte_tmle
#' @export
print.tte_tmle <- function(x, ...) {
  cat("Target Trial Emulation via Time-to-Event TMLE (TTE-TMLE)\n")
  cat("--------------------------------------------------------\n")
  cat("Follow-up Horizon (tau):", x$params$tau, "\n")
  cat("Grace Period:", x$params$grace_period, "\n")
  cat("Targeting Mode:", if (x$params$pure_gcomp) "Pure G-computation (Untargeted)" else "TMLE (Targeted)", "\n")
  cat("Multiplier Bootstrap Critical Value (95% Band):", round(x$crit_val_simul, 3), "\n\n")
  
  times_to_show <- c(min(x$params$grace_period, x$params$tau), 
                     round(x$params$tau / 2), 
                     x$params$tau)
  times_to_show <- unique(pmin(times_to_show, x$params$tau))
  
  cat("Primary Estimand: Marginal Risk Difference (RD = Treat - Control)\n")
  sub_rd <- x$res_tmle %>% 
    dplyr::filter(t %in% times_to_show) %>%
    dplyr::select(t, control_no_tx, treat_tx, rd, se_RD, rd_ci_lower, rd_ci_upper, rd_band_lower, rd_band_upper)
  print(as.data.frame(sub_rd), digits = 4, row.names = FALSE)
  
  cat("\nSecondary Derivative: Cumulative Hazard Ratio (Delta Method)\n")
  sub_hr <- x$res_tmle %>%
    dplyr::filter(t %in% times_to_show) %>%
    dplyr::select(t, hr_cum, hr_ci_lower, hr_ci_upper, hr_band_lower, hr_band_upper, hr_daily)
  print(as.data.frame(sub_hr), digits = 4, row.names = FALSE)
  
  invisible(x)
}

#' Plot method for tte_tmle
#'
#' @description Produces discrete-time step visualizations for survival curves,
#'   marginal risk differences with simultaneous confidence bands, or cumulative hazard ratios.
#'   Defaults to modern \code{ggplot2} step plots, with fallback to base R graphics.
#'
#' @param x An object of class \code{"tte_tmle"}.
#' @param which Plot type: \code{"survival"}, \code{"rd"} (risk difference), or \code{"hr"} (hazard ratio).
#' @param use_ggplot Logical. If \code{TRUE} (default), renders a publication-ready \code{ggplot2} step plot.
#'   If \code{FALSE}, renders a base R graphic.
#' @param step Logical. If \code{TRUE} (default), renders step functions matching discrete time intervals.
#'   If \code{FALSE}, renders continuous linear interpolation.
#' @param ... Additional arguments passed to base \code{plot} or \code{matplot} when \code{use_ggplot = FALSE}.
#' @return A \code{ggplot} object if \code{use_ggplot = TRUE}, or invisible \code{NULL} if base R plotting is used.
#' @export
plot.tte_tmle <- function(x, which = c("survival", "rd", "hr"), use_ggplot = TRUE, step = TRUE, ...) {
  which <- match.arg(which)
  
  if (use_ggplot && requireNamespace("ggplot2", quietly = TRUE)) {
    return(plot_tte_tmle_ggplot(x, which = which, step = step))
  }
  
  plot_type <- if (step) "s" else "l"
  
  if (which == "survival") {
    matplot(x$res_tmle$t, cbind(x$res_tmle$treat_tx, x$res_tmle$control_no_tx),
            type = plot_type, lty = c(1, 2), lwd = 2, col = c("blue", "red"),
            xlab = "Time (days)", ylab = "Marginal Survival Probability",
            main = "TTE-TMLE Marginal Survival Curves", ...)
    legend("bottomleft", legend = c("Treat Arm (Initiate <= Grace)", "Control Arm (Never Treat)"),
           col = c("blue", "red"), lty = c(1, 2), lwd = 2)
    grid()
  } else if (which == "rd") {
    plot(x$res_tmle$t, x$res_tmle$rd, type = plot_type, lwd = 2, col = "purple",
         ylim = range(c(x$res_tmle$rd_band_lower, x$res_tmle$rd_band_upper), na.rm = TRUE),
         xlab = "Time (days)", ylab = "Marginal Risk Difference (Treat - Control)",
         main = "TTE-TMLE Marginal Risk Difference & Simultaneous Bands", ...)
    lines(x$res_tmle$t, x$res_tmle$rd_ci_lower, type = plot_type, lty = 2, col = "darkgray")
    lines(x$res_tmle$t, x$res_tmle$rd_ci_upper, type = plot_type, lty = 2, col = "darkgray")
    lines(x$res_tmle$t, x$res_tmle$rd_band_lower, type = plot_type, lty = 3, col = "black", lwd = 1.5)
    lines(x$res_tmle$t, x$res_tmle$rd_band_upper, type = plot_type, lty = 3, col = "black", lwd = 1.5)
    abline(h = 0, lty = 2, col = "red")
    legend("topright", legend = c("Risk Difference", "95% Pointwise CI", "95% Simultaneous Band"),
           col = c("purple", "darkgray", "black"), lty = c(1, 2, 3), lwd = c(2, 1, 1.5))
    grid()
  } else if (which == "hr") {
    y_vals <- c(x$res_tmle$hr_band_lower, x$res_tmle$hr_band_upper, x$res_tmle$hr_cum)
    y_vals <- y_vals[is.finite(y_vals) & y_vals > 0]
    y_lim <- if (length(y_vals) > 0) range(pmin(pmax(y_vals, 0.05), 20)) else c(0.1, 5)
    
    plot(x$res_tmle$t, x$res_tmle$hr_cum, type = plot_type, lwd = 2, col = "darkgreen", log = "y",
         ylim = y_lim,
         xlab = "Time (days)", ylab = "Cumulative Hazard Ratio (Log Scale)",
         main = "TTE-TMLE Cumulative Hazard Ratio & Simultaneous Bands", ...)
    lines(x$res_tmle$t, x$res_tmle$hr_ci_lower, type = plot_type, lty = 2, col = "darkgray")
    lines(x$res_tmle$t, x$res_tmle$hr_ci_upper, type = plot_type, lty = 2, col = "darkgray")
    lines(x$res_tmle$t, x$res_tmle$hr_band_lower, type = plot_type, lty = 3, col = "black", lwd = 1.5)
    lines(x$res_tmle$t, x$res_tmle$hr_band_upper, type = plot_type, lty = 3, col = "black", lwd = 1.5)
    abline(h = 1, lty = 2, col = "red")
    legend("topright", legend = c("Cumulative HR", "95% Pointwise CI", "95% Simultaneous Band", "Null HR = 1"),
           col = c("darkgreen", "darkgray", "black", "red"), lty = c(1, 2, 3, 2), lwd = c(2, 1, 1.5, 1))
    grid()
  }
}

#' Autoplot method for tte_tmle
#'
#' @description Produces a \code{ggplot2} step visualization for \code{tte_tmle} objects.
#'
#' @param object An object of class \code{"tte_tmle"}.
#' @param which Plot type: \code{"survival"}, \code{"rd"}, or \code{"hr"}.
#' @param step Logical. If \code{TRUE} (default), renders step functions matching discrete time intervals.
#' @param ... Additional arguments (currently unused).
#' @return A \code{ggplot} object.
#' @importFrom ggplot2 autoplot
#' @exportS3Method ggplot2::autoplot
autoplot.tte_tmle <- function(object, which = c("survival", "rd", "hr"), step = TRUE, ...) {
  which <- match.arg(which)
  plot_tte_tmle_ggplot(object, which = which, step = step)
}

#' Plot TTE-TMLE Results using ggplot2
#'
#' @description Generates publication-ready ggplot2 step visualizations for survival curves,
#'   marginal risk differences with simultaneous confidence bands, or cumulative hazard ratios.
#'
#' @param x An object of class \code{"tte_tmle"}.
#' @param which Character string specifying the plot type: \code{"survival"}, \code{"rd"}, or \code{"hr"}.
#' @param step Logical. If \code{TRUE} (default), uses native \code{geom_step()} for discrete-time steps.
#'   If \code{FALSE}, uses \code{geom_line()} for linear interpolation.
#'
#' @return A \code{ggplot} object.
#' @export
plot_tte_tmle_ggplot <- function(x, which = c("survival", "rd", "hr"), step = TRUE) {
  which <- match.arg(which)
  res <- x$res_tmle
  gp <- x$params$grace_period
  tau <- x$params$tau
  step_dir <- "hv"
  
  if (which == "survival") {
    df0 <- data.frame(t = 0, treat_tx = 1.0, control_no_tx = 1.0)
    df <- rbind(df0, as.data.frame(res[, c("t", "treat_tx", "control_no_tx")]))
    df_long <- tidyr::pivot_longer(df, cols = c("treat_tx", "control_no_tx"),
                                   names_to = "arm", values_to = "surv")
    df_long$arm <- factor(df_long$arm,
                          levels = c("treat_tx", "control_no_tx"),
                          labels = c("Treatment (Initiate <= Grace)", "Control (Never Treat)"))
    
    p <- ggplot2::ggplot(df_long, ggplot2::aes(x = t, y = surv, color = arm, linetype = arm))
    if (step) {
      p <- p + ggplot2::geom_step(direction = step_dir, linewidth = 0.7)
    } else {
      p <- p + ggplot2::geom_line(linewidth = 0.7)
    }
    p <- p +
      ggplot2::geom_vline(xintercept = gp, linetype = "dashed", color = "gray50", alpha = 0.7) +
      ggplot2::scale_color_manual(values = c("Treatment (Initiate <= Grace)" = "#1f77b4", "Control (Never Treat)" = "#d62728")) +
      ggplot2::scale_linetype_manual(values = c("Treatment (Initiate <= Grace)" = "solid", "Control (Never Treat)" = "dashed")) +
      ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
      ggplot2::scale_x_continuous(breaks = seq(0, tau, by = if (tau <= 30) 5 else 10)) +
      ggplot2::labs(
        title = "TTE-TMLE Marginal Counterfactual Survival Curves",
        subtitle = paste0("Protocol Grace Period = ", gp, " days | Discrete-time step functions"),
        x = "Follow-up Time (days)",
        y = "Marginal Survival Probability S(t)",
        color = "Regime",
        linetype = "Regime"
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        legend.position = "bottom",
        panel.grid.minor = ggplot2::element_blank()
      )
    return(p)
    
  } else if (which == "rd") {
    df0 <- data.frame(t = 0, rd = 0, rd_ci_lower = 0, rd_ci_upper = 0,
                      rd_band_lower = 0, rd_band_upper = 0)
    cols <- c("t", "rd", "rd_ci_lower", "rd_ci_upper", "rd_band_lower", "rd_band_upper")
    df <- rbind(df0, as.data.frame(res[, cols]))
    
    p <- ggplot2::ggplot(df, ggplot2::aes(x = t)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "firebrick", alpha = 0.8) +
      ggplot2::geom_vline(xintercept = gp, linetype = "dotted", color = "gray50", alpha = 0.8)
    
    if (step) {
      p <- p +
        ggplot2::geom_step(ggplot2::aes(y = rd_band_upper, linetype = "95% Simultaneous Band"), direction = step_dir, color = "gray30", linewidth = 0.5) +
        ggplot2::geom_step(ggplot2::aes(y = rd_band_lower, linetype = "95% Simultaneous Band"), direction = step_dir, color = "gray30", linewidth = 0.5) +
        ggplot2::geom_step(ggplot2::aes(y = rd_ci_upper, linetype = "95% Pointwise CI"), direction = step_dir, color = "#1f77b4", linewidth = 0.5) +
        ggplot2::geom_step(ggplot2::aes(y = rd_ci_lower, linetype = "95% Pointwise CI"), direction = step_dir, color = "#1f77b4", linewidth = 0.5) +
        ggplot2::geom_step(ggplot2::aes(y = rd, color = "Risk Difference (TMLE)"), direction = step_dir, linewidth = 0.7)
    } else {
      p <- p +
        ggplot2::geom_line(ggplot2::aes(y = rd_band_upper, linetype = "95% Simultaneous Band"), color = "gray30", linewidth = 0.5) +
        ggplot2::geom_line(ggplot2::aes(y = rd_band_lower, linetype = "95% Simultaneous Band"), color = "gray30", linewidth = 0.5) +
        ggplot2::geom_line(ggplot2::aes(y = rd_ci_upper, linetype = "95% Pointwise CI"), color = "#1f77b4", linewidth = 0.5) +
        ggplot2::geom_line(ggplot2::aes(y = rd_ci_lower, linetype = "95% Pointwise CI"), color = "#1f77b4", linewidth = 0.5) +
        ggplot2::geom_line(ggplot2::aes(y = rd, color = "Risk Difference (TMLE)"), linewidth = 0.7)
    }
    
    p <- p +
      ggplot2::scale_color_manual(name = "", values = c("Risk Difference (TMLE)" = "#2ca02c")) +
      ggplot2::scale_linetype_manual(name = "Uncertainty Limits",
                                    values = c("95% Simultaneous Band" = "dotted", "95% Pointwise CI" = "dashed")) +
      ggplot2::scale_x_continuous(breaks = seq(0, tau, by = if (tau <= 30) 5 else 10)) +
      ggplot2::labs(
        title = "TTE-TMLE Marginal Risk Difference & Simultaneous Confidence Bands",
        subtitle = paste0("Protocol Grace Period = ", gp, " days | Simultaneous multiplier critical value = ", round(x$crit_val_simul, 3)),
        x = "Follow-up Time (days)",
        y = "Risk Difference RD(t) = S1(t) - S0(t)"
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        legend.position = "bottom",
        legend.box = "horizontal",
        panel.grid.minor = ggplot2::element_blank()
      )
    return(p)
    
  } else if (which == "hr") {
    cols <- c("t", "hr_cum", "hr_ci_lower", "hr_ci_upper", "hr_band_lower", "hr_band_upper")
    df <- as.data.frame(res[, cols])
    df <- df[is.finite(df$hr_cum) & df$hr_cum > 0, ]
    
    p <- ggplot2::ggplot(df, ggplot2::aes(x = t)) +
      ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "firebrick", alpha = 0.8) +
      ggplot2::geom_vline(xintercept = gp, linetype = "dotted", color = "gray50", alpha = 0.8)
    
    if (step) {
      p <- p +
        ggplot2::geom_step(ggplot2::aes(y = hr_band_upper, linetype = "95% Simultaneous Band"), direction = step_dir, color = "gray30", linewidth = 0.5) +
        ggplot2::geom_step(ggplot2::aes(y = hr_band_lower, linetype = "95% Simultaneous Band"), direction = step_dir, color = "gray30", linewidth = 0.5) +
        ggplot2::geom_step(ggplot2::aes(y = hr_ci_upper, linetype = "95% Pointwise CI"), direction = step_dir, color = "#1f77b4", linewidth = 0.5) +
        ggplot2::geom_step(ggplot2::aes(y = hr_ci_lower, linetype = "95% Pointwise CI"), direction = step_dir, color = "#1f77b4", linewidth = 0.5) +
        ggplot2::geom_step(ggplot2::aes(y = hr_cum, color = "Cumulative Hazard Ratio (TMLE)"), direction = step_dir, linewidth = 0.7)
    } else {
      p <- p +
        ggplot2::geom_line(ggplot2::aes(y = hr_band_upper, linetype = "95% Simultaneous Band"), color = "gray30", linewidth = 0.5) +
        ggplot2::geom_line(ggplot2::aes(y = hr_band_lower, linetype = "95% Simultaneous Band"), color = "gray30", linewidth = 0.5) +
        ggplot2::geom_line(ggplot2::aes(y = hr_ci_upper, linetype = "95% Pointwise CI"), color = "#1f77b4", linewidth = 0.5) +
        ggplot2::geom_line(ggplot2::aes(y = hr_ci_lower, linetype = "95% Pointwise CI"), color = "#1f77b4", linewidth = 0.5) +
        ggplot2::geom_line(ggplot2::aes(y = hr_cum, color = "Cumulative Hazard Ratio (TMLE)"), linewidth = 0.7)
    }
    
    p <- p +
      ggplot2::scale_y_log10() +
      ggplot2::scale_color_manual(name = "", values = c("Cumulative Hazard Ratio (TMLE)" = "#2ca02c")) +
      ggplot2::scale_linetype_manual(name = "Uncertainty Limits",
                                    values = c("95% Simultaneous Band" = "dotted", "95% Pointwise CI" = "dashed")) +
      ggplot2::scale_x_continuous(breaks = seq(0, tau, by = if (tau <= 30) 5 else 10)) +
      ggplot2::labs(
        title = "TTE-TMLE Cumulative Hazard Ratio & Simultaneous Confidence Bands",
        subtitle = paste0("Protocol Grace Period = ", gp, " days | Delta-method log-cumulative hazard ratio"),
        x = "Follow-up Time (days)",
        y = "Cumulative Hazard Ratio (Log Scale)"
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        legend.position = "bottom",
        legend.box = "horizontal",
        panel.grid.minor = ggplot2::element_blank()
      )
    return(p)
  }
}

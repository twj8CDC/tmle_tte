#' Modular Single-Dataset Simulation Runner Across All Estimators and Specification Regimes
#'
#' @description Generates a single Monte Carlo dataset for a given scenario (1 to 8)
#'   and evaluates:
#'   1. Non-Parametric Estimators ONCE per dataset:
#'      - G-COMP (BART) [Python bart-survival v0.2.5]
#'      - TMLE + BART (BART for both g and Q) [Python bart-survival v0.2.5 + R TMLE/EIC]
#'      - TTE-TMLE (SuperLearner)
#'   2. Parametric Benchmarks under Correct and Misspecified functional forms (never omitting variables):
#'      - Naive Cox (Corr Q vs Mis Q)
#'      - Landmark Cox (Corr Q vs Mis Q; RMSTD = NA)
#'      - Time-Dependent Cox (Corr Q vs Mis Q)
#'      - Time-Matched Cox (Corr Q + g vs Mis Q + g)
#'      - CCW-CPHM (Corr Q + g vs Mis Q + g)
#'      - CCW-KM (Corr g vs Mis g)
#'   3. TTE-TMLE (GLM) across all 4 regimes:
#'      - Regime 1 (Corr Q + Corr g)
#'      - Regime 2 (Mis Q + Corr g)
#'      - Regime 3 (Corr Q + Mis g)
#'      - Regime 4 (Mis Q + Mis g)

library(dplyr)
library(tibble)
library(survival)

# Source core library
source("R/00_utils.R")
source("R/01_cloning.R")
source("R/02_ipcw_models.R")
source("R/03_q_models.R")
source("R/04_targeting.R")
source("R/05_gcomp.R")
source("R/06_inference.R")
source("R/07_tte_tmle.R")
source("R/08_benchmarks.R")
source("R/09_bart_gcomp.R")

# Source DGPs and Truth Engine
source("sim/dgp/01_dgp_structural_baseline.R")
source("sim/dgp/02_dgp_misspecification.R")
source("sim/dgp/03_dgp_waning_nph.R")
source("sim/dgp/04_dgp_prognostic_efficiency.R")
source("sim/dgp/05_dgp_positivity_strain.R")
source("sim/dgp/06_dgp_time_varying_feedback.R")
source("sim/dgp/07_dgp_positivity_tail_depletion.R")
source("sim/dgp/08_dgp_depletion_of_susceptibles.R")
source("sim/true_values.R")

get_scenario_spec <- function(scenario_num, n = 1000, tau = 30, grace_period = 10, seed = 42, truth_precomp = NULL) {
  nl_w12 <- "I(exp(0.7 * (w1 - 3) / 1.5)) + I(((w1 - 3) / 1.5)^2) + w2 + I(w1 > 3.8 & w2 == 1)"
  lin_w12 <- "w1 + w2"

  if (scenario_num == 1) {
    df_raw <- generate_scenario1_data(n = n, tau = tau, grace_period = grace_period, seed = seed)
    truth <- if (!is.null(truth_precomp)) truth_precomp else compute_scenario1_truth(tau = tau, grace_period = grace_period)
    covars_all <- c("x1")
    p_form_corr <- "x1"
    p_form_mis  <- "I(x1 > 0)"
    q_form_corr <- "x1"
    q_form_mis  <- "I(x1 > 0)"
  } else if (scenario_num == 2) {
    df_raw <- generate_scenario2_data(n = n, tau = tau, grace_period = grace_period, seed = seed)
    truth <- if (!is.null(truth_precomp)) truth_precomp else compute_scenario2_truth(tau = tau, grace_period = grace_period)
    covars_all <- c("w1", "w2")
    p_form_corr <- nl_w12
    p_form_mis  <- lin_w12
    q_form_corr <- nl_w12
    q_form_mis  <- lin_w12
  } else if (scenario_num == 3) {
    df_raw <- generate_scenario3_data(n = n, tau = tau, grace_period = grace_period, seed = seed)
    truth <- if (!is.null(truth_precomp)) truth_precomp else compute_scenario3_truth(tau = tau, grace_period = grace_period)
    covars_all <- c("w1", "w2")
    p_form_corr <- nl_w12
    p_form_mis  <- lin_w12
    q_form_corr <- nl_w12
    q_form_mis  <- lin_w12
  } else if (scenario_num == 4) {
    df_raw <- generate_scenario4_data(n = n, tau = tau, grace_period = grace_period, seed = seed)
    truth <- if (!is.null(truth_precomp)) truth_precomp else compute_scenario4_truth(tau = tau, grace_period = grace_period)
    covars_all <- c("w1", "w2", "w_prog")
    p_form_corr <- paste(nl_w12, "+ w_prog")
    p_form_mis  <- "w1 + w2 + w_prog"
    q_form_corr <- paste(nl_w12, "+ I(exp(0.55 * w_prog)) + I(w_prog^2)")
    q_form_mis  <- "w1 + w2 + w_prog"
  } else if (scenario_num == 5) {
    df_raw <- generate_scenario5_data(n = n, tau = tau, grace_period = grace_period, seed = seed)
    truth <- if (!is.null(truth_precomp)) truth_precomp else compute_scenario5_truth(tau = tau, grace_period = grace_period)
    covars_all <- c("w1", "w2")
    p_form_corr <- nl_w12
    p_form_mis  <- lin_w12
    q_form_corr <- nl_w12
    q_form_mis  <- lin_w12
  } else if (scenario_num == 6) {
    df_raw <- generate_scenario6_data(n = n, tau = tau, grace_period = grace_period, seed = seed)
    truth <- if (!is.null(truth_precomp)) truth_precomp else compute_scenario6_truth(tau = tau, grace_period = grace_period)
    covars_all <- c("w1", "w2")
    p_form_corr <- nl_w12
    p_form_mis  <- lin_w12
    q_form_corr <- nl_w12
    q_form_mis  <- lin_w12
  } else if (scenario_num == 7) {
    df_raw <- generate_scenario7_data(n = n, tau = tau, grace_period = grace_period, seed = seed)
    truth <- if (!is.null(truth_precomp)) truth_precomp else compute_scenario7_truth(tau = tau, grace_period = grace_period)
    covars_all <- c("w1", "w2")
    p_form_corr <- nl_w12
    p_form_mis  <- lin_w12
    q_form_corr <- nl_w12
    q_form_mis  <- lin_w12
  } else if (scenario_num == 8) {
    df_raw <- generate_scenario8_data(n = n, tau = tau, grace_period = grace_period, seed = seed)
    truth <- if (!is.null(truth_precomp)) truth_precomp else compute_scenario8_truth(tau = tau, grace_period = grace_period)
    covars_all <- c("w1", "w2")
    p_form_corr <- nl_w12
    p_form_mis  <- lin_w12
    q_form_corr <- nl_w12
    q_form_mis  <- lin_w12
  } else {
    stop("Unsupported scenario_num: ", scenario_num)
  }

  list(
    df_raw = df_raw,
    truth = truth,
    covars_all = covars_all,
    p_form_corr = p_form_corr,
    p_form_mis = p_form_mis,
    q_form_corr = q_form_corr,
    q_form_mis = q_form_mis
  )
}

run_simulation_dataset <- function(scenario_num = 1,
                                   n = 1000,
                                   tau = 30,
                                   grace_period = 10,
                                   seed = 42,
                                   bart_trees = 50,
                                   bart_ndpost = 250,
                                   bart_nskip = 250,
                                   bart_thin = 4,
                                   bart_chains = 4L,
                                   fast_only = FALSE,
                                   cox_only = FALSE,
                                   truth_precomp = NULL) {
  if (!is.null(seed)) set.seed(seed)

  spec <- get_scenario_spec(
    scenario_num = scenario_num,
    n = n,
    tau = tau,
    grace_period = grace_period,
    seed = seed,
    truth_precomp = truth_precomp
  )
  df_raw <- spec$df_raw
  truth <- spec$truth
  covars_all <- spec$covars_all

  df_std <- prepare_tte_data(df_raw)
  clones <- if (!cox_only) create_clones(df_std, covars = covars_all, grace_period = grace_period, tau = tau) else NULL

  # True scalar estimands
  rd_true_15 <- truth$rd_true[min(15, length(truth$rd_true))]
  rd_true_30 <- truth$rd_true[min(30, length(truth$rd_true))]
  rmstd_true_30 <- sum(truth$rd_true[1:tau], na.rm = TRUE)
  hr_marg_true_30 <- if ("hr_true" %in% names(truth) && length(truth$hr_true) >= tau) {
    truth$hr_true[tau]
  } else {
    log(pmax(pmin(truth$S1_true[tau], 1 - 1e-6), 1e-6)) / log(pmax(pmin(truth$S0_true[tau], 1 - 1e-6), 1e-6))
  }

  # Hall-Wellner / Nair simultaneous band start day: >= 5 active treated events
  ev_trt_days <- sort(ceiling(df_std$time[df_std$tx_ever == 1L & df_std$t_switch <= grace_period & df_std$event == 1L & df_std$time >= df_std$t_switch]))
  t_band_start <- if (length(ev_trt_days) >= 5L) max(1L, as.integer(ev_trt_days[5L])) else as.integer(grace_period)

  eval_band_cov <- function(band_low, band_upp, true_curve, t_start = t_band_start) {
    T_len <- min(length(band_low), length(true_curve))
    idx <- seq(from = min(t_start, T_len), to = T_len)
    is_cov <- (band_low[idx] <= true_curve[idx]) & (true_curve[idx] <= band_upp[idx])
    list(
      unif = as.integer(all(is_cov, na.rm = TRUE)),
      prop = mean(is_cov, na.rm = TRUE),
      width = mean(band_upp[idx] - band_low[idx], na.rm = TRUE)
    )
  }

  make_row <- function(estimator,
                       spec_group,
                       regime_id,
                       elapsed_sec,
                       rd15 = NA_real_, se15 = NA_real_,
                       ci15_l = if (is.na(rd15) || is.na(se15)) NA_real_ else rd15 - 1.96 * se15,
                       ci15_u = if (is.na(rd15) || is.na(se15)) NA_real_ else rd15 + 1.96 * se15,
                       rd30 = NA_real_, se30 = NA_real_,
                       ci30_l = if (is.na(rd30) || is.na(se30)) NA_real_ else rd30 - 1.96 * se30,
                       ci30_u = if (is.na(rd30) || is.na(se30)) NA_real_ else rd30 + 1.96 * se30,
                       band30_l = NA_real_, band30_u = NA_real_,
                       rmstd = NA_real_, se_rmstd = NA_real_,
                       rmstd_l = if (is.na(rmstd) || is.na(se_rmstd)) NA_real_ else rmstd - 1.96 * se_rmstd,
                       rmstd_u = if (is.na(rmstd) || is.na(se_rmstd)) NA_real_ else rmstd + 1.96 * se_rmstd,
                       hr = NA_real_, hr_l = NA_real_, hr_u = NA_real_,
                       hr_cond = NA_real_, hr_cond_l = NA_real_, hr_cond_u = NA_real_,
                       band_cov = list(unif = NA_integer_, prop = NA_real_, width = NA_real_)) {
    z30 <- if (is.na(rd30) || is.na(se30)) NA_real_ else rd30 / pmax(se30, 1e-4)
    p30 <- if (is.na(z30)) NA_real_ else 2 * (1 - stats::pnorm(abs(z30)))
    tibble(
      estimator = estimator,
      spec_group = spec_group,
      regime_id = regime_id,
      elapsed_sec = elapsed_sec,
      # RD(15)
      rd_15_est = rd15,
      rd_15_true = if (is.na(rd15)) NA_real_ else rd_true_15,
      bias_15 = if (is.na(rd15)) NA_real_ else rd15 - rd_true_15,
      se_15 = se15,
      ci_lower_15 = ci15_l,
      ci_upper_15 = ci15_u,
      covered_15 = if (is.na(rd15) || is.na(ci15_l) || is.na(ci15_u)) NA_integer_ else as.integer(ci15_l <= rd_true_15 & rd_true_15 <= ci15_u),
      ci_width_15 = ci15_u - ci15_l,
      # RD(30)
      rd_30_est = rd30,
      rd_30_true = if (is.na(rd30)) NA_real_ else rd_true_30,
      bias_30 = if (is.na(rd30)) NA_real_ else rd30 - rd_true_30,
      rel_bias_30 = if (is.na(rd30)) NA_real_ else (rd30 - rd_true_30) / rd_true_30 * 100,
      se_30 = se30,
      ci_lower_30 = ci30_l,
      ci_upper_30 = ci30_u,
      covered_30 = if (is.na(rd30) || is.na(ci30_l) || is.na(ci30_u)) NA_integer_ else as.integer(ci30_l <= rd_true_30 & rd_true_30 <= ci30_u),
      ci_width_30 = ci30_u - ci30_l,
      band_lower_30 = band30_l,
      band_upper_30 = band30_u,
      # RMSTD(30)
      rmstd_est = rmstd,
      rmstd_true = if (is.na(rmstd)) NA_real_ else rmstd_true_30,
      rmstd_bias = if (is.na(rmstd)) NA_real_ else rmstd - rmstd_true_30,
      se_rmstd = se_rmstd,
      rmstd_ci_lower = rmstd_l,
      rmstd_ci_upper = rmstd_u,
      rmstd_covered = if (is.na(rmstd) || is.na(rmstd_l) || is.na(rmstd_u)) NA_integer_ else as.integer(rmstd_l <= rmstd_true_30 & rmstd_true_30 <= rmstd_u),
      rmstd_ci_width = rmstd_u - rmstd_l,
      # Marginal HR(30)
      hr_marginal_est = hr,
      hr_marginal_true = if (is.na(hr)) NA_real_ else hr_marg_true_30,
      hr_marginal_bias = if (is.na(hr)) NA_real_ else hr - hr_marg_true_30,
      hr_ci_lower_30 = hr_l,
      hr_ci_upper_30 = hr_u,
      hr_covered_30 = if (is.na(hr) || is.na(hr_l) || is.na(hr_u)) NA_integer_ else as.integer(hr_l <= hr_marg_true_30 & hr_marg_true_30 <= hr_u),
      hr_ci_width_30 = hr_u - hr_l,
      # Conditional Cox HR (where applicable; NA for marginal-only estimators)
      hr_cond_est = hr_cond,
      hr_cond_ci_lower = hr_cond_l,
      hr_cond_ci_upper = hr_cond_u,
      hr_cond_ci_width = hr_cond_u - hr_cond_l,
      # Hypothesis testing & Simultaneous band diagnostics
      z_stat_30 = z30,
      p_val_30 = p30,
      reject_H0 = if (is.na(p30)) NA_integer_ else as.integer(p30 < 0.05),
      uniform_band_covered = band_cov$unif,
      prop_time_points_covered = band_cov$prop,
      mean_band_width = band_cov$width
    )
  }

  rows <- list()

  # =========================================================================
  # PART A: Parametric Benchmarks under Correct & Misspecified Regimes
  # =========================================================================
  param_specs <- list(
    list(
      spec_group = "Parametric (Correct)",
      regime_id = "Regime 1 (Corr Q + Corr g)",
      p_form = spec$p_form_corr,
      q_form = spec$q_form_corr
    ),
    list(
      spec_group = "Parametric (Misspecified)",
      regime_id = "Regime 4 (Mis Q + Mis g)",
      p_form = spec$p_form_mis,
      q_form = spec$q_form_mis
    )
  )

  ipcw_cache <- list()

  for (ps in param_specs) {
    sg <- ps$spec_group
    rid <- ps$regime_id
    pf <- ps$p_form
    qf <- ps$q_form

    # 1. Landmark Cox (Conditional post-landmark HR only; baseline marginal RD/RMSTD/HR undefined)
    t0 <- proc.time()["elapsed"]
    fit_lm <- fit_landmark_cox(df_std, landmark_time = grace_period, covars = covars_all, tau = tau, q_form = qf)
    t_lm <- as.numeric(proc.time()["elapsed"] - t0)

    coef_lm <- stats::coef(fit_lm)
    hr_lm_cond <- if ("tx_landmark" %in% names(coef_lm)) unname(exp(coef_lm["tx_landmark"])) else NA_real_
    ci_hr_lm_cond <- if ("tx_landmark" %in% names(coef_lm)) {
      unname(exp(suppressWarnings(stats::confint(fit_lm))["tx_landmark", ]))
    } else {
      c(NA_real_, NA_real_)
    }

    rows[[length(rows) + 1L]] <- make_row(
      estimator = "Landmark Cox",
      spec_group = sg,
      regime_id = rid,
      elapsed_sec = t_lm,
      hr_cond = hr_lm_cond,
      hr_cond_l = ci_hr_lm_cond[1],
      hr_cond_u = ci_hr_lm_cond[2]
    )

    # 2. Time-Dependent Cox (Conditional instantaneous HR only; static grace-period RD/RMSTD/HR undefined)
    t0 <- proc.time()["elapsed"]
    fit_td <- fit_td_cox(df_std, covars = covars_all, tau = tau, grace_period = grace_period, q_form = qf)
    t_td <- as.numeric(proc.time()["elapsed"] - t0)

    rows[[length(rows) + 1L]] <- make_row(
      estimator = "Time-Dependent Cox",
      spec_group = sg,
      regime_id = rid,
      elapsed_sec = t_td,
      hr_cond = unname(fit_td$hr_cond),
      hr_cond_l = unname(fit_td$hr_cond_ci[1]),
      hr_cond_u = unname(fit_td$hr_cond_ci[2])
    )

    # 3. Time-Matched Cox (Baseline-standardized curves + patient-level nonparametric bootstrap)
    t0 <- proc.time()["elapsed"]
    fit_tm <- fit_time_matched_cox(
      data = df_std,
      covars = covars_all,
      p_form = pf,
      q_form = qf,
      tau = tau,
      grace_period = grace_period,
      n_boot = 200L,
      seed = seed
    )
    t_tm <- as.numeric(proc.time()["elapsed"] - t0)

    rows[[length(rows) + 1L]] <- make_row(
      estimator = "Time-Matched Cox",
      spec_group = sg,
      regime_id = rid,
      elapsed_sec = t_tm,
      rd15 = fit_tm$rd_15, se15 = fit_tm$se_rd_15,
      rd30 = fit_tm$rd_tau, se30 = fit_tm$se_rd_tau,
      rmstd = fit_tm$rmstd, se_rmstd = fit_tm$se_rmstd,
      hr = unname(fit_tm$hr), hr_l = unname(fit_tm$hr_ci[1]), hr_u = unname(fit_tm$hr_ci[2]),
      hr_cond = unname(fit_tm$hr_cond), hr_cond_l = unname(fit_tm$hr_cond_ci[1]), hr_cond_u = unname(fit_tm$hr_cond_ci[2])
    )

    if (!cox_only) {
      # Precompute IPCW weights for pf
      t0_ipcw <- proc.time()["elapsed"]
      ipcw_res <- compute_ipcw_weights(
        clones = clones,
        covars = covars_all,
        grace_period = grace_period,
        dt = 1,
        p_form = pf,
        use_sl = FALSE
      )
      t_ipcw <- as.numeric(proc.time()["elapsed"] - t0_ipcw)
      ipcw_cache[[sg]] <- ipcw_res
      panel_w <- ipcw_res$panel_w

      # 4. CCW-KM
      t0 <- proc.time()["elapsed"]
      fit_km <- fit_ccw_km(panel_w, tau = tau, n_boot = 2000, seed = seed)
      t_km <- t_ipcw + as.numeric(proc.time()["elapsed"] - t0)
      res_km <- fit_km$results
      band_cov_km <- eval_band_cov(res_km$rd_band_lower, res_km$rd_band_upper, truth$rd_true)

      rows[[length(rows) + 1L]] <- make_row(
        estimator = "CCW-KM",
        spec_group = sg,
        regime_id = rid,
        elapsed_sec = t_km,
        rd15 = res_km$rd[15], se15 = res_km$se_RD[15],
        ci15_l = res_km$rd_ci_lower[15], ci15_u = res_km$rd_ci_upper[15],
        rd30 = res_km$rd[tau], se30 = res_km$se_RD[tau],
        ci30_l = res_km$rd_ci_lower[tau], ci30_u = res_km$rd_ci_upper[tau],
        band30_l = res_km$rd_band_lower[tau], band30_u = res_km$rd_band_upper[tau],
        rmstd = fit_km$rmstd_tau, se_rmstd = fit_km$se_rmstd,
        hr = fit_km$hr_marginal_tau, hr_l = fit_km$hr_ci_lower_tau, hr_u = fit_km$hr_ci_upper_tau,
        band_cov = band_cov_km
      )
    }
  }

  # =========================================================================
  # PART B: TTE-TMLE (GLM) Across All 4 Regimes (Regimes 1, 2, 3, 4)
  # =========================================================================
  if (!cox_only) {
    glm_regimes <- list(
      list(r = 1L, regime_id = "Regime 1 (Corr Q + Corr g)", spec_group = "Parametric (Correct)",      pf = spec$p_form_corr, qf = spec$q_form_corr),
      list(r = 2L, regime_id = "Regime 2 (Mis Q + Corr g)",  spec_group = "Parametric (Mis Q + Corr g)", pf = spec$p_form_corr, qf = spec$q_form_mis),
      list(r = 3L, regime_id = "Regime 3 (Corr Q + Mis g)",  spec_group = "Parametric (Corr Q + Mis g)", pf = spec$p_form_mis,  qf = spec$q_form_corr),
      list(r = 4L, regime_id = "Regime 4 (Mis Q + Mis g)",   spec_group = "Parametric (Misspecified)", pf = spec$p_form_mis,  qf = spec$q_form_mis)
    )

    for (gr in glm_regimes) {
      t0 <- proc.time()["elapsed"]
      fit_glm <- tte_tmle(
        data = df_std,
        grace_period = grace_period,
        tau = tau,
        covars = covars_all,
        p_form = gr$pf,
        q_form = gr$qf,
        use_sl = FALSE,
        n_boot = 2000,
        seed = seed
      )
      t_glm <- as.numeric(proc.time()["elapsed"] - t0)
      res_glm <- fit_glm$res_tmle
      band_cov_glm <- eval_band_cov(res_glm$rd_band_lower, res_glm$rd_band_upper, truth$rd_true)

      rows[[length(rows) + 1L]] <- make_row(
        estimator = "TTE-TMLE (GLM)",
        spec_group = gr$spec_group,
        regime_id = gr$regime_id,
        elapsed_sec = t_glm,
        rd15 = res_glm$rd[15], se15 = res_glm$se_RD[15],
        ci15_l = res_glm$rd_ci_lower[15], ci15_u = res_glm$rd_ci_upper[15],
        rd30 = res_glm$rd[tau], se30 = res_glm$se_RD[tau],
        ci30_l = res_glm$rd_ci_lower[tau], ci30_u = res_glm$rd_ci_upper[tau],
        band30_l = res_glm$rd_band_lower[tau], band30_u = res_glm$rd_band_upper[tau],
        rmstd = sum(res_glm$rd[1:tau], na.rm = TRUE), se_rmstd = fit_glm$se_rmstd,
        hr = res_glm$hr_cum[tau], hr_l = res_glm$hr_ci_lower[tau], hr_u = res_glm$hr_ci_upper[tau],
        band_cov = band_cov_glm
      )
    }
  }

  # =========================================================================
  # PART C: Non-Parametric Estimators (Run ONCE per dataset!)
  # =========================================================================
  if (!fast_only && !cox_only) {
    # 7 & 8. G-COMP (BART) and TMLE + BART (BART for both g and Q via Python bart-survival v0.2.5)
    bart_pair <- fit_bart_pair(
      data = df_std,
      clones = clones,
      covars = covars_all,
      grace_period = grace_period,
      tau = tau,
      ntree = bart_trees,
      ndpost = bart_ndpost,
      nskip = bart_nskip,
      thin = bart_thin,
      chains = bart_chains,
      n_boot = 2000,
      seed = seed
    )

    # 7. G-COMP (BART)
    fit_bg <- bart_pair$gcomp
    res_bg <- fit_bg$results
    band_cov_bg <- eval_band_cov(res_bg$band_lower, res_bg$band_upper, truth$rd_true)

    rows[[length(rows) + 1L]] <- make_row(
      estimator = "G-COMP (BART)",
      spec_group = "Non-Parametric",
      regime_id = "Non-Parametric",
      elapsed_sec = fit_bg$elapsed_sec,
      rd15 = res_bg$rd[15], se15 = res_bg$se_RD[15],
      ci15_l = res_bg$cri_lower[15], ci15_u = res_bg$cri_upper[15],
      rd30 = res_bg$rd[tau], se30 = res_bg$se_RD[tau],
      ci30_l = res_bg$cri_lower[tau], ci30_u = res_bg$cri_upper[tau],
      band30_l = res_bg$band_lower[tau], band30_u = res_bg$band_upper[tau],
      rmstd = fit_bg$rmstd, se_rmstd = fit_bg$se_rmstd,
      rmstd_l = fit_bg$rmstd_ci_lower, rmstd_u = fit_bg$rmstd_ci_upper,
      hr = res_bg$hr_cum[tau], hr_l = res_bg$hr_ci_lower[tau], hr_u = res_bg$hr_ci_upper[tau],
      band_cov = band_cov_bg
    )

    # 8. TTE-TMLE (BART)
    fit_tb <- bart_pair$tmle
    res_tb <- fit_tb$results
    band_cov_tb <- eval_band_cov(res_tb$rd_band_lower, res_tb$rd_band_upper, truth$rd_true)

    rows[[length(rows) + 1L]] <- make_row(
      estimator = "TTE-TMLE (BART)",
      spec_group = "Non-Parametric",
      regime_id = "Non-Parametric",
      elapsed_sec = fit_tb$elapsed_sec,
      rd15 = res_tb$rd[15], se15 = res_tb$se_RD[15],
      ci15_l = res_tb$rd_ci_lower[15], ci15_u = res_tb$rd_ci_upper[15],
      rd30 = res_tb$rd[tau], se30 = res_tb$se_RD[tau],
      ci30_l = res_tb$rd_ci_lower[tau], ci30_u = res_tb$rd_ci_upper[tau],
      band30_l = res_tb$rd_band_lower[tau], band30_u = res_tb$rd_band_upper[tau],
      rmstd = sum(res_tb$rd[1:tau], na.rm = TRUE), se_rmstd = fit_tb$se_rmstd,
      hr = res_tb$hr_cum[tau], hr_l = res_tb$hr_ci_lower[tau], hr_u = res_tb$hr_ci_upper[tau],
      band_cov = band_cov_tb
    )

    # 9. TTE-TMLE (SuperLearner)
    t0 <- proc.time()["elapsed"]
    fit_sl <- tte_tmle(
      data = df_std,
      grace_period = grace_period,
      tau = tau,
      covars = covars_all,
      use_sl = TRUE,
      sl_lib = c("SL.mean", "SL.glm.time.int", "SL.hal.time", "SL.earth.int"),
      trim_quantile = c(0.001, 0.999),
      n_boot = 2000,
      seed = seed
    )
    t_sl <- as.numeric(proc.time()["elapsed"] - t0)
    res_sl <- fit_sl$res_tmle
    band_cov_sl <- eval_band_cov(res_sl$rd_band_lower, res_sl$rd_band_upper, truth$rd_true)

    rows[[length(rows) + 1L]] <- make_row(
      estimator = "TTE-TMLE (SuperLearner)",
      spec_group = "Non-Parametric",
      regime_id = "Non-Parametric",
      elapsed_sec = t_sl,
      rd15 = res_sl$rd[15], se15 = res_sl$se_RD[15],
      ci15_l = res_sl$rd_ci_lower[15], ci15_u = res_sl$rd_ci_upper[15],
      rd30 = res_sl$rd[tau], se30 = res_sl$se_RD[tau],
      ci30_l = res_sl$rd_ci_lower[tau], ci30_u = res_sl$rd_ci_upper[tau],
      band30_l = res_sl$rd_band_lower[tau], band30_u = res_sl$rd_band_upper[tau],
      rmstd = sum(res_sl$rd[1:tau], na.rm = TRUE), se_rmstd = fit_sl$se_rmstd,
      hr = res_sl$hr_cum[tau], hr_l = res_sl$hr_ci_lower[tau], hr_u = res_sl$hr_ci_upper[tau],
      band_cov = band_cov_sl
    )
  }

  out_df <- dplyr::bind_rows(rows) %>%
    dplyr::mutate(
      scenario_id = paste0("Scenario ", scenario_num),
      n = n,
      .before = 1
    )

  return(out_df)
}

# Backwards-compatible wrapper
run_simulation_condition <- function(scenario_num = 1,
                                     regime_num = NULL,
                                     n = 1000,
                                     tau = 30,
                                     grace_period = 10,
                                     seed = 42,
                                     mc_cores = 1,
                                     bart_trees = 50,
                                     bart_ndpost = 250,
                                     bart_nskip = 250,
                                     bart_thin = 4,
                                     bart_chains = 4L,
                                     fast_only = FALSE,
                                     truth_precomp = NULL) {
  res <- run_simulation_dataset(
    scenario_num = scenario_num,
    n = n,
    tau = tau,
    grace_period = grace_period,
    seed = seed,
    bart_trees = bart_trees,
    bart_ndpost = bart_ndpost,
    bart_nskip = bart_nskip,
    bart_thin = bart_thin,
    bart_chains = bart_chains,
    fast_only = fast_only,
    truth_precomp = truth_precomp
  )
  if (!is.null(regime_num)) {
    target_reg <- switch(
      as.character(regime_num),
      "1" = c("Regime 1 (Corr Q + Corr g)", "Non-Parametric"),
      "2" = c("Regime 2 (Mis Q + Corr g)", "Non-Parametric"),
      "3" = c("Regime 3 (Corr Q + Mis g)", "Non-Parametric"),
      "4" = c("Regime 4 (Mis Q + Mis g)", "Non-Parametric"),
      unique(res$regime_id)
    )
    res <- res %>% dplyr::filter(regime_id %in% target_reg)
  }
  res
}

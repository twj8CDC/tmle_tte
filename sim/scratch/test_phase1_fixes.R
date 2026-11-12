suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(splines)
  library(survival)
})

source("R/00_utils.R")
source("R/01_cloning.R")
source("R/02_ipcw_models.R")
source("R/03_q_models.R")
source("R/04_targeting.R")
source("R/05_gcomp.R")
source("R/06_inference.R")
source("R/08_benchmarks.R")
source("sim/dgp/01_dgp_structural_baseline.R")
source("sim/dgp/02_dgp_misspecification.R")
source("sim/dgp/04_dgp_prognostic_efficiency.R")
source("sim/dgp/05_dgp_positivity_strain.R")
source("sim/dgp/07_dgp_positivity_tail_depletion.R")
source("sim/true_values.R")

run_test_scen <- function(scen_num, seeds = 1:15, ctrl_mode = "ns3_int", align_treat = TRUE) {
  nl_w12 <- "I(exp(0.7 * (w1 - 3) / 1.5)) + I(((w1 - 3) / 1.5)^2) + w2 + I(w1 > 3.8 & w2 == 1)"
  if (scen_num == 1) {
    gen_fn <- generate_scenario1_data; truth <- compute_scenario1_truth()
    covs <- "x1"; pf <- "x1"; qf_corr <- "x1"; qf_mis <- "I(x1 > 0)"
  } else if (scen_num == 2) {
    gen_fn <- generate_scenario2_data; truth <- compute_scenario2_truth()
    covs <- c("w1", "w2"); pf <- nl_w12; qf_corr <- nl_w12; qf_mis <- "w1 + w2"
  } else if (scen_num == 4) {
    gen_fn <- generate_scenario4_data; truth <- compute_scenario4_truth()
    covs <- c("w1", "w2", "w_prog")
    pf <- paste(nl_w12, "+ w_prog")
    qf_corr <- paste(nl_w12, "+ I(exp(0.55 * w_prog)) + I(w_prog^2)")
    qf_mis <- "w1 + w2 + w_prog"
  } else if (scen_num == 5) {
    gen_fn <- generate_scenario5_data; truth <- compute_scenario5_truth()
    covs <- c("w1", "w2"); pf <- nl_w12; qf_corr <- nl_w12; qf_mis <- "w1 + w2"
  } else if (scen_num == 7) {
    gen_fn <- generate_scenario7_data; truth <- compute_scenario7_truth()
    covs <- c("w1", "w2"); pf <- nl_w12; qf_corr <- nl_w12; qf_mis <- "w1 + w2"
  }

  grace_period <- 10
  tau <- 30
  rows <- list()

  for (s in seeds) {
    dat <- gen_fn(n = 1000, seed = s)
    df <- prepare_tte_data(dat)
    clones <- create_clones(df, covars = covs, grace_period = grace_period, tau = tau)

    # 1. Fit obs_panel for m_denom_ctrl
    obs_panel <- clones %>%
      filter(clone_arm == 1L) %>%
      mutate(
        tmax = pmin(if_else(tx_ever == 1L, t_switch, time), time, grace_period),
        n_int = pmax(ceiling(tmax), 1L)
      ) %>%
      uncount(weights = n_int, .id = "t") %>%
      mutate(
        t_sq = t^2,
        event_tx = if_else(tx_ever == 1L & t == t_switch & time >= t_switch & t <= grace_period, 1L, 0L)
      )

    ctrl_rhs <- switch(
      ctrl_mode,
      "ns3_int" = paste0("splines::ns(t, df = 3) * (", pf, ")"),
      "ns4_lin_int" = paste0("splines::ns(t, df = 4) + (", pf, ") + t:(", pf, ")"),
      "quad_int" = paste0("(t + t_sq) * (", pf, ")"),
      "ns4_int" = paste0("splines::ns(t, df = 4) * (", pf, ")")
    )
    num_rhs <- switch(
      ctrl_mode,
      "ns3_int" = "splines::ns(t, df = 3)",
      "ns4_lin_int" = "splines::ns(t, df = 4)",
      "quad_int" = "t + t_sq",
      "ns4_int" = "splines::ns(t, df = 4)"
    )

    m_denom_ctrl <- suppressWarnings(glm(as.formula(paste("event_tx ~", ctrl_rhs)), family = binomial(), data = obs_panel))
    m_num_ctrl <- suppressWarnings(glm(as.formula(paste("event_tx ~", num_rhs)), family = binomial(), data = obs_panel))

    # Predict S_{0,A}(t | W) for t = 1..grace_period for each patient
    g_grid <- expand_grid(id = df$id, t = seq_len(grace_period)) %>%
      left_join(df %>% select(id, all_of(covs)), by = "id") %>%
      mutate(t_sq = t^2)

    g_grid$lam_denom <- bound_precision(suppressWarnings(predict(m_denom_ctrl, newdata = g_grid, type = "response")))
    g_grid$lam_num   <- bound_precision(suppressWarnings(predict(m_num_ctrl, newdata = g_grid, type = "response")))

    g_grid <- g_grid %>%
      arrange(id, t) %>%
      group_by(id) %>%
      mutate(
        S_A_denom = cumprod(1 - lam_denom),
        S_A_num   = cumprod(1 - lam_num)
      ) %>%
      ungroup()

    S_m_lookup <- g_grid %>%
      filter(t == grace_period) %>%
      select(id, S_m_denom = S_A_denom, S_m_num = S_A_num)

    # Also treat_day10_data for fallback / comparison
    treat_d10 <- clones %>% filter(arm == "treat_tx", time > grace_period) %>%
      mutate(took_tx = if_else(tx_ever == 1L & t_switch <= grace_period, 1L, 0L))
    m_denom_treat <- suppressWarnings(glm(as.formula(paste("took_tx ~", pf)), family = binomial(), data = treat_d10))
    m_num_treat <- suppressWarnings(glm(took_tx ~ 1, family = binomial(), data = treat_d10))

    # Expand clones to panel
    panel <- clones %>%
      mutate(n_int = pmax(ceiling(time_cc), 0L)) %>%
      uncount(weights = n_int, .id = "t") %>%
      mutate(
        t0 = t - 1L,
        t1 = t,
        t_sq = t^2,
        event_cc_day = if_else(event_cc == 1L & time_cc > t0 & time_cc <= t1, 1L, 0L),
        is_treated = if_else(tx_ever == 1L & t >= t_switch, 1, 0)
      ) %>%
      left_join(g_grid %>% select(id, t, S_A_denom, S_A_num), by = c("id", "t")) %>%
      left_join(S_m_lookup, by = "id") %>%
      mutate(
        S_A_denom = if_else(t > grace_period, S_m_denom, S_A_denom),
        S_A_num   = if_else(t > grace_period, S_m_num, S_A_num)
      )

    if (align_treat) {
      panel <- panel %>%
        mutate(
          w_ipcw_time = case_when(
            arm == "control_no_tx" ~ S_A_num / bound_precision(S_A_denom),
            arm == "treat_tx" & is_treated == 0 & t <= grace_period ~
              pmax(S_A_denom - S_m_denom, 0) / (bound_precision(S_A_denom) * bound_precision(1 - S_m_denom)),
            arm == "treat_tx" & is_treated == 1 & tx_ever == 1L & t_switch <= grace_period ~
              1 / bound_precision(1 - S_m_denom),
            TRUE ~ 0
          ),
          w_unstab_time = case_when(
            arm == "control_no_tx" ~ 1 / bound_precision(S_A_denom),
            arm == "treat_tx" ~ w_ipcw_time,
            TRUE ~ 0
          )
        )
    }

    # Arm-stratified trimming on positive weights (keep structural zeros as 0)
    panel_w <- panel %>%
      group_by(arm) %>%
      mutate(
        q_lo = quantile(w_ipcw_time[w_ipcw_time > 0], 0.001, na.rm = TRUE),
        q_hi = quantile(w_ipcw_time[w_ipcw_time > 0], 0.999, na.rm = TRUE),
        qu_lo = quantile(w_unstab_time[w_unstab_time > 0], 0.001, na.rm = TRUE),
        qu_hi = quantile(w_unstab_time[w_unstab_time > 0], 0.999, na.rm = TRUE),
        w_trim = if_else(w_ipcw_time > 0, pmin(pmax(w_ipcw_time, q_lo), q_hi), 0),
        w_unstab_trim = if_else(w_unstab_time > 0, pmin(pmax(w_unstab_time, qu_lo), qu_hi), 0),
        post_grace = if_else(t > grace_period, 1, 0)
      ) %>%
      ungroup()

    # Check CCW-KM
    km <- compute_ccw_km_inference(panel_w, tau = tau, n_boot = 500, seed = s)$results

    # Check CCW-CPHM (Regime 1: qf_corr) filtering zero weights
    cphm <- fit_ccw_cphm(panel_w %>% filter(w_unstab_trim > 0), covars = covs, tau = tau, q_form = qf_corr)

    # Check TTE-TMLE Regime 1 & Regime 2
    run_tmle_reg <- function(qf) {
      qm <- fit_q_models(panel_w, covars = covs, q_form = qf)
      tg <- target_hazard(panel_w, qm$pred_Q_init, grace_period = grace_period)

      set.seed(s)
      comp_t <- df %>% filter(tx_ever == 1L, t_switch <= grace_period) %>% pull(t_switch)
      if (length(comp_t) == 0) comp_t <- grace_period

      cf_base <- df %>% select(id, tx_ever, t_switch, all_of(covs)) %>%
        mutate(sampled_t = sample(comp_t, n(), replace = TRUE))

      cf_data <- cf_base %>%
        uncount(weights = tau, .id = "t") %>%
        mutate(t_sq = t^2, post_grace = if_else(t > grace_period, 1, 0)) %>%
        uncount(weights = 2, .id = "clone_arm") %>%
        mutate(
          arm = if_else(clone_arm == 1L, "control_no_tx", "treat_tx"),
          cf_t_switch = case_when(
            arm == "control_no_tx" ~ 999,
            arm == "treat_tx" & tx_ever == 1L & t_switch <= grace_period ~ t_switch,
            arm == "treat_tx" ~ sampled_t
          ),
          is_treated = if_else(t >= cf_t_switch, 1, 0)
        ) %>%
        left_join(g_grid %>% select(id, t, S_A_denom, S_A_num), by = c("id", "t")) %>%
        left_join(S_m_lookup, by = "id") %>%
        mutate(
          S_A_denom = if_else(t > grace_period, S_m_denom, S_A_denom),
          S_A_num   = if_else(t > grace_period, S_m_num, S_A_num)
        )

      q0_ctrl <- predict(qm$m_Q_ctrl, newdata = cf_data, type = "response")
      q0_trt  <- predict(qm$m_Q_treat, newdata = cf_data, type = "response")
      cf_data$cf_Q_init_hazard <- bound_precision(if_else(cf_data$arm == "treat_tx", q0_trt, q0_ctrl))
      cf_data$cf_Q_init_logit <- qlogis(cf_data$cf_Q_init_hazard)

      q_arm <- panel_w %>% group_by(arm) %>% summarise(lo = first(q_lo), hi = first(q_hi), .groups = "drop")
      lo_c <- q_arm$lo[q_arm$arm == "control_no_tx"]; hi_c <- q_arm$hi[q_arm$arm == "control_no_tx"]
      lo_t <- q_arm$lo[q_arm$arm == "treat_tx"];      hi_t <- q_arm$hi[q_arm$arm == "treat_tx"]

      cf_data <- cf_data %>%
        mutate(
          cf_H_raw = case_when(
            arm == "control_no_tx" ~ S_A_num / bound_precision(S_A_denom),
            arm == "treat_tx" & is_treated == 0 & t <= grace_period ~
              pmax(S_A_denom - S_m_denom, 0) / (bound_precision(S_A_denom) * bound_precision(1 - S_m_denom)),
            arm == "treat_tx" & is_treated == 1 ~ 1 / bound_precision(1 - S_m_denom),
            TRUE ~ 0
          ),
          cf_H = if_else(arm == "control_no_tx", pmin(pmax(cf_H_raw, lo_c), hi_c),
                         if_else(cf_H_raw > 0, pmin(pmax(cf_H_raw, lo_t), hi_t), 0)),
          cf_H_ctrl = if_else(arm == "control_no_tx", cf_H, 0),
          cf_H_treat_pre = if_else(arm == "treat_tx" & t <= grace_period, cf_H, 0),
          cf_H_treat_post = if_else(arm == "treat_tx" & t > grace_period, cf_H, 0),
          Q_star_hazard = plogis(cf_Q_init_logit + tg$eps_ctrl * cf_H_ctrl + tg$eps_treat_pre * cf_H_treat_pre + tg$eps_treat_post * cf_H_treat_post)
        ) %>%
        arrange(id, arm, t) %>%
        group_by(id, arm) %>%
        mutate(indiv_survival = cumprod(1 - Q_star_hazard)) %>%
        ungroup()

      tmle_surv <- cf_data %>% group_by(arm, t) %>% summarise(survival = mean(indiv_survival), .groups = "drop")
      res_tmle <- tmle_surv %>% select(t, arm, survival) %>% pivot_wider(names_from = arm, values_from = survival) %>% mutate(rd = treat_tx - control_no_tx)
      inf <- compute_eic_inference(tg$panel_w, cf_data, tmle_surv, res_tmle, tg$eps_ctrl, tg$eps_treat_pre, tg$eps_treat_post, n_boot = 500, seed = s)
      inf$res_tmle
    }

    tmle_r1 <- run_tmle_reg(qf_corr)
    tmle_r2 <- run_tmle_reg(qf_mis)

    rows[[length(rows) + 1]] <- tibble(
      scen = scen_num,
      km_bias15 = km$rd[15] - truth$rd_true[15],
      km_bias30 = km$rd[30] - truth$rd_true[30],
      km_cov30  = as.numeric(km$rd_ci_lower[30] <= truth$rd_true[30] & truth$rd_true[30] <= km$rd_ci_upper[30]),
      cphm_bias30 = cphm$rd_tau - truth$rd_true[30],
      r1_bias30 = tmle_r1$rd[30] - truth$rd_true[30],
      r1_cov30  = as.numeric(tmle_r1$rd_ci_lower[30] <= truth$rd_true[30] & truth$rd_true[30] <= tmle_r1$rd_ci_upper[30]),
      r2_bias30 = tmle_r2$rd[30] - truth$rd_true[30],
      r2_cov30  = as.numeric(tmle_r2$rd_ci_lower[30] <= truth$rd_true[30] & truth$rd_true[30] <= tmle_r2$rd_ci_upper[30])
    )
  }
  bind_rows(rows) %>% summarise(across(-scen, mean)) %>% mutate(scen = scen_num, mode = ctrl_mode, .before = 1)
}

cat("=== Testing ns3_int + align_treat=TRUE across Scenarios 1, 2, 4, 5, 7 (15 seeds each) ===\n")
res_ns3 <- bind_rows(lapply(c(1, 2, 4, 5, 7), run_test_scen, seeds = 1:15, ctrl_mode = "ns3_int", align_treat = TRUE))
print(as.data.frame(res_ns3), digits = 4)

cat("\n=== Testing ns4_lin_int + align_treat=TRUE across Scenarios 1, 2, 4, 5, 7 (15 seeds each) ===\n")
res_ns4lin <- bind_rows(lapply(c(1, 2, 4, 5, 7), run_test_scen, seeds = 1:15, ctrl_mode = "ns4_lin_int", align_treat = TRUE))
print(as.data.frame(res_ns4lin), digits = 4)

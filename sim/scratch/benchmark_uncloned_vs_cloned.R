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

# Function to run Uncloned Stochastic TMLE
run_uncloned_tmle <- function(data, 
                              grace_period = 10, 
                              tau = 30, 
                              covars = "x1", 
                              q_form = NULL, 
                              p_form = NULL,
                              method = c("conditional", "marginal"),
                              integration = c("analytical", "mc")) {
  
  method <- match.arg(method)
  integration <- match.arg(integration)
  
  df_std <- prepare_tte_data(data)
  
  # 1. Uncloned person-time panel
  panel <- df_std %>%
    mutate(n_int = pmin(time, tau)) %>%
    uncount(weights = n_int, .id = "t") %>%
    mutate(
      t_sq = t^2,
      event_day = if_else(event == 1L & t == time, 1L, 0L),
      is_tx = if_else(tx_ever == 1L & t >= t_switch, 1L, 0L),
      init_day = if_else(tx_ever == 1L & t == t_switch, 1L, 0L)
    )
  
  # 2. Discrete initiation hazard model
  init_data <- panel %>%
    filter(t <= grace_period, (tx_ever == 0L | t <= t_switch))
  
  if (is.null(p_form)) {
    p_form_str <- paste(covars, collapse = " + ")
  } else {
    p_form_str <- p_form
  }
  
  init_formula <- as.formula(paste("init_day ~", p_form_str, "+ t + t_sq"))
  m_init <- suppressWarnings(glm(init_formula, family = binomial(), data = init_data))
  
  # Predict initiation probabilities across k in 1..grace_period
  grid_k <- expand_grid(
    id = unique(df_std$id),
    k = 1:grace_period
  ) %>%
    left_join(df_std %>% select(id, all_of(covars)), by = "id") %>%
    mutate(t = k, t_sq = k^2)
  
  grid_k$pred_lambda <- bound_precision(suppressWarnings(predict(m_init, newdata = grid_k, type = "response")), tol = 1e-4)
  
  # Compute S_A(t | W) and P(T_A = k | W)
  patient_k_probs <- grid_k %>%
    arrange(id, k) %>%
    group_by(id) %>%
    mutate(
      S_A = cumprod(1 - pred_lambda),
      S_prev = lag(S_A, default = 1),
      P_k = S_prev * pred_lambda
    ) %>%
    ungroup()
  
  patient_S_m <- patient_k_probs %>%
    filter(k == grace_period) %>%
    select(id, S_m = S_A)
  
  patient_t_switch_prob <- df_std %>%
    filter(tx_ever == 1L, t_switch <= grace_period) %>%
    select(id, t_switch) %>%
    left_join(patient_k_probs %>% select(id, k, P_switch = P_k), by = c("id", "t_switch" = "k")) %>%
    select(id, P_switch)
  
  # Marginal distribution of complier times
  complier_times <- df_std %>%
    filter(tx_ever == 1L, t_switch <= grace_period) %>%
    pull(t_switch)
  if (length(complier_times) == 0) complier_times <- grace_period
  p_star_marg <- table(factor(complier_times, levels = 1:grace_period)) / length(complier_times)
  
  # 3. Compute Clever Covariates
  panel <- panel %>%
    left_join(patient_k_probs %>% select(id, k, S_A), by = c("id", "t" = "k")) %>%
    left_join(patient_S_m, by = "id") %>%
    left_join(patient_t_switch_prob, by = "id") %>%
    arrange(id, t) %>%
    group_by(id) %>%
    fill(S_A, .direction = "down") %>%
    ungroup() %>%
    mutate(
      S_A = if_else(is.na(S_A), S_m, S_A),
      p_star_switch = if_else(tx_ever == 1L & t_switch <= grace_period, as.numeric(p_star_marg[as.character(t_switch)]), 0),
      # H0
      H0 = if_else(is_tx == 0L, 1 / bound_precision(S_A, tol = 1e-4), 0),
      # H1 Conditional
      H1_cond = case_when(
        is_tx == 0L & t <= grace_period ~ (S_A - S_m) / (bound_precision(S_A, tol=1e-4) * bound_precision(1 - S_m, tol=1e-4)),
        is_tx == 1L & tx_ever == 1L & t_switch <= grace_period ~ 1 / bound_precision(1 - S_m, tol = 1e-4),
        TRUE ~ 0
      ),
      # H1 Marginal
      H1_marg = case_when(
        is_tx == 0L & t <= grace_period ~ {
          vapply(t, function(cur_t) sum(p_star_marg[(cur_t+1):grace_period]), numeric(1)) / bound_precision(S_A, tol=1e-4)
        },
        is_tx == 1L & tx_ever == 1L & t_switch <= grace_period ~ {
          p_star_switch / bound_precision(P_switch, tol = 1e-4)
        },
        TRUE ~ 0
      )
    )
  
  # Choose H1
  panel$H1 <- if (method == "conditional") panel$H1_cond else panel$H1_marg
  
  # 4. Outcome Hazard Model (Q-model)
  if (is.null(q_form)) {
    q_form_str <- paste(covars, collapse = " + ")
  } else {
    q_form_str <- q_form
  }
  
  panel <- panel %>% mutate(post_grace = if_else(t > grace_period, 1L, 0L))
  q_formula <- as.formula(paste("event_day ~", q_form_str, "+ t + t_sq + post_grace + is_tx"))
  m_Q <- suppressWarnings(glm(q_formula, family = binomial(), data = panel))
  
  panel$Q0_hazard <- bound_precision(suppressWarnings(predict(m_Q, type = "response")), tol = 1e-4)
  panel$Q0_logit <- qlogis(panel$Q0_hazard)
  
  # 5. Separate Targeting Fluctuation
  m_fluct_0 <- suppressWarnings(glm(event_day ~ -1 + H0 + offset(Q0_logit), family = binomial(), data = panel))
  m_fluct_1 <- suppressWarnings(glm(event_day ~ -1 + H1 + offset(Q0_logit), family = binomial(), data = panel))
  
  eps0 <- unname(coef(m_fluct_0)["H0"])
  eps1 <- unname(coef(m_fluct_1)["H1"])
  if (is.na(eps0)) eps0 <- 0
  if (is.na(eps1)) eps1 <- 0
  
  # 6. Counterfactual G-Computation Marginalization
  # Target 0: Never Treat
  cf_ctrl <- expand_grid(
    id = unique(df_std$id),
    t = 1:tau
  ) %>%
    left_join(df_std %>% select(id, all_of(covars)), by = "id") %>%
    left_join(patient_k_probs %>% select(id, k, S_A), by = c("id", "t" = "k")) %>%
    left_join(patient_S_m, by = "id") %>%
    arrange(id, t) %>%
    group_by(id) %>%
    fill(S_A, .direction = "down") %>%
    ungroup() %>%
    mutate(
      S_A = if_else(is.na(S_A), S_m, S_A),
      t_sq = t^2,
      post_grace = if_else(t > grace_period, 1L, 0L),
      is_tx = 0L,
      H0 = 1 / bound_precision(S_A, tol = 1e-4)
    )
  
  cf_ctrl$Q0_hazard <- bound_precision(suppressWarnings(predict(m_Q, newdata = cf_ctrl, type = "response")), tol = 1e-4)
  cf_ctrl$Q0_logit <- qlogis(cf_ctrl$Q0_hazard)
  cf_ctrl$Q_star_0 <- plogis(cf_ctrl$Q0_logit + eps0 * cf_ctrl$H0)
  
  cf_ctrl <- cf_ctrl %>%
    arrange(id, t) %>%
    group_by(id) %>%
    mutate(S0_indiv = cumprod(1 - Q_star_0)) %>%
    ungroup()
  
  S0_pop <- cf_ctrl %>% group_by(t) %>% summarise(S0 = mean(S0_indiv), .groups = "drop")
  
  # Target 1: Treated Arm
  if (integration == "analytical") {
    cf_treat_grid <- expand_grid(
      id = unique(df_std$id),
      k_init = 1:grace_period,
      t = 1:tau
    ) %>%
      left_join(df_std %>% select(id, all_of(covars)), by = "id") %>%
      left_join(patient_k_probs %>% select(id, k, S_A), by = c("id", "t" = "k")) %>%
      left_join(patient_S_m, by = "id") %>%
      left_join(patient_k_probs %>% select(id, k, P_k), by = c("id", "k_init" = "k")) %>%
      arrange(id, k_init, t) %>%
      group_by(id, k_init) %>%
      fill(S_A, .direction = "down") %>%
      ungroup() %>%
      mutate(
        S_A = if_else(is.na(S_A), S_m, S_A),
        t_sq = t^2,
        post_grace = if_else(t > grace_period, 1L, 0L),
        is_tx = if_else(t >= k_init, 1L, 0L),
        H1 = if (method == "conditional") {
          case_when(
            t < k_init & t <= grace_period ~ (S_A - S_m) / (bound_precision(S_A, tol=1e-4) * bound_precision(1 - S_m, tol=1e-4)),
            t >= k_init ~ 1 / bound_precision(1 - S_m, tol = 1e-4),
            TRUE ~ 0
          )
        } else {
          case_when(
            t < k_init & t <= grace_period ~ vapply(t, function(cur_t) sum(p_star_marg[(cur_t+1):grace_period]), numeric(1)) / bound_precision(S_A, tol=1e-4),
            t >= k_init ~ as.numeric(p_star_marg[as.character(k_init)]) / bound_precision(P_k, tol = 1e-4),
            TRUE ~ 0
          )
        }
      )
    
    cf_treat_grid$Q0_hazard <- bound_precision(suppressWarnings(predict(m_Q, newdata = cf_treat_grid, type = "response")), tol = 1e-4)
    cf_treat_grid$Q0_logit <- qlogis(cf_treat_grid$Q0_hazard)
    cf_treat_grid$Q_star_1 <- plogis(cf_treat_grid$Q0_logit + eps1 * cf_treat_grid$H1)
    
    cf_treat_grid <- cf_treat_grid %>%
      arrange(id, k_init, t) %>%
      group_by(id, k_init) %>%
      mutate(S1_indiv_k = cumprod(1 - Q_star_1)) %>%
      ungroup()
    
    weights_k <- patient_k_probs %>%
      left_join(patient_S_m, by = "id") %>%
      mutate(
        w = if (method == "conditional") P_k / bound_precision(1 - S_m, tol = 1e-4) else as.numeric(p_star_marg[as.character(k)])
      ) %>%
      select(id, k_init = k, w)
    
    S1_pop <- cf_treat_grid %>%
      left_join(weights_k, by = c("id", "k_init")) %>%
      group_by(id, t) %>%
      summarise(S1_indiv = sum(w * S1_indiv_k), .groups = "drop") %>%
      group_by(t) %>%
      summarise(S1 = mean(S1_indiv), .groups = "drop")
    
  } else {
    # Monte Carlo integration
    cf_base <- df_std %>%
      mutate(sampled_t = sample(complier_times, n(), replace = TRUE)) %>%
      select(id, sampled_t, all_of(covars))
    
    cf_treat_mc <- cf_base %>%
      uncount(weights = tau, .id = "t") %>%
      left_join(patient_k_probs %>% select(id, k, S_A), by = c("id", "t" = "k")) %>%
      left_join(patient_S_m, by = "id") %>%
      left_join(patient_k_probs %>% select(id, k, P_k), by = c("id", "sampled_t" = "k")) %>%
      arrange(id, t) %>%
      group_by(id) %>%
      fill(S_A, .direction = "down") %>%
      ungroup() %>%
      mutate(
        S_A = if_else(is.na(S_A), S_m, S_A),
        t_sq = t^2,
        post_grace = if_else(t > grace_period, 1L, 0L),
        is_tx = if_else(t >= sampled_t, 1L, 0L),
        H1 = if (method == "conditional") {
          case_when(
            t < sampled_t & t <= grace_period ~ (S_A - S_m) / (bound_precision(S_A, tol=1e-4) * bound_precision(1 - S_m, tol=1e-4)),
            t >= sampled_t ~ 1 / bound_precision(1 - S_m, tol = 1e-4),
            TRUE ~ 0
          )
        } else {
          case_when(
            t < sampled_t & t <= grace_period ~ vapply(t, function(cur_t) sum(p_star_marg[(cur_t+1):grace_period]), numeric(1)) / bound_precision(S_A, tol=1e-4),
            t >= sampled_t ~ as.numeric(p_star_marg[as.character(sampled_t)]) / bound_precision(P_k, tol = 1e-4),
            TRUE ~ 0
          )
        }
      )
    
    cf_treat_mc$Q0_hazard <- bound_precision(suppressWarnings(predict(m_Q, newdata = cf_treat_mc, type = "response")), tol = 1e-4)
    cf_treat_mc$Q0_logit <- qlogis(cf_treat_mc$Q0_hazard)
    cf_treat_mc$Q_star_1 <- plogis(cf_treat_mc$Q0_logit + eps1 * cf_treat_mc$H1)
    
    cf_treat_mc <- cf_treat_mc %>%
      arrange(id, t) %>%
      group_by(id) %>%
      mutate(S1_indiv = cumprod(1 - Q_star_1)) %>%
      ungroup()
    
    S1_pop <- cf_treat_mc %>%
      group_by(t) %>%
      summarise(S1 = mean(S1_indiv), .groups = "drop")
  }
  
  res <- S1_pop %>%
    left_join(S0_pop, by = "t") %>%
    mutate(rd = S1 - S0)
  
  return(list(
    res = res,
    eps0 = eps0,
    eps1 = eps1
  ))
}

cat("Testing run_uncloned_tmle across configurations on single dataset...\n")
dat <- generate_scenario1_data(n = 1000, seed = 100)
truth <- compute_scenario1_truth()

t1 <- system.time(m_cond_ana <- run_uncloned_tmle(dat, method = "conditional", integration = "analytical"))
t2 <- system.time(m_cond_mc  <- run_uncloned_tmle(dat, method = "conditional", integration = "mc"))
t3 <- system.time(m_marg_ana <- run_uncloned_tmle(dat, method = "marginal", integration = "analytical"))
t4 <- system.time(m_cloned   <- tte_tmle(dat, covars = "x1", p_form = "x1", q_form = "x1", use_sl = FALSE))

cat(sprintf("Runtimes:\n  Uncloned (Cond + Analytical): %.2fs\n  Uncloned (Cond + MC):         %.2fs\n  Uncloned (Marg + Analytical): %.2fs\n  Cloned TTE-TMLE:              %.2fs\n",
            t1["elapsed"], t2["elapsed"], t3["elapsed"], t4["elapsed"]))

cat("\nEstimates at Day 15:\n")
cat(sprintf("  Ground Truth:                  %.4f\n", truth$rd_true[15]))
cat(sprintf("  Cloned TTE-TMLE:               %.4f (Bias: %+.4f)\n", m_cloned$res_tmle$rd[15], m_cloned$res_tmle$rd[15] - truth$rd_true[15]))
cat(sprintf("  Uncloned (Cond + Analytical):  %.4f (Bias: %+.4f)\n", m_cond_ana$res$rd[15], m_cond_ana$res$rd[15] - truth$rd_true[15]))
cat(sprintf("  Uncloned (Cond + MC):          %.4f (Bias: %+.4f)\n", m_cond_mc$res$rd[15], m_cond_mc$res$rd[15] - truth$rd_true[15]))
cat(sprintf("  Uncloned (Marg + Analytical):  %.4f (Bias: %+.4f)\n", m_marg_ana$res$rd[15], m_marg_ana$res$rd[15] - truth$rd_true[15]))

cat("\nEstimates at Day 30:\n")
cat(sprintf("  Ground Truth:                  %.4f\n", truth$rd_true[30]))
cat(sprintf("  Cloned TTE-TMLE:               %.4f (Bias: %+.4f)\n", m_cloned$res_tmle$rd[30], m_cloned$res_tmle$rd[30] - truth$rd_true[30]))
cat(sprintf("  Uncloned (Cond + Analytical):  %.4f (Bias: %+.4f)\n", m_cond_ana$res$rd[30], m_cond_ana$res$rd[30] - truth$rd_true[30]))
cat(sprintf("  Uncloned (Cond + MC):          %.4f (Bias: %+.4f)\n", m_cond_mc$res$rd[30], m_cond_mc$res$rd[30] - truth$rd_true[30]))
cat(sprintf("  Uncloned (Marg + Analytical):  %.4f (Bias: %+.4f)\n", m_marg_ana$res$rd[30], m_marg_ana$res$rd[30] - truth$rd_true[30]))

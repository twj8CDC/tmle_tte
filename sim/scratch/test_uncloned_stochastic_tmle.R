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

cat("=== Generating Scenario 1 Data (N=1000, seed=42) ===\n")
dat <- generate_scenario1_data(n = 1000, tau = 30, grace_period = 10, seed = 42)
truth <- compute_scenario1_truth(tau = 30, grace_period = 10)
rd15_true <- truth$rd_true[15]
rd30_true <- truth$rd_true[30]
cat(sprintf("Ground Truth: RD(15) = %.4f, RD(30) = %.4f\n", rd15_true, rd30_true))

# Standardize columns
df_std <- prepare_tte_data(dat)

cat("\n=== Running Existing Cloned TTE-TMLE (GLM) ===\n")
fit_cloned <- tte_tmle(
  data = dat,
  grace_period = 10,
  tau = 30,
  covars = "x1",
  p_form = "x1",
  q_form = "x1",
  use_sl = FALSE
)
res_cloned <- fit_cloned$res_tmle
rd15_cloned <- res_cloned$rd[res_cloned$t == 15]
rd30_cloned <- res_cloned$rd[res_cloned$t == 30]
cat(sprintf("Cloned TTE-TMLE: RD(15) = %.4f (Bias: %.4f), RD(30) = %.4f (Bias: %.4f)\n",
            rd15_cloned, rd15_cloned - rd15_true,
            rd30_cloned, rd30_cloned - rd30_true))
cat("Cloned eps:\n")
print(fit_cloned$eps)

cat("\n=== Building Uncloned Stochastic Survival TMLE ===\n")

# 1. Expand raw data to uncloned person-time panel
tau <- 30
grace_period <- 10

uncloned_panel <- df_std %>%
  mutate(n_int = pmin(time, tau)) %>%
  uncount(weights = n_int, .id = "t") %>%
  mutate(
    event_day = if_else(event == 1L & t == time, 1L, 0L),
    # Treatment indicator: 1 if patient initiated on or before day t
    is_tx = if_else(tx_ever == 1L & t >= t_switch, 1L, 0L),
    # Event of initiation at day t:
    init_day = if_else(tx_ever == 1L & t == t_switch, 1L, 0L)
  )

cat(sprintf("Uncloned panel rows: %d (vs Cloned panel rows: %d)\n",
            nrow(uncloned_panel), nrow(fit_cloned$panel_w)))

# 2. Fit treatment initiation hazard model on untreated person-time during grace period
init_data <- uncloned_panel %>%
  filter(t <= grace_period, (tx_ever == 0L | t <= t_switch)) %>%
  mutate(t_sq = t^2)

m_init <- glm(init_day ~ x1 + t + t_sq, family = binomial(), data = init_data)

# Predict initiation hazard for each person at each day k in 1..m
grid_k <- expand_grid(
  id = unique(df_std$id),
  k = 1:grace_period
) %>%
  left_join(df_std %>% select(id, x1), by = "id") %>%
  mutate(t = k, t_sq = k^2)

grid_k$pred_lambda <- bound_precision(predict(m_init, newdata = grid_k, type = "response"), tol = 1e-4)

# Calculate S_0,A(t | W) = prod_{s=1}^t (1 - lambda(s | W))
grid_k <- grid_k %>%
  arrange(id, k) %>%
  group_by(id) %>%
  mutate(
    S_A = cumprod(1 - pred_lambda),
    # P_0(T_A = k | W) = S_A(k-1) * lambda(k)
    S_prev = lag(S_A, default = 1),
    P_k = S_prev * pred_lambda
  ) %>%
  ungroup()

# Extract S_A at each day and at day m
S_A_wide <- grid_k %>%
  select(id, k, S_A, P_k)

S_A_m <- grid_k %>%
  filter(k == grace_period) %>%
  select(id, S_m = S_A)

# Marginal empirical distribution of complier times
complier_times <- df_std %>%
  filter(tx_ever == 1L, t_switch <= grace_period) %>%
  pull(t_switch)

p_star_marg <- table(factor(complier_times, levels = 1:grace_period)) / length(complier_times)

# 3. Attach weights / clever covariates to uncloned panel
uncloned_panel <- uncloned_panel %>%
  left_join(S_A_wide %>% rename(t = k), by = c("id", "t")) %>%
  left_join(S_A_m, by = "id") %>%
  arrange(id, t) %>%
  group_by(id) %>%
  fill(S_A, .direction = "down") %>%
  ungroup() %>%
  mutate(
    S_A = if_else(is.na(S_A), S_m, S_A), # for t > grace_period, S_A stays S_m
    # Untreated Clever Covariate:
    # H0(t) = 1 / S_A(t) if untreated (is_tx == 0), else 0
    H0 = if_else(is_tx == 0L, 1 / bound_precision(S_A, tol = 1e-4), 0),
    # Treated Clever Covariate (Conditional):
    # If untreated at t <= m: (S_A(t) - S_m) / (S_A(t) * (1 - S_m))
    # If initiated at t_switch <= m and t >= t_switch: 1 / (1 - S_m)
    # If untreated at t > m: 0
    H1_cond = case_when(
      is_tx == 0L & t <= grace_period ~ (S_A - S_m) / (bound_precision(S_A, tol=1e-4) * bound_precision(1 - S_m, tol=1e-4)),
      is_tx == 1L & tx_ever == 1L & t_switch <= grace_period ~ 1 / bound_precision(1 - S_m, tol = 1e-4),
      TRUE ~ 0
    ),
    # Treated Clever Covariate (Marginal):
    # Post-switch: p_star(t_switch) / P_0(t_switch | W)
    p_star_val = if_else(tx_ever == 1L & t_switch <= grace_period, as.numeric(p_star_marg[as.character(t_switch)]), 0),
    H1_marg = case_when(
      is_tx == 0L & t <= grace_period ~ {
        # sum_{j=t+1}^m p*(j) / S_A(t)
        vapply(t, function(cur_t) sum(p_star_marg[(cur_t+1):grace_period]), numeric(1)) / bound_precision(S_A, tol=1e-4)
      },
      is_tx == 1L & tx_ever == 1L & t_switch <= grace_period ~ {
        p_star_val / bound_precision(P_k, tol = 1e-4)
      },
      TRUE ~ 0
    )
  )

cat("Summary of H0:\n")
print(summary(uncloned_panel$H0))
cat("Summary of H1_cond:\n")
print(summary(uncloned_panel$H1_cond))
cat("Summary of H1_marg:\n")
print(summary(uncloned_panel$H1_marg))

# 4. Fit single pooled outcome hazard model Q0(t, is_tx, W)
uncloned_panel <- uncloned_panel %>%
  mutate(t_sq = t^2, post_grace = if_else(t > grace_period, 1L, 0L))

m_Q <- glm(event_day ~ x1 + t + t_sq + post_grace + is_tx, family = binomial(), data = uncloned_panel)
uncloned_panel$Q0_pred <- bound_precision(predict(m_Q, type = "response"), tol = 1e-4)
uncloned_panel$Q0_logit <- qlogis(uncloned_panel$Q0_pred)

# 5. Fluctuation Model
cat("\n--- Testing Fluctuation: Separate 1-param glm vs Joint 2-param glm ---\n")
# Option A: Joint
m_fluct_joint <- glm(event_day ~ -1 + H0 + H1_cond + offset(Q0_logit), family = binomial(), data = uncloned_panel)
cat("Joint fluctuation coefs (Conditional):\n")
print(coef(m_fluct_joint))

# Option B: Separate
m_fluct_0 <- glm(event_day ~ -1 + H0 + offset(Q0_logit), family = binomial(), data = uncloned_panel)
m_fluct_1_cond <- glm(event_day ~ -1 + H1_cond + offset(Q0_logit), family = binomial(), data = uncloned_panel)
m_fluct_1_marg <- glm(event_day ~ -1 + H1_marg + offset(Q0_logit), family = binomial(), data = uncloned_panel)
cat(sprintf("Separate coefs: eps0 = %.5f, eps1_cond = %.5f, eps1_marg = %.5f\n",
            coef(m_fluct_0)["H0"], coef(m_fluct_1_cond)["H1_cond"], coef(m_fluct_1_marg)["H1_marg"]))

# 6. G-Computation / Marginalization
# Counterfactual Control (Never Treat): A_t = 0 for all t
cf_ctrl <- expand_grid(
  id = unique(df_std$id),
  t = 1:tau
) %>%
  left_join(df_std %>% select(id, x1), by = "id") %>%
  left_join(S_A_wide %>% rename(t = k), by = c("id", "t")) %>%
  left_join(S_A_m, by = "id") %>%
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

cf_ctrl$Q0_hazard <- bound_precision(predict(m_Q, newdata = cf_ctrl, type = "response"), tol = 1e-4)
cf_ctrl$Q0_logit <- qlogis(cf_ctrl$Q0_hazard)

# Fluctuate control
eps0 <- coef(m_fluct_0)["H0"]
cf_ctrl <- cf_ctrl %>%
  mutate(
    # Targeted hazard
    Q_star_0 = plogis(Q0_logit + eps0 * H0),
    # Untargeted pure gcomp hazard
    Q_pure_0 = Q0_hazard
  ) %>%
  arrange(id, t) %>%
  group_by(id) %>%
  mutate(
    S_indiv_star_0 = cumprod(1 - Q_star_0),
    S_indiv_pure_0 = cumprod(1 - Q_pure_0)
  ) %>%
  ungroup()

S0_star_pop <- cf_ctrl %>% group_by(t) %>% summarise(S0 = mean(S_indiv_star_0), .groups="drop")
S0_pure_pop <- cf_ctrl %>% group_by(t) %>% summarise(S0 = mean(S_indiv_pure_0), .groups="drop")

# Counterfactual Treated: Analytical expectation over initiation times k in 1..m
# For each person i and each initiation day k:
cf_treat_grid <- expand_grid(
  id = unique(df_std$id),
  k_init = 1:grace_period,
  t = 1:tau
) %>%
  left_join(df_std %>% select(id, x1), by = "id") %>%
  left_join(S_A_wide %>% rename(t = k), by = c("id", "t")) %>%
  left_join(S_A_m, by = "id") %>%
  arrange(id, k_init, t) %>%
  group_by(id, k_init) %>%
  fill(S_A, .direction = "down") %>%
  ungroup() %>%
  mutate(
    S_A = if_else(is.na(S_A), S_m, S_A),
    t_sq = t^2,
    post_grace = if_else(t > grace_period, 1L, 0L),
    is_tx = if_else(t >= k_init, 1L, 0L),
    # Pre-switch: (S_A(t) - S_m) / (S_A(t) * (1 - S_m))
    # Post-switch: 1 / (1 - S_m)
    H1_cond = case_when(
      t < k_init & t <= grace_period ~ (S_A - S_m) / (bound_precision(S_A, tol=1e-4) * bound_precision(1 - S_m, tol=1e-4)),
      t >= k_init ~ 1 / bound_precision(1 - S_m, tol = 1e-4),
      TRUE ~ 0
    )
  )

cf_treat_grid$Q0_hazard <- bound_precision(predict(m_Q, newdata = cf_treat_grid, type = "response"), tol = 1e-4)
cf_treat_grid$Q0_logit <- qlogis(cf_treat_grid$Q0_hazard)

eps1_cond <- coef(m_fluct_1_cond)["H1_cond"]
cf_treat_grid <- cf_treat_grid %>%
  mutate(
    Q_star_1 = plogis(Q0_logit + eps1_cond * H1_cond),
    Q_pure_1 = Q0_hazard
  ) %>%
  arrange(id, k_init, t) %>%
  group_by(id, k_init) %>%
  mutate(
    S_indiv_star_k = cumprod(1 - Q_star_1),
    S_indiv_pure_k = cumprod(1 - Q_pure_1)
  ) %>%
  ungroup()

# Weights w_k(W_i) = P_0(T_A = k | W_i) / (1 - S_m(W_i)) for conditional
# Or w_k = p_star_marg[k] for marginal
weights_k <- grid_k %>%
  left_join(S_A_m, by = "id") %>%
  mutate(
    w_cond = P_k / bound_precision(1 - S_m, tol = 1e-4),
    w_marg = as.numeric(p_star_marg[as.character(k)])
  ) %>%
  select(id, k_init = k, w_cond, w_marg)

# Integrate over k
cf_treat_indiv <- cf_treat_grid %>%
  left_join(weights_k, by = c("id", "k_init")) %>%
  group_by(id, t) %>%
  summarise(
    S_indiv_star_cond = sum(w_cond * S_indiv_star_k),
    S_indiv_pure_cond = sum(w_cond * S_indiv_pure_k),
    S_indiv_star_marg = sum(w_marg * S_indiv_star_k),
    S_indiv_pure_marg = sum(w_marg * S_indiv_pure_k),
    .groups = "drop"
  )

# Marginalize over population
res_uncloned <- cf_treat_indiv %>%
  group_by(t) %>%
  summarise(
    S1_star_cond = mean(S_indiv_star_cond),
    S1_pure_cond = mean(S_indiv_pure_cond),
    S1_star_marg = mean(S_indiv_star_marg),
    S1_pure_marg = mean(S_indiv_pure_marg),
    .groups = "drop"
  ) %>%
  left_join(S0_star_pop %>% rename(S0_star = S0), by = "t") %>%
  left_join(S0_pure_pop %>% rename(S0_pure = S0), by = "t") %>%
  left_join(truth %>% select(t, S0_true, S1_true, rd_true), by = "t") %>%
  mutate(
    rd_star_cond = S1_star_cond - S0_star,
    rd_pure_cond = S1_pure_cond - S0_pure,
    rd_star_marg = S1_star_marg - S0_star,
    rd_pure_marg = S1_pure_marg - S0_pure
  )

cat("\n=== ESTIMATES COMPARISON AT DAY 15 AND DAY 30 ===\n")
t15 <- res_uncloned %>% filter(t == 15)
t30 <- res_uncloned %>% filter(t == 30)

cat(sprintf("Day 15 Truth: RD = %.4f\n", t15$rd_true))
cat(sprintf("  Cloned TTE-TMLE:       RD = %.4f (Bias = %+.4f)\n", rd15_cloned, rd15_cloned - t15$rd_true))
cat(sprintf("  Uncloned TMLE (Cond):  RD = %.4f (Bias = %+.4f)\n", t15$rd_star_cond, t15$rd_star_cond - t15$rd_true))
cat(sprintf("  Uncloned TMLE (Marg):  RD = %.4f (Bias = %+.4f)\n", t15$rd_star_marg, t15$rd_star_marg - t15$rd_true))
cat(sprintf("  Uncloned Pure G-Comp:  RD = %.4f (Bias = %+.4f)\n", t15$rd_pure_cond, t15$rd_pure_cond - t15$rd_true))

cat(sprintf("\nDay 30 Truth: RD = %.4f\n", t30$rd_true))
cat(sprintf("  Cloned TTE-TMLE:       RD = %.4f (Bias = %+.4f)\n", rd30_cloned, rd30_cloned - t30$rd_true))
cat(sprintf("  Uncloned TMLE (Cond):  RD = %.4f (Bias = %+.4f)\n", t30$rd_star_cond, t30$rd_star_cond - t30$rd_true))
cat(sprintf("  Uncloned TMLE (Marg):  RD = %.4f (Bias = %+.4f)\n", t30$rd_star_marg, t30$rd_star_marg - t30$rd_true))
cat(sprintf("  Uncloned Pure G-Comp:  RD = %.4f (Bias = %+.4f)\n", t30$rd_pure_cond, t30$rd_pure_cond - t30$rd_true))

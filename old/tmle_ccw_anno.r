library(dplyr)
library(tidyr)

run_ccw_tmle_analysis2 <- function(data,
                                   grace_period = 10,
                                   tau = 30,
                                   covars = c("x1"),
                                   dt = 1,
                                   event,
                                   event_time, 
                                   tx_obs, 
                                   tx_time_obs,
                                   sl_lib = c("SL.mean", "SL.glm", "SL.earth.int", "SL.ranger")) {
    
    if (!require("SuperLearner", quietly = TRUE)) stop("Package 'SuperLearner' is required.")

    # Create the 2-way interaction MARS wrapper for the SuperLearner library
    if ("SL.earth.int" %in% sl_lib && !exists("SL.earth.int", mode = "function")) {
        assign("SL.earth.int", function(..., degree = 2) {
            SuperLearner::SL.earth(..., degree = degree)
        }, envir = .GlobalEnv)
    }

    # Helper function to clean up bounding logic for probabilities and weights
    ptrim <- function(x, lower = 1e-4, upper = 1 - lower) {
        pmin(pmax(x, lower), upper)
    }

    # =====================================================================
    # STEP 1: THE ARCHITECTURE (Clone & Censor)
    # Goal: Eliminate immortal time bias by enforcing a strict Day 0.
    # Theory: We duplicate the cohort to evaluate two competing strategies.
    # If a clone violates their assigned strategy, they are artificially
    # censored. Death acts as a terminal event, not a censoring event.
    # =====================================================================

    data <- data %>%
        rename(event = event, time = event_time, tx_ever = tx_obs, t_switch = tx_time_obs)

    clones <- data %>%
        select(id, tx_ever, t_switch, event, time, all_of(covars)) %>%
        uncount(weights = 2, .id = "clone_arm") %>%
        mutate(
            arm = if_else(clone_arm == 1, "control_no_tx", "treat_tx"),
            clone_id = paste0(id, "_", clone_arm)
        )

    clones <- clones %>%
        mutate(
            # dev time is deviation time, but a better way to think of it is the last day they are fully compliant
            # nontreatment obs breaks the rule on day 5, then their last compliant day (dev_time) is day 4
            # treatment obs doesn't take the drug by end of the grace period day (ie day 10), 
            # then they are non-compliant on day 11 and their last compliant day (dev_time) is day 10 (last day of grace period) 
            # NOTE: treatment is assumed to always occur prior to event if they co-occur on the same day
            # so if non-treatment obs even on day 5 and breaks rule also on day 5. Then they will always have broken the rule
            # before the event is recorded (that is the logical assumption). So that obs last compliant day (dev_time) is day 4.
            dev_time = case_when(
                arm == "control_no_tx" & tx_ever == 1 ~ t_switch - 1,
                arm == "treat_tx" & !(tx_ever == 1 & t_switch <= grace_period) & time > grace_period ~ grace_period,
                TRUE ~ NA_real_
            ),
            # time_cc is then always the time of the event or censoring 
            time_cc = if_else(!is.na(dev_time), pmin(time, dev_time), time) %>% pmin(tau),
            # event_cc is true when the obs had an event and it occurs <= to the artificial censoring
            event_cc = if_else(event == 1 & time <= time_cc, 1, 0)
        )

    # =====================================================================
    # STEP 2: THE MECHANISM OF MISSINGNESS (IPCW / g-models)
    # Goal: Solve confounding by modeling protocol deviance (censoring).
    # Theory: We predict the probability of breaking the rules (event_tx).
    # H (the clever covariate) is the cumulative inverse probability of
    # remaining uncensored: H = 1 / P(missing).
    # =====================================================================

    # Create the control cencsoring model
    
    X_ctrl <- clones %>%
        filter(clone_arm == 1) %>%
        mutate(
            tmax = if_else(tx_ever == 1, t_switch, time),
            n_int = pmax(ceiling(tmax / dt), 1)
        ) %>%
        uncount(weights = n_int, .id = "t") %>%
        mutate(
            t_sq = t^2,
            event_tx = if_else(tx_ever == 1 & t == t_switch, 1, 0)
        )

    # NOTE: The treatment arm can only deviate if no treatment occurs <= grace_day. 
    # We only need to examine for tx after that point
    X_trt <- clones %>%
        filter(arm == "treat_tx", time > grace_period) %>%
        mutate(took_tx = if_else(tx_ever == 1 & t_switch <= grace_period, 1, 0))

    # Train Denominator g-models (SuperLearner)
    m_denom_ctrl <- SuperLearner(
        Y = X_ctrl$event_tx, 
        X = X_ctrl %>% select(all_of(covars), t, t_sq) %>% as.data.frame(), 
        family = binomial(), 
        SL.library = sl_lib, 
        cvControl = list(V = 2))

    m_denom_treat <- SuperLearner(
        Y = X_trt$took_tx, 
        X = X_trt %>% select(all_of(covars)) %>% as.data.frame(), 
        family = binomial(), 
        SL.library = sl_lib, 
        cvControl = list(V = 2))

    # Train Numerator models (Stabilization)
    m_num_ctrl <- glm(event_tx ~ t + t_sq, family = binomial(), data = X_ctrl)
    m_num_treat <- glm(took_tx ~ 1, family = binomial(), data = X_trt)

    # Expand Clones & Apply Weights
    panel <- clones %>%
        mutate(n_int = pmax(ceiling(time_cc / dt), 0)) %>%
        uncount(weights = n_int, .id = "t") %>%
        mutate(
            t0 = (t - 1) * dt,
            t1 = t * dt,
            t_sq = t^2,
            event_cc_day = if_else(event_cc == 1 & time_cc > t0 & time_cc <= t1, 1, 0),
            is_treated = if_else(tx_ever == 1 & t >= t_switch, 1, 0)
        )

    # predict probability of censor
    pred_ctrl_denom <- predict(
        m_denom_ctrl, 
        newdata = panel %>% select(all_of(covars), t, t_sq) %>% as.data.frame(),  # ctrl prediction
        onlySL = TRUE)$pred[, 1]

    pred_treat_denom <- predict(
        m_denom_treat, 
        newdata = panel %>% select(all_of(covars)) %>% as.data.frame() , 
        onlySL = TRUE)$pred[, 1]

    panel <- panel %>%
        mutate(
            p_ctrl_denom = pred_ctrl_denom %>% ptrim(),
            p_ctrl_num = predict(m_num_ctrl, newdata = ., type = "response") %>% ptrim(),
            p_treat_denom = pred_treat_denom %>% ptrim(),
            p_treat_num = predict(m_num_treat, newdata = ., type = "response") %>% ptrim(),
            w_inc = case_when(
                arm == "control_no_tx" ~ (1 - p_ctrl_num) / (1 - p_ctrl_denom),
                arm == "treat_tx" ~ case_when(t == grace_period + 1 ~ p_treat_num / p_treat_denom, TRUE ~ 1)
            )
        ) %>%
        arrange(clone_id, t) %>%
        group_by(clone_id) %>%
        mutate(w_ipcw_time = cumprod(w_inc)) %>%
        ungroup()

    q_lim <- quantile(panel$w_ipcw_time, probs = c(0.01, 0.99), na.rm = TRUE)

    panel_w <- panel %>%
        mutate(
            w_trim = w_ipcw_time %>% ptrim(lower = q_lim[1], upper = q_lim[2]),
            post_grace = if_else(t > grace_period, 1, 0)
        )

    fit_ccw <- survival::survfit(Surv(t0, t1, event_cc_day) ~ arm, data = panel_w, weights = w_trim, id = clone_id, robust = TRUE)

    # =====================================================================
    # STEP 3: THE LAZY IMPUTERS (Initial Q-models)
    # Goal: Estimate baseline survival under strict compliance.
    # Theory: We switch the outcome to the true clinical event (Death).
    # These models only see person-time where clones perfectly obeyed the
    # rules. They are naive to the missing data bias caused by censoring.
    # =====================================================================

    panel_ctrl <- panel_w %>% filter(arm == "control_no_tx")
    panel_treat <- panel_w %>% filter(arm == "treat_tx")

    X_Q_ctrl <- panel_ctrl %>%
        select(all_of(covars), t, t_sq, post_grace, is_treated) %>%
        as.data.frame()
    X_Q_treat <- panel_treat %>%
        select(all_of(covars), t, t_sq, post_grace, is_treated) %>%
        as.data.frame()

    m_Q_ctrl <- SuperLearner(Y = panel_ctrl$event_cc_day, X = X_Q_ctrl, family = binomial(), SL.library = sl_lib, cvControl = list(V = 2))
    m_Q_treat <- SuperLearner(Y = panel_treat$event_cc_day, X = X_Q_treat, family = binomial(), SL.library = sl_lib, cvControl = list(V = 2))

    X_Q_full <- panel_w %>%
        select(all_of(covars), t, t_sq, post_grace, is_treated) %>%
        as.data.frame()
    pred_Q_ctrl <- predict(m_Q_ctrl, newdata = X_Q_full, onlySL = TRUE)$pred[, 1]
    pred_Q_treat <- predict(m_Q_treat, newdata = X_Q_full, onlySL = TRUE)$pred[, 1]

    pred_Q_init <- if_else(panel_w$arm == "treat_tx", pred_Q_treat, pred_Q_ctrl)

    # =====================================================================
    # STEP 4: THE TARGETING STEP (Double Robust Fluctuation)
    # Goal: Fix the Q-model's extrapolation errors.
    # Theory: The logistic regression calculates how much (eps) to shift
    # the initial guesses based on the mechanism of missingness (H) across
    # three distinct biological/protocol states.
    # =====================================================================

    panel_w <- panel_w %>%
        mutate(
            Q_init_hazard = pred_Q_init %>% ptrim(),
            Q_init_logit = qlogis(Q_init_hazard),
            H_ctrl = if_else(arm == "control_no_tx", w_trim, 0),
            H_treat_pre = if_else(arm == "treat_tx" & t <= grace_period, w_trim, 0),
            H_treat_post = if_else(arm == "treat_tx" & t > grace_period, w_trim, 0)
        )

    m_fluctuate <- glm(event_cc_day ~ -1 + H_ctrl + H_treat_pre + H_treat_post + offset(Q_init_logit), family = binomial(), data = panel_w)

    eps_ctrl <- coef(m_fluctuate)["H_ctrl"]
    eps_treat_pre <- coef(m_fluctuate)["H_treat_pre"]
    eps_treat_post <- coef(m_fluctuate)["H_treat_post"]

    eps_ctrl <- ifelse(is.na(eps_ctrl), 0, eps_ctrl)
    eps_treat_pre <- ifelse(is.na(eps_treat_pre), 0, eps_treat_pre)
    eps_treat_post <- ifelse(is.na(eps_treat_post), 0, eps_treat_post)

    # =====================================================================
    # STEP 5: THE PERFECT UNIVERSE (Targeted G-computation)
    # Goal: Simulate the counterfactual where nobody breaks the rules.
    # Theory: We construct an artificial grid where censoring does not exist.
    # We predict initial hazards, calculate counterfactual H weights, and
    # apply the targeting update (Q*) to retrieve the unbiased hazard.
    # =====================================================================

    # complier_times <- clones %>%
    #     filter(arm == "treat_tx" & tx_ever == 1 & t_switch <= grace_period) %>%
    #     pull(t_switch)
    # if (length(complier_times) == 0) complier_times <- grace_period

    # cf_base <- data %>%
    #     select(id, tx_ever, t_switch, all_of(covars)) %>%
    #     mutate(sampled_t = sample(complier_times, n(), replace = TRUE))

    cf_data <- data %>%
        select(id, tx_ever, t_switch, all_of(covars)) %>%
        uncount(weights = tau, .id = "t") %>%
        uncount(weights = 2, .id = "clone_arm") %>%
        mutate(
            t_sq = t^2,
            post_grace = if_else(t > grace_period, 1, 0),
            arm = if_else(clone_arm == 1, "control_no_tx", "treat_tx"),
            cf_t_switch = case_when(
                arm == "control_no_tx" ~ 999,
                arm == "treat_tx" & tx_ever == 1 & t_switch <= grace_period ~ t_switch,
                # arm == "treat_tx" ~ sampled_t
                arm == "treat_tx" ~ grace_period # using grace_period over the sample distributions
            ),
            is_treated = if_else(t >= cf_t_switch, 1, 0)
        )

    X_cf <- cf_data %>%
        select(all_of(covars), t, t_sq, post_grace, is_treated) %>%
        as.data.frame()
    X_cf_treat_only <- cf_data %>%
        select(all_of(covars)) %>%
        as.data.frame()

    cf_pred_Q_ctrl <- predict(m_Q_ctrl, newdata = X_cf, onlySL = TRUE)$pred[, 1]
    cf_pred_Q_treat <- predict(m_Q_treat, newdata = X_cf, onlySL = TRUE)$pred[, 1]
    cf_pred_Q_init <- if_else(cf_data$arm == "treat_tx", cf_pred_Q_treat, cf_pred_Q_ctrl)

    cf_pred_ctrl_denom <- predict(m_denom_ctrl, newdata = X_cf, onlySL = TRUE)$pred[, 1]
    cf_pred_treat_denom <- predict(m_denom_treat, newdata = X_cf_treat_only, onlySL = TRUE)$pred[, 1]

    # =====================================================================
    # STEP 6: MARGINALIZATION
    # Goal: Calculate the final Population Average Treatment Effect (ATE).
    # Theory: Convert the targeted daily hazards (Q*) into individual
    # cumulative survival probabilities, then average across the population.
    # =====================================================================

    cf_data <- cf_data %>%
        mutate(
            cf_Q_init_hazard = cf_pred_Q_init %>% ptrim(),
            cf_Q_init_logit = qlogis(cf_Q_init_hazard),
            cf_p_ctrl_denom = cf_pred_ctrl_denom %>% ptrim(),
            cf_p_ctrl_num = predict(m_num_ctrl, newdata = ., type = "response") %>% ptrim(),
            cf_p_treat_denom = cf_pred_treat_denom %>% ptrim(),
            cf_p_treat_num = predict(m_num_treat, newdata = ., type = "response") %>% ptrim(),
            w_inc = case_when(
                arm == "control_no_tx" ~ (1 - cf_p_ctrl_num) / (1 - cf_p_ctrl_denom),
                arm == "treat_tx" ~ case_when(t == grace_period + 1 ~ cf_p_treat_num / cf_p_treat_denom, TRUE ~ 1)
            )
        ) %>%
        arrange(id, arm, t) %>%
        group_by(id, arm) %>%
        mutate(cf_H_raw = cumprod(w_inc)) %>%
        ungroup() %>%
        mutate(
            cf_H = cf_H_raw %>% ptrim(lower = q_lim[1], upper = q_lim[2]),
            cf_H_ctrl = if_else(arm == "control_no_tx", cf_H, 0),
            cf_H_treat_pre = if_else(arm == "treat_tx" & t <= grace_period, cf_H, 0),
            cf_H_treat_post = if_else(arm == "treat_tx" & t > grace_period, cf_H, 0),
            Q_star_hazard = plogis(cf_Q_init_logit + eps_ctrl * cf_H_ctrl + eps_treat_pre * cf_H_treat_pre + eps_treat_post * cf_H_treat_post)
        ) %>%
        group_by(id, arm) %>%
        mutate(indiv_survival = cumprod(1 - Q_star_hazard)) %>%
        ungroup()

    tmle_surv <- cf_data %>%
        group_by(arm, t) %>%
        summarise(survival = mean(indiv_survival), .groups = "drop")

    res_tmle <- tmle_surv %>%
        select(t, arm, survival) %>%
        pivot_wider(names_from = arm, values_from = survival) %>%
        mutate(rd = treat_tx - control_no_tx)

    # =====================================================================
    # STEP 7: THE EMPIRICAL INFLUENCE CURVE (Variance Estimation)
    # =====================================================================

    panel_w <- panel_w %>%
        mutate(
            Q_star_hazard = plogis(Q_init_logit + eps_ctrl * H_ctrl + eps_treat_pre * H_treat_pre + eps_treat_post * H_treat_post),
            D_martingale_step = w_trim * (event_cc_day - Q_star_hazard) / pmax(1 - Q_star_hazard, 1e-6)
        ) %>%
        arrange(clone_id, t) %>%
        group_by(clone_id) %>%
        mutate(cum_martingale = cumsum(D_martingale_step)) %>%
        ungroup()

    emp_martingale <- panel_w %>% select(id, arm, t, cum_martingale)

    cf_eic <- cf_data %>%
        select(id, arm, t, indiv_survival) %>%
        left_join(emp_martingale, by = c("id", "arm", "t")) %>%
        arrange(id, arm, t) %>%
        group_by(id, arm) %>%
        fill(cum_martingale, .direction = "down") %>%
        mutate(cum_martingale = replace_na(cum_martingale, 0)) %>%
        ungroup()

    marginal_surv <- tmle_surv %>% select(arm, t, psi = survival)

    eic_rd <- cf_eic %>%
        left_join(marginal_surv, by = c("arm", "t")) %>%
        mutate(IC_a = (indiv_survival - psi) - (indiv_survival * cum_martingale)) %>%
        select(id, arm, t, IC_a) %>%
        pivot_wider(names_from = arm, values_from = IC_a) %>%
        mutate(IC_RD = treat_tx - control_no_tx) %>%
        group_by(t) %>%
        summarise(
            var_IC_RD = var(IC_RD) / n(),
            se_RD = sqrt(var_IC_RD),
            .groups = "drop"
        )

    res_tmle <- res_tmle %>%
        left_join(eic_rd %>% select(t, se_RD), by = "t") %>%
        mutate(
            rd_ci_lower = rd - 1.96 * se_RD,
            rd_ci_upper = rd + 1.96 * se_RD
        )

    # =====================================================================
    # STEP 8: EXTRACT CCW RESULTS (For Comparison)
    # =====================================================================

    sum_log <- summary(fit_ccw, times = 1:tau, extend = TRUE)
    res_log <- tibble(
        t = sum_log$time,
        tx = sum_log$strata,
        sv = sum_log$surv,
        se = sum_log$std.err
    ) %>%
        mutate(tx = str_replace(tx, "arm=", "")) %>%
        pivot_wider(names_from = tx, values_from = c(sv, se)) %>%
        mutate(
            rd = sv_treat_tx - sv_control_no_tx,
            rd_se = sqrt(se_control_no_tx^2 + se_treat_tx^2),
        ) %>%
        mutate(
            ci_lower = rd - 1.96 * rd_se,
            ci_upper = rd + 1.96 * rd_se,
        )

    return(list(fit = fit_ccw, data = panel_w, res_log = res_log, res_tmle = res_tmle))
}

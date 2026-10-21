run_ccw_tmle_analysis <- function(data,
                             grace_period = 10,
                             tau = 30,
                             covars = c("x1"),
                             p_form = NULL, 
                             q_form = NULL, 
                             dt = 1,
                             use_sl = FALSE,
                             sl_lib = c("SL.mean", "SL.glm", "SL.earth.int", "SL.ranger")) {
    
    if (use_sl) {
        if (!require("SuperLearner", quietly = TRUE)) stop("Package 'SuperLearner' is required.")
        
        # Create the 2-way interaction MARS wrapper if it doesn't exist
        if ("SL.earth.int" %in% sl_lib && !exists("SL.earth.int", mode = "function")) {
            assign("SL.earth.int", function(..., degree = 2) {
                SuperLearner::SL.earth(..., degree = degree)
            }, envir = .GlobalEnv)
        }
    }

    if (is.null(p_form)) p_form <- paste(covars, collapse = " + ")
    if (is.null(q_form)) q_form <- paste(covars, collapse = " + ")

    # 1. Clone Creation
    data = data %>% 
        rename(time = event_time, tx_ever = tx_obs, t_switch = tx_time_obs)
    clones <- data %>%
        dplyr::select(id, tx_ever, t_switch, event, time, all_of(covars)) %>%
        tidyr::uncount(weights = 2, .id = "clone_arm") %>%
        mutate(
            arm = if_else(clone_arm == 1L, "control_no_tx", "treat_tx"),
            clone_id = paste0(id, "_", clone_arm)
        )

    # 2. STRICT Interval Alignment
    clones <- clones %>%
        mutate(
            dev_time = case_when(
                arm == "control_no_tx" & tx_ever == 1L ~ t_switch - 1,
                arm == "treat_tx" & !(tx_ever == 1L & t_switch <= grace_period) & time > grace_period ~ grace_period,
                TRUE ~ NA_real_
            ),
            time_cc = if_else(!is.na(dev_time), pmin(time, dev_time), time),
            time_cc = pmin(time_cc, tau),
            event_cc = if_else(event == 1L & time <= time_cc, 1L, 0L)
        )

    # 3. Fit IPCW Models 
    obs_panel <- clones %>% filter(clone_arm == 1L) %>% 
        mutate(
            tmax = if_else(tx_ever == 1L, t_switch, time),
            n_int = pmax(ceiling(tmax / dt), 1L)
        ) %>%
        tidyr::uncount(weights = n_int, .id = "t") %>%
        mutate(
            t_sq = t^2,
            event_tx = if_else(tx_ever == 1L & t == t_switch, 1L, 0L)
        )

    treat_day10_data <- clones %>%
        filter(arm == "treat_tx", time > grace_period) %>%
        mutate(took_tx = if_else(tx_ever == 1L & t_switch <= grace_period, 1L, 0L))

    if (use_sl) {
        X_ctrl <- obs_panel %>% dplyr::select(all_of(covars), t, t_sq) %>% as.data.frame()
        m_denom_ctrl <- SuperLearner(Y = obs_panel$event_tx, X = X_ctrl, family = binomial(), SL.library = sl_lib, cvControl = list(V = 2))
        
        X_treat <- treat_day10_data %>% dplyr::select(all_of(covars)) %>% as.data.frame()
        m_denom_treat <- SuperLearner(Y = treat_day10_data$took_tx, X = X_treat, family = binomial(), SL.library = sl_lib, cvControl = list(V = 2))
    } else {
        m_denom_ctrl <- suppressWarnings(glm(as.formula(paste("event_tx ~", p_form, "+ t + t_sq")), family = binomial(), data = obs_panel))
        m_denom_treat <- suppressWarnings(glm(as.formula(paste("took_tx ~", p_form)), family = binomial(), data = treat_day10_data))
    }

    m_num_ctrl   <- suppressWarnings(glm(event_tx ~ t + t_sq, family = binomial(), data = obs_panel))
    m_num_treat  <- suppressWarnings(glm(took_tx ~ 1, family = binomial(), data = treat_day10_data))

    # 4. Expand Clones & Apply Weights
    panel <- clones %>%
        mutate(n_int = pmax(ceiling(time_cc / dt), 0L)) %>% 
        tidyr::uncount(weights = n_int, .id = "t") %>%
        mutate(
            t0 = (t - 1L) * dt,
            t1 = t * dt,
            t_sq = t^2,
            event_cc_day = if_else(event_cc == 1L & time_cc > t0 & time_cc <= t1, 1L, 0L),
            is_treated = if_else(tx_ever == 1L & t >= t_switch, 1, 0)
        )
        
    if (use_sl) {
        X_pred_ctrl <- panel %>% dplyr::select(all_of(covars), t, t_sq) %>% as.data.frame()
        X_pred_treat <- panel %>% dplyr::select(all_of(covars)) %>% as.data.frame()
        
        pred_ctrl_denom <- predict(m_denom_ctrl, newdata = X_pred_ctrl, onlySL = TRUE)$pred[, 1]
        pred_treat_denom <- predict(m_denom_treat, newdata = X_pred_treat, onlySL = TRUE)$pred[, 1]
    } else {
        pred_ctrl_denom <- predict(m_denom_ctrl, newdata = panel, type = "response")
        pred_treat_denom <- predict(m_denom_treat, newdata = panel, type = "response")
    }

    panel <- panel %>%
        mutate(
            p_ctrl_denom = pmin(pmax(pred_ctrl_denom, 1e-4), 1 - 1e-4),
            p_ctrl_num   = pmin(pmax(predict(m_num_ctrl, newdata = ., type = "response"), 1e-4), 1 - 1e-4),
            p_treat_denom = pmin(pmax(pred_treat_denom, 1e-4), 1 - 1e-4),
            p_treat_num   = pmin(pmax(predict(m_num_treat, newdata = ., type = "response"), 1e-4), 1 - 1e-4)
        ) %>%
        mutate(
            w_inc = case_when(
                arm == "control_no_tx" ~ (1 - p_ctrl_num) / (1 - p_ctrl_denom),
                arm == "treat_tx" ~ case_when(
                    t == grace_period + 1L ~ p_treat_num / p_treat_denom,
                    TRUE ~ 1
                )
            )
        ) %>%
        arrange(clone_id, t) %>%
        group_by(clone_id) %>%
        mutate(w_ipcw_time = cumprod(w_inc)) %>%
        ungroup()

    q_lim <- quantile(panel$w_ipcw_time, probs = c(0.01, 0.99), na.rm = TRUE)
    
    panel_w <- panel %>% 
        mutate(
            w_trim = pmin(pmax(w_ipcw_time, q_lim[1]), q_lim[2]),
            post_grace = if_else(t > grace_period, 1, 0) 
        )

    fit_ccw <- survival::survfit(Surv(t0, t1, event_cc_day) ~ arm, data = panel_w, weights = w_trim, id = clone_id, robust=TRUE)

    # 5. TMLE STEP 1 & 2: Stratified Q-model and Targeting
    panel_ctrl <- panel_w %>% filter(arm == "control_no_tx")
    panel_treat <- panel_w %>% filter(arm == "treat_tx")

    if (use_sl) {
        X_Q_ctrl <- panel_ctrl %>% dplyr::select(all_of(covars), t, t_sq, post_grace, is_treated) %>% as.data.frame()
        X_Q_treat <- panel_treat %>% dplyr::select(all_of(covars), t, t_sq, post_grace, is_treated) %>% as.data.frame()
            
        m_Q_ctrl <- SuperLearner(Y = panel_ctrl$event_cc_day, X = X_Q_ctrl, family = binomial(), SL.library = sl_lib, cvControl = list(V = 2))
        m_Q_treat <- SuperLearner(Y = panel_treat$event_cc_day, X = X_Q_treat, family = binomial(), SL.library = sl_lib, cvControl = list(V = 2))
        
        X_Q_full <- panel_w %>% dplyr::select(all_of(covars), t, t_sq, post_grace, is_treated) %>% as.data.frame()
        pred_Q_ctrl <- predict(m_Q_ctrl, newdata = X_Q_full, onlySL = TRUE)$pred[, 1]
        pred_Q_treat <- predict(m_Q_treat, newdata = X_Q_full, onlySL = TRUE)$pred[, 1]
        
        pred_Q_init <- if_else(panel_w$arm == "treat_tx", pred_Q_treat, pred_Q_ctrl)
    } else {
        Q_formula <- as.formula(paste("event_cc_day ~", q_form, "+ t + t_sq + post_grace + is_treated"))
        
        m_Q_ctrl <- suppressWarnings(glm(Q_formula, family = binomial(), data = panel_ctrl))
        m_Q_treat <- suppressWarnings(glm(Q_formula, family = binomial(), data = panel_treat))
        
        pred_Q_ctrl <- predict(m_Q_ctrl, newdata = panel_w, type = "response")
        pred_Q_treat <- predict(m_Q_treat, newdata = panel_w, type = "response")
        
        pred_Q_init <- if_else(panel_w$arm == "treat_tx", pred_Q_treat, pred_Q_ctrl)
    }

    panel_w <- panel_w %>%
        mutate(
            Q_init_hazard = pmin(pmax(pred_Q_init, 1e-4), 1 - 1e-4),
            Q_init_logit = qlogis(Q_init_hazard),
            H_ctrl = if_else(arm == "control_no_tx", w_trim, 0),
            H_treat_pre = if_else(arm == "treat_tx" & t <= grace_period, w_trim, 0),
            H_treat_post = if_else(arm == "treat_tx" & t > grace_period, w_trim, 0)
        )

    m_fluctuate <- suppressWarnings(glm(event_cc_day ~ -1 + H_ctrl + H_treat_pre + H_treat_post + offset(Q_init_logit), family = binomial(), data = panel_w))
    eps_ctrl <- coef(m_fluctuate)["H_ctrl"]
    eps_treat_pre <- coef(m_fluctuate)["H_treat_pre"]
    eps_treat_post <- coef(m_fluctuate)["H_treat_post"]
    
    eps_ctrl <- ifelse(is.na(eps_ctrl), 0, eps_ctrl)
    eps_treat_pre <- ifelse(is.na(eps_treat_pre), 0, eps_treat_pre)
    eps_treat_post <- ifelse(is.na(eps_treat_post), 0, eps_treat_post)

    # 6. TMLE STEP 3: G-Computation / Marginalization
    # Get empirical distribution of complier switch times to mirror a random dynamic regime
    complier_times <- clones %>% 
        filter(arm == "treat_tx" & tx_ever == 1L & t_switch <= grace_period) %>% 
        pull(t_switch)
    if(length(complier_times) == 0) complier_times <- grace_period # Safety fallback
    
    cf_base <- data %>%
        dplyr::select(id, tx_ever, t_switch, all_of(covars)) %>%
        mutate(sampled_t = sample(complier_times, n(), replace = TRUE))
        
    cf_data <- cf_base %>%
        tidyr::uncount(weights = tau, .id = "t") %>%
        mutate(t_sq = t^2, post_grace = if_else(t > grace_period, 1, 0)) %>%
        tidyr::uncount(weights = 2, .id = "clone_arm") %>%
        mutate(
            arm = if_else(clone_arm == 1L, "control_no_tx", "treat_tx"),
            cf_t_switch = case_when(
                arm == "control_no_tx" ~ 999,
                arm == "treat_tx" & tx_ever == 1L & t_switch <= grace_period ~ t_switch,
                arm == "treat_tx" ~ sampled_t # THE FIX: Sample from empirical distribution!
            ),
            is_treated = if_else(t >= cf_t_switch, 1, 0)
        )

    if (use_sl) {
        X_cf <- cf_data %>% dplyr::select(all_of(covars), t, t_sq, post_grace, is_treated) %>% as.data.frame()
        X_cf_treat_only <- cf_data %>% dplyr::select(all_of(covars)) %>% as.data.frame()

        cf_pred_Q_ctrl <- predict(m_Q_ctrl, newdata = X_cf, onlySL = TRUE)$pred[, 1]
        cf_pred_Q_treat <- predict(m_Q_treat, newdata = X_cf, onlySL = TRUE)$pred[, 1]
        cf_pred_Q_init <- if_else(cf_data$arm == "treat_tx", cf_pred_Q_treat, cf_pred_Q_ctrl)

        cf_pred_ctrl_denom <- predict(m_denom_ctrl, newdata = X_cf, onlySL = TRUE)$pred[, 1]
        cf_pred_treat_denom <- predict(m_denom_treat, newdata = X_cf_treat_only, onlySL = TRUE)$pred[, 1]
    } else {
        cf_pred_Q_ctrl <- predict(m_Q_ctrl, newdata = cf_data, type = "response")
        cf_pred_Q_treat <- predict(m_Q_treat, newdata = cf_data, type = "response")
        cf_pred_Q_init <- if_else(cf_data$arm == "treat_tx", cf_pred_Q_treat, cf_pred_Q_ctrl)
        
        cf_pred_ctrl_denom <- predict(m_denom_ctrl, newdata = cf_data, type = "response")
        cf_pred_treat_denom <- predict(m_denom_treat, newdata = cf_data, type = "response")
    }

    cf_data <- cf_data %>%
        mutate(
            cf_Q_init_hazard = pmin(pmax(cf_pred_Q_init, 1e-4), 1 - 1e-4),
            cf_Q_init_logit  = qlogis(cf_Q_init_hazard),
            cf_p_ctrl_denom  = pmin(pmax(cf_pred_ctrl_denom, 1e-4), 1 - 1e-4),
            cf_p_ctrl_num    = pmin(pmax(predict(m_num_ctrl, newdata = ., type = "response"), 1e-4), 1 - 1e-4),
            cf_p_treat_denom = pmin(pmax(cf_pred_treat_denom, 1e-4), 1 - 1e-4),
            cf_p_treat_num   = pmin(pmax(predict(m_num_treat, newdata = ., type = "response"), 1e-4), 1 - 1e-4)
        ) %>%
        mutate(
            w_inc = case_when(
                arm == "control_no_tx" ~ (1 - cf_p_ctrl_num) / (1 - cf_p_ctrl_denom),
                arm == "treat_tx" ~ case_when(t == grace_period + 1L ~ cf_p_treat_num / cf_p_treat_denom, TRUE ~ 1)
            )
        ) %>%
        arrange(id, arm, t) %>%
        group_by(id, arm) %>%
        mutate(cf_H_raw = cumprod(w_inc)) %>%
        ungroup() %>%
        mutate(
            cf_H = pmin(pmax(cf_H_raw, q_lim[1]), q_lim[2]),
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
        dplyr::select(t, arm, survival) %>%
        pivot_wider(names_from = arm, values_from = survival) %>%
        mutate(rd = treat_tx - control_no_tx)

    # 7. TMLE STEP 4: The Empirical Influence Curve
    panel_w <- panel_w %>%
        mutate(
            Q_star_hazard = plogis(Q_init_logit + eps_ctrl * H_ctrl + eps_treat_pre * H_treat_pre + eps_treat_post * H_treat_post),
            D_martingale_step = w_trim * (event_cc_day - Q_star_hazard) / pmax(1 - Q_star_hazard, 1e-6)
        ) %>%
        arrange(clone_id, t) %>%
        group_by(clone_id) %>%
        mutate(cum_martingale = cumsum(D_martingale_step)) %>%
        ungroup()

    emp_martingale <- panel_w %>% dplyr::select(id, arm, t, cum_martingale)

    cf_eic <- cf_data %>%
        dplyr::select(id, arm, t, indiv_survival) %>%
        left_join(emp_martingale, by = c("id", "arm", "t")) %>%
        arrange(id, arm, t) %>%
        group_by(id, arm) %>%
        tidyr::fill(cum_martingale, .direction = "down") %>%
        mutate(cum_martingale = tidyr::replace_na(cum_martingale, 0)) %>%
        ungroup()

    marginal_surv <- tmle_surv %>% dplyr::select(arm, t, psi = survival)

    eic_rd <- cf_eic %>%
        left_join(marginal_surv, by = c("arm", "t")) %>%
        mutate(IC_a = (indiv_survival - psi) - (indiv_survival * cum_martingale)) %>%
        dplyr::select(id, arm, t, IC_a) %>%
        pivot_wider(names_from = arm, values_from = IC_a) %>%
        mutate(IC_RD = treat_tx - control_no_tx) %>%
        group_by(t) %>%
        summarise(
            var_IC_RD = var(IC_RD) / n(), 
            se_RD = sqrt(var_IC_RD),
            .groups = "drop"
        )

    res_tmle <- res_tmle %>%
        left_join(eic_rd %>% dplyr::select(t, se_RD), by = "t") %>%
        mutate(
            rd_ci_lower = rd - 1.96 * se_RD,
            rd_ci_upper = rd + 1.96 * se_RD
        )

    # 8. Extract CCW Results
    sum_log <- summary(fit_ccw, times = 1:tau, extend = TRUE)
    res_log <- tibble(
        t = sum_log$time, 
        tx = sum_log$strata, 
        sv = sum_log$surv,
        se = sum_log$std.err) %>%
        mutate(tx = str_replace(tx, "arm=", ""))  %>% 
        pivot_wider(names_from = tx, values_from = c(sv,se)) %>% 
        mutate(
            rd = sv_treat_tx - sv_control_no_tx,
            rd_se = sqrt(se_control_no_tx^2 + se_treat_tx^2),
        ) %>% 
        mutate(      # 95% Confidence Intervals
            ci_lower = rd - 1.96 * rd_se,
            ci_upper = rd + 1.96 * rd_se,
      )

    return(list(fit = fit_ccw, data = panel_w, res_log = res_log, res_tmle = res_tmle))
}

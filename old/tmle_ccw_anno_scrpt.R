library(dplyr)
library(tidyr)
library(readr)
library(stringr)
set.seed(100)
{
    TAU    <- 365.24   # 1-year horizon (in days): the time we report risks at
    GRACE  <- 182      # ~6-month grace period for the "surgery within 6 months" strategy
    EPS    <- 1e-6     # tiny number to keep probabilities away from exactly 0 or 1
    covars     <- c("age", "sex", "charlson", "perf", "stage", "emergency", "deprivation")
    cat_covars <- c("sex", "charlson", "perf", "stage", "emergency")


    dat = read_csv("ije-2019-08-1035-File008 1.csv")
    dat = dat %>% 
        mutate(
        id            = as.character(id),
        surgery       = as.integer(surgery),
        death         = as.integer(death),
        fup_obs       = as.numeric(fup_obs),
        timetosurgery = as.numeric(timetosurgery)
    )

    # --- 
    # Force categoricals to one-hots to ensure Superlearner can use them
    # --- 
    # install.packages("fastDummies")
    dat2 = fastDummies::dummy_cols(dat, cat_covars) %>% 
        select(-all_of(cat_covars)) %>% 
        select(-all_of(c("sex_1", "charlson_0", "perf_0","stage_1", "emergency_0"))) %>% 
        mutate(timetosurgery = ifelse(is.na(timetosurgery), Inf, timetosurgery)) 

    # ---
    # extract updated names
    # ---
    covars2 = c()
    for(n  in names(dat2)){
        if(str_replace(n, "_[0-9]", "") %in% covars){
            covars2 = append(covars2, n)
        }
    }

    # ---
    # set variables
    # ---
    data = dat2
    grace_period = GRACE
    tau = TAU
    covars= covars2
    dt = 1                          # time interval length, 7 = weeks, 30 = months ...
    event = "death"
    event_time = "fup_obs"
    tx_obs = "surgery"
    tx_time_obs = "timetosurgery"
    sl_lib = c("SL.mean", "SL.glm", "SL.earth.int", "SL.ranger")
        
    if (!require("SuperLearner", quietly = TRUE)) stop("Package 'SuperLearner' is required.")

    # ======================================================================
    # Set-up helper functions
    # Create the 2-way interaction MARS wrapper for the SuperLearner library
    # =====================================================================
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
    options(width=400)
    data <- data %>%
        rename(event = event, time = event_time, tx_ever = tx_obs, t_switch = tx_time_obs)

    clones <- data %>%
        select(id, tx_ever, t_switch, event, time, all_of(covars)) %>%
        uncount(weights = 2, .id = "clone_arm") %>%
        mutate(
            arm = if_else(clone_arm == 1, "control_no_tx", "treat_tx"),
            clone_id = paste0(id, "_", clone_arm)
        ) %>% 
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

    # ---
    # Isolate training datasets
    # ---
    # Create the control censoring model
    # Note that this is a ragged long-form dataset (different number of observed days per person)
    X_ctrl <- clones %>%
        filter(arm == "control_no_tx") %>%
        mutate(
            tmax = if_else(tx_ever == 1, t_switch, time),
            n_int = pmax(ceiling(tmax / dt), 1)
        ) %>%
        uncount(weights = n_int, .id = "t") %>%
        mutate(
            t_sq = t^2,
            event_tx = if_else(tx_ever == 1 & t == t_switch, 1, 0)
        ) 
        


    # NOTE: The treatment arm noncompliance can only occur at the grace period day.
    # We only need to observe risk of non-compliance at that single day.
    X_trt <- clones %>%
        filter(arm == "treat_tx", time > grace_period) %>%
        mutate(took_tx = if_else(tx_ever == 1 & t_switch <= grace_period, 1, 0))  


    # ---
    # Train Denominator g-models (SuperLearner)
    # ---
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

    # ---
    # Train Numerator models (Stabilization)
    # ---
    m_num_ctrl <- glm(event_tx ~ t + t_sq, family = binomial(), data = X_ctrl)
    m_num_treat <- glm(took_tx ~ 1, family = binomial(), data = X_trt)

    # ---
    # Expand Clones & Apply Weights
    # ---

    # ---
    # create frame to predict weights from
    # ---
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

    # ---
    # predict probability of censor
    # censoring weight is contructed of a scaling numerator and probability of censor denominator
    # ---
    pred_ctrl_denom <- predict(
        m_denom_ctrl, 
        newdata = panel %>% select(all_of(covars), t, t_sq) %>% as.data.frame(),  # ctrl prediction
        onlySL = TRUE)$pred[, 1] %>% 
        ptrim()

    pred_treat_denom <- predict(
        m_denom_treat, 
        newdata = panel %>% select(all_of(covars)) %>% as.data.frame() , 
        onlySL = TRUE)$pred[, 1] %>% 
        ptrim()

    p_ctrl_num = predict(m_num_ctrl, newdata = panel, type = "response") %>% ptrim()
    p_treat_num = predict(m_num_treat, newdata = panel, type = "response") %>% ptrim()

    # ---
    # combine all predictions and construct stable probabilities of non-compliance
    # ---
    panel <- panel %>%
        mutate(
            p_ctrl_denom = pred_ctrl_denom,
            p_ctrl_num = p_ctrl_num,
            p_treat_denom = pred_treat_denom,
            p_treat_num = p_treat_num,
            w_inc = case_when(
                arm == "control_no_tx" ~ (1 - p_ctrl_num) / (1 - p_ctrl_denom),
                arm == "treat_tx" ~ case_when(t == grace_period + 1 ~ p_treat_num / p_treat_denom, TRUE ~ 1)
            )
        ) %>%
        arrange(clone_id, t) %>%
        group_by(clone_id) %>%
        mutate(w_ipcw_time = cumprod(w_inc)) %>%
        ungroup()

    # ---
    # trim quantile to promote stability
    # ---
    q_lim <- quantile(panel$w_ipcw_time, probs = c(0.01, 0.99), na.rm = TRUE)
    panel_w <- panel %>%
        mutate(
            w_trim = w_ipcw_time %>% ptrim(lower = q_lim[1], upper = q_lim[2]),
            post_grace = if_else(t > grace_period, 1, 0)
        )

    # =====================================================================
    # STEP 3: Initial Q-models
    # Goal: Estimate baseline survival under strict compliance.
    # Theory: We switch the outcome to the true clinical event (Death).
    # These models only see person-time where clones perfectly obeyed the
    # rules. They are naive to the missing data bias caused by censoring.
    # =====================================================================

    panel_ctrl <- panel_w %>% filter(arm == "control_no_tx")
    panel_treat <- panel_w %>% filter(arm == "treat_tx") # outcome model is over the full follow-up-period

    X_Q_ctrl <- panel_ctrl %>%
        select(all_of(covars), t, t_sq, post_grace, is_treated) %>%
        as.data.frame()
    X_Q_treat <- panel_treat %>%
        select(all_of(covars), t, t_sq, post_grace, is_treated) %>%
        as.data.frame()

    # ---
    # train  two LS models
    # ---
    m_Q_ctrl <- SuperLearner(Y = panel_ctrl$event_cc_day, X = X_Q_ctrl, family = binomial(), SL.library = sl_lib, cvControl = list(V = 2))
    m_Q_treat <- SuperLearner(Y = panel_treat$event_cc_day, X = X_Q_treat, family = binomial(), SL.library = sl_lib, cvControl = list(V = 2))

    # ---
    # predict treated and control for each
    # ---
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

    m_fluctuate <- glm(
        event_cc_day ~ -1 + H_ctrl + H_treat_pre + H_treat_post + offset(Q_init_logit), 
        family = binomial(), 
        data = panel_w)

    eps_ctrl <- coef(m_fluctuate)["H_ctrl"]
    eps_treat_pre <- coef(m_fluctuate)["H_treast_pre"]
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

    # ---
    # full long form data
    # notice that the every patient is represented by the max number of days in the follow-up period
    # ---
    cf_data <- data %>%
        select(id, tx_ever, t_switch, all_of(covars)) %>%
        uncount(weights = floor(tau), .id = "t") %>%
        uncount(weights = 2, .id = "clone_arm") %>%
        mutate(
            t_sq = t^2,
            post_grace = if_else(t > grace_period, 1, 0),
            arm = if_else(clone_arm == 1, "control_no_tx", "treat_tx"),
            cf_t_switch = case_when(
                arm == "control_no_tx" ~ 999,
                arm == "treat_tx" & tx_ever == 1 & t_switch <= grace_period ~ t_switch,
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

    # ---
    # generate all predictions, for all possible data points
    # ----
    cf_pred_Q_ctrl <- predict(m_Q_ctrl, newdata = X_cf, onlySL = TRUE)$pred[, 1] %>% ptrim()
    cf_pred_Q_treat <- predict(m_Q_treat, newdata = X_cf, onlySL = TRUE)$pred[, 1] %>% ptrim()
    cf_pred_Q_init <- if_else(cf_data$arm == "treat_tx", cf_pred_Q_treat, cf_pred_Q_ctrl)

    cf_p_ctrl_denom <- predict(m_denom_ctrl, newdata = X_cf, onlySL = TRUE)$pred[, 1] %>% ptrim()
    cf_p_treat_denom <- predict(m_denom_treat, newdata = X_cf_treat_only, onlySL = TRUE)$pred[, 1] %>% ptrim()

    cf_p_ctrl_num = predict(m_num_ctrl, newdata = cf_data, type = "response") %>% ptrim()
    cf_p_treat_num = predict(m_num_treat, newdata = cf_data, type = "response") %>% ptrim()
    # =====================================================================
    # STEP 6: MARGINALIZATION
    # Goal: Calculate the final Population Average Treatment Effect (ATE).
    # Theory: Convert the targeted daily hazards (Q*) into individual
    # cumulative survival probabilities, then average across the population.
    # =====================================================================
    cf_data <- cf_data %>%
        mutate(
            cf_Q_init_hazard = cf_pred_Q_init,
            # ---
            # transform to logit for targeting step
            # ---
            cf_Q_init_logit = qlogis(cf_Q_init_hazard),
            cf_p_ctrl_denom,
            cf_p_ctrl_num,
            cf_p_treat_denom,
            cf_p_treat_num,
            w_inc = case_when(
                arm == "control_no_tx" ~ (1 - cf_p_ctrl_num) / (1 - cf_p_ctrl_denom),
                arm == "treat_tx" ~ case_when(t == grace_period + 1 ~ cf_p_treat_num / cf_p_treat_denom, TRUE ~ 1)
            )
        ) %>%
        arrange(id, arm, t) %>%
        group_by(id, arm) %>%
        mutate(cf_H_raw = cumprod(w_inc)) %>%
        ungroup() %>%
        # ---
        # updated predictions
        # ---
        mutate(
            cf_H = cf_H_raw %>% ptrim(lower = q_lim[1], upper = q_lim[2]),
            cf_H_ctrl = if_else(arm == "control_no_tx", cf_H, 0),
            cf_H_treat_pre = if_else(arm == "treat_tx" & t <= grace_period, cf_H, 0),
            cf_H_treat_post = if_else(arm == "treat_tx" & t > grace_period, cf_H, 0),
            Q_star_hazard = plogis(cf_Q_init_logit + eps_ctrl * cf_H_ctrl + eps_treat_pre * cf_H_treat_pre + eps_treat_post * cf_H_treat_post)
        ) %>%
        # ---
        # calculate survival per individual
        # ---
        group_by(id, arm) %>%
        mutate(indiv_survival = cumprod(1 - Q_star_hazard)) %>%
        ungroup()

    # ---
    # calculate marginal survival per arm
    # ---
    tmle_surv <- cf_data %>%
        group_by(arm, t) %>%
        summarise(survival = mean(indiv_survival), .groups = "drop")

    # ---
    # calculate marginal risk difference 
    # ---
    res_tmle <- tmle_surv %>%
        select(t, arm, survival) %>%
        pivot_wider(names_from = arm, values_from = survival) %>%
        mutate(rd = treat_tx - control_no_tx)



    # =====================================================================
    # STEP 7: THE EMPIRICAL INFLUENCE CURVE (Variance Estimation)
    # =====================================================================

    # ---
    # Generate g-model residual (martingale residuals)
    # calculate for observed data
    # ---
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

    # ---
    # attach martingale residuals to the full-length cf datasets
    # ---
    cf_eic <- cf_data %>%
        select(id, arm, t, indiv_survival) %>%
        left_join(emp_martingale, by = c("id", "arm", "t")) %>%
        arrange(id, arm, t) %>%
        group_by(id, arm) %>%
        fill(cum_martingale, .direction = "down") %>%
        mutate(cum_martingale = replace_na(cum_martingale, 0)) %>% # catch for immediate non-compliers
        ungroup()

    marginal_surv <- tmle_surv %>% select(arm, t, psi = survival)

    # --- 
    # Extract variance and Standardd Error from treatment - no treatment
    # Variance and Standard Error are calculated for each time independently 
    # ---
    eic_rd <- cf_eic %>%
        left_join(marginal_surv, by = c("arm", "t")) %>%
        # ---
        # calculate marginal survival influence curve value, per pat-arm-time
        # ---
        mutate(IC_a = (indiv_survival - psi) - (indiv_survival * cum_martingale)) %>%
        select(id, arm, t, IC_a) %>%
        pivot_wider(names_from = arm, values_from = IC_a) %>%
        # ---
        # calculate risk difference IC (simple difference of arm IC)
        # ---
        mutate(IC_RD = treat_tx - control_no_tx) %>%
        # variance and se calculated per each time
        group_by(t) %>%
        summarise(
            var_IC_RD = var(IC_RD) / n(),
            se_RD = sqrt(var_IC_RD),
            .groups = "drop"
        )



    # Generate Wald-Type CIs
    res_tmle <- res_tmle %>%
        left_join(eic_rd %>% select(t, se_RD), by = "t") %>%
        mutate(
            rd_ci_lower = rd - 1.96 * se_RD,
            rd_ci_upper = rd + 1.96 * se_RD
        )


    # =====================================================================
    # 1 year point estimate
    # =====================================================================
    cat("CCW + TMLE: 1 Year Risk Difference")
    print(res_tmle %>% filter(t == 365))
}

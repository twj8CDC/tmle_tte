library(tidyverse)
library(survival)
library(BART)
library(randomForest)
library(glmnet)
library(glue)

library(tidyverse)
library(stringr)


get_weibull_trajectory <- function(id, shape, s0, st, st2, t_switch, t_switch2, t_max, prefix) {
    u <- runif(1)
    target_h <- -log(pmax(u, 1e-10))

    # 1. Evaluate Event Time
    h_at_switch <- (t_switch / s0)^shape
    h_at_switch2 <- h_at_switch + (pmax(t_switch2^shape - t_switch^shape, 0)) / st^shape

    if (target_h <= h_at_switch) {
        event_time <- s0 * (target_h)^(1 / shape)
    } else if (target_h <= h_at_switch2) {
        rem_h <- target_h - h_at_switch
        event_time <- (rem_h * (st^shape) + (t_switch^shape))^(1 / shape)
    } else {
        rem_h <- target_h - h_at_switch2
        event_time <- (rem_h * (st2^shape) + (t_switch2^shape))^(1 / shape)
    }

    # 2. Vectorized Trajectory Calculation (NO MORE FOR LOOPS)
    days <- 1:t_max
    sv <- numeric(t_max)
    hz <- numeric(t_max)

    # Create boolean masks for each phase
    idx_p1 <- days < t_switch
    idx_p2 <- days >= t_switch & days < t_switch2
    idx_p3 <- days >= t_switch2

    # Vectorized Math for Phase 1
    if (any(idx_p1)) {
        d1 <- days[idx_p1]
        sv[idx_p1] <- exp(-(d1 / s0)^shape)
        hz[idx_p1] <- (shape / s0) * (d1 / s0)^(shape - 1)
    }

    # Vectorized Math for Phase 2
    if (any(idx_p2)) {
        d2 <- days[idx_p2]
        sv[idx_p2] <- exp(-((t_switch / s0)^shape + (pmax(d2^shape - t_switch^shape, 0) / st^shape)))
        hz[idx_p2] <- (shape / st) * (d2 / st)^(shape - 1)
    }

    # Vectorized Math for Phase 3
    if (any(idx_p3)) {
        d3 <- days[idx_p3]
        sv[idx_p3] <- exp(-((t_switch / s0)^shape +
            (pmax(t_switch2^shape - t_switch^shape, 0) / st^shape) +
            (pmax(d3^shape - t_switch2^shape, 0) / st2^shape)))
        hz[idx_p3] <- (shape / st2) * (d3 / st2)^(shape - 1)
    }

    # 3. Dynamic Column Naming using base R setNames
    list(
        summary = setNames(
            data.frame(
                id,
                min(ceiling(event_time), t_max),
                as.numeric(event_time <= t_max)
            ),
            c("id", paste0(prefix, "_event_time"), paste0(prefix, "_event"))
        ),
        trajectory = setNames(
            data.frame(id, days, sv, hz),
            c("id", "day", paste0(prefix, "_survival"), paste0(prefix, "_hazard"))
        )
    )
}
library(tidyr)
library(dplyr)

generate_survival_sim <- function(
  df_covs,
  t_max = 30,
  shape = 1,
  eq_tx_prob = "qnorm(0.6) + x1",
  eq_scale_0 = "t_max * exp(-(x1 * log(2)) / shape)",
  eq_scale_t = "t_max * exp(-(tx * log(0.5) + x1 * log(2)) / shape)",
  eq_scale_t2 = "t_max * exp(-(tx * log(.9) + x1 * log(2)) / shape)",
  alpha_beta = c(5, 3), # if alpha>beta then mass is right shifted
  max_tx_day = 10,
  min_supp_quant = 0.1,
  t2_time = 15
) {
    n <- nrow(df_covs)

    df_out = df_covs
    lp_tx <- with(df_covs, eval(parse(text = eq_tx_prob)))
    df_out["id"] = seq(1,nrow(df_covs))
    df_out["tx"] = rbinom(n, 1, pnorm(lp_tx))


    df_covs["tx"] <- 1 # just copy into df_cov to make it easier
    t_switch <- ceiling(rbeta(n, alpha_beta[1], alpha_beta[2]) * max_tx_day)
    t_switch2 <- t_switch + t2_time

    # set scale params
    s0 <- with(df_covs, eval(parse(text = eq_scale_0)))
    st <- with(df_covs, eval(parse(text = eq_scale_t)))
    st2 <- with(df_covs, eval(parse(text = eq_scale_t2)))

    # 2. Run Individual Trajectories (The Factual/Observed Data)
    # results <- map(1:n, ~ get_weibull_trajectory(
    #     .x, shape, s0[.x], st[.x], st2[.x],
    #     t_switch[.x], t_switch2[.x], t_max,
    #     "ntx"
    # ))

    # no one takes the drug ever
    traj_0 <- map(1:n, ~ get_weibull_trajectory(
        .x, shape, s0[.x], st[.x], st2[.x],
        t_switch = Inf, t_switch2 = Inf, t_max,
        "ntx"
    ))

    traj_dyn <- map(1:n, ~ get_weibull_trajectory(
        .x, shape, s0[.x], st[.x], st2[.x],
        t_switch[.x], t_switch2[.x], t_max,
        "dtx"
    ))

    # establish the minnimum time of tx that there is sufficient number of observations
    min_supp <- quantile(t_switch, min_supp_quant)
    traj_stat <- map(1:n, ~ get_weibull_trajectory(
        .x, shape, s0[.x], st[.x], st2[.x],
        min_supp, min_supp + t2_time, t_max,
        "stx"
    ))
    # print(bind_rows(traj_dyn))
    
    # ---------------------------------------------------------
    # 3. Consolidate Data (Extracting Summaries and Trajectories)
    # ---------------------------------------------------------
    
    # A. Extract and bind all daily trajectories
    traj_daily <- map_dfr(traj_0, "trajectory") %>%
        left_join(map_dfr(traj_dyn, "trajectory"), by = c("id", "day")) %>%
        left_join(map_dfr(traj_stat, "trajectory"), by = c("id", "day"))
        
    # B. Extract and bind all event summaries
    traj_summary <- map_dfr(traj_0, "summary") %>%
        left_join(map_dfr(traj_dyn, "summary"), by = "id") %>%
        left_join(map_dfr(traj_stat, "summary"), by = "id")
        
    # C. Combine into the "God Mode" dataset (true1)
    true1 <- traj_daily %>%
        left_join(traj_summary, by = "id") %>%
        left_join(df_out %>% select(id, tx), by = "id") %>% # Bring in assignment
        mutate(
            event = ifelse(tx == 1, dtx_event, ntx_event),
            event_time = ifelse(tx == 1, dtx_event_time, ntx_event_time)
        )
    
    # ---------------------------------------------------------
    # 4. Construct the Observational Dataset (df_out)
    # ---------------------------------------------------------
    df_out <- df_out %>% 
        left_join(
            true1 %>% distinct(id, event, event_time),
            by = "id"
        ) %>%
        mutate(tx_time_potential = t_switch) %>% 
        
        # The ITB Generator (Observational Missingness)
        mutate(tx_obs = case_when(
            tx == 0 ~ 0,
            tx == 1 & event_time < tx_time_potential ~ 0, 
            TRUE ~ 1
        )) %>% 
        mutate(tx_time_obs = ifelse(tx_obs == 1, tx_time_potential, Inf)) %>%
        select(-tx_time_potential)

    # observed data
    # true1 = bind_rows(traj_0) %>% 
    #     left_join(bind_rows(traj_dyn), c("id","day")) %>% 
    #     left_join(bind_rows(traj_stat) %>% select(-c(stx_event_time,stx_event)), c("id","day")) %>% 
    #     left_join(df_out, by="id") %>% 
    #     mutate(event = ifelse(tx == 1, dtx_event, ntx_event)) %>% 
    #     mutate(event_time = ifelse(tx == 1, dtx_event_time, ntx_event_time))
    
    # df_out = df_out %>% 
    #     left_join(
    #         true1 %>% 
    #         distinct(id, event, event_time),
    #         "id"
    #     ) %>%
    #     mutate(tx_time = t_switch) %>% 
    #     mutate(tx_time = ifelse(tx==0, Inf, tx_time)) %>% 
    #     mutate(tx_obs = case_when(
    #         tx == 0 ~ 0,                                  # Never assigned to drug
    #         tx == 1 & event_time < tx_time ~ 0,         # Died before they could take it
    #         TRUE ~ 1                                      # Survived and took it
    #     )) %>% 
    #     mutate(tx_time_obs = ifelse(tx_obs == 1, tx_time, Inf))
        
        

    true2 = true1 %>%
        group_by(day) %>%
        summarise(
        # --------------------------------------------------------
        # 1. The True Marginal Hazards (Weighted by Survival)
        # This perfectly tracks the derivative of the Marginal Survival Curve
        # --------------------------------------------------------
        marg_ntx_hz = sum(ntx_hazard * ntx_survival) / sum(ntx_survival),
        marg_dtx_hz = sum(dtx_hazard * dtx_survival) / sum(dtx_survival),
        marg_stx_hz = sum(stx_hazard * stx_survival) / sum(stx_survival),
        
        # The Target Estimands for CCW (Dynamic) and BART (Static)
        true_marg_dyn_hr  = marg_dtx_hz / marg_ntx_hz,
        true_marg_stat_hr = marg_stx_hz / marg_ntx_hz,

        # --------------------------------------------------------
        # 2. The True Conditional (Biological) Hazard Ratios
        # This is immune to survivor bias (Geometric mean of ratios)
        # --------------------------------------------------------
        cond_dyn_hr  = exp(mean(log(dtx_hazard / ntx_hazard))),
        cond_stat_hr = exp(mean(log(stx_hazard / ntx_hazard))),
        
        # --------------------------------------------------------
        # 3. Marginal Survival (Risk Difference Benchmarks)
        # --------------------------------------------------------
        marg_ntx_surv = mean(ntx_survival),
        marg_dtx_surv = mean(dtx_survival),
        marg_stx_surv = mean(stx_survival)
        ) %>%
        mutate(
            drd_surv = marg_dtx_surv-marg_ntx_surv, 
            srd_surv = marg_stx_surv-marg_ntx_surv
            )  %>% 
        mutate(min_supp = ifelse(day >= min_supp, 1,0)) %>% 
        ungroup()


    return(list(obs_data = df_out,sv_data=true1, est=true2, min_supp = min_supp))
}


#######################################
# Estimatorsd
cph_pred_survival <- function(fit, coh) {
    n <- nrow(coh)
    # L = min(coh$time)
    # 1. Create Counterfactuals (No tstart/tstop/id/evt needed!)
    pop_control <- coh %>% mutate(tx_ever = 0)
    pop_treated <- coh %>% mutate(tx_ever = 1)
    new_dat_std <- bind_rows(pop_control, pop_treated)
    # 2. Predict the Curves
    # survfit will naturally return a clean matrix here because there are no overlapping intervals
    cph_sv <- survfit(fit, newdata = new_dat_std)
    out <- tibble(
        time = cph_sv$time,
        ct = rowMeans(cph_sv$surv[,1:n]),
        tx = rowMeans(cph_sv$surv[,(n+1):(n*2)])
    ) %>%
        mutate(sd_diff = tx - ct)
    return(out)
}


cph_pred_survival_tdc <- function(dat, fit, min_support) {
    # Assuming 'fit' is your trained tdc model and 'dat' is sim_result$data
    n <- nrow(dat)
    t_max <- max(dat$time)
    min_support <- min_support

    # 1. Build the Control Counterfactual (tx_td = 0 always)
    # One interval per person: Day 0 to t_max
    pop_control <- dat %>%
        mutate(
            clone_id = id, # Keep original IDs for controls
            tstart = 0,
            tstop = t_max,
            tx_td = 0
        )

    # 2. Build the Treated Counterfactual (Respecting Positivity)
    # Two intervals per person: Day 0 to min_support (tx=0), then min_support to t_max (tx=1)
    pop_treated_p1 <- dat %>%
        mutate(
            clone_id = id + n, # Shift IDs so they don't overlap with controls
            tstart = 0,
            tstop = min_support,
            tx_td = 0 # Grace period: Untreated
        )

    pop_treated_p2 <- dat %>%
        mutate(
            clone_id = id + n, # Same shifted ID links this to p1
            tstart = min_support,
            tstop = t_max,
            tx_td = 1 # Post-positivity threshold: Treated
        )

    # Bind the treated intervals and sort them so survfit processes them chronologically
    pop_treated <- bind_rows(pop_treated_p1, pop_treated_p2) %>%
        arrange(clone_id, tstart)


    # 3. Combine both Counterfactual Arms
    new_dat_dynamic <- bind_rows(pop_control, pop_treated) %>%
        mutate(evt = 0)

    # 4. Predict the Dynamic Curves
    # The 'id = clone_id' argument is CRITICAL. It tells survfit to stitch pop_treated_p1 and p2 together.

    cph_sv <- survfit(fit, newdata = new_dat_dynamic, id = clone_id)
    out <- tibble(
        sv = cph_sv$surv,
        time = cph_sv$time,
        grp = rep(1:length(cph_sv$strata), cph_sv$strata)
    ) %>%
        filter(floor(time) == time) %>%
        mutate(tx = ifelse(grp > n, "tx", "ct")) %>%
        group_by(tx, time) %>%
        summarise(sv = mean(sv)) %>%
        pivot_wider(names_from = tx, values_from = sv) %>%
        mutate(sd_diff = tx - ct)
    return(out)
}


cph_time_match <- function(data, eq = "tx_ever + x1", caliper = .05) {
    # data = sim1$obs_data %>% 
    #     rename(time= event_time, tx_ever = tx_obs, t_switch = tx_time_obs)

    cph_tx <- data %>% filter(tx_ever == 1)
    cph_cntrl <- data %>% filter(tx == 0)
    cph_match_df <- tibble()
    for (i in 1:nrow(cph_tx)) {
        tx_obs <- cph_tx[i, ]
        # 1. Find the controls alive at the exact time the treated person took the drug
        cntrl_match <- cph_cntrl %>%
            filter(
                time >= tx_obs$t_switch &
                abs(x1 - tx_obs$x1) < caliper
            ) 
            
        if (nrow(cntrl_match) > 0) {
            cph_match_df <- cph_match_df %>%
                # 2. Shift the Treated patient's time
                bind_rows(
                    tx_obs %>%
                        mutate(
                            match_id = i,
                            time = time - t_switch  # Post-treatment survival only
                        )
                ) %>%
                # 3. Shift the Matched Controls' time
                bind_rows(
                    cntrl_match %>%
                        sample_n(pmin(nrow(cntrl_match), 10), replace = FALSE) %>%
                        mutate(
                            match_id = i,
                            time = time - tx_obs$t_switch # Shift by the treated person's t_switch!
                        )
                )
        }
    }

    cph_match_df <- cph_match_df %>%
        group_by(match_id) %>%
        mutate(weight = ifelse(
            tx_ever == 1,
            1,
            1 / sum(tx_ever == 0)
        ))
    # print(cph_match_df %>% ungroup %>% count(tx_ever))
    # print(
    #     cph_match_df %>% group_by(tx_ever) %>% reframe(quantile(x1,c(.1,.5,.9)))   # matched
    # )
    fit <- coxph(
        as.formula(paste0("Surv(time, event) ~", eq, "+ strata(match_id)")),
        data = cph_match_df,
        weights = weight,
        cluster = id, # CRITICAL: id is the original person's ID, not the match_id
        robust = TRUE,
        x=TRUE,
        y=TRUE
    )
    fit2 <- coxph(
        as.formula(paste0("Surv(time, event) ~", eq)),
        data = cph_match_df,
        weights = weight,
        cluster = id, # CRITICAL: id is the original person's ID, not the match_id
        robust = TRUE,
        x=TRUE,
        y=TRUE
    )
    return(
        list(result = exp(coef(fit)), coh = cph_match_df, fit = fit, fit2 = fit2)
    )
}

# --- 1. cox proportional hazards suite ---
# this function tests the "standard" ways people analyze this data
run_cph_suite <- function(sim_result, eq_cph = "tx_ever + x1", landmark_day = NULL, caliper=0.05) {
    dat <- sim_result$obs_data %>% 
        rename(time= event_time, tx_ever = tx_obs, t_switch = tx_time_obs)
    # use the simulation's min_support as the default landmark
    lm_day <- ifelse(is.null(landmark_day), sim_result$min_supp, landmark_day)
    
    # a. naive model: ignores timing, uses "ever treated" as a baseline flag
    outcome <- "Surv(time,event)~"
    fit_naive <- coxph(as.formula(paste0(outcome, eq_cph)), data = dat)
    hr_naive <- exp(coef(fit_naive))["tx_ever"]
    # print(glue("hr_naive {hr_naive}"))

    # b. landmark model: only looks at people who survived to 'lm_day'
    dat_lm <- dat %>% filter(time >= lm_day)
    outcome <- "Surv(time,event)~"
    fit_lm <- coxph(as.formula(paste0(outcome, eq_cph)), data = dat_lm)
    hr_lm <- exp(coef(fit_lm))["tx_ever"]
    
    # print(glue("hr_landmark {hr_lm}"))

    # c. time-dependent covariate (tdc): the "gold standard" for cox
    # re-structure data into start/stop intervals
    dat_long <- tmerge(dat, dat, id = id, tstop = time, evt = event(time, event))
    dat_long <- tmerge(dat_long, dat, id = id, tx_td = tdc(t_switch - .01)) # the eps here makes sure that if tt_tx =5 and tt_event=5, that the patient is counted as tx when the event occurs
    #   fit_tdc = coxph(Surv(tstart, tstop, evt) ~ tx_td + x1, data = dat_long)
    outcome <- "Surv(tstart, tstop, evt) ~ "
    eq_cph_tv <- str_replace_all(eq_cph, "tx_ever", "tx_td")
    fit_tdc <- coxph(as.formula(paste0(outcome, eq_cph_tv)), data = dat_long)
    hr_tdc <- exp(coef(fit_tdc))["tx_td"]

    # d. time-match
    cph_tm = cph_time_match(dat, eq = eq_cph, caliper =caliper)
    hr_tm = cph_tm$result["tx_ever"]
    



    # return a summary table of the estimated hazard ratios
    return(
        list(
            result = tibble(
                method = c("cph_naive", "cph_landmark", "cph_tdc", "cph_tm"),
                hr_estimate = c(hr_naive, hr_lm, hr_tdc, hr_tm)
            ),
            lnd_coh = dat_lm,
            tm_coh = cph_tm$coh,
            fit_lm = fit_lm,
            fit_tdc = fit_tdc,
            fit_tm = cph_tm$fit,
            fit_tm2 = cph_tm$fit2, 
            data = dat
        )
    )
}


# general functio nfor detting the rd
cph_sv_diff = function(cph, data, min_support){
    sd_df = tibble()
    cph_sv_lm = cph_pred_survival(cph$fit_lm,cph$lnd_coh)
    cph_lm_sd= cph_sv_lm %>% filter(time %in% c(5,10,15,20,25,30))
    sd_df = bind_rows(sd_df,
        tibble(method = "cph_lm", sd=cph_lm_sd)
    )

    cph_sv_lm = cph_pred_survival(cph$fit_tm2,cph$tm_coh)
    cph_lm_sd= cph_sv_lm %>% filter(time %in% c(5,10,15,20,25,30))
    sd_df = bind_rows(sd_df,
        tibble(method = "cph_tm", sd=cph_lm_sd)
    )

    cph_sv = cph_pred_survival_tdc(data, cph$fit_tdc, min_support)
    cph_tdc_sd = cph_sv %>% filter(time %in% c(5,10,15,20,25,30))
    sd_df = bind_rows(sd_df,
        tibble(method = "cph_tdc", sd=cph_tdc_sd)
    )
    return(sd_df)
}


# --- 2. bart survival function ---
# uses machine learning to estimate s(t | x, tx_td)
run_bart_survival <- function(sim_result, covars, n_trees = 50, n_cores = 1, trunc = 1, k_p = 1, keepevery = 10, nskip = 500, grace=10) {
    dat <- sim_result$obs_data %>%
        rename(time= event_time, tx_ever = tx_obs, t_switch = tx_time_obs) %>% 
        mutate(time = ceiling(time / trunc) * trunc)
        

    t_max <- max(dat$time)
    n_obs <- nrow(dat)
    times <- sort(unique(dat$time))
    k <- length(times)

    # 1. Build Person-Period Training Data
    # FIXED: Replaced 'n' with 'n_obs' and ordered explicitly by id, then t
    long <- tibble(id = rep(1:n_obs, each = k), t = rep(times, n_obs)) %>%
        left_join(dat, by = "id") %>%
        mutate(tx_td = ifelse(t_switch <= t, 1, 0)) %>%
        select(id, t, all_of(covars), tx_td, time, event)
    return(long)

    long_train <- long %>%
        filter(t <= time) %>%
        mutate(event_td = ifelse(time == t, event, 0))

    # 2. Train Model
    post <- BART::mc.surv.bart(
        x.train = data.frame(long_train) %>% select(t, tx_td, all_of(covars)),
        y.train = long_train$event_td,
        k = k_p,
        ntree = n_trees,
        keepevery = keepevery,
        nskip = nskip,
        mc.cores = n_cores
    )

    # 3. Create Counterfactual Datasets
    # min_support <- sim_result$min_support
    # z0_set <- times[times < min_support]
    # z1_set <- times[times >= min_support]


    # df_pred1 <- long %>%
    #     select(t, tx_td, all_of(covars)) %>%
    #     mutate(tx_td = rep(c(rep(0, length(z0_set)), rep(1, length(z1_set))), n_obs))

    # df_pred0 <- long %>%
    #     select(t, tx_td, all_of(covars)) %>%
    #     mutate(tx_td = 0)

    # 3. updated counterfactuals to match CCW dynamic tx
    
    
    # 4. Predict Hazards (prob.test)
    b_pred1 <- predict(post, data.frame(df_pred1))
    b_pred0 <- predict(post, data.frame(df_pred0))

    # 5. Extract Hazards and execute Manual G-Computation
    # Matrices are ndpost x (N*K). Transpose to (N*K) x ndpost
    haz1 <- t(b_pred1$prob.test)
    haz0 <- t(b_pred0$prob.test)

    draw_cols <- paste0("d_", 1:ncol(haz1))
    colnames(haz1) <- draw_cols
    colnames(haz0) <- draw_cols

    # FIXED: Manually calculate S(t) = cumprod(1 - h(t)) per person, per draw
    surv1_ind <- bind_cols(id = long$id, t = long$t, as_tibble(haz1)) %>%
        group_by(id) %>%
        mutate(across(starts_with("d_"), ~ cumprod(1 - .x))) %>%
        ungroup()

    surv0_ind <- bind_cols(id = long$id, t = long$t, as_tibble(haz0)) %>%
        group_by(id) %>%
        mutate(across(starts_with("d_"), ~ cumprod(1 - .x))) %>%
        ungroup()

    # Marginal Survival (Average across individuals for each time and draw)
    marg_s1 <- surv1_ind %>%
        select(-id) %>%
        group_by(t) %>%
        summarise(across(starts_with("d_"), mean))
    marg_s0 <- surv0_ind %>%
        select(-id) %>%
        group_by(t) %>%
        summarise(across(starts_with("d_"), mean))

    # Risk/Survival Difference
    sd_diff <- marg_s1 %>%
        select(starts_with("d_")) - marg_s0 %>%
        select(starts_with("d_"))
    
    m1 = marg_s1 %>% 
        pivot_longer(starts_with("d_"), names_to = "draw", values_to = "m1") %>% 
        group_by(t) %>% 
        summarise(m1_est = mean(m1))  
    m0 = marg_s0 %>% 
        pivot_longer(starts_with("d_"), names_to = "draw", values_to = "m0") %>% 
        group_by(t) %>% 
        summarise(m0_est = mean(m0))  
    

    sd_res <- bind_cols(t = marg_s1$t, sd_diff) %>%
        pivot_longer(starts_with("d_"), names_to = "draw", values_to = "sd") %>%
        group_by(t) %>%
        summarise(
            sd_estimate = mean(sd),
            conf_low = quantile(sd, 0.025),
            conf_high = quantile(sd, 0.975)
        ) %>%
        mutate(
            m1=m1$m1_est, 
            m0 = m0$m0_est
        )

    # FIXED: Marginal Hazard Ratio (Ratio of the Marginal Hazards to prevent convexity bias)
    marg_h1 <- bind_cols(id = long$id, t = long$t, as_tibble(haz1)) %>%
        select(-id) %>%
        group_by(t) %>%
        summarise(across(starts_with("d_"), mean))
    marg_h0 <- bind_cols(id = long$id, t = long$t, as_tibble(haz0)) %>%
        select(-id) %>%
        group_by(t) %>%
        summarise(across(starts_with("d_"), mean))

    hr_ratio <- marg_h1 %>% select(starts_with("d_")) / marg_h0 %>% select(starts_with("d_"))
    hr_res <- bind_cols(t = marg_h1$t, hr_ratio) %>%
        filter(t >= min_support) %>%
        pivot_longer(starts_with("d_"), names_to = "draw", values_to = "hr") %>%
        group_by(t) %>%
        summarise(
            hr_estimate = mean(hr),
            conf_low = quantile(hr, 0.025),
            conf_high = quantile(hr, 0.975)
        )


    # alt
    marg_h1 <- bind_cols(id = long$id, t = long$t, as_tibble(haz1)) %>%
        select(-id) #%>%
        # group_by(t) 
        # summarise(across(starts_with("d_"), mean))
    marg_h0 <- bind_cols(id = long$id, t = long$t, as_tibble(haz0)) %>%
        select(-id) #%>%
        # group_by(t) 
        # summarise(across(starts_with("d_"), mean))
    
    hr_ratio <- log((marg_h1 %>% ungroup %>%  select(starts_with("d_")) / marg_h0 %>% ungroup %>% select(starts_with("d_")))) %>% 
        bind_cols(t = long$t)   %>% 
        group_by(t) %>% 
        summarise(across(starts_with("d_"), mean)) 
        

    hr_res2 <- hr_ratio %>% 
        filter(t >= min_support) %>%
        pivot_longer(starts_with("d_"), names_to = "draw", values_to = "hr") %>%
        group_by(t) %>%
        summarise(
            hr_mean= mean(exp(hr)),
            hr_median = quantile(exp(hr), .5),
            conf_low = quantile(exp(hr), 0.025),
            conf_high = quantile(exp(hr), 0.975)
        )

    return(list(hr = hr_res, hr2 = hr_res2, rd = sd_res, model = post))
}




cox_ipcw_by_arm <- function(d_arm, covars, grace_period, EPS = 1e-6) {
    # Identify deviations occurring within the grace period for this specific arm
    d_arm <- d_arm %>%
        mutate(
            t_limit = pmin(grace_period, time_cc),
            # A deviation event occurred if dev was 1 and happened at or before t_limit
            dev_event = if_else(dev == 1L & !is.na(dev_time) & dev_time <= t_limit, 1L, 0L)
        )
    # Fit Cox model to the probability of STAYING in the protocol (not deviating)
    # We use the baseline covariates to model the "hazard" of deviating
    f <- as.formula(paste0("Surv(t_limit, dev_event) ~ ", paste(covars, collapse = " + ")))
    fit <- coxph(f, data = d_arm, ties = "efron", x = TRUE)
    # Baseline cumulative hazard H0(t) for deviation
    bh <- basehaz(fit, centered = FALSE) %>%
        as_tibble() %>%
        rename(H0 = hazard, time = time)
    # Linear predictor (risk score) for each individual's deviation risk
    lp <- predict(fit, type = "lp")
    idx <- findInterval(d_arm$t_limit, bh$time)
    H0_t <- ifelse(idx == 0, 0, bh$H0[idx])
    # Probability of remaining uncensored (adherent) = exp(-Cumulative Hazard)
    # This uses the relationship S(t) = exp(-H(t))
    P_uncens <- exp(-H0_t * exp(lp))
    # The IPCW weight is the inverse of the probability of remaining in the protocol
    w_cox <- 1 / pmax(P_uncens, EPS)
    d_arm %>%
        mutate(
            w_cox = w_cox,
            P_uncens_cox = P_uncens
        )
}



run_tmle_survival <- function(cloned_w,
                              covars = c("x1"),
                              tau = 30) {
    # 1. Prep data for TMLE at a specific time point (TAU)
    # We filter to those who either had the event or reached the end of follow-up
    dat_tmle <- cloned_w %>%
        filter(event_cc == 1 | time_cc >= tau) %>%
        mutate(
            Y = event_cc,
            A = if_else(arm == "treat_tx", 1, 0)
        )

    # 2. Fit Ensemble (Outcome Model Q)
    # Basic GLM
    f_glm <- as.formula(paste("Y ~ A +", paste(covars, collapse = " + ")))
    m_glm <- glm(f_glm, data = dat_tmle, weights = w_trim, family = binomial)

    # Random Forest
    print(tibble(dat_tmle))
    print(unique(dat_tmle$Y))
    m_rf <- randomForest::randomForest(
        x = dat_tmle[, c("A", covars)],
        y = as.factor(dat_tmle$Y),
        ntree = 500
    )

    # Elastic Net
    X_mat_full <- model.matrix(as.formula(paste("~ A +", paste(covars, collapse = "+"))), data = dat_tmle)[, -1]
    m_enet <- glmnet::cv.glmnet(X_mat_full, dat_tmle$Y, weights = dat_tmle$w_trim, family = "binomial", alpha = 0.5)

    # 3. Initial Predictions (Q0)
    dat_A1 <- dat_tmle %>% mutate(A = 1)
    dat_A0 <- dat_tmle %>% mutate(A = 0)

    X_new_1 <- model.matrix(as.formula(paste("~ A +", paste(covars, collapse = "+"))), data = dat_A1)[, -1]
    X_new_0 <- model.matrix(as.formula(paste("~ A +", paste(covars, collapse = "+"))), data = dat_A0)[, -1]

    p_glm_obs <- predict(m_glm, newdata = dat_tmle, type = "response")
    p_glm_1 <- predict(m_glm, newdata = dat_A1, type = "response")
    p_glm_0 <- predict(m_glm, newdata = dat_A0, type = "response")

    p_rf_obs <- predict(m_rf, newdata = dat_tmle[, c("A", covars)], type = "prob")[, 2]
    p_rf_1 <- predict(m_rf, newdata = dat_A1[, c("A", covars)], type = "prob")[, 2]
    p_rf_0 <- predict(m_rf, newdata = dat_A0[, c("A", covars)], type = "prob")[, 2]

    p_enet_obs <- predict(m_enet, newx = X_mat_full, type = "response", s = "lambda.min")[, 1]
    p_enet_1 <- predict(m_enet, newx = X_new_1, type = "response", s = "lambda.min")[, 1]
    p_enet_0 <- predict(m_enet, newx = X_new_0, type = "response", s = "lambda.min")[, 1]

    # Ensemble Average & Bounding
    Q0_obs <- pmin(pmax((p_glm_obs + p_rf_obs + p_enet_obs) / 3, 0.001), 0.999)
    Q0_1 <- pmin(pmax((p_glm_1 + p_rf_1 + p_enet_1) / 3, 0.001), 0.999)
    Q0_0 <- pmin(pmax((p_glm_0 + p_rf_0 + p_enet_0) / 3, 0.001), 0.999)

    # 4. Targeted Update (Propensity g = 0.5 due to cloning)
    g <- 0.5
    dat_tmle$H <- ifelse(dat_tmle$A == 1, 1 / g, -1 / (1 - g))

    ##### potential update
    # Instead of g <- 0.5, estimate the weighted propensity score
    # m_g <- glm(A ~ x1, data = dat_tmle, weights = w_trim, family = binomial)
    # g1_obs <- predict(m_g, type = "response")
    # Update H to use the individual-specific probabilities
    # dat_tmle$H <- ifelse(dat_tmle$A == 1, 1 / g1_obs, -1 / (1 - g1_obs))

    m_update <- glm(Y ~ -1 + H, offset = qlogis(Q0_obs), weights = w_trim, family = binomial, data = dat_tmle)
    epsilon <- coef(m_update)["H"]

    # Final Estimates (Star)
    Q1_star <- plogis(qlogis(Q0_1) + epsilon * (1 / g))
    Q0_star <- plogis(qlogis(Q0_0) + epsilon * (-1 / (1 - g)))

    tmle_risk_1 <- weighted.mean(Q1_star, w = dat_tmle$w_trim)
    tmle_risk_0 <- weighted.mean(Q0_star, w = dat_tmle$w_trim)
    rd_estimate <- tmle_risk_1 - tmle_risk_0

    # 5. Inference (Influence Curve)
    # Note: Standard TMLE variance uses H * (Y - Q_star)
    Q_star_obs <- ifelse(dat_tmle$A == 1, Q1_star, Q0_star)
    ic_val <- (dat_tmle$H * (dat_tmle$Y - Q_star_obs)) + (Q1_star - Q0_star - rd_estimate)

    n <- nrow(dat_tmle)
    var_ic <- var(ic_val * dat_tmle$w_trim)
    se_rd <- sqrt(var_ic / n)

    return(tibble::tibble(
        Statistic = c("RD", "SE", "LL", "UL", "p-val", "sd"),
        Value = c(
            rd_estimate,
            se_rd,
            rd_estimate - 1.96 * se_rd,
            rd_estimate + 1.96 * se_rd,
            2 * (1 - pnorm(abs(rd_estimate / se_rd))),
            -1 * rd_estimate
        )
    ))
}

cox_ipcw_by_arm <- function(d_arm, covars, grace_period, EPS = 1e-6) {
    # 1. Identify deviations
    d_arm <- d_arm %>%
        mutate(
            t_limit = pmin(grace_period, time_cc),
            dev_event = if_else(dev == 1L & !is.na(dev_time) & dev_time <= t_limit, 1L, 0L)
        )
    
    # 2. SAFEGUARD: If no deviations occurred, everyone has a 100% chance of staying in protocol
    if (sum(d_arm$dev_event) == 0) {
        return(d_arm %>% mutate(w_cox = 1, P_uncens_cox = 1))
    }

    # 3. Fit the Cox Model
    f <- as.formula(paste0("Surv(t_limit, dev_event) ~ ", paste(covars, collapse = " + ")))
    fit <- coxph(f, data = d_arm, ties = "efron", x = TRUE)
    
    # 4. Extract Expected Cumulative Hazard
    # type = "expected" automatically evaluates the cumulative hazard for each person up to their t_limit
    cum_haz <- predict(fit, type = "expected")
    
    # 5. Convert to Survival Probability S(t) = exp(-H(t))
    P_uncens <- exp(-cum_haz)
    
    # 6. Apply Weights
    d_arm %>%
        mutate(
            w_cox = 1 / pmax(P_uncens, EPS),
            P_uncens_cox = P_uncens
        )
}

run_ccw_analysis <- function(data,
                             grace_period = 10,
                             tau = 30,
                             covars = c("x1"),
                             dt = 1) {
    
    # 1. Clone Creation
    clones <- data %>%
        rename(time = event_time, tx_ever = tx_obs, t_switch = tx_time_obs) %>%
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

    # ---------------------------------------------------------
    # 4. Fit IPCW Models (The Paradigm Shift)
    # ---------------------------------------------------------
    
    # MODEL A: Control Arm (Daily Hazard of Deviation)
    obs_panel <- clones %>% filter(clone_arm == 1L) %>% 
        mutate(
            tmax = if_else(tx_ever == 1L, t_switch, time),
            tmax = pmin(tmax, grace_period),
            n_int = pmax(ceiling(tmax / dt), 1L)
        ) %>%
        tidyr::uncount(weights = n_int, .id = "t") %>%
        mutate(
            t_sq = t^2,
            event_tx = if_else(tx_ever == 1L & t == t_switch, 1L, 0L)
        )

    m_denom_ctrl <- glm(as.formula(paste("event_tx ~", paste(covars, collapse = " + "), "+ t + t_sq")), family = binomial(), data = obs_panel)
    m_num_ctrl   <- glm(event_tx ~ t + t_sq, family = binomial(), data = obs_panel)

    # MODEL B: Treat Arm (Single Probability of Complying by Day 10)
    # We only calculate this for clones who actually survive past the grace period to face the censoring event!
    treat_day10_data <- clones %>%
        filter(arm == "treat_tx", time > grace_period) %>%
        mutate(took_tx = if_else(tx_ever == 1L & t_switch <= grace_period, 1L, 0L))

    m_denom_treat <- glm(as.formula(paste("took_tx ~", paste(covars, collapse = " + "))), family = binomial(), data = treat_day10_data)
    m_num_treat   <- glm(took_tx ~ 1, family = binomial(), data = treat_day10_data)

    # ---------------------------------------------------------
    # 5. Expand Clones & Apply Mathematical Compliance Paths
    # ---------------------------------------------------------
    panel <- clones %>%
        mutate(n_int = pmax(ceiling(time_cc / dt), 0L)) %>% 
        tidyr::uncount(weights = n_int, .id = "t") %>%
        mutate(
            t0 = (t - 1L) * dt,
            t1 = t * dt,
            event_cc_day = if_else(event_cc == 1L & time_cc > t0 & time_cc <= t1, 1L, 0L),
            t_sq = t^2
        ) %>%
        mutate(
            # Predict daily probabilities for Control
            p_ctrl_denom = pmin(pmax(predict(m_denom_ctrl, newdata = ., type = "response"), 1e-4), 1 - 1e-4),
            p_ctrl_num   = pmin(pmax(predict(m_num_ctrl, newdata = ., type = "response"), 1e-4), 1 - 1e-4),
            
            # Predict single baseline probability for Treat
            p_treat_denom = pmin(pmax(predict(m_denom_treat, newdata = ., type = "response"), 1e-4), 1 - 1e-4),
            p_treat_num   = pmin(pmax(predict(m_num_treat, newdata = ., type = "response"), 1e-4), 1 - 1e-4)
        ) %>%
        mutate(
            w_inc = case_when(
                arm == "control_no_tx" ~ case_when(
                    t > grace_period ~ 1, # Protocol locked
                    TRUE ~ (1 - p_ctrl_num) / (1 - p_ctrl_denom)
                ),
                arm == "treat_tx" ~ case_when(
                    # THE FIX: Apply the weight update EXACTLY once, on the first day they are missing the non-takers
                    t == grace_period + 1L ~ p_treat_num / p_treat_denom,
                    TRUE ~ 1
                )
            )
        ) %>%
        arrange(clone_id, t) %>%
        group_by(clone_id) %>%
        mutate(w_ipcw_time = cumprod(w_inc)) %>%
        ungroup()

    # 6. Trim and Fit
    q_lim <- quantile(panel$w_ipcw_time, probs = c(0.01, 0.99), na.rm = TRUE)
    panel_w <- panel %>% mutate(w_trim = pmin(pmax(w_ipcw_time, q_lim[1]), q_lim[2]))

    fit_ccw <- survival::survfit(
        Surv(t0, t1, event_cc_day) ~ arm,
        data = panel_w,
        weights = w_trim,
        id = clone_id
    )

    # 7. Extract Results
    sum_log <- summary(fit_ccw, times = 1:tau, extend = TRUE)
    res_log <- tibble(
        t = sum_log$time,
        tx = sum_log$strata,
        sv = sum_log$surv
    ) %>%
        pivot_wider(names_from = tx, values_from = sv) %>%
        mutate(rd = `arm=treat_tx` - `arm=control_no_tx`)

    return(list(fit = fit_ccw, data = panel_w, res_log = res_log))
}

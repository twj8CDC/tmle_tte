library(tidyverse)
library(survival)
library(BART)
library(randomForest)
library(glmnet)
library(glue)

library(tidyverse)
library(stringr)

get_weibull_trajectory <- function(id, shape, s0, st, st2, t_switch, t_switch2, t_max) {
    u <- runif(1)
    target_h <- -log(pmax(u, 1e-10))

    # FIXED: Added pmax protection for the interval differences
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

    days <- 1:t_max

    # FIXED: S(t) calculation replaced append() with pre-allocated vector for stability
    sv <- numeric(t_max)
    for (d in days) {
        if (d < t_switch) {
            sv[d] <- exp(-(d / s0)^shape)
        } else if (d < t_switch2) {
            sv[d] <- exp(-((t_switch / s0)^shape + (pmax(d^shape - t_switch^shape, 0) / st^shape)))
        } else {
            sv[d] <- exp(-((t_switch / s0)^shape +
                (pmax(t_switch2^shape - t_switch^shape, 0) / st^shape) +
                (pmax(d^shape - t_switch2^shape, 0) / st2^shape)))
        }
    }

    list(
        summary = data.frame(
            id = id,
            time = min(ceiling(event_time), t_max),
            event = as.numeric(event_time <= t_max)
        ),
        trajectory = data.frame(id = id, day = days, survival = sv)
    )
}

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

    # 1. Evaluate Assignment and Scales
    lp_tx <- with(df_covs, eval(parse(text = eq_tx_prob)))
    tx_status <- rbinom(n, 1, pnorm(lp_tx))
    df_covs["tx"] <- tx_status # just copy into df_cov to make it easier
    t_switch <- ceiling(rbeta(n, alpha_beta[1], alpha_beta[2]) * max_tx_day)
    t_switch[tx_status == 0] <- Inf
    t_switch2 <- t_switch + t2_time

    s0 <- with(df_covs, eval(parse(text = eq_scale_0)))
    st <- with(df_covs, eval(parse(text = eq_scale_t)))
    st2 <- with(df_covs, eval(parse(text = eq_scale_t2)))

    # FIXED: Replaced str_replace with str_replace_all and word boundaries (\\b)
    # This prevents 'tx_coef' from turning into '1_coef' if passed in the string
    st_cf <- with(df_covs, eval(parse(text = str_replace_all(eq_scale_t, "\\btx\\b", "1"))))
    st_cf2 <- with(df_covs, eval(parse(text = str_replace_all(eq_scale_t2, "\\btx\\b", "1"))))

    # 2. Run Individual Trajectories (The Factual/Observed Data)
    results <- map(1:n, ~ get_weibull_trajectory(
        .x, shape, s0[.x], st[.x], st2[.x],
        t_switch[.x], t_switch2[.x], t_max
    ))

    # 3. Calculate Benchmarks (The "Truth" for Comparison)
    min_supp <- quantile(t_switch[t_switch != Inf], min_supp_quant)

    truth <- map_dfr(1:t_max, function(curr_t) {
        # A. E(S | T = t_obs)
        s_obs_all <- map_dbl(results, ~ .x$trajectory$survival[curr_t])
        s_obs_tx <- mean(s_obs_all[tx_status == 1])
        s_obs_ct <- mean(s_obs_all[tx_status == 0])

        # B. E(S | T = T^1) and E(S | T = T^0)
        s_cf0 <- exp(-(curr_t / s0)^shape)
        hz_ct <- (shape / s0) * (curr_t / s0)^(shape - 1)

        if (curr_t < min_supp) {
            s_cf1 <- exp(-(curr_t / s0)^shape)
            hz_tx <- hz_ct
        } else if (curr_t < min_supp + t2_time) {
            s_cf1 <- exp(-((min_supp / s0)^shape + (curr_t^shape - min_supp^shape) / st_cf^shape))
            hz_tx <- (shape / st_cf) * (curr_t / st_cf)^(shape - 1)
        } else {
            s_cf1 <- exp(-(
                (min_supp / s0)^shape +
                    ((min_supp + t2_time)^shape - min_supp^shape) / st_cf^shape +
                    (curr_t^shape - (min_supp + t2_time)^shape) / st_cf2^shape
            ))
            hz_tx <- (shape / st_cf2) * (curr_t / st_cf2)^(shape - 1)
        }

        tibble(
            t = curr_t,
            s_observed = mean(s_obs_all),
            s_observed_tx = s_obs_tx,
            s_observed_ct = s_obs_ct,
            s_cf_tx = mean(s_cf1),
            s_cf_ct = mean(s_cf0),
            true_hr_t = mean(hz_tx) / mean(hz_ct),
            true_hr_t2 = exp(mean(log(hz_tx / hz_ct)))
        )
    })

    # 4. Final Output
    return(list(
        data = bind_rows(map(results, "summary")) %>%
            bind_cols(df_covs) %>%
            mutate(t_switch_sim = t_switch) %>%
            mutate(tx_status_sim = tx_status) %>%
            mutate(t_switch = ifelse(t_switch_sim <= time, t_switch_sim, Inf)) %>%
            mutate(tx_ever = as.numeric(t_switch <= time)) %>%
            mutate(h0 = 1 / s0, h1 = 1 / st, h21 = 1 / st2),
        truth = truth,
        min_support = min_supp
    ))
}

# --- 1. Core Weibull Logic (Internal) ---
# get_weibull_trajectory <- function(id, shape, s0, st, st2, t_switch, t_switch2, t_max) {
#     u <- runif(1)
#     target_h <- -log(pmax(u, 1e-10))
#     h_at_switch <- (t_switch / s0)^shape
#     h_at_switch2 <- h_at_switch + (t_switch2^shape - t_switch^shape) / st^shape

#     if (target_h <= h_at_switch) {
#         event_time <- s0 * (target_h)^(1 / shape)
#     } else if (target_h <= h_at_switch2) {
#         rem_h <- target_h - h_at_switch
#         event_time <- (rem_h * (st^shape) + (t_switch^shape))^(1 / shape)
#     } else {
#         rem_h <- target_h - h_at_switch2
#         event_time <- (rem_h * (st2^shape) + (t_switch2^shape))^(1 / shape)
#     }

#     days <- 1:t_max
#     # S(t) calculation
#     sv <- c()
#     for (d in days) {
#         sv <- if (d < t_switch) {
#             sv <- append(
#                 sv,
#                 exp(-(d / s0)^shape)
#             )
#         } else if (d < t_switch2) {
#             sv <- append(
#                 sv,
#                 exp(-((t_switch / s0)^shape + (pmax(d^shape - t_switch^shape, 0) / st^shape)))
#             )
#         } else {
#             sv <- append(
#                 sv,
#                 exp(-((t_switch / s0)^shape +
#                     (pmax(t_switch2^shape - t_switch^shape, 0) / st^shape) +
#                     (pmax(d^shape - t_switch2^shape, 0) / st2^shape)
#                 ))
#             )
#         }
#     }

#     list(
#         summary = data.frame(
#             id = id, time = min(ceiling(event_time), t_max),
#             event = as.numeric(event_time <= t_max)
#         ),
#         trajectory = data.frame(id = id, day = days, survival = sv)
#     )
# }

# generate_survival_sim <- function(
#   df_covs,
#   t_max = 30,
#   shape = 1,
#   eq_tx_prob = "qnorm(0.6) + x1",
#   eq_scale_0 = "t_max * exp(-(x1 * log(2)) / shape)",
#   eq_scale_t = "t_max * exp(-(tx * log(0.5) + x1 * log(2)) / shape)",
#   eq_scale_t2 = "t_max * exp(-(tx * log(.9) + x1 * log(2)) / shape)",
#   alpha_beta = c(5, 3), # if alpha>beta then mass is right shifted
#   max_tx_day = 10,
#   min_supp_quant = 0.1,
#   t2_time = 15
# ) {
#     n <- nrow(df_covs)

#     # 1. Evaluate Assignment and Scales
#     lp_tx <- with(df_covs, eval(parse(text = eq_tx_prob)))
#     tx_status <- rbinom(n, 1, pnorm(lp_tx))
#     df_covs["tx"] <- tx_status # just copy into df_cov to make it easier
#     t_switch <- ceiling(rbeta(n, alpha_beta[1], alpha_beta[2])  * max_tx_day)
#     # t_switch <- ceiling(rbeta(n, alpha_beta[1], alpha_beta[2]) * df_covs$x1/max(df_covs$x1) * max_tx_day)
#     t_switch[tx_status == 0] <- Inf
#     t_switch2 <- t_switch + t2_time


#     s0 <- with(df_covs, eval(parse(text = eq_scale_0)))
#     st <- with(df_covs, eval(parse(text = eq_scale_t)))
#     st2 <- with(df_covs, eval(parse(text = eq_scale_t2)))
#     st_cf <- with(df_covs, eval(parse(text = str_replace(eq_scale_t, "tx", "1"))))
#     st_cf2 <- with(df_covs, eval(parse(text = str_replace(eq_scale_t2, "tx", "1"))))

#     # 2. Run Individual Trajectories (The Factual/Observed Data)
#     results <- map(1:n, ~ get_weibull_trajectory(
#         .x, shape, s0[.x], st[.x], st2[.x],
#         t_switch[.x], t_switch2[.x], t_max
#     ))

#     # 3. Calculate Benchmarks (The "Truth" for Comparison)
#     min_supp <- quantile(t_switch[t_switch != Inf], min_supp_quant)

#     truth <- map_dfr(1:t_max, function(curr_t) {
#         # A. E(S | T = t_obs)
#         # The average survival of what actually happened in the simulation
#         s_obs_all <- map_dbl(results, ~ .x$trajectory$survival[curr_t])
#         # 2. calculate e(s | t = t_obs) split by their actual assignment
#         # (tx_status is the vector of 1s and 0s we generated earlier)
#         s_obs_tx <- mean(s_obs_all[tx_status == 1])
#         s_obs_ct <- mean(s_obs_all[tx_status == 0])

#         # B. E(S | T = T^1) and E(S | T = T^0)
#         # The Marginal Counterfactuals (Target Trial)
#         s_cf0 <- exp(-(curr_t / s0)^shape)
#         hz_ct <- (shape / s0) * (curr_t / s0)^(shape - 1)
#         if (curr_t < min_supp) {
#             s_cf1 <- exp(-(curr_t / s0)^shape)
#             hz_tx <- hz_ct
#         } else if (curr_t < min_supp + t2_time) {
#             s_cf1 <- exp(-((min_supp / s0)^shape + (curr_t^shape - min_supp^shape) / st_cf^shape))
#             hz_tx <- (shape / st_cf) * (curr_t / st_cf)^(shape - 1)
#         } else {
#             s_cf1 <- exp(-(
#                 (min_supp / s0)^shape +
#                     ((min_supp + t2_time)^shape - min_supp^shape) / st_cf^shape +
#                     (curr_t^shape - (min_supp + t2_time)^shape) / st_cf2^shape
#             ))
#             hz_tx <- (shape / st_cf2) * (curr_t / st_cf2)^(shape - 1)
#         }


#         tibble(
#             t = curr_t,
#             s_observed = mean(s_obs_all), # E(S | T = t_obs)
#             s_observed_tx = s_obs_tx, # factual survival of the treated
#             s_observed_ct = s_obs_ct, # factual survival of the controls
#             s_cf_tx = mean(s_cf1), # E(S | T = T^1)
#             s_cf_ct = mean(s_cf0), # E(S | T = T^0)
#             # hz_tx = hz_tx,
#             # hz_ct = hz_ct,
#             true_hr_t = mean(hz_tx) / mean(hz_ct),
#             true_hr_t2 = exp(mean(log(hz_tx / hz_ct)))
#         )
#     })


#     # 4. Final Output
#     return(list(
#         data = bind_rows(map(results, "summary")) %>%
#             bind_cols(df_covs) %>%
#             mutate(t_switch_sim = t_switch) %>%
#             mutate(tx_status_sim = tx_status) %>%
#             mutate(t_switch = ifelse(t_switch_sim <= time, t_switch_sim, Inf)) %>%
#             mutate(tx_ever = as.numeric(t_switch <= time)) %>%
#             mutate(h0 = 1/s0, h1 = 1/st, h21=1/st2) # note these are approx hazard (not real)
#             ,             # tx_ever is the real addition
#         truth = truth,
#         min_support = min_supp
#     ))
# }

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


cph_time_match <- function(sim1, eq = "tx_ever + x1", caliper = .05) {
    cph_tx <- sim1$data %>% filter(tx_ever == 1)
    cph_cntrl <- sim1$data %>% filter(tx == 0)
    cph_match_df <- tibble()
    for (i in 1:nrow(cph_tx)) {
        tx_obs <- cph_tx[i, ]
        cntrl_match <- cph_cntrl %>%
            filter(
                time >= tx_obs$t_switch &
                    abs(x1 - tx_obs$x1) < caliper
            )
        if (nrow(cntrl_match) > 0) {
            cph_match_df <- cph_match_df %>%
                bind_rows(
                    cph_tx[i, ] %>%
                        mutate(match_id = i)
                ) %>%
                bind_rows(
                    cntrl_match %>%
                        sample_n(pmin(nrow(cntrl_match), 10), replace = FALSE) %>%
                        mutate(match_id = i)
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
        robust = TRUE
    )
    return(
        list(result = exp(coef(fit)), coh = cph_match_df, fit = fit)
    )
}

# --- 1. cox proportional hazards suite ---
# this function tests the "standard" ways people analyze this data
run_cph_suite <- function(sim_result, eq_cph = "tx_ever + x1", landmark_day = null) {
    dat <- sim_result$data
    # use the simulation's min_support as the default landmark
    lm_day <- ifelse(is.null(landmark_day), sim_result$min_support, landmark_day)

    # a. naive model: ignores timing, uses "ever treated" as a baseline flag
    # (this is usually heavily biased by immortal time)
    #   fit_naive = coxph(Surv(time, event) ~ tx_ever + x1, data = dat)
    # fit_naive <- coxph(Surv(time, event) ~ formula(eq_cph), data = dat)
    outcome <- "Surv(time,event)~"
    fit_naive <- coxph(as.formula(paste0(outcome, eq_cph)), data = dat)
    hr_naive <- exp(coef(fit_naive))["tx_ever"]

    # b. landmark model: only looks at people who survived to 'lm_day'
    dat_lm <- dat %>% filter(time >= lm_day)
    #   fit_lm = coxph(Surv(time, event) ~ tx_ever + x1, data = dat_lm)
    # fit_lm <- coxph(Surv(time, event) ~ formulat(eq_cph), data = dat_lm)
    outcome <- "Surv(time,event)~"
    fit_lm <- coxph(as.formula(paste0(outcome, eq_cph)), data = dat_lm)
    hr_lm <- exp(coef(fit_lm))["tx_ever"]

    # c. time-dependent covariate (tdc): the "gold standard" for cox
    # re-structure data into start/stop intervals
    dat_long <- tmerge(dat, dat, id = id, tstop = time, evt = event(time, event))
    dat_long <- tmerge(dat_long, dat, id = id, tx_td = tdc(t_switch - .01)) # the eps here makes sure that if tt_tx =5 and tt_event=5, that the patient is counted as tx when the event occurs

    #   fit_tdc = coxph(Surv(tstart, tstop, evt) ~ tx_td + x1, data = dat_long)
    outcome <- "Surv(tstart, tstop, evt) ~ "
    eq_cph <- str_replace_all(eq_cph, "tx_ever", "tx_td")
    fit_tdc <- coxph(as.formula(paste0(outcome, eq_cph)), data = dat_long)
    hr_tdc <- exp(coef(fit_tdc))["tx_td"]

    # return a summary table of the estimated hazard ratios
    return(
        list(
            result = tibble(
                method = c("cph_naive", "cph_landmark", "cph_tdc"),
                hr_estimate = c(hr_naive, hr_lm, hr_tdc)
            ),
            lnd_coh = dat_lm,
            fit_lm = fit_lm,
            fit_tdc = fit_tdc
        )
    )
}

# --- 2. bart survival function ---
# uses machine learning to estimate s(t | x, tx_td)
run_bart_survival <- function(sim_result, covars, n_trees = 50, n_cores = 1, trunc = 1, k_p = 1, keepevery = 10, nskip = 500) {
    dat <- sim_result$data %>%
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
    min_support <- sim_result$min_support
    z0_set <- times[times < min_support]
    z1_set <- times[times >= min_support]

    df_pred1 <- long %>%
        select(t, tx_td, all_of(covars)) %>%
        mutate(tx_td = rep(c(rep(0, length(z0_set)), rep(1, length(z1_set))), n_obs))

    df_pred0 <- long %>%
        select(t, tx_td, all_of(covars)) %>%
        mutate(tx_td = 0)

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
                             dt = 1, 
                             min_support
                             ) {
    # 1. Clone creation: each person becomes 2 rows, one per strategy
    clones <- data %>%
        dplyr::select(id, tx, t_switch, event, time, all_of(covars)) %>%
        tidyr::uncount(weights = 2, .id = "clone_arm") %>%
        mutate(
            arm = if_else(clone_arm == 1L, "control_no_tx", "treat_tx")
        )

    # 2. Define protocol deviation time (dev_time)
    # Control arm deviates if they get tx within grace
    # Treat arm deviates if they do NOT get tx within grace
    clones <- clones %>%
        mutate(
            # did they have tx within the grace period?
            tx_in_grace = (tx == 1L & !is.na(t_switch) & t_switch <= grace_period),
            dev_time = case_when(
                arm == "control_no_tx" & tx_in_grace ~ t_switch,
                arm == "treat_tx" & !tx_in_grace & time >= grace_period ~ grace_period,
                TRUE ~ NA_real_
            ),
            dev = if_else(!is.na(dev_time), 1L, 0L)
        )

    # 3. Protocol-censored outcomes for trial emulation
    # time_cc = follow-up time after protocol censoring + admin censoring at tau
    clones <- clones %>%
        mutate(
            time_cc = if_else(!is.na(dev_time), pmin(time, dev_time), time),
            time_cc = pmin(time_cc, tau),
            event_cc = if_else(event == 1L & time <= time_cc, 1L, 0L)
        )

    # 4. IPCW weight creation (pooled logistic, during grace period)
    panel <- clones %>%
        mutate(
            tmax  = pmin(grace_period, time_cc),
            n_int = pmax(ceiling(tmax / dt), 0L)
        ) %>%
        tidyr::uncount(weights = n_int, .id = "k") %>%
        mutate(
            t0 = (k - 1L) * dt, # interval start
            t1 = k * dt, # interval end
            t = t1, # simple time term
            t_sq = t^2
        ) %>%
        mutate(
            cens = if_else(
                dev == 1L & !is.na(dev_time) & dev_time > t0 & dev_time <= t1,
                1L, 0L
            )
        ) %>%
        arrange(id, clone_arm, t1) %>%
        group_by(id, clone_arm) %>%
        filter(cumsum(cens) <= 1) %>%
        filter(t >= min_support) %>% 
        ungroup()
    return(panel)

    # 5. Fit pooled logistic models separately by arm
    denom_formula <- as.formula(paste("cens ~", paste(covars, collapse = " + "), "+ t + t_sq"))
    num_formula <- cens ~ t + t_sq

    # Helper function to calculate weights
    get_arm_weights <- function(arm_data) {
        if (nrow(arm_data) == 0) {
            return(tibble())
        }
        m_denom <- glm(denom_formula, data = arm_data, family = binomial())
        m_num <- glm(num_formula, data = arm_data, family = binomial())

        arm_data %>%
            mutate(
                p_unc_denom = pmax(1 - predict(m_denom, type = "response"), 1e-6),
                p_unc_num   = pmax(1 - predict(m_num, type = "response"), 1e-6),
                w_inc       = p_unc_num / p_unc_denom
            ) %>%
            group_by(id, clone_arm) %>%
            summarise(w_ipcw = prod(w_inc), .groups = "drop")
    }

    weights_all <- bind_rows(
        get_arm_weights(panel %>% filter(arm == "control_no_tx")),
        get_arm_weights(panel %>% filter(arm == "treat_tx"))
    )

    # 6. Apply weights and trim
    clones_w <- clones %>%
        left_join(weights_all, by = c("id", "clone_arm")) %>%
        mutate(w_ipcw = tidyr::replace_na(w_ipcw, 1))

    # q_lim <- quantile(clones_w$w_ipcw, probs = c(0.05, 0.95), na.rm = TRUE)
    q_lim <- quantile(clones_w$w_ipcw, probs = c(0.01, 0.99), na.rm = TRUE)
    clones_w <- clones_w %>%
        mutate(w_trim = pmin(pmax(w_ipcw, q_lim[1]), q_lim[2]))

    # 7. Fit weighted KM (FIXED: Added id = id)
    fit_ccw <- survival::survfit(
        Surv(time_cc, event_cc) ~ arm,
        data = clones_w,
        weights = w_trim,
        id = id 
    )


    # Compute Cox weights separately by arm using simulation names
    clones_cox <- bind_rows(
        cox_ipcw_by_arm(clones_w %>% filter(arm == "control_no_tx"), covars, grace_period = grace_period),
        cox_ipcw_by_arm(clones_w %>% filter(arm == "treat_tx"), covars, grace_period = grace_period)
    )
    q_cox <- quantile(clones_cox$w_cox, probs = c(0.01, 0.99), na.rm = TRUE)
    clones_cox <- clones_cox %>%
        mutate(
            w_cox_trim = pmin(pmax(w_cox, q_cox[[1]]), q_cox[[2]])
        )
        
    # Fit Cox weighted KM (FIXED: Added id = id)
    fit_ccw_cox <- survival::survfit(
        Surv(time_cc, event_cc) ~ arm,
        data = clones_cox,
        weights = w_cox_trim,
        id = id
    )

    # TMLE
    # tmle <- run_tmle_survival(clones_w, covars = covars, tau = tau)
    tmle <- NULL

    # 8. Extract Results (FIXED: used 'times' and 'extend = TRUE')
    sum_log <- summary(fit_ccw, times = 1:tau, extend = TRUE)
    res_log <- tibble(
        t = sum_log$time,
        tx = sum_log$strata,
        sv = sum_log$surv
    ) %>%
        pivot_wider(names_from = tx, values_from = sv) %>%
        mutate(rd = `arm=treat_tx` - `arm=control_no_tx`)

    sum_cox <- summary(fit_ccw_cox, times = 1:tau, extend = TRUE)
    res_cox <- tibble(
        t = sum_cox$time,
        tx = sum_cox$strata,
        sv = sum_cox$surv
    ) %>%
        pivot_wider(names_from = tx, values_from = sv) %>%
        mutate(rd = `arm=treat_tx` - `arm=control_no_tx`)

    return(
        list(
            fit = fit_ccw,
            data = clones_w,
            res_log = res_log,
            res_cox = res_cox,
            res_tmle = tmle
        )
    )
}



# run_ccw_analysis <- function(data,
#                              grace_period = 10,
#                              tau = 30,
#                              covars = c("x1"),
#                              dt = 1) {
#     # 1. Clone creation: each person becomes 2 rows, one per strategy
#     clones <- data %>%
#         dplyr::select(id, tx, t_switch, event, time, all_of(covars)) %>%
#         tidyr::uncount(weights = 2, .id = "clone_arm") %>%
#         mutate(
#             arm = if_else(clone_arm == 1L, "control_no_tx", "treat_tx")
#         )

#     # 2. Define protocol deviation time (dev_time)
#     # Control arm deviates if they get tx within grace
#     # Treat arm deviates if they do NOT get tx within grace
#     clones <- clones %>%
#         mutate(
#             # did they have tx within the grace period?
#             tx_in_grace = (tx == 1L & !is.na(t_switch) & t_switch <= grace_period),
#             dev_time = case_when(
#                 arm == "control_no_tx" & tx_in_grace ~ t_switch,
#                 arm == "treat_tx" & !tx_in_grace & time >= grace_period ~ grace_period,
#                 TRUE ~ NA_real_
#             ),
#             dev = if_else(!is.na(dev_time), 1L, 0L)
#         )

#     # 3. Protocol-censored outcomes for trial emulation
#     # time_cc = follow-up time after protocol censoring + admin censoring at tau
#     clones <- clones %>%
#         mutate(
#             time_cc = if_else(!is.na(dev_time), pmin(time, dev_time), time),
#             time_cc = pmin(time_cc, tau),
#             event_cc = if_else(event == 1L & time <= time_cc, 1L, 0L)
#         )

#     # 4. IPCW weight creation (pooled logistic, during grace period)
#     panel <- clones %>%
#         mutate(
#             tmax  = pmin(grace_period, time_cc),
#             n_int = pmax(ceiling(tmax / dt), 0L)
#         ) %>%
#         tidyr::uncount(weights = n_int, .id = "k") %>%
#         mutate(
#             t0 = (k - 1L) * dt, # interval start
#             t1 = k * dt, # interval end
#             t = t1, # simple time term
#             t_sq = t^2
#         ) %>%
#         mutate(
#             cens = if_else(
#                 dev == 1L & !is.na(dev_time) & dev_time > t0 & dev_time <= t1,
#                 1L, 0L
#             )
#         ) %>%
#         arrange(id, clone_arm, t1) %>%
#         group_by(id, clone_arm) %>%
#         filter(cumsum(cens) <= 1) %>%
#         ungroup()

#     # 5. Fit pooled logistic models separately by arm
#     denom_formula <- as.formula(paste("cens ~", paste(covars, collapse = " + "), "+ t + t_sq"))
#     num_formula <- cens ~ t + t_sq

#     # Helper function to calculate weights
#     get_arm_weights <- function(arm_data) {
#         if (nrow(arm_data) == 0) {
#             return(tibble())
#         }
#         m_denom <- glm(denom_formula, data = arm_data, family = binomial())
#         m_num <- glm(num_formula, data = arm_data, family = binomial())

#         arm_data %>%
#             mutate(
#                 p_unc_denom = pmax(1 - predict(m_denom, type = "response"), 1e-6),
#                 p_unc_num   = pmax(1 - predict(m_num, type = "response"), 1e-6),
#                 w_inc       = p_unc_num / p_unc_denom
#             ) %>%
#             group_by(id, clone_arm) %>%
#             summarise(w_ipcw = prod(w_inc), .groups = "drop")
#     }

#     weights_all <- bind_rows(
#         get_arm_weights(panel %>% filter(arm == "control_no_tx")),
#         get_arm_weights(panel %>% filter(arm == "treat_tx"))
#     )

#     # 6. Apply weights and trim
#     clones_w <- clones %>%
#         left_join(weights_all, by = c("id", "clone_arm")) %>%
#         mutate(w_ipcw = tidyr::replace_na(w_ipcw, 1))

#     q_lim <- quantile(clones_w$w_ipcw, probs = c(0.01, 0.99), na.rm = TRUE)
#     clones_w <- clones_w %>%
#         mutate(w_trim = pmin(pmax(w_ipcw, q_lim[1]), q_lim[2]))

#     # 7. Fit weighted KM
#     fit_ccw <- survival::survfit(
#         Surv(time_cc, event_cc) ~ arm,
#         data = clones_w,
#         weights = w_trim
#     )


#     # Compute Cox weights separately by arm using simulation names
#     clones_cox <- bind_rows(
#         cox_ipcw_by_arm(clones_w %>% filter(arm == "control_no_tx"), covars, grace_period = 10),
#         cox_ipcw_by_arm(clones_w %>% filter(arm == "treat_tx"), covars, grace_period = 10)
#     )
#     q_cox <- quantile(clones_cox$w_cox, probs = c(0.01, 0.99), na.rm = TRUE)
#     clones_cox <- clones_cox %>%
#         mutate(
#             w_cox_trim = pmin(pmax(w_cox, q_cox[[1]]), q_cox[[2]])
#         )
#     fit_ccw_cox <- survfit(
#         Surv(time_cc, event_cc) ~ arm,
#         data = clones_cox,
#         weights = w_cox_trim
#     )

#     # TMLE
#     # tmle <- run_tmle_survival(clones_w, covars = covars, tau = tau)
#     tmle <- NULL

#     res_log <- tibble(
#         t = summary(fit_ccw, time = 1:tau)$time,
#         tx = summary(fit_ccw, time = 1:tau)$strata,
#         sv = summary(fit_ccw, time = 1:tau)$surv
#     ) %>%
#         pivot_wider(names_from = tx, values_from = sv) %>%
#         mutate(rd = `arm=treat_tx` - `arm=control_no_tx`)
#     res_cox <- tibble(
#         t = summary(fit_ccw_cox, time = 1:tau)$time,
#         tx = summary(fit_ccw_cox, time = 1:tau)$strata,
#         sv = summary(fit_ccw_cox, time = 1:tau)$surv
#     ) %>%
#         pivot_wider(names_from = tx, values_from = sv) %>%
#         mutate(rd = `arm=treat_tx` - `arm=control_no_tx`)

#     return(
#         list(
#             fit = fit_ccw,
#             data = clones_w,
#             res_log = res_log,
#             res_cox = res_cox,
#             res_tmle = tmle
#         )
#     )
# }

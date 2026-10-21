library(tidyverse)
library(survival)
library(BART)
library(randomForest)
library(glmnet)
library(glue)

# --- 1. Core Piecewise Weibull Logic ---
get_weibull_trajectory <- function(id, shape, s0, st, st2, t_switch, t_switch2, t_max) {
  u <- runif(1)
  target_h <- -log(pmax(u, 1e-10))
  
  # Calculate cumulative hazards at the switch points
  h_at_switch <- (t_switch / s0)^shape
  h_at_switch2 <- h_at_switch + (max(0, t_switch2^shape - t_switch^shape)) / st^shape

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
  sv <- map_dbl(days, function(d) {
    if (d < t_switch) {
      exp(-(d / s0)^shape)
    } else if (d < t_switch2) {
      exp(-((t_switch / s0)^shape + (max(0, d^shape - t_switch^shape) / st^shape)))
    } else {
      exp(-((t_switch / s0)^shape + 
            (max(0, t_switch2^shape - t_switch^shape) / st^shape) + 
            (max(0, d^shape - t_switch2^shape) / st2^shape)))
    }
  })

  list(
    summary = tibble(id = id, time = min(ceiling(event_time), t_max), 
                     event = as.numeric(event_time <= t_max)),
    trajectory = tibble(id = id, day = days, survival = sv)
  )
}

generate_survival_sim <- function(df_covs, t_max = 30, shape = 1, 
                                  eq_tx_prob = "qnorm(0.4) + 0.5*x1",
                                  eq_scale_0 = "t_max * exp(-(x1 * log(2)) / shape)",
                                  eq_scale_t = "t_max * exp(-(tx * log(0.5) + x1 * log(2)) / shape)",
                                  eq_scale_t2 = NULL,
                                  max_tx_day = 10, t2_offset = 15) {
  n <- nrow(df_covs)
  if(is.null(eq_scale_t2)) eq_scale_t2 <- eq_scale_t
  
  lp_tx <- with(df_covs, eval(parse(text = eq_tx_prob)))
  tx_status <- rbinom(n, 1, pnorm(lp_tx))
  df_covs$tx_assigned <- tx_status
  
  t_switch <- ceiling(runif(n, 1, max_tx_day))
  t_switch[tx_status == 0] <- Inf
  t_switch2 <- t_switch + t2_offset

  s0  <- with(df_covs, eval(parse(text = eq_scale_0)))
  st  <- with(df_covs %>% mutate(tx=1), eval(parse(text = eq_scale_t)))
  st2 <- with(df_covs %>% mutate(tx=1), eval(parse(text = eq_scale_t2)))
  
  results <- map(1:n, ~ get_weibull_trajectory(.x, shape, s0[.x], st[.x], st2[.x], 
                                                t_switch[.x], t_switch2[.x], t_max))
  
  truth <- map_dfr(c(15, 20, 30), function(curr_t) {
    s_cf0 <- mean(exp(-(curr_t / s0)^shape))
    # Truth assumes immediate treatment for target trial comparison
    s_cf1 <- mean(exp(-(1 / s0)^shape + (max(0, curr_t^shape - 1^shape) / st^shape)))
    tibble(t = curr_t, true_s0 = s_cf0, true_s1 = s_cf1, true_sd = s_cf1 - s_cf0)
  })

  list(
    data = bind_rows(map(results, "summary")) %>% bind_cols(df_covs) %>%
           mutate(t_switch = ifelse(tx_assigned == 1 & t_switch <= time, t_switch, Inf),
                  tx_ever = as.numeric(!is.infinite(t_switch))),
    truth = truth
  )
}

# --- 2. Analytics Suite ---

run_cph_suite <- function(sim_res, landmarks = c(15, 20, 30)) {
  dat <- sim_res$data
  dat_long <- tmerge(dat, dat, id = id, tstop = time, evt = event(time, event))
  dat_long <- tmerge(dat_long, dat, id = id, tx_td = tdc(t_switch - 0.001))
  
  fit <- coxph(Surv(tstart, tstop, evt) ~ tx_td + x1, data = dat_long)
  
  est <- map_dfr(landmarks, function(l) {
    new_dat <- expand.grid(tx_td = c(0, 1), x1 = dat$x1)
    s_curve <- survfit(fit, newdata = new_dat)
    # The survfit object returns a matrix where rows are time and cols are individuals
    t_idx <- which.min(abs(s_curve$time - l))
    s_vals <- s_curve$surv[t_idx, ]
    s0_avg <- mean(s_vals[1:nrow(dat)])
    s1_avg <- mean(s_vals[(nrow(dat)+1):(nrow(dat)*2)])
    tibble(t = l, method = "cph_tdc", sd_estimate = s1_avg - s0_avg)
  })
  
  list(estimates = est, fit = fit)
}

run_bart_suite <- function(sim_res, landmarks = c(15, 20, 30), n_trees = 50) {
  dat <- sim_res$data
  n <- nrow(dat)
  times <- sort(unique(dat$time))
  
  long_train <- tibble(t = rep(times, n), id = rep(dat$id, each = length(times))) %>%
    left_join(dat, by = "id") %>%
    mutate(tx_td = ifelse(t_switch - 0.001 <= t, 1, 0)) %>%
    filter(t <= time) %>%
    mutate(y = ifelse(t == time & event == 1, 1, 0))

  post <- BART::mc.surv.bart(x.train = long_train %>% select(t, tx_td, x1), 
                             y.train = long_train$y, ntree = n_trees)
  
  est <- map_dfr(landmarks, function(l) {
    # Generate predictions for every person under both conditions at the landmark
    pred0 <- tibble(t = l, tx_td = 0, x1 = dat$x1)
    pred1 <- tibble(t = l, tx_td = 1, x1 = dat$x1)
    
    p0 <- predict(post, pred0)
    p1 <- predict(post, pred1)
    
    s0_avg <- mean(p0$surv.test)
    s1_avg <- mean(p1$surv.test)
    tibble(t = l, method = "bart", sd_estimate = s1_avg - s0_avg)
  })
  
  list(estimates = est, model = post)
}

run_ccw_tmle_suite <- function(sim_res, landmarks = c(15, 20, 30)) {
  dat <- sim_res$data
  clones <- dat %>% uncount(2, .id = "clone_arm") %>%
    mutate(arm = if_else(clone_arm == 1, "control", "treated"),
           dev = case_when(arm == "control" & tx_ever == 1 ~ 1,
                           arm == "treated" & tx_ever == 0 ~ 1, TRUE ~ 0))
  
  m_w <- glm(dev ~ x1, data = clones, family = binomial)
  clones$w <- 1 / (1 - predict(m_w, type = "response"))
  
  est <- map_dfr(landmarks, function(l) {
    c_dat <- clones %>% mutate(time_cc = pmin(time, l), 
                               event_cc = if_else(event == 1 & time <= l, 1, 0))
    fit <- survfit(Surv(time_cc, event_cc) ~ arm, data = c_dat, weights = w)
    s <- summary(fit, times = l)
    ccw_sd <- s$surv[2] - s$surv[1]
    
    # TMLE update
    tmle_dat <- c_dat %>% filter(time_cc >= l | event_cc == 1) %>% 
      mutate(Y = event_cc, A = if_else(arm == "treated", 1, 0))
    m_q <- glm(Y ~ A + x1, data = tmle_dat, family = binomial, weights = w)
    q1 <- predict(m_q, mutate(tmle_dat, A=1), type="response")
    q0 <- predict(m_q, mutate(tmle_dat, A=0), type="response")
    
    g <- 0.5 
    h <- ifelse(tmle_dat$A == 1, 1/g, -1/(1-g))
    m_eps <- glm(Y ~ -1 + h, offset = qlogis(pmin(pmax(ifelse(tmle_dat$A==1, q1, q0), .001), .999)), 
                 family=binomial, data=tmle_dat, weights=w)
    eps <- coef(m_eps)
    
    q1_star <- plogis(qlogis(q1) + eps*(1/g))
    q0_star <- plogis(qlogis(q0) + eps*(-1/(1-g)))
    tmle_sd <- (mean(q1_star) - mean(q0_star)) * -1 
    
    tibble(t = l, ccw_sd = ccw_sd, tmle_sd = tmle_sd)
  }) %>% 
    pivot_longer(cols = contains("sd"), names_to = "method", values_to = "sd_estimate")
  
  list(estimates = est)
}
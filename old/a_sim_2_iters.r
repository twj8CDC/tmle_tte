# set-up
setwd("survival_sim")
source("sim/a_sim_func2.r")
set.seed(99)

# Sim1
# This simulation includes a binary tx and a
# single continuous covariate. The outcome-tx relationship 
# # is confounded# by the single covariate.
# # In the simulation the tx coefficient varies from 0 to 5.
# # Survival Difference and Hazard Ratios are the primary estimands
# f_hr <- data.frame()
# f_rd = data.frame()
# for (i in 1:2) {
#     i = i*10
#     set.seed(i)
    
#     n <- 300
#     tx_coef = .5
#     alpha_coef = 7
#     beta_coef = 3
#     t2_time_coef = 5

#     df_covs <- data.frame(x1 = rnorm(n, 0, 1))
        
#     sim1 <- generate_survival_sim(
#         df_covs,
#         t_max = 30,
#         shape = 1.5,
#         eq_tx_prob = "qnorm(0.5) + x1",
#         # eq_tx_prob = "0",
#         eq_scale_0 = "t_max * exp(-(x1 * log(1.2)) / shape)", # no treatment effect before they are treated
#         eq_scale_t = glue("t_max * exp(-(tx * log({tx_coef}) + x1 * log(1.2)) / shape)"),
#         eq_scale_t2 = glue("t_max * exp(-(tx * log({tx_coef}) + x1 * log(1.2)) / shape)"),
#         alpha_beta = c(alpha_coef, beta_coef),
#         max_tx_day = 10,
#         min_supp_quant = 0.1,
#         t2_time = t2_time_coef
#     )


#     # run models
#     cph_df <- run_cph_suite(sim1, landmark_day = 11)
#     cph_tm = cph_time_match(sim1, caliper = 0.05)
    
#     b <- run_bart_survival(
#         sim1,
#         covars = c("x1"),
#         n_tree = 50,
#         n_cores = 8,
#         trunc = 1,
#         k_p = 2 # turn down k to keep the unknowns drifting back to 1
#     )
#     ccw <- run_ccw_analysis(sim1$data) #standard
    
#     # collect results
#     hr_df <- bind_rows(
#         cph_df$result,
#         data.frame(
#             method = "cph_match",
#             hr_estimate = cph_tm$result["tx_ever"]
#         ),
#         data.frame(
#             method = c("bart_mean", "bart_median"),
#             hr_estimate = c(b$hr$`mean(value)`, b$hr$`quantile(value, 0.5)`)
#         ) 
#     ) %>% 
#         mutate(
#             i = i,
#             tx_coef, alpha_coef, beta_coef, t2_time_coef,
#             min_support = sim1$min_support
#             )
    
#     f_hr <- bind_rows(f_hr, hr_df) 
#     b_max_time = length(b$rd$`mean(value)`)
#     b_max_time
#     t_max
#     rd_df = tibble(
#         bart_sd = b$rd$`mean(value)`[b_max_time], 
#         ccw_log_sd = ccw$res_log$rd[30],
#         ccw_cox_sd = ccw$res_cox$rd[30],
#         ccw_tmle_sd = ccw$res_tmle %>% filter(Statistic == "sd") %>% pull(Value)
#         ) %>% 
#         pivot_longer(cols=everything(), names_to="method", values_to="sd_estimate") %>% 
#         mutate(
#             i = i,
#             tx_coef, alpha_coef, beta_coef, t2_time_coef,
#             min_support = sim1$min_support,
#             sd_truth = sim1$truth %>% filter(t==30) %>% mutate(sd = s_cf_tx - s_cf_ct) %>% pull(sd)
#             ) 
#     rd_df
#     f_rd = bind_rows(f_rd, rd_df)
# }
# options(width=400)
# b$hr
# g = sim1$truth %>%
#     pivot_longer(cols = c(s_observed_tx, s_observed_ct, s_cf_tx, s_cf_ct)) %>%
#     mutate(tx = ifelse(str_detect(name, "tx"), "tx", "ct")) %>%
#     mutate(type = ifelse(str_detect(name, "cf"), "counterfact", "observed")) %>%
#     ggplot(aes(x = t, y = value, col = type, shape = tx)) +
#     geom_point() +
#     ylim(c(0, 1)) +
#     scale_color_discrete() +
#     theme_minimal()

# saveRDS(list(
#     rd = f_rd, 
#     hr = f_hr,
#     sim1 = sim1,
#     plt = g),
#     "res/sim1.rds"
#     )

# b$hr
# mean(
#     rowMeans(
#         b$bres %>% select(contains("hz_tx_"))/
#         b$bres %>% select(contains("hz_ct_"))
#         )[sim1$min_support:nrow(b$bres)]
# )

# sim1r = readRDS("res/sim1.rds")

# sim1r$hr %>%
#     filter(method != "bart_mean") %>% # drop because known not correct
#     group_by(method) %>%
#     mutate(hr_diff = log(hr_estimate) - log(tx_coef)) %>% 
#     summarise(
#         bias = mean(hr_diff),
#         rmse = sqrt(mean((hr_diff)^2))
#     ) %>% 
#     arrange(rmse)
# sim1r$rd %>%
#     group_by(method) %>%
#     mutate(rd_diff = sd_estimate - sd_truth) %>% 
#     summarise(
#         bias = mean(rd_diff),
#         rmse = sqrt(mean((rd_diff)^2))
#     ) %>% 
#     arrange(rmse)



########################################
# simulation 2
# break by simple non-linearity

f_hr <- data.frame()
f_sd = data.frame()
for (i in 1:10) {
    i = i*2033
    set.seed(i)
    
    n <- 300
    # tx_coef = runif(1,0,5)
    tx_coef = 1.2
    interact_coef = 2
    alpha_coef = 10
    beta_coef = 1
    t2_time_coef = 5
    t_max = 30
    shape = 1.5

    df_covs <- data.frame(x1 = rnorm(n, 3, 2))
    sim1 <- generate_survival_sim(
        df_covs,
        t_max = t_max,
        shape = shape,
        eq_tx_prob = "qnorm(0.5) + x1 - mean(x1)",
        # eq_tx_prob = "0",
        # eq_scale_0 = "t_max*3 * exp(-(
        #         x1  * log(1.3)
        #     ) / shape)", # no treatment effect before they are treated
        # eq_scale_t = glue("t_max*3 * exp(-(
        #     tx*log({tx_coef}) + x1*log(1.3) + x1*tx*log({interact_coef})
        #     ) / shape)"),
        # eq_scale_t2 = glue("t_max*3 * exp(-(
        #     tx*log({tx_coef}) + x1*log(1.3) + x1*tx*log({interact_coef})
        #     ) / shape)"),

        eq_scale_0 = "t_max*3 * exp(-(
                x1  * log(1.2)
            ) / shape)", # no treatment effect before they are treated
        eq_scale_t = glue("t_max*3 * exp(-(
            tx*log({tx_coef}) + x1*log(1.2)
            ) / shape)"),
        eq_scale_t2 = glue("t_max*3 * exp(-(
            tx*log({tx_coef}) + x1*log(1.2)
            ) / shape)"),
        # alpha_beta = c(alpha_coef, beta_coef),
        alpha_beta = c(alpha_coef,beta_coef),
        max_tx_day = 10,
        min_supp_quant = 0.1,
        t2_time = t2_time_coef
    )


        
    hr_df = tibble()
    sd_df = tibble()
    cph_df <- run_cph_suite(sim1, eq_cph = "tx_ever + x1 + tx_ever*x1", landmark_day = 11)
    # landmark
    fit = cph_df$fit_lm
    beta_tx <- coef(fit)["tx_ever"]
    beta_int <- coef(fit)["tx_ever:x1"]
    individual_log_hrs <- beta_tx + (beta_int * cph_df$lnd_coh$x1)
    individual_log_hrs <- beta_tx + (beta_int * sim1$data$x1)
    approx_marginal_hr <- exp(mean(individual_log_hrs))
    hr_df = bind_rows(hr_df,
        tibble(method = "cph_lm", hr = approx_marginal_hr)
    )
    cph_sv_lm = cph_pred_survival(cph_df$fit_lm, cph_df$lnd_coh)
    cph_lm_sd= cph_sv_lm %>% tail(1) %>% pull(sd_diff)
    sd_df = bind_rows(sd_df,
        tibble(method = "cph_lm", sd=cph_lm_sd)
    )

    # tdc
    fit = cph_df$fit_tdc
    beta_tx <- coef(fit)["tx_td"]
    beta_int <- coef(fit)["tx_td:x1"]
    individual_log_hrs <- beta_tx + (beta_int * sim1$data$x1)
    approx_marginal_hr <- exp(mean(individual_log_hrs))
    hr_df = bind_rows(hr_df,
        tibble(method = "cph_tdc", hr = approx_marginal_hr)
    )
    cph_sv = cph_pred_survival_tdc(sim1$data, fit, sim1$min_support)
    cph_tdc_sd = cph_sv %>% tail(1) %>% pull(sd_diff)
    sd_df = bind_rows(sd_df,
        tibble(method = "cph_tdc", sd=cph_tdc_sd)
    )


    cph_tm = cph_time_match(sim1, eq="tx_ever+x1+tx_ever*x1", caliper = 0.05)
    fit = cph_tm$fit
    beta_tx <- coef(fit)["tx_ever"]
    beta_int <- coef(fit)["tx_ever:x1"]
    individual_log_hrs <- beta_tx + (beta_int * cph_tm$coh$x1)    
    approx_marginal_hr <- exp(mean(individual_log_hrs))
    hr_df = bind_rows(hr_df,
        tibble(method = "cph_tm", hr = approx_marginal_hr)
    )



    # b <- run_bart_survival(
    #     sim1,
    #     covars = c("x1"),
    #     n_tree = 50,
    #     n_cores = 8,
    #     trunc = 1,
    #     k_p = 1.5, # turn down k to keep the unknowns drifting back to 1
    #     nskip=100,
    #     keepevery = 5
    # )

    # b_hr = mean(b$hr$hr_estimate)
    # b_hr2 =mean(b$hr2$hr_mean)
    # b_hr3 =mean(b$hr2$hr_median)
    # hr_df = bind_rows(hr_df,
    #     tibble(
    #         method = c("bart_hom", "bart_moh_mean", "bart_moh_median"), 
    #         hr = c(b_hr, b_hr2, b_hr3),
    #         )
    # )
    # # hr_df    
    # # sim1$truth
    # # Survival Difference
    # b_sd = b$rd %>% tail(1) %>% pull(sd_estimate)
    # sd_df = bind_rows(sd_df,
    #     tibble(method = "bart", sd = b_sd)
    # )
    
    
    ccw <- run_ccw_analysis(sim1$data, tau=t_max, min_support = sim1$min_support) #standard
    
    ccw_log_sd = ccw$res_log %>% tail(1) %>% pull(rd)
    ccw_cox_sd = ccw$res_cox %>% tail(1) %>% pull(rd)
    # ccw_log_sd
    # ccw_cox_sd

    sd_df = bind_rows(sd_df,
        tibble(method = "ccw_log", sd = ccw_log_sd),
        tibble(method = "ccw_cox", sd = ccw_cox_sd )
    )
    plot(sim1$truth$t, sim1$truth$s_cf_ct, col="red", ylim=c(0,1))    
    points(sim1$truth$t, sim1$truth$s_cf_tx, col="red")

    points(1:30, ccw$res_log$`arm=control_no_tx`, col="green")
    points(1:30, ccw$res_log$`arm=treat_tx`, col="green")
    
    points(1:30, ccw$res_cox$`arm=control_no_tx`, col="purple")
    points(1:30, ccw$res_cox$`arm=treat_tx`, col="purple")
    
    
    points(b$rd$t, b$rd$m1, col="black")
    points(b$rd$t, b$rd$m0, col="black")


    sd_df = sd_df %>% 
        mutate(sd_truth = sim1$truth %>% tail(1) %>% mutate(sd = s_cf_tx - s_cf_ct) %>% pull(sd)) %>% 
        mutate(sd_diff = sd - sd_truth) %>% 
        mutate(
            i = i,
            tx_coef, alpha_coef, beta_coef, t2_time_coef,
            min_support = sim1$min_support
            ) 

    hr_df = hr_df %>% 
        mutate(hr_truth = sim1$truth %>% filter(t >= sim1$min_support) %>% summarise(mean(true_hr_t2)) %>% pull()) %>% 
        mutate(hr_diff = hr - hr_truth) %>% 
        mutate(
            i = i,
            tx_coef, alpha_coef, beta_coef, t2_time_coef,
            min_support = sim1$min_support
            ) 
    sd_df
    f_hr = bind_rows(f_hr, hr_df)
    f_sd = bind_rows(f_sd, sd_df)
}

ccw$data %>% select(arm, x1, w_trim) %>% 
    ggplot(aes(x=x1, y=w_trim, col=arm)) +
    geom_point()



ccw$data %>%
    group_by(arm) %>%
    summarize(
        raw_mean_x1 = mean(x1),
        weighted_mean_x1 = weighted.mean(x1, w = w_trim)
    )



ccw$data %>% filter(x1 > 6 | x1 < -2)
sim1$data %>% 
    group_by(tx_ever) %>% 
    summarise(mean(x1))


f_sd 
f_sd %>% 
    group_by(method) %>% 
    summarise(bias = mean(sd_diff), rmse = sqrt(mean(sd_diff ^2)))

f_hr %>% 
    group_by(method) %>% 
    summarise(bias = mean(hr_diff), rmse = sqrt(mean(hr_diff ^2)))

sim1$data %>% group_by(tx_ever) %>% summarise(mean(time), mean(x1), n())
cph_df$lnd_coh %>% group_by(tx_ever) %>% summarise(mean(time), mean(x1), n())
cph_tm$coh %>% group_by(tx_ever) %>% summarise(mean(time), mean(x1), n())




plot(sim1$truth$t, sim1$truth$s_cf_tx, ylim = c(0,1))
points(sim1$truth$t, sim1$truth$s_cf_ct)

sim1$truth
# saveRDS(list(
#     rd = f_rd, 
#     hr = f_hr,
#     sim1 = sim1,
#     plt = g),
#     "res/sim2.rds"
#     )
f_hr
sim1$truth
f_rd
f_hr



b$hr2
f_hr %>%
    filter(method != "bart_mean") %>% # drop because known not correct
    group_by(method) %>%
    mutate(hr_diff = hr_estimate -hr_truth) %>% 
    summarise(
        bias = mean(hr_diff),
        rmse = sqrt(mean((hr_diff)^2))
    ) %>% 
    arrange(rmse)
f_hrd
tibble(f_hr) %>% print(n=100)
f_rd %>%
    group_by(method) %>%
    mutate(rd_diff = sd_estimate - sd_truth) %>% 
    summarise(
        bias = mean(rd_diff, na.rm=TRUE),
        rmse = sqrt(mean((rd_diff)^2, na.rm=TRUE))
    ) %>% 
    arrange(rmse)
f_rd
f_hr
f_hr

sim1$truth %>%
    pivot_longer(cols = c(s_observed_tx, s_observed_ct, s_cf_tx, s_cf_ct)) %>%
    mutate(tx = ifelse(str_detect(name, "tx"), "tx", "ct")) %>%
    mutate(type = ifelse(str_detect(name, "cf"), "counterfact", "observed")) %>%
    ggplot(aes(x = t, y = value, col = type, shape = tx)) +
    geom_point() +
    ylim(c(0, 1)) +
    scale_color_discrete() +
    theme_minimal()

dev.off()

########################################
##################sim 3
###############
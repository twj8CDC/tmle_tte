# Need to test out the sim-iters
options(width=400)


# need to establish the true Dynamic Treatment Rule HR
source("sim/a_sim_func2.r")


i = i*2033
set.seed(i)

{
    n <- 300
    # tx_coef = runif(1,0,5)
    tx_coef = 1.8
    interact_coef = 2
    alpha_coef = 8
    beta_coef = 2
    t2_time_coef = 5
    t_max = 30
    shape = 2
    base_scale_mult = 1
    df_covs <- data.frame(x1 = rnorm(n, 3, 2))
}

source("sim/a_sim_func2.r")
sim1 <- generate_survival_sim(
        df_covs,
        t_max = t_max,
        shape = shape,
        eq_tx_prob = "qnorm(0.5) + .3 * (x1 - mean(x1))",
        eq_scale_0 = glue("
        t_max*{base_scale_mult} * 
        exp(-(
        x1  * log(1.8)
            ) / shape)"), # no treatment effect before they are treated
        eq_scale_t = glue("t_max*{base_scale_mult} * 
        exp(-(
        tx*log({tx_coef}) + x1*log(1.8)
            ) / shape)"),
        eq_scale_t2 = glue("t_max*{base_scale_mult} * 
        exp(-(
        tx*log({tx_coef}) + x1*log(1.8)
            ) / shape)"),
        # alpha_beta = c(alpha_coef, beta_coef),
        alpha_beta = c(alpha_coef,beta_coef),
        max_tx_day = 10,
        min_supp_quant = 0.1,
        t2_time = t2_time_coef
    )
sim1$obs_data

# true_marg_dyn_hr: This is the exact number your CCW Cox model is trying to target. It natively bakes in the 10-day delay penalty, AND it accounts for the changing population composition over time.
# true_marg_stat_hr: This is the exact number your Time-Varying BART model is trying to target when enforcing a strict Day 0/Min-Support rule.
# cond_stat_hr: This will perfectly output the 0.5 to 0.9 biological shifts you hard-coded into the scales, completely ignoring survivor bias.
# epi convention is that tx starts at the beginning of the day and death occurs at the end of the day.
# so tx=5, event=5 (iow. same-day tx and event) implies that the patient took the tx before the event occurs

source("sim/a_sim_func2.r")
source("sim/tmle_ccw.r")
hr_lst = tibble()
rd_lst = tibble()
rdtmle_lst = tibble()
for(i in 1:2){
    i=i*6
    {
        n <- 300
        # tx_coef = runif(1,0,5)
        tx_coef = .5
        interact_coef = 2
        alpha_coef = 5
        beta_coef = 3
        t2_time_coef = 15
        t_max = 30
        shape = 1.2
        base_scale_mult = 3
        # df_covs <- data.frame(x1 = rnorm(n, 3, 2))
    }
    df_covs <- data.frame(x1 = rnorm(n, 3, 2))
    sim1 <- generate_survival_sim(
        df_covs,
        t_max = t_max,
        shape = shape,
        # eq_tx_prob = "qnorm(0.5) + .3* (x1 - mean(x1))",
        # eq_tx_prob = "qnorm(0.4) + 0.3 * (x1 - 3)^2 - 1.2",
        eq_tx_prob = "qnorm(0.4) + 0.2 * (x1 - 3)",
        # eq_scale_0 = glue("t_max*{base_scale_mult} * exp(-((x1-3)^2  * log(1.8)) / shape)"), # no treatment effect before they are treated
        eq_scale_0 = glue("t_max*{base_scale_mult} * exp(-(x1 * log(1.2)) / shape)"),
        # eq_scale_t = glue("t_max*{base_scale_mult} * exp(-(tx*log({tx_coef}) + (x1-3)^2*log(1.8) + x1 * tx * log(.2)) / shape)"),
        eq_scale_t = glue("t_max*{base_scale_mult} * exp(-(tx * (x1 > 4.5) * log(0.2) + x1 * log(1.2)) / shape)"),
        eq_scale_t2 = glue("t_max*{base_scale_mult} * exp(-(tx * (x1 > 4.5) * log(0.2) + x1 * log(1.2)) / shape)"),
        # eq_scale_t2 = glue("t_max*{base_scale_mult} * exp(-(tx*log({tx_coef}) + (x1-3)^2*log(1.8) + x1 * tx * log(1)) / shape)"),
        # alpha_beta = c(alpha_coef, beta_coef),
        alpha_beta = c(alpha_coef,beta_coef),
        max_tx_day = 10,
        min_supp_quant = 0.1,
        t2_time = t2_time_coef
    )
    cph1 = run_cph_suite(sim1, eq_cph="tx_ever + x1", landmark_day = 11)
    cph_hr = cph1$result
    hr_lst = rbind(hr_lst, cph_hr)

    
    #################################################
    # CCW
    #################################################
    dtruth = sim1$est %>%
        select(day, marg_ntx_surv, marg_dtx_surv, drd_surv)
    ccw = run_ccw_analysis(sim1$obs_data)
    ccw_res = ccw$res_log
    
    #################################################
    # CCW + TMLE
    #################################################
    ccwtmle = run_ccw_tmle_analysis(
        sim1$obs_data, 
        use_sl = TRUE,
        sl_lib = c("SL.mean", "SL.glm", "SL.earth", "SL.xgboost")) #"SL.ranger",
    ccwtmle_res = ccwtmle$res_tmle 
    
    #################################################
    # Metrics
    #################################################
    rd_out = dtruth %>% 
        left_join(ccw_res, by= c("day" = "t")) %>% 
        mutate(dd = drd_surv-rd) %>% 
        select(day, dd)
    rd_lst = bind_rows(rd_lst, rd_out)

    rd_tmle_out = dtruth %>% 
        left_join(ccwtmle_res, by= c("day" = "t")) %>% 
        mutate(dd = drd_surv-rd) %>% 
        select(day, dd)
    rdtmle_lst = bind_rows(rdtmle_lst, rd_tmle_out)        
}

plot(ccw$res_log$t, ccw$res_log[["arm=control_no_tx"]], ylim=c(0,1))
points(ccw$res_log$t, ccw$res_log[["arm=treat_tx"]])
points(sim1$est$day, sim1$est$marg_dtx_surv,col="red")
points(sim1$est$day, sim1$est$marg_ntx_surv, col="red")

points(ccwtmle$res_tmle$t, ccwtmle$res_tmle$control_no_tx, ylim=c(0,1), col="blue")
points(ccwtmle$res_tmle$t, ccwtmle$res_tmle$treat_tx, col='blue')
# points(sim1$est$day, sim1$est$marg_dtx_surv,col="red")
# points(sim1$est$day, sim1$est$marg_ntx_surv, col="red")



rd_lst %>% 
    group_by(day) %>% 
    summarise(round(mean(dd),3), round(var(dd),3)) %>% print(n=30)



rdtmle_lst %>% 
    group_by(day) %>% 
    summarise(round(mean(dd),3), round(var(dd),3)) %>% 
    print(n=30)





# hr_lst %>% 
#     group_by(method) %>% 
#     summarise(mean(hr_estimate), var(hr_estimate))

# dtruth = sim1$est %>% filter(day %in% c(5,10,15,20,25,30))  %>% 
#     select(day, marg_ntx_surv, marg_dtx_surv, drd_surv)
    

# ccwtmle_res$rd_ci_lower
# ccwtmle_res$rd_ci_upper



# run_bart_survival(sim1, )
# source("sim/a_sim_func2.r")
# run_bart_survival(
#         sim1,
#         covars = c("x1"),
#         n_tree = 50,
#         n_cores = 8,
#         trunc = 1,
#         k_p = 1.5, # turn down k to keep the unknowns drifting back to 1
#         nskip=100,
#         keepevery = 5,
#         grace=10
#     )

# ccw = run_ccw_analysis(sim1$obs_data)
# ccw_res = ccw$res_log %>% filter(t %in% c(5,10,15,20,25,30))
# dtruth %>% 
#     left_join(ccw_res, by= c("day" = "t")) %>% 
#     mutate(dd = drd_surv-rd) %>% 
#     select(day, dd)


# source("sim/tmle_ccw.r")
# ccwtmle = run_ccw_tmle_analysis(sim1$obs_data)
# ccwtmle$res_log %>% filter(t %in% c(5,10,15,20,25,30))

# source("sim/a_sim_func2.r")
# ccw = run_ccw_analysis(sim1$obs_data)
# ccw$res_log %>% filter(t %in% c(5,10,15,20,25,30))


source("sim/a_sim_func2.r")
source("sim/tmle_ccw.r")
hr_lst = tibble()
rd_lst = tibble()
rdtmle_lst = tibble()
for(i in 1:2){
    n <- 3000
    tx_coef <- 0.2     
    alpha_coef <- 5
    beta_coef <- 3
    t2_time_coef <- 15 
    t_max <- 30
    shape <- 1.2
    base_scale_mult <- 3

    # Add x2 as a binary comorbidity (e.g., Diabetes)
    df_covs <- data.frame(
        x1 = rnorm(n, 3, 2),    # Age proxy
        x2 = rbinom(n, 1, 0.5)  # Comorbidity
    )

    sim1 <- generate_survival_sim(
        df_covs,
        t_max = t_max,
        shape = shape,
        
        # THE TRAP: Doctors ONLY heavily prescribe to Frail AND Diabetic patients.
        # (x1 > 4) * x2 evaluates to 1 ONLY if both conditions are met.
        # Baseline tx prob is low (qnorm(0.1)), but jumps massively for this subgroup.
        eq_tx_prob = "qnorm(0.1) + 2.0 * (x1 > 4) * x2",
        
        # THE TRAP: That exact same subgroup has a massive baseline mortality risk.
        # We add a tiny bit of main effect (0.1 * x1) just to give the GLM some bait.
        eq_scale_0 = glue("
        t_max*{base_scale_mult} * exp(-(
            1.5 * (x1 > 4) * x2 + 0.1 * x1
        ) / shape)"), 
        
        # Early Phase (Days 1-15): Strong Treatment Effect
        eq_scale_t = glue("
        t_max*{base_scale_mult} * exp(-(
            tx*log({tx_coef}) + 1.5 * (x1 > 4) * x2 + 0.1 * x1
        ) / shape)"),
        
        # Late Phase (Days 16-30): Waning Effect
        eq_scale_t2 = glue("
        t_max*{base_scale_mult} * exp(-(
            tx*log(1.0) + 1.5 * (x1 > 4) * x2 + 0.1 * x1
        ) / shape)"),
        
        alpha_beta = c(alpha_coef, beta_coef),
        max_tx_day = 10,
        min_supp_quant = 0.1,
        t2_time = t2_time_coef
    )
    
     #################################################
    # CCW
    #################################################
    dtruth = sim1$est %>%
        select(day, marg_ntx_surv, marg_dtx_surv, drd_surv)
    ccw = run_ccw_analysis(sim1$obs_data, covars=c("x1","x2"))
    ccw_res = ccw$res_log
    
    #################################################
    # CCW + TMLE
    #################################################
    ccwtmle = run_ccw_tmle_analysis(
        sim1$obs_data, 
        use_sl = TRUE,
        covars=c("x1","x2"),
        sl_lib = c("SL.mean", "SL.glm", "SL.earth", "SL.xgboost")) #"SL.ranger",
    ccwtmle_res = ccwtmle$res_tmle 
    
    #################################################
    # Metrics
    #################################################
    rd_out = dtruth %>% 
        left_join(ccw_res, by= c("day" = "t")) %>% 
        mutate(dd = drd_surv-rd) %>% 
        select(day, dd)
    rd_lst = bind_rows(rd_lst, rd_out)

    rd_tmle_out = dtruth %>% 
        left_join(ccwtmle_res, by= c("day" = "t")) %>% 
        mutate(dd = drd_surv-rd) %>% 
        select(day, dd)
    rdtmle_lst = bind_rows(rdtmle_lst, rd_tmle_out)     
}
plot(ccw$res_log$t, ccw$res_log[["arm=control_no_tx"]], ylim=c(0,1))
points(ccw$res_log$t, ccw$res_log[["arm=treat_tx"]])
points(sim1$est$day, sim1$est$marg_dtx_surv,col="red")
points(sim1$est$day, sim1$est$marg_ntx_surv, col="red")

points(ccwtmle$res_tmle$t, ccwtmle$res_tmle$control_no_tx, ylim=c(0,1), col="blue")
points(ccwtmle$res_tmle$t, ccwtmle$res_tmle$treat_tx, col='blue')
# points(sim1$est$day, sim1$est$marg_dtx_surv,col="red")
# points(sim1$est$day, sim1$est$marg_ntx_surv, col="red")



rd_lst %>% 
    group_by(day) %>% 
    summarise(round(mean(dd),3), round(var(dd),3)) %>% print(n=30)



rdtmle_lst %>% 
    group_by(day) %>% 
    summarise(round(mean(dd),3), round(var(dd),3)) %>% 
    print(n=30)




##########################################
source("sim/a_sim_func2.r")
source("sim/tmle_ccw_anno.r")
hr_lst = tibble()
rd_lst = tibble()
rdtmle_lst = tibble()
for(i in 1:1){
    n <- 2000
    
    # tx_coef = 0.03 acts as the "Miracle Cure". log(0.03) is ~ -3.5. 
    # This perfectly cancels out the massive +3.5 lethal penalty of the disease.
    tx_coef <- 0.03     
    alpha_coef <- 5
    beta_coef <- 3
    t2_time_coef <- 10
    t_max <- 30
    shape <- 1.2
    
    # Bumped base scale slightly so the healthy people survive nicely
    base_scale_mult <- 4 

    # x1 is Age-proxy, x2 is Diabetes
    df_covs <- data.frame(
        x1 = rnorm(n, 3, 2),    
        x2 = rbinom(n, 1, 0.5)  
    )

    sim1 <- generate_survival_sim(
        df_covs,
        t_max = t_max,
        shape = shape,
        
        # THE TRAP: Base rate is 5%. But if Frail AND Diabetic, rate is 97%.
        eq_tx_prob = "qnorm(0.05) + 3.5 * (x1 > 4) * x2",
        
        # THE LETHAL BASELINE: +3.5 penalty for the subgroup. They die fast.
        eq_scale_0 = glue("
        t_max*{base_scale_mult} * exp(-(
            3.5 * (x1 > 4) * x2 - 0.1 * x1
        ) / shape)"), 
        
        # EARLY PHASE: Miracle Cure. tx_coef neutralizes the 3.5 penalty.
        eq_scale_t = glue("
        t_max*{base_scale_mult} * exp(-(
            tx*log({tx_coef}) + 3.5 * (x1 > 4) * x2 - 0.1 * x1
        ) / shape)"),
        
        # LATE PHASE: The Wane. Drug stops working, lethal penalty returns.
        eq_scale_t2 = glue("
        t_max*{base_scale_mult} * exp(-(
            tx*log(1.0) + 3.5 * (x1 > 4) * x2 - 0.1 * x1
        ) / shape)"),
        
        alpha_beta = c(alpha_coef, beta_coef),
        max_tx_day = 10,
        min_supp_quant = 0.1,
        t2_time = t2_time_coef
    )
     #################################################
    # CCW
    #################################################
    dtruth = sim1$est %>%
        select(day, marg_ntx_surv, marg_dtx_surv, drd_surv)
    ccw = run_ccw_analysis(sim1$obs_data, covars=c("x1","x2"))
    ccw_res = ccw$res_log
    
    #################################################
    # CCW + TMLE
    #################################################
    sim1$obs_data %>% head()
    source("sim/tmle_ccw_anno.r")
    ccwtmle = run_ccw_tmle_analysis2(
        sim1$obs_data, 
        # use_sl = TRUE,
        covars=c("x1","x2"),
        sl_lib = c("SL.mean", "SL.glm", "SL.earth", "SL.xgboost")) #"SL.ranger",
    ccwtmle_res = ccwtmle$res_tmle 
    
    #################################################
    # Metrics
    #################################################
    rd_out = dtruth %>% 
        left_join(ccw_res, by= c("day" = "t")) %>% 
        mutate(dd = drd_surv-rd) %>% 
        select(day, dd)
    rd_lst = bind_rows(rd_lst, rd_out)

    rd_tmle_out = dtruth %>% 
        left_join(ccwtmle_res, by= c("day" = "t")) %>% 
        mutate(dd = drd_surv-rd) %>% 
        select(day, dd)
    rdtmle_lst = bind_rows(rdtmle_lst, rd_tmle_out)     
}
hr_lst
ccwtmle$res_tmle









library(ggplot2)
library(dplyr)
library(tidyr)
sim1$obs_data %>% head()
f1 = glm(tx ~ x1 + x2 , famil=binomial, sim1$obs_data)
p1 = predict(f1, type="response")

p1 %>% head()
sum(round(sim1$obs_data$tx-p1,2) < 0)
plot(p1, jitter(sim1$obs_data$tx))

# 1. Reconstruct the Average Estimated Curves
# We take the True Risk Difference from one of your simulation objects
true_rd <- sim1$est %>%
  select(day, true_rd = drd_surv)

# Calculate the average GLM estimate across your simulation loop
glm_summary <- rd_lst %>%
  group_by(day) %>%
  summarise(glm_bias = mean(dd, na.rm = TRUE)) %>%
  left_join(true_rd, by = "day") %>%
  mutate(glm_est = true_rd + glm_bias)

# Calculate the average TMLE estimate across your simulation loop
tmle_summary <- rdtmle_lst %>%
  group_by(day) %>%
  summarise(tmle_bias = mean(dd, na.rm = TRUE)) %>%
  left_join(true_rd, by = "day") %>%
  mutate(tmle_est = true_rd + tmle_bias)

# 2. Combine and format the data for ggplot
plot_data <- true_rd %>%
  left_join(glm_summary %>% select(day, glm_est), by = "day") %>%
  left_join(tmle_summary %>% select(day, tmle_est), by = "day") %>%
  pivot_longer(
    cols = c(true_rd, glm_est, tmle_est),
    names_to = "Model",
    values_to = "Risk_Difference"
  ) %>%
  mutate(
    Model = factor(Model,
                   levels = c("true_rd", "tmle_est", "glm_est"),
                   labels = c("True Marginal Effect", "TMLE + SuperLearner", "Parametric IPW-GLM"))
  )

# 3. Generate the Presentation-Ready Plot
p <- ggplot(plot_data, aes(x = day, y = Risk_Difference, color = Model, linetype = Model)) +
  
  # A. The Structural Break Annotation
#   geom_vline(xintercept = 15, linetype = "dotted", color = "gray50", size = 1) +
#   annotate("text", x = 15.5, y = max(plot_data$Risk_Difference, na.rm=TRUE) * 0.9,
        #    label = "Treatment Effect Wanes", hjust = 0, color = "gray50", fontface = "italic") +

  # B. The Curves
  geom_line(size = 1.2) +

  # C. Colors matching your description (Red = True, Blue = TMLE, Black = GLM)
  scale_color_manual(values = c("True Marginal Effect" = "#e41a1c",
                                "TMLE + SuperLearner" = "#377eb8",
                                "Parametric IPW-GLM" = "#000000")) +
  
  # Making GLM dashed helps visual separation when lines cross
  scale_linetype_manual(values = c("True Marginal Effect" = "solid",
                                   "TMLE + SuperLearner" = "solid",
                                   "Parametric IPW-GLM" = "dashed")) +
  
  # D. Labels and Formatting
  labs(
    # title = "Causal Survival Analysis: The Masked Miracle Cure",
    # subtitle = "TMLE dynamically corrects for masked synergistic confounding while GLM fails after structural break.",
    x = "Days of Follow-up",
    y = "Marginal Risk Difference (Treated - Control)",
    color = NULL, # Hide legend title for cleaner look
    linetype = NULL
  ) +
  
  # E. Minimalist Theme (Perfect for Slides)
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 18),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(), # Cleans up the x-axis
    legend.text = element_text(size = 14)
  )

# Display the plot
print(p)

# Optional: Save for your presentation
# ggsave("tmle_vs_glm_presentation.png", plot = p, width = 10, height = 6, dpi = 300)











# scenario 2
###############################################################################
hr_lst = tibble()
rd_lst = tibble()
rdtmle_lst = tibble()
for(i in 1:1){
    n <- 800 # Small sample size keeps statistical power on the edge!
    
    # THE TRUE SIGNAL: The drug really works (40% hazard reduction).
    tx_coef <- 0.6     
    alpha_coef <- 5
    beta_coef <- 3
    t_max <- 30
    shape <- 1.2
    base_scale_mult <- 3 

    # x1 is a massive biological prognostic factor (e.g., tumor grade)
    df_covs <- data.frame(
        x1 = rnorm(n, 0, 1)    
    )

    sim1 <- generate_survival_sim(
        df_covs,
        t_max = t_max,
        shape = shape,
        
        # MILD LINEAR CONFOUNDING: GLM gets the weights perfectly correct. 
        # The GLM will NOT be biased in this scenario.
        eq_tx_prob = "qnorm(0.5) + 0.2 * x1",
        
        # MASSIVE PROGNOSTIC NOISE: x1 heavily dictates baseline survival.
        # This creates the giant variance that will blind the GLM.
        eq_scale_0 = glue("
        t_max*{base_scale_mult} * exp(-(
            3.0 * x1
        ) / shape)"), 
        
        # TRUE TREATMENT EFFECT: tx_coef = 0.6 reduces the hazard.
        eq_scale_t = glue("
        t_max*{base_scale_mult} * exp(-(
            tx*log({tx_coef}) + 3.0 * x1
        ) / shape)"),
        
        # NO WANING EFFECT: Keep it simple to isolate precision
        eq_scale_t2 = glue("
        t_max*{base_scale_mult} * exp(-(
            tx*log({tx_coef}) + 3.0 * x1
        ) / shape)"),
        
        alpha_beta = c(alpha_coef, beta_coef),
        max_tx_day = 10,
        min_supp_quant = 0.1,
        t2_time = 15 # Doesn't matter here since tx_coef is the same in both phases
    )
     #################################################
    # CCW
    #################################################
    dtruth = sim1$est %>%
        select(day, marg_ntx_surv, marg_dtx_surv, drd_surv)
    ccw = run_ccw_analysis(sim1$obs_data, covars=c("x1"))
    ccw_res = ccw$res_log
    
    #################################################
    # CCW + TMLE
    #################################################
    ccwtmle = run_ccw_tmle_analysis2(
        sim1$obs_data, 
        # use_sl = TRUE,
        covars=c("x1"),
        sl_lib = c("SL.mean", "SL.glm", "SL.earth"))#, "SL.xgboost")) #"SL.ranger",
    ccwtmle_res = ccwtmle$res_tmle 
    ccwtmle$res_log %>% tail()
    ccwtmle$res_tmle %>% tail()
     #################################################
    # Metrics
    #################################################
    rd_out = dtruth %>% 
        left_join(ccw_res, by= c("day" = "t")) %>% 
        mutate(dd = drd_surv-rd) %>% 
        select(day, dd)
    rd_lst = bind_rows(rd_lst, rd_out)

    rd_tmle_out = dtruth %>% 
        left_join(ccwtmle_res, by= c("day" = "t")) %>% 
        mutate(dd = drd_surv-rd) %>% 
        select(day, dd)
    rdtmle_lst = bind_rows(rdtmle_lst, rd_tmle_out)     
}
dtruth %>% tail()



ccw
ccwtmle$res_tmle %>% tail()
rd_lst %>% 
    group_by(day) %>% 
    summarise(round(mean(dd),3), round(var(dd),3)) %>% print(n=30)



rdtmle_lst %>% 
    group_by(day) %>% 
    summarise(round(mean(dd),3), round(var(dd),3)) %>% 
    print(n=30)


library(ggplot2)
library(dplyr)
library(tidyr)

# 1. Prepare the Data
# We extract the True Risk Difference (drd_surv) from the simulation object
true_rd <- sim1$est %>%
  select(day, true_rd = drd_surv)

# Format the GLM results (Assuming 'glm_results' is the output from 'get_glm_ccw_ci')
# If you are using a loop, this would be the averaged results
glm_plot_data <- ccwtmle$res_log %>%
    select(day=t, rd, ci_lower, ci_upper) %>%
    mutate(Model = "Parametric IPW-GLM")

# Format the TMLE results (Assuming 'ccwtmle$res_tmle' is the output)
# Note: Ensure the column names match your TMLE output structure
tmle_plot_data <- ccwtmle$res_tmle %>%
  select(day = t, rd, ci_lower = rd_ci_lower, ci_upper = rd_ci_upper) %>%
  mutate(Model = "TMLE + SuperLearner")

# Combine the GLM and TMLE data, then join with the True RD
plot_data <- bind_rows(glm_plot_data, tmle_plot_data) %>%
  left_join(true_rd, by = "day") %>%
  # Reorder factors so TMLE plots on top of GLM if they overlap
  mutate(Model = factor(Model, levels = c("TMLE + SuperLearner", "Parametric IPW-GLM")))

# 2. Generate the Presentation-Ready Plot
p <- ggplot(plot_data, aes(x = day)) +
  
  # A. The Null Hypothesis Line (y = 0)
  # This is crucial for showing that the GLM CI crosses zero while TMLE does not.
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30", size = 1) +
  annotate("text", x = 15, y = -0.005, label = "Null Hypothesis (No Effect)", 
           color = "gray30", fontface = "italic", hjust = 0) +

  # B. The Confidence Interval Ribbons
  # We use alpha for transparency so you can see both models
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = Model), alpha = 0.1) +

  # C. The Point Estimates (Lines)
  geom_line(aes(y = rd, color = Model, linetype = Model), size = 1.2) +

  # D. The True Signal (Red Line)
  geom_line(aes(y = true_rd), color = "#e41a1c", size = 1, linetype = "solid") +
  
  # E. Colors matching your description
  scale_color_manual(values = c("TMLE + SuperLearner" = "#377eb8",
                                "Parametric IPW-GLM" = "#000000")) +
  scale_fill_manual(values = c("TMLE + SuperLearner" = "#377eb8",
                               "Parametric IPW-GLM" = "#000000")) +
  scale_linetype_manual(values = c("TMLE + SuperLearner" = "solid",
                                   "Parametric IPW-GLM" = "dashed")) +
  
  # F. Labels and Formatting
  labs(
    # title = "Rescuing the Signal: Gaining Statistical Power",
    # subtitle = "TMLE learns prognostic factors, shrinking the CI and preventing a false negative.",
    x = "Days of Follow-up",
    y = "Marginal Risk Difference (Treated - Control)",
    color = NULL,
    fill = NULL,
    linetype = NULL
  ) +
  
  # G. Minimalist Theme (Perfect for Slides)
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 18),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.text = element_text(size = 14)
  )

# Display the plot
print(p)

# Optional: Save for your presentation
# ggsave("tmle_efficiency_gain.png", plot = p, width = 10, height = 6, dpi = 300)
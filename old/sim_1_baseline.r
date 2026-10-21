
library(survival)
library(tidyverse)
library(survminer)
library(BART)
# library(BART) # if bart3 not installed


# sv, SV, Sv all equal survival (probability/function)
# hz, HZ, Hz all equal hazard (function)

# Set-Up
# Generate the event time, along with the true Survival and Hazard using
#   a weibull function that includes:
#   a conditional switch of the underlying hazard as determined by the 
#   time varying/delayed start covariate.
#   Only incoproartes 1 covariate, 1 timepooint switch at the moment

# - scale parameter is the time in which 62.3% of the population in that group will have an event
#   larger scale == longer time to event
#   The parameterization focuses on switching the scale value at the time of tx to reflect the 
#   hazard under the new model.

get_weibull_full_trajectory = function(
    shape,          # shape param
    scale_0,        # baseline scale param (sets baseline hazard)
    scale_t,        # adjusted scale parame (sets hazard after the time-delayed covariate is introduced)
    t_switch,       # time that the switch occurs
    t_max,          # maximum number of days to generate true sv and hz
    id
    ) { 
        # draw a random hzard value (which will correspond to event time)
        # If the target hazard is <= the hazard at switch, 
        # we can say that the event occurs before the switch happens.
        # Otherwise, we say the event occurs in correspondence to
        # the combined hazard up to the switch and the hazard occuring 
        # with the new scale value.
        u <- runif(1)
        target_hazard <- -log(u)
        h_at_switch <- (t_switch / scale_0)^shape
        if (target_hazard <= h_at_switch) {
            event_time <- scale_0 * (target_hazard)^(1/shape)
        } else {
            rem_h <- target_hazard - h_at_switch
            event_time <- (rem_h * (scale_t^shape) + (t_switch^shape))^(1/shape)
        }
        # Generate the true Hazard and Survival distributions over the time range
        days <- 1:t_max
        hz <- numeric(t_max)
        sv <- numeric(t_max)
        for (t in days) {
            # pre switch
            if (t < t_switch) {
                hz[t] <- (shape / scale_0) * (t / scale_0)^(shape - 1)
                sv[t] <- exp(-(t / scale_0)^shape)
            } else {                
                # post switch
                hz[t] <- (shape / scale_t) * (t / scale_t)^(shape - 1) 
                sv[t] <- exp(- ( (t_switch/scale_0)^shape + ((t^shape - t_switch^shape)/scale_t^shape) ) )
            }
        }
        return(list(
            summary = data.frame(
                time = min(ceiling(event_time), t_max), 
                event = if(event_time <= t_max) 1 else 0 # onlyy adminstrative censoring applied
                ),
            trajectory = data.frame(day = days, hazard = hz, survival = sv, id=id)
        ))
}


# To evaluate the our estimates correctly, we need to compare 
#   against the counterfactual SV distributions.
#   These include:
#   - SV if everyone had been treated from the first day of acceptable Positivity 
#   - SV if everyone had never received treated
#   NOTE: acceptable Positivity means that it is the earliest day we have strong 
#       enough representation of the tx group
#       ie. if no one, or only 1 out of 100 people receives treatment before day 10 in our sample, 
#           then we really can't say we have support to estimate what tx would like from before day 10
#       the estimators discussed should all converge to this counterfactual scenario

# Function to get counterfactual survival vectors
get_true_survival = function(shape, scale_0, scale_t, t_min_support, t_max) {
    t_seq <- 1:t_max
    # for the counterfactual 1
    surv_treated <- ifelse(
        t_seq < t_min_support,
        exp(-(t_seq / scale_0)^shape),
        exp(-( (t_min_support/scale_0)^shape + (t_seq^shape - t_min_support^shape)/scale_t^shape ))
    )
    surv_control <- exp(-(t_seq / scale_0)^shape)
    # hazards
    hz_treated <- ifelse(
        t_seq < t_min_support,
        (shape / scale_0) * (t_seq / scale_0)^(shape - 1),
        (shape / scale_t) * (t_seq / scale_t)^(shape - 1)
    )
    hz_control = (shape / scale_0) * (t_seq / scale_0)^(shape - 1)
    return(list(
        sv_tx = surv_treated, 
        sv_ct = surv_control,
        hz_tx = hz_treated,
        hz_ct = hz_control)
        )
}



####################################
# Simulation 1
# Primary estimand - Survival Risk Difference as the Average Treatmetn Effect at Day 30
# Secondary estimand - Hazard Ratio 
####################
seed=102
n <- 300
shape <- 1              # keep the hazard constant over time
t_max=30                # keep the 




# create base data
x1 <- rnorm(n,mean=0,sd=1)
tx <- rbinom(n, 1,  pnorm(qnorm(.6) + x1) ) 
tt_tx = ceiling(rbeta(n,3,5) * 10)/tx
hist(tt_tx)
table(tx, tt_tx)
table(tx, tt_tx<Inf)

# overlap  of the tx and control x1 distributions
hist(x1[tx==1], xlim=c(-3,3))
hist(x1[tx==0], add=TRUE, col="red")

# scale parameters (bigger == lower hazard == lower risk of event)
s0 <- 30 * exp(-(x1 * log(2)) / shape)              #NOTE: tx is not included in the base hazard s0
st <- 30 * exp(-(tx * log(.5) + x1 * log(2)) / shape)

# initalize result lists
results_list <- list() 
sv_hz = list()
for(i in 1:n) {
  res = get_weibull_full_trajectory(shape, s0[i], st[i], tt_tx[i], t_max, id=i)
  # Store metadata
  res$id <- i
  res$x1 <- x1[i]
  res$t_switch <- tt_tx[i]
  sv_hz[[i]] = list(res$trajectory)
  results_list[[i]] <- data.frame(res$summary, id=res$id, x1=res$x1, t_switch=res$t_switch)
}

# reformat sim_data and true sv_hz
sim_data = tibble(bind_rows(results_list))
sv_hz_c = tibble(bind_rows(sv_hz)) %>% 
    left_join(sim_data)

# get counterfactuals
# counterfactuals shape only varies on the x1 variable
cf_s0 =  30 * exp(-(x1 * log(2)) / shape)              
cf_st = 30 * exp(-(1 * log(.5) + x1 * log(2)) / shape)
cntfct_df = tibble()
(min_support= quantile(sim_data$t_switch, .1))
for(i in 1:n){
  cntfct_df = bind_rows(
    cntfct_df, 
    as_tibble(
      get_true_survival(shape, scale_0=cf_s0[i], scale_t=cf_st[i], t_min_support = min_support, t_max=t_max) 
     ) %>% 
     mutate(id=i, t= row_number(), x1 = x1[i], tt_tx = tt_tx[i])
  )
}


#########################
# Review obs and couunterfactual data
all_rl = cntfct_df %>% 
  left_join(sv_hz_c %>% select(hazard, survival, id, day), by=c("id"="id", "t"="day")) %>% 
  mutate(tx_td = as.numeric(t>=tt_tx)) %>% 
  mutate(tx_obs = as.numeric(tt_tx != Inf))


# Preparre plot data
tmp = all_rl %>% 
  group_by(t, tx_obs)%>% 
  summarise(
    sv_tx =mean(sv_tx), 
    sv_ct = mean(sv_ct),
    sv_obs = mean(survival)
  )
# marginal counterfactual estimates
tmp2 = cntfct_df%>% 
  group_by(t) %>% 
  summarize(svt = mean(sv_tx),svc = mean(sv_ct))

tmp
msk = tmp$tx_obs==1
# add in true and counterfactuals conditioned on observed tx status
# tx are initiated tx at the minimum support date
plot(tmp[msk,]$t, tmp[msk,]$sv_obs, ylim=c(0,1))  
points(tmp[!msk,]$t, tmp[!msk,]$sv_obs, ylim=c(0,1))  
points(tmp[msk,]$t, tmp[msk,]$sv_tx, ylim=c(0,1), col="red")  
points(tmp[!msk,]$t, tmp[!msk,]$sv_ct, ylim=c(0,1), col="red")  
# add in marginal counterfactual estimates
points(tmp2$t, tmp2$svt, col="blue")
points(tmp2$t, tmp2$svc, col="blue")


#######################
# write about theseplots
#########

####################
# NOTE: min-support is obvious in the counterfactuals
# prior to time 3, both tx and ct are aligned and after day 3 they diverge
##################

# counterfactual ATE and HR
cntfct_sm = cntfct_df %>% 
  mutate(
    diff = sv_tx-sv_ct,
    quo = hz_tx/hz_ct
  ) %>% 
  group_by(t) %>% 
  summarize(mean(diff), mean(quo), sv_tx = mean(sv_tx), sv_ct=mean(sv_ct)) %>% 
  print(n=30)


obs_sm = sv_hz_c %>% 
  group_by(id) %>%
  mutate(tx = as.numeric(any(t_switch!=Inf))) %>% 
  # print(n=400)
  group_by(day, tx) %>% 
  summarize(sv = mean(survival)) %>% 
  pivot_wider(names_from=tx, values_from=sv) %>%
  rename(tx=`1`, ct=`0`)



#####################################################################################
############
# Cox Models
###########
# Strategy 1, t0 tx assignment
fit1 <- coxph(
  Surv(time, event) ~ tx + x1, 
  data = sim_data %>% mutate(tx = as.numeric(t_switch != Inf))
  )
summary(fit1)
###########################
# Warning !!!!!!!!!! DATA LEAK
sim_data %>% 
  filter(time < t_switch) %>% 
  filter(t_switch != Inf)
# this result includes patients who died before the treated
# while this scenario wouldn't be possible in un-simulated data
# it does highlight the issues a simple cox modeling approach
# because there are many events at day 1, this improper 
# flagging of tx individuals creates a false association of risk 
# closer to the day 1 time (all individuals are being said to have been tx
# at day 0). Increase association would increase the hazard from the mitigating
# effect we see overall in this tx scenario
###########################


# straegy 2, tx only if t_switch >= event
# this is how we would observe treated
sim_data %>% mutate(tx = as.numeric(t_switch != Inf)) %>% 
  group_by(event) %>% count(tx)
# create the observed tx data
obs_sim_data = sim_data %>% mutate(tx = ifelse(t_switch <= time, 1,0))
obs_sim_data %>% group_by(event) %>% count(tx)
fit2 <- coxph(
  Surv(time, event) ~ tx + x1, 
  data = obs_sim_data
  )
summary(fit2)
# this is the classic immortal time bias case
# overly confident in the effect that there are few people to die
# only people who are treated are those who didn't already die 
# and as we can see most of the people who were counterfactually
# treated, were counterfactually dead by the time tx occured.


##################################
# strategy 3, time-based matching
  # doing an easy caliper match because only 1 covariate
cph_tx = obs_sim_data %>% filter(tx==1) 
cph_cntrl = obs_sim_data %>% filter(tx==0)
cph_match_df = tibble()
for(i in 1:nrow(cph_tx)) {
  tx_obs = cph_tx[i,]
  cntrl_match = cph_cntrl %>% 
    filter(
      time >= tx_obs$t_switch &
      abs(x1 - tx_obs$x1) < 0.3
    ) 
  if(nrow(cntrl_match) > 0) {
    cph_match_df = cph_match_df %>% 
      bind_rows(cph_tx[i,]) %>% 
      bind_rows(cntrl_match %>% sample_n(10, replace=TRUE))
  } 
}


# Check the selection bias (if any on the observed variables)
obs_sim_data %>% group_by(tx) %>% reframe(quantile(x1, c(.1,.5,.9)))  # original
cph_match_df %>% group_by(tx) %>% reframe(quantile(x1,c(.1,.5,.9)))   # matched
fit3 = coxph(Surv(time, event) ~ tx + x1, data= cph_match_df)
summary(fit3)
# So even with matching on time and the single covariate we observe
# a very strong inhibitory effect estimate of the tx, that is closer to 
# the immortal time biased effect, than the true result.
# The true result is still 0.3, since in this example it hazards are time-invariant
# and there is no effect heterogeneity.
obs_sim_data %>% group_by(tx) %>% count(event) %>% mutate(perc = n/sum(n))
cph_match_df %>% group_by(tx) %>% count(event) %>% mutate(perc = n/sum(n))
# I think what we should have seen was all of the early events
# from the cntrl group be removed, which would have decreased the effective
# hazard rate of the cntrl group, making it closer to that of the tx group
# and moving the HR closer to the true HR 0.3 and away from the naive estimate.

##################################
# strategy 4, landmark
# we choose a landmark day of 11 because we know that all individuals have been treated
cph_lndmk = obs_sim_data %>% 
  filter(time >= 11)
cph_lndmk %>% group_by(tx) %>% count(event) %>% mutate(perc = n/sum(n))
fit4 = coxph(Surv(time, event) ~ tx + x1, data= cph_lndmk)
summary(fit4)
# removes immortal time  bia by remove individuals by ensuring we only observe patients
# who have the event after the full tx window (1-10)


######################################
# strategy 5, time-varying
# Setup the time-varying dataset
sim_long <- tmerge(obs_sim_data, obs_sim_data, id=id, tstop=time, evt=event(time, event))
sim_long <- tmerge(sim_long, obs_sim_data, id=id, tx_td = tdc(t_switch))
# tibble(sim_long)
# Fit the model
fit5 <- coxph(Surv(tstart, tstop, evt) ~ tx_td + x1, data = sim_long)
print(summary(fit5))


#################################################
# BART
options(width=400)
# library(BART3)
trunc = 1

# create base bart sim data
# this is a long format dataset, where each time-point (day) is its own observation/person
# this format is for all OBSERVED EVENT DAYs, ie events occur on day 1,3,4,5,8,9
# the dataset will have 6 ovservations for each individual
times = sort(unique(obs_sim_data$time))
k = length(times)
b_sim_data = tibble(t = rep(times, n), id=rep(1:n, each=k)) %>% 
  left_join(obs_sim_data, by="id") %>% 
  mutate(tx_td = ifelse(t_switch<=times,1,0)) 

# create bart training data
# from the long format we keep only the observation days upto and including the event(censor) time
# all times prior to the event time
b_train = b_sim_data %>% 
  filter(t<=time) %>% 
  mutate(event_td = ifelse(time == t, event, 0)) 


# Train model
post = BART::mc.surv.bart(
    x.train = data.frame(b_train %>% select(t, x1, tx_td)),
    y.train = b_train$event_td,
    k=1,
    ntree=25, 
    keepevery=10,
    mc.cores=8
)

# for windows you can't use parallel
# post = BART::surv.bart(
#     x.train = data.frame(b_train %>% select(t, x1, tx_td)),
#     y.train = b_train$event_td,
#     k=1,
#     ntree=50, 
#     keepevery=10,
# )
# generate predictions for the observed 
# predictions are generated for all observed days (even those after the event)
p1 = predict(post, data.frame(b_sim_data %>% select(t, x1, tx_td)))

# Lets look at these estimates first, see how they align with the expected observed
res1 = b_sim_data %>% 
  mutate(b_preds = p1$surv.test.mean)  %>% 
  # mutate(c_preds=as.vector(fits$surv)) %>% 
  left_join(
    all_rl %>% select(sv_tx, sv_ct, hz_tx, hz_ct, sv_obs=survival, hz_obs=hazard, id,t)
    , by=c("id","t")
    )
res2 = res1 %>% 
  group_by(tx, t) %>% 
  summarise(b = mean(b_preds), o = mean(sv_obs))

# Plot Predicted vs Obseved (under the observed x1 values)
res2 %>%
  pivot_longer(cols = c(b, o), names_to = "type", values_to = "value") %>%
  ggplot(aes(x = t, y = value, color = type)) +
  geom_point() +
  scale_color_manual(values = c("b" = "blue", "o" = "red"), 
                     labels = c("b" = "BART", "o" = "Observations")) +
  ylim(0, 1) +
  theme_minimal()


# now we wantt to predict under the counterfactuals,
# this is how we want to generate our HR and RD estimates
z0_set = times[times<min_support]
z1_set = times[times>=min_support]
# Predicty for treated
# Remeber as in the simulated counterfactuals, the tx is expected to start at the
# minimum support value
b_pred1 = predict(post, data.frame(
  b_sim_data %>% 
    select(t, x1, tx_td) %>% 
    mutate(tx_td=rep(c(rep(0, length(z0_set)), rep(1, length(z1_set))),n))
))
b_pred2 = predict(
  post, 
  data.frame(
    b_sim_data %>% 
      select(t, x1, tx_td) %>% 
      mutate(tx_td=0)
    )
)

res3 = res1 %>% 
  mutate(
    b_sv_tx = b_pred1$surv.test.mean,
    b_sv_ct = b_pred2$surv.test.mean
    )
res4 = res3 %>% 
  group_by(t) %>%
  summarize(ot=mean(sv_tx), oc=mean(sv_ct), bt=mean(b_sv_tx), bc=mean(b_sv_ct)) %>% 
  mutate(ord = ot-oc) %>% 
  mutate(brd = bt-bc)
res4_long <- res4 %>%
  pivot_longer(
    cols = c(ot, oc, bt, bc), 
    names_to = "group", 
    values_to = "survival"
  ) %>%
  mutate(
    # Create helper columns for the legend logic
    Source = ifelse(grepl("^o", group), "Marginal-Survival", "BART"),
    Treatment = ifelse(grepl("t$", group), "Tx", "Control")
  )

ggplot(res4_long, aes(x = t, y = survival, color = Source, shape = Treatment)) +
  geom_point(size = 2.5, alpha = 0.8) +
  # scale_color_manual(values = c("Observed" = "red", "Predicted" = "blue")) +
  scale_color_discrete() +
  ylim(0, 1) +
  labs(
    title = "Survival: Observed vs Predicted",
    y = "Mean Survival",
    x = "Time (t)"
  ) +
  theme_minimal()


#############################
# get estimates
br1 = res1 %>% 
  select(t, id) %>% 
  mutate(
    as_tibble(t(b_pred1$prob.test))%>% 
    rename_with(~paste0("hz_tx_",.x))
    ) %>%
  mutate(
    as_tibble(t(b_pred2$prob.test)) %>% 
    rename_with(~paste0("hz_ct_",.x)) 
    ) %>%
  mutate(
    as_tibble(t(b_pred1$surv.test)) %>% 
    rename_with(~paste0("sv_tx_",.x)) 
    ) %>%
  mutate(
    as_tibble(t(b_pred2$surv.test)) %>% 
    rename_with(~paste0("sv_ct_",.x)) 
    ) %>%
  mutate( 
    as_tibble(t(b_pred1$prob.test/b_pred2$prob.test)) %>% 
    rename_with(~paste0("hr",.x)) 
    ) %>%
  group_by(t) %>% 
  summarise(across(contains("V"), mean))
  

# #  HR
# quo1 = br1 %>% 
#   select(contains("hr")) %>% 
#   mutate(t=br1$t) %>% 
#   pivot_longer(cols=contains("hr")) %>% 
#   filter(t >=min_support) #%>% 
#   summarise(mean(value), quantile(value,.025), quantile(value, .975)) 

quo1 = br1 %>% 
  select(contains("hr")) %>% 
  mutate(t=br1$t) %>% 
  pivot_longer(cols=contains("hr")) %>% 
  filter(t >=min_support) %>% 
  # trim for extra stability (outliers can cause extremelhy unstable estimates)
  filter(value < quantile(value, .999)) %>% 
  filter(value > quantile(value, .001)) %>% 
  summarise(mean(value), quantile(value,.025), quantile(value,.5), quantile(value, .975)) 
quo1



# Risk Difference (RD)
bdiff = (br1 %>% ungroup() %>% select(contains("sv_tx")) -
  br1 %>% ungroup() %>% select(contains("sv_ct"))) %>% 
  mutate(t=br1$t) %>% 
  pivot_longer(cols=contains("v")) %>% 
  group_by(t)%>%  
  summarise(mean(value), quantile(value,.025), quantile(value, .975)) 
# show estimate for the last few time points
# goal with time 30
bdiff %>% filter(t>25) %>% print(n=30)

# check the expected
cntfct_sm %>% filter(t>25) %>% print(n=30)


#############################
# Clone-Censor Weights
# recast sim data to match expected structure
TAU = t_max
DT=1 # interval length set to days
EPS = 0.000001
covars = c("x1")
# GRACE=min_support
GRACE=11 # GRACE is set at 11 because we know the max time of tx is day 10 in this simulation

# rename variables to match desired
dat0 = sim_data %>% 
    mutate(
        id= as.character((id)),
        x1 = as.double(x1),
        tt_tx = ifelse(t_switch == Inf, NA, t_switch),
        tx = ifelse(t_switch != Inf, 1,0),
        event = as.double(event),
        tt_event = as.numeric(time)
    ) %>% 
    select(-c(t_switch, time))

# create base_clone data
clones <- dat0 %>%    
  uncount(weights = 2, .id = "clone_arm") %>%
  mutate(
    arm = if_else(clone_arm == 1L, "cntrl", "tx")
  ) %>%
  mutate(
      tx_in_grace= (tx == 1 & !is.na(tt_tx) & tt_tx <= GRACE)
  ) %>% 
  mutate(   
      dev_time = case_when(
          arm == "cntrl" & tx_in_grace ~ tt_tx, #cntrl arm but has tx in grace
          arm == "tx" & !tx_in_grace & tt_event >= GRACE ~ GRACE, #tx arm but does not have tx in grace and has an even after grace
          TRUE ~ NA # no deviation from the arm
          )
    ) %>% 
    mutate(
        dev = if_else(!is.na(dev_time), 1, 0)
    ) %>%
  mutate(
    time_cc = if_else(!is.na(dev_time), pmin(tt_event, dev_time), tt_event), # takes min of dev_time and tt_event, so if the event occured before deviation, it won't appear
    time_cc = pmin(time_cc, TAU),   # EOS truncation
    event_cc = if_else(event == 1 & tt_event <= time_cc, 1, 0) # identifies if event occured or censored
    # event_cc is our updated event status
  )


panel0 <- clones3 %>%
  mutate(
    # artificial censoring only matters during the grace period
    tmax  = pmin(GRACE, time_cc),
    n_int = pmax(ceiling(tmax / DT), 0L) # make time into nice days, no fractions
  )

# long form
panel1 <- panel0 %>%
    uncount(weights = n_int, .id = "k") %>%
    mutate(
        t0 = (k - 1L) * DT,   # interval start
        t1 = k * DT,          # interval end
        t  = t1,              # simple time term
        t_sq = t^2             # squared time
    ) %>% 
    # create the daily censoring outcome
    mutate(
        cens = if_else(
            dev == 1 & !is.na(dev_time) & dev_time > t0 & dev_time <= t1,
        1, 0    
        )
    ) %>% 
    arrange(id, clone_arm, t1) %>%
    group_by(id, clone_arm) %>%
    mutate(cens_cum = cumsum(cens)) %>%
    ungroup() %>%
    filter(cens_cum <= 1) %>%
    select(-cens_cum)



# Prop Censor 
denom_formula <- as.formula(
  paste("cens ~", paste(covars, collapse = " + "), "+ t + t_sq") 
) # primary IPCW 
num_formula <- cens ~ t + t_sq # stablizing component

panel_control <- panel1 %>% filter(arm == "cntrl")
panel_treat   <- panel1 %>% filter(arm == "tx")

m_denom_control <- glm(denom_formula, data = panel_control, family = binomial())
m_num_control   <- glm(num_formula,   data = panel_control, family = binomial())

m_denom_treat <- glm(denom_formula, data = panel_treat, family = binomial())
m_num_treat   <- glm(num_formula,   data = panel_treat, family = binomial())



panel_control <- panel_control %>%
  mutate(
    p_cens_denom = predict(m_denom_control, type = "response"),
    p_cens_num   = predict(m_num_control,   type = "response"),
    p_unc_denom  = pmax(1 - p_cens_denom, EPS),
    p_unc_num    = pmax(1 - p_cens_num,   EPS),
    w_inc        = p_unc_num / p_unc_denom
  )


panel_treat <- panel_treat %>%
  mutate(
    p_cens_denom = predict(m_denom_treat, type = "response"),
    p_cens_num   = predict(m_num_treat,   type = "response"),
    p_unc_denom  = pmax(1 - p_cens_denom, EPS),
    p_unc_num    = pmax(1 - p_cens_num,   EPS),
    w_inc        = p_unc_num / p_unc_denom
  )



panel_control <- panel_control %>%
  arrange(id, clone_arm, t1) %>%
  group_by(id, clone_arm) %>%
  mutate(w_ipcw = cumprod(w_inc)) %>%
  ungroup()

panel_treat <- panel_treat %>%
  arrange(id, clone_arm, t1) %>%
  group_by(id, clone_arm) %>%
  mutate(w_ipcw = cumprod(w_inc)) %>%
  ungroup()
panel_w <- bind_rows(panel_control, panel_treat)

clone_weights <- panel_w %>%
  arrange(id, clone_arm, t1) %>%
  group_by(id, clone_arm) %>%
  summarise(w_ipcw_final = last(w_ipcw), .groups = "drop")

clones_w <- clones %>%
  left_join(clone_weights, by = c("id","clone_arm")) %>%
  mutate(
    # if a clone had zero intervals (rare), set weight = 1
    w_ipcw_final = replace_na(w_ipcw_final, 1)
  )

# trim weights
q <- quantile(clones_w$w_ipcw_final, probs = c(0.01, 0.99), na.rm = TRUE)


clones_w <- clones_w %>%
  mutate(
    w_ipcw_trim = pmin(pmax(w_ipcw_final, q[[1]]), q[[2]])
  )


 # FIt
fit_ccw <- survfit(
  Surv(time_cc, event_cc) ~ arm,
  data = clones_w,
  weights = w_ipcw_trim
)


# (15) 1-year (TAU) death risk and risk difference from weighted KM
s <- summary(fit_ccw, times = times)

km_tab <- tibble(
  t = rep(times,2),
  arm = s$strata,
  surv_tau = s$surv,
  risk_death_tau = 1 - s$surv
)

# Combine with BART Results and Observed CF estimates
km_long = km_tab %>% 
  rename(survival = surv_tau)%>% 
  mutate(Treatment = ifelse(arm == "arm=cntrl", "Control", "Tx")) %>% 
  mutate(Source = "CCW-KM") %>% 
  select(t, survival, Source, Treatment)

cond_survival_long = tmp %>% 
  rename(survival=sv_tx) %>% 
  mutate(Treatment = ifelse(tx_obs==0, "Control", "Tx")) %>%
  mutate(Source="Obs-Cond-Survival") %>% 
  select(t, survival, Source, Treatment)


res4_long %>%
  bind_rows(km_long) %>% 
  bind_rows(cond_survival_long) %>% 
  ggplot(aes(x = t, y = survival, color = Source, shape = Treatment)) +
  geom_point(size = 2.5, alpha = 0.8) +
  scale_color_discrete()+
  # scale_color_manual(values = c("Observed" = "red", "Predicted" = "blue")) +
  ylim(0, 1) +
  labs(
    title = "Survival: Observed vs Predicted",
    y = "Mean Survival",
    x = "Time (t)"
  ) +
  theme_minimal()



km_diff = km_tab %>% filter(t==30) %>% pull(surv_tau)
km_diff[2]-km_diff[1]
cntfct_sm %>% filter(t>25) %>% print(n=30)
bdiff %>% filter(t>25) %>% print(n=30)







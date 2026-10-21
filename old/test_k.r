library(tidyverse)




TAU    <- 365.24   # 1-year horizon (in days): the time we report risks at
GRACE  <- 182      # ~6-month grace period for the "surgery within 6 months" strategy
EPS    <- 1e-6     # tiny number to keep probabilities away from exactly 0 or 1
B      <- 500      # bootstrap replicates for confidence intervals
                   #   (bump this up to 2000+ for a final run; 500 keeps knitting fast)
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

# fix the categoricals to one-hots
# install.packages("fastDummies")
dat2 = fastDummies::dummy_cols(dat, cat_covars) %>% 
    select(-all_of(cat_covars)) %>% 
    select(-all_of(c("sex_1", "charlson_0", "perf_0","stage_1", "emergency_0"))) %>% 
    mutate(timetosurgery = ifelse(is.na(timetosurgery), Inf, timetosurgery)) %>% 
    glimpse

# extract updated names
covars2 = c()
for(n  in names(dat2)){
    if(str_replace(n, "_[0-9]", "") %in% covars){
        covars2 = append(covars2, n)
    }
}


glm(surgery ~ death + fup_obs dat2

source("sim/tmle_ccw_anno.r")
ccwtmle = run_ccw_tmle_analysis2(
        dat2,
        grace_period = GRACE,
        tau = floor(TAU),
        event = "death",
        event_time = "fup_obs",
        tx_obs = "surgery",
        tx_time_obs = "timetosurgery",
        covars=covars2,
        sl_lib = c("SL.mean", "SL.glm", "SL.earth", "SL.xgboost")) 

ccwtmle_res = ccwtmle$res_tmle 

ccwtmle_res %>% filter(row_number() > 360)

dat2 %>% head() %>% glimpse
unique(dat2$surgery)

pmin(c(20,30,40), 30)
c(20,30,40) %>% pmin(30)

source("sim/tmle_ccw_anno_scrpt.R")


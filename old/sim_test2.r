library(dplyr)
library(survival)

set.seed(42)
n_patients <- 500

# ---------------------------------------------------------
# 1. GENERATE THE TRUE DATA
# ---------------------------------------------------------
# Initialize population
df <- tibble(
  id = 1:n_patients,
  A1 = rbinom(n_patients, 1, 0.5), # 50% start Tx on Day 1
  Y1 = rbinom(n_patients, 1, ifelse(A1 == 1, 0.05, 0.10))
)

# Day 2 (Only for those who survived Day 1)
df <- df %>%
  mutate(
    A2 = ifelse(Y1 == 1, NA, ifelse(A1 == 1, 1, rbinom(n(), 1, 0.5))),
    Dur2 = ifelse(Y1 == 1, NA, ifelse(A1 == 1, 2, ifelse(A2 == 1, 1, 0))),
    Y2 = ifelse(Y1 == 1, NA, rbinom(n(), 1, case_when(
      Dur2 == 0 ~ 0.10,
      Dur2 == 1 ~ 0.05,
      Dur2 == 2 ~ 0.02
    )))
  )

# Day 3 (Only for those who survived Day 2)
df <- df %>%
  mutate(
    A3 = ifelse(Y1 == 1 | Y2 == 1, NA, ifelse(A2 == 1, 1, rbinom(n(), 1, 0.5))),
    Dur3 = ifelse(Y1 == 1 | Y2 == 1, NA, case_when(
      A1 == 1 ~ 3,
      A2 == 1 & A1 == 0 ~ 2,
      A3 == 1 & A2 == 0 ~ 1,
      TRUE ~ 0
    )),
    Y3 = ifelse(Y1 == 1 | Y2 == 1, NA, rbinom(n(), 1, case_when(
      Dur3 == 0 ~ 0.10,
      Dur3 == 1 ~ 0.05,
      Dur3 == 2 ~ 0.02,
      Dur3 == 3 ~ 0.01
    )))
  )

# Convert to Long Format (Person-Time)
df_long <- bind_rows(
  df %>% select(id, A=A1, Y=Y1) %>% mutate(time=1, dur=ifelse(A==1,1,0)),
  df %>% filter(Y1 == 0) %>% select(id, A=A2, Y=Y2, dur=Dur2) %>% mutate(time=2),
  df %>% filter(Y1 == 0 & Y2 == 0) %>% select(id, A=A3, Y=Y3, dur=Dur3) %>% mutate(time=3)
)

# ---------------------------------------------------------
# 2. CLONE-CENSOR-WEIGHT (CCW) ESTIMATE
# ---------------------------------------------------------
# Create the "Always Treat" Clones
ccw_always_treat <- df_long %>%
  group_by(id) %>%
  # Artificially censor the moment they deviate from A=1
  mutate(censor = ifelse(A == 0, 1, 0),
         cum_censor = cumsum(censor)) %>%
  filter(cum_censor == 0 | (cum_censor == 1 & censor == 1)) %>%
  # Only keep rows before they were censored
  filter(censor == 0)

ccw_always_treat %>% arrange(id, time)

# Calculate the discrete daily hazard directly from the long data
ccw_survival <- ccw_always_treat %>%
  group_by(time) %>%
  summarize(
    events = sum(Y == 1, na.rm = TRUE),
    at_risk = n(),
    hazard = events / at_risk,
    .groups = 'drop'
  ) %>%
  mutate(surv = cumprod(1 - hazard))
ccw_survival


# Calculate simple KM for the "Always Treat" clones
# (Weights omitted because without covariates, they cancel to 1)
ccw_survival <- summary(survfit(Surv(time-1, time, Y) ~ 1, data = ccw_always_treat, id=id))$surv
cat("CCW 'Always Treat' Survival (Day 1, 2, 3):\n", ccw_survival, "\n\n")

# ---------------------------------------------------------
# 3. NAIVE G-COMP (Memoryless Trap)
# ---------------------------------------------------------
# Fit model WITHOUT duration
fit_naive <- glm(Y ~ factor(time) + A, data = df_long, family = binomial)

# Simulate "Always Treat" Counterfactual (A=1 for all days)
cf_naive <- data.frame(time = 1:3, A = 1)
cf_naive$risk <- predict(fit_naive, newdata = cf_naive, type = "response")
cf_naive_survival <- cumprod(1 - cf_naive$risk)

cat("Naive G-Comp Survival (Misses the truth):\n", cf_naive_survival, "\n\n")

# ---------------------------------------------------------
# 4. CORRECTED G-COMP (With Duration)
# ---------------------------------------------------------
# Fit model WITH duration
fit_correct <- glm(Y ~ factor(time) + A + factor(dur), data = df_long, family = binomial)
predict(fit_correct, type="response")[1:10]


# Simulate "Always Treat" Counterfactual (A=1, dur increments)
cf_correct <- data.frame(time = 1:3, A = 1, dur = 1:3)
cf_correct$risk <- predict(fit_correct, newdata = cf_correct, type = "response")
cf_correct_survival <- cumprod(1 - cf_correct$risk)

cat("Corrected G-Comp Survival (Perfect Alignment):\n", cf_correct_survival, "\n\n")
cat("Oracle Truth for Day 3 is: 0.92169\n")
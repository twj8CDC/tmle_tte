suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(splines)
  library(glmnet)
  library(SuperLearner)
})

source("R/00_utils.R")
source("R/01_cloning.R")
source("R/02_ipcw_models.R")
source("R/03_q_models.R")
source("sim/dgp/02_dgp_misspecification.R")

cat("Generating Scenario 2 data (n = 200, seed = 42)...\n")
d2 <- generate_scenario2_data(n = 200, seed = 42)
d2_clean <- prepare_tte_data(d2)

# Build a small longitudinal panel to test SL.hal.time directly
clones <- create_clones(d2_clean, covars = c("w1", "w2"), grace_period = 10, tau = 30)
panel <- clones %>%
  dplyr::mutate(n_int = pmax(ceiling(time_cc), 0L)) %>%
  tidyr::uncount(weights = n_int, .id = "t") %>%
  dplyr::mutate(
    event_cc_day = dplyr::if_else(event_cc == 1L & time_cc > (t - 1L) & time_cc <= t, 1L, 0L),
    is_treated = dplyr::if_else(tx_ever == 1L & t >= t_switch, 1L, 0L)
  )

X_train <- panel %>% dplyr::select(t, w1, w2, is_treated)
Y_train <- panel$event_cc_day
X_test  <- X_train[1:50, ]
ids     <- panel$id

cat("Testing upgraded SL.hal.time directly from R/00_utils.R...\n")
res <- SL.hal.time(Y = Y_train, X = X_train, newX = X_test, family = stats::binomial(),
                   obsWeights = rep(1, length(Y_train)), id = ids)

cat("Direct pred length:", length(res$pred), "| Range:", range(res$pred), "\n")

# Verify predict method matches
pred_test <- predict.SL.hal.time(res$fit, newdata = X_test)
stopifnot(all.equal(res$pred, pred_test, tolerance = 1e-10) == TRUE)
cat("SUCCESS: Direct prediction and predict method are strictly identical!\n")

# Test in SuperLearner
register_sl_earth_int()

cat("Fitting SuperLearner with SL.mean and updated SL.hal.time...\n")
sl_fit <- SuperLearner(
  Y = Y_train,
  X = X_train,
  newX = X_test,
  family = binomial(),
  SL.library = c("SL.mean", "SL.hal.time"),
  cvControl = list(V = 2L),
  id = ids
)
print(sl_fit)
cat("SUCCESS: SuperLearner with updated SL.hal.time completed cleanly!\n")

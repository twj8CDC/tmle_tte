# ttetmle: Target Trial Emulation via Time-to-Event TMLE

[![R-CMD-check](https://img.shields.io/badge/R--package-0.1.0-blue.svg)](DESCRIPTION)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

`ttetmle` provides a fast, principled implementation of **Target Trial Emulation** combined with **Survival Targeted Minimum Loss-Based Estimation (TMLE)**. It is designed for observational time-to-event studies with dynamic treatment regimes and initiation grace periods, eliminating immortal time bias and confounding.

---

## Key Features

- **Clone-Censor-Weighting (CCW)**: Emulates target trials with dynamic grace periods while avoiding artificial censoring bias via time-varying IPCW.
- **Stratified Hazard Estimation & Decoupled Clever Covariates**: Efficiently updates discrete-time hazard models without pooling across incompatible regimes.
- **Machine Learning Integration**: Supports GLM and ensemble learning via `SuperLearner` (with built-in wrappers for splines, Earth/MARS, Ranger random forests, and BART).
- **Rigorous Asymptotic & Simultaneous Inference**:
  - Clustered patient-level Empirical Influence Curve (EIC) standard errors and pointwise 95% confidence intervals.
  - Multiplier bootstrap simultaneous confidence bands across follow-up time horizons.
  - Delta-method cumulative hazard ratios and pointwise/simultaneous inference.
- **Benchmarking Suite**: Out-of-the-box comparative models including Naive Cox (raw immortal time bias), Landmarking Cox, Time-Dependent Cox (`td-Cox`), Clone-Censor-Weighting Cox (`ccw-cphm`), CCW Kaplan-Meier, and BART G-computation.

---

## Installation

Install the development version directly from GitHub via `remotes` or `pak`:

```r
# If not already installed:
install.packages("remotes")

# Install ttetmle from GitHub:
remotes::install_github("<your-github-username>/tmle_tte")
```

Or using `pak`:

```r
# install.packages("pak")
pak::pkg_install("<your-github-username>/tmle_tte")
```

---

## Quick Start Example

Here is a minimal reproducible example demonstrating how to run an end-to-end TTE-TMLE analysis on observational survival data.

```r
library(ttetmle)
library(dplyr)
library(ggplot2)

# 1. Simulate a synthetic observational cohort
set.seed(42)
n <- 500
x1 <- rnorm(n, mean = 0, sd = 1)
p_tx <- plogis(-0.2 + 0.4 * x1)
tx_obs <- rbinom(n, 1, p_tx)

# Treatment initiation time during a 10-day grace period
tx_time_obs <- ifelse(
  tx_obs == 1,
  pmax(1, pmin(10, round(rbeta(n, 4, 2) * 10))),
  Inf
)

# Event time with treatment benefit
haz_base <- 0.03 * exp(0.3 * x1)
evt_time_raw <- rexp(n, rate = haz_base)
event_time <- pmax(1, pmin(30, ceiling(evt_time_raw)))
event <- ifelse(evt_time_raw <= 30, 1L, 0L)

df <- tibble(
  id = seq_len(n),
  x1 = x1,
  tx_obs = tx_obs,
  tx_time_obs = tx_time_obs,
  event_time = event_time,
  event = event
)

# 2. Fit TTE-TMLE
fit <- tte_tmle(
  data = df,
  grace_period = 10,       # Grace period duration for treatment strategy
  tau = 30,                # Evaluation horizon (e.g., 30 days)
  covars = c("x1"),        # Baseline confounders
  use_sl = FALSE,          # Set TRUE to use SuperLearner ensemble
  n_boot = 5000,           # Multiplier bootstrap iterations for simultaneous bands
  seed = 42
)

# 3. Print summary results
print(fit)

# 4. Access marginal survival, risk differences, and hazard ratios
head(fit$res_tmle)
# Columns include:
# - control_no_tx, treat_tx: Counterfactual marginal survival curves S(t)
# - rd: Risk Difference S_1(t) - S_0(t)
# - rd_ci_lower, rd_ci_upper: Pointwise 95% CIs
# - rd_band_lower, rd_band_upper: Multiplier bootstrap simultaneous 95% bands
# - hr_cum, hr_ci_lower, hr_ci_upper: Cumulative hazard ratios and CIs

# 5. Visualize
# Survival curves with pointwise CIs
plot(fit, which = "survival")

# Marginal Risk Difference with simultaneous confidence bands
plot(fit, which = "rd")

# Cumulative Hazard Ratio trajectory
plot(fit, which = "hr")
```

---

## Data Requirements

Input datasets can be supplied as a standard `data.frame` or `tibble`. The required variables are:

| Parameter in `tte_tmle` | Default Column Name | Description |
|:------------------------|:--------------------|:------------|
| `id_col` | `"id"` | Unique patient/individual identifier |
| `time_col` | `"event_time"` | Observed time to event or censoring |
| `event_col` | `"event"` | Binary event indicator (`1` = event, `0` = censored) |
| `tx_col` | `"tx_obs"` | Binary treatment indicator (`1` = ever received treatment, `0` = never) |
| `tx_time_col` | `"tx_time_obs"` | Time of treatment initiation (use `Inf` for never treated) |
| `covars` | `c("x1")` | Character vector of baseline covariate column names |

Custom column names can be passed directly to `tte_tmle()`:
```r
fit <- tte_tmle(
  data = my_data,
  id_col = "patient_id",
  time_col = "followup_days",
  event_col = "death",
  tx_col = "drug_initiated",
  tx_time_col = "day_of_initiation",
  covars = c("age", "biomarker", "sex")
)
```

---

## Comparing with Standard Benchmarks

The package provides direct wrappers to evaluate competing estimators side-by-side:

```r
# Naive Cox (subject to immortal time bias)
fit_naive <- fit_naive_cox(df, covars = c("x1"))

# Landmarking Cox (evaluated at grace period cutoff)
fit_lm <- fit_landmark_cox(df, landmark_time = 10, covars = c("x1"), tau = 30)

# Time-Dependent Cox (td-Cox)
fit_td <- fit_td_cox(df, covars = c("x1"), tau = 30, grace_period = 10)

# CCW Kaplan-Meier & CCW Cox Proportional Hazards
clones <- create_clones(df, grace_period = 10)
panel_w <- compute_ipcw_weights(clones, grace_period = 10, tau = 30, covars = c("x1"))
fit_km <- fit_ccw_km(panel_w, tau = 30)
fit_cphm <- fit_ccw_cphm(panel_w, covars = c("x1"), tau = 30)
```

---

## Package Architecture

```
R/
├── 00_utils.R         # Helper utilities, discrete grids, custom SuperLearner wrappers
├── 01_cloning.R       # Clone-Censor expansion for dynamic regimes
├── 02_ipcw_models.R   # Time-varying IPCW estimation and weight stabilization
├── 03_q_models.R      # Initial outcome conditional hazard models Q(t, A, W)
├── 04_targeting.R     # Clever covariate fluctuation and submodel targeting
├── 05_gcomp.R         # Marginal survival curve integration via G-computation
├── 06_inference.R     # Clustered EIC influence curves, Delta-method HR, Multiplier bands
├── 07_tte_tmle.R      # High-level user interface tte_tmle(), print(), and autoplot()
├── 08_benchmarks.R    # Naive Cox, Landmarking, td-Cox, CCW-KM, CCW-CPHM
└── 09_bart_gcomp.R    # Non-parametric Bayesian Additive Regression Trees (BART)
```

---

## License

This package is licensed under the [MIT License](LICENSE).


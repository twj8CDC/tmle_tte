# TTE-TMLE Roadmap, Bug Fixes & Package Architecture TODO

**Created:** 2026-11-11  
**Scope:** Two-phase roadmap covering (1) immediate IPCW bug fixes & benchmark estimand alignment, and (2) general R package refactoring to support 4 explicit treatment-arm execution modes.

---

## Phase 1: Immediate Bug Fixes & Estimand Alignment *(Completed 2026-11-11)*

### 1.1 Fix Missing $t \times W$ Interaction & Pooled Trimming in `R/02_ipcw_models.R`
- [x] **Fix `m_denom_ctrl` model specification (`R/02_ipcw_models.R`):**
  - Updated `m_denom_ctrl` to `splines::ns(t, df = 3) * (<p_form>)` (outperformed `df = 4` and linear `t:` under positivity strain) and `m_num_ctrl` to `splines::ns(t, df = 3)` with `compute_grace_survival_grid()` so splines are never extrapolated past `grace_period`.
- [x] **Fix arm-pooled weight trimming & tail separation (`R/02_ipcw_models.R`, `R/05_gcomp.R`):**
  - Implemented arm-stratified trimming (`compute_arm_trim_bounds()`, `apply_arm_weight_trimming()`) preserving structural `0`s, plus explicit `max_weight = 10` (`2 * max_weight = 20` unstabilized) ceiling to prevent spline tail separation (`eps_ctrl`).
- [x] **Verify `R/09_bart_gcomp.R` numerator consistency *(Completed with Phase 2)*:**
  - Aligned `m_num_ctrl` to `splines::ns(t, df = min(3, n_grace_int - 1))` and switched BART IPCW trimming to arm-stratified quantile + `max_weight = 10` capping via shared `R/02_ipcw_models.R` helpers alongside Phase 2 BART splicing.

### 1.2 Align `CCW-KM` & `CCW-CPHM` Grace-Period Weighting to Complier Timing Target (`S1_beta`)
- [x] **Resolve Days 1–10 estimand mismatch (`R/02_ipcw_models.R`, `R/05_gcomp.R`, `R/07_tte_tmle.R`, `R/08_benchmarks.R`):**
  - Added `align_grace_weights = TRUE` (default) applying conditional stochastic grace-period weights $H^{(1)}_{\text{cond}}(t)$ on `treat_tx` during Days $1 \dots m$, synced counterfactual weights in `marginalize_gcomp()`, and filtered structural zero weights before `survival::coxph()` in `fit_ccw_cphm()`.
- [x] **Re-run parametric benchmarks across all 800 datasets (`Scenarios 1–8` $\times$ `100` iters):**
  - Re-evaluated all 12,800 parametric rows (`fast_only = TRUE` on exact benchmark seeds), merged with preserved 2,400 Non-Parametric rows into `sim/simulation_results_100iter.rds` & `.csv`, and regenerated `sim/simulation_results_100iter.xlsx`, all 8 `.png` plots, and `sim/benchmark_100iter_visual_report.md`.

---

## Phase 2: Package Refactor — Treatment-Arm Execution Modes for `tte_tmle()` *(Completed 2026-11-12)*

Design goal: Make treatment-arm initiation policies explicit options inside [`tte_tmle()`](file:///Users/jacobtiegs/projects/tmle_tte/R/07_tte_tmle.R) via `tx_strategy = c("conditional_complier", "empirical_complier", "natural_grace")` (with `"conditional_complier"` as default; defer `"custom"` for a later pass) and strict **$Q$-counterfactual and $g$-clever-covariate pairing**.

### 2.1 Mode Definitions & All-Observation Assignment Rules

| Mode | Policy Name | Counterfactual Initiation Rule (All Patients from Baseline $W_i$) | Matching `g` / Grace-Period Clever Covariate | Notes |
| :--- | :--- | :--- | :--- | :--- |
| **1 (Default)** | **Conditional Complier** (`"conditional_complier"`) | Assign **all** observations across $k \in \{1 \dots m\}$ weighted by $w_k(W_i) = \frac{\hat{P}_0(T_A = k \mid W_i)}{1 - \hat{S}_{0,A}(m \mid W_i)}$ derived from `m_denom_ctrl`. | **Conditional Dynamic Weights:** <br>Pre-switch: $\frac{S_{0,A}(t \mid W) - S_{0,A}(m \mid W)}{S_{0,A}(t \mid W)(1 - S_{0,A}(m \mid W))}$ <br>Post-switch: $\frac{1}{1 - S_{0,A}(m \mid W)}$ | Derives $w_k(W_i)$ directly from `m_denom_ctrl` on `obs_panel` (properly accounting for competing mortality and avoiding positivity violations when $T_A$ depends on $W$). |
| **2** | **Marginal Empirical Complier** (`"empirical_complier"`) | Assign **all** observations (both observed compliers and non-compliers) across $k \in \{1 \dots m\}$ weighted by empirical complier probabilities $p^*(k)$. | **Marginal Dynamic Weights:** <br>Pre-switch: $\frac{\sum_{j=t+1}^m p^*(j)}{S_{0,A}(t \mid W)}$ <br>Post-switch: $\frac{p^*(T_A)}{P_0(T_A \mid W)}$ | Applying $p^*(k)$ to *all* obs avoids post-baseline survivorship filtering ($T \ge t_{\text{switch}}$) and hybrid mixtures when $T_A$ depends on $W$. |
| **3** | **Hernán Natural Grace Period** (`"natural_grace"`) | Keep natural initiation hazard $\hat{P}_0(T_A = k \mid W_i)$ for $k \in \{1 \dots m\}$, and assign remaining probability $\hat{S}_{0,A}(m \mid W_i)$ to **`grace_period + 1`** (or for MC/observed: `true_tx = true_tx_time`, `true_untx = grace_period + 1`). | **Standard Hernán CCW Weights:** <br>$w = 1$ for $t \le m$, and $\frac{1}{P(T_A \le m \mid T > m, W)}$ for $t > m$. | Matches `R/01_cloning.R` censoring at `dev_time = grace_period` (`is_treated = 0` on Day `m`, `1` starting Day `m+1`). Note prevalence dependence in docs. |
| **4 (Later)** | **Custom Strategy** (`"custom"`) *(Deferred)* | User-supplied deterministic switch day (e.g., `1` for immediate baseline treatment), custom vector of switch times, or discrete probability vector / function `p_star_fn(data, m)`. | **General Radon–Nikodym Ratio:** <br>Cumulative $\prod_{k=1}^t \frac{g^*_k}{g_{0,k}}$ matching user policy $g^*$. | Enables protocols such as "treat within 3 days", uniform initiation, or immediate $t=1$ initiation. |

### 2.2 Counterfactual Integration Engineering (Exact Analytical Expectation)
- [x] **Replace single `sample()` draw with analytical integration in `R/05_gcomp.R`:**
  - Compute exact expected individual survival curve:
    $$\hat{S}_i^{(1)}(t) = \sum_{k=1}^{m_{\max}} w_k(W_i)\,\hat{S}^*(t \mid T_A = k, W_i)$$
  - Eliminates Monte Carlo seed jitter and avoids Jensen's inequality error (do **not** plug mean/RMST switch day $\bar{k}(W_i)$ into a single survival curve, since $S(t \mid T_A=k, W)$ is non-linear in $k$).
  - For **all-observation assignment** in Modes 1 & 2, apply weights $w_k$ uniformly to every patient $i = 1 \dots N$ based purely on baseline $W_i$, removing post-baseline survivorship bias (`tx_ever == 1`).
- [x] **Zero-overhead BART analytical integration (`R/09_bart_gcomp.R` & `inst/python/run_bart_survival.py`):**
  - Because `SurvBART` currently takes only `(t, covars, is_treated)` (no `time_since_treatment` feature), predicting `is_treated = 0` and `is_treated = 1` across the two $N \times \tau$ grids ($2 \times N \times \tau$ total rows—identical to current memory footprint!) supplies every hazard needed to construct **all** $k = 1 \dots m$ survival curves via array slicing (`q0` for $t < k$, `q1` for $t \ge k$).
  - If `time_since_treatment` is added to `Q` later, implement chunked/in-Python aggregation over $k = 1 \dots m$ inside `run_bart_survival.py` so only the marginalized $(M \times N \times \tau)$ survival array is returned to R.
- [x] **Align EIC inference (`R/06_inference.R`):**
  - Verify individual influence curves use the analytically marginalized $\hat{S}_i^{(1)}(t)$ and corresponding $H^{(1)}(t)$ clever covariates for valid variance and multiplier bootstrap bands.

---

## Phase 3: Secondary Thesis & Benchmark Improvements (Backlog)
- [ ] **Upgrade All BART Defaults & Benchmarks to `SurvBART.preset_standard()` (`R/09_bart_gcomp.R`, `sim/run_100iter_benchmark.R`, `sim/run_simulation_condition.R`):**
  - Switch `fit_bart_pair()`, `fit_bart_gcomp()`, `fit_bart_tmle()`, and simulation runners from fast (`trees = 30`, `ndpost = 100`, `nskip = 100`, `thin = 2`, `chains = 1`) to `SurvBART.preset_standard()` (`trees = 50`, `draws = 1000` via `ndpost = 250, thin = 4`, `tune = 250`, `chains = 4`, yielding 1,000 thinned posterior draws). Verified on Scenarios 1–3 (`sim/scratch/test_gcomp_bart_standard_s1_s2.R`, `sim/scratch/test_bart_preset_standard_s3_20iter.R`): fixes fast-mode interval/band under-coverage and cuts Scenario 3 Day-15 bias (`-0.0286 -> -0.0071` for `TTE-TMLE (BART)`, `90%` uniform band coverage).
- [ ] **Upgrade `SL.hal.time()` & SuperLearner CV (`R/00_utils.R`, `R/03_q_models.R`):**
  - Add nonlinear covariate basis expansion (`splines::ns(w, df = 3)` + pairwise interactions) alongside the 1-day cumulative time steps `t_ge_k` / `trt_t_ge_k` in `SL.hal.time()` (which currently passes linear-only covariates), test `alpha = 1` (pure LASSO without ridge smearing across adjacent step days), and move `SuperLearner` from `cvControl = list(V = 2)` unclustered to `V = 5` patient-clustered CV (`id`).
- [ ] **New Benchmark Scenario (Covariate-Dependent Treatment Timing):** Create a DGP & oracle truth generator where initiation day $T_A \in \{1 \dots m\}$ among compliers explicitly depends on baseline covariates $W$ (e.g., high-risk patients start earlier around Days 1–3, low-risk patients later around Days 7–10), separating Mode 1 (`"conditional_complier"`), Mode 2 (`"empirical_complier"`), and Mode 3 (`"natural_grace"`), and testing day-level positivity behavior.
- [ ] **Scenario 3 Non-Proportional Hazards (`R/03_q_models.R`):** Add time-by-treatment interaction (`is_treated * splines::ns(t, df=3)` or `is_treated * t`) in parametric Q-models (or document why constant-HR `is_treated` causes Day-15 bias and ~50% uniform band coverage under sharp Day-15 efficacy waning).
- [ ] **Scenario 6 Boundary Documentation:** Document explicitly in paper/vignette that baseline TTE-TMLE requires Longitudinal TMLE (`ltmle`) when post-baseline time-varying confounders suffer treatment-confounder feedback.

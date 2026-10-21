# TTE-TMLE Project Master Plan & Validation Roadmap

**Project Title:** TTE-TMLE: Target Trial Emulation via Time-to-Event Targeted Minimum Loss-Based Estimation for Immortal Time Bias  
**R Package Name:** `ttetmle` (v0.1.0)  
**Author / Lead:** Jacob Tiegs  
**Last Updated:** October 2026  
**Status:** Package Skeleton & Core Modules Implemented; Advancing to Simulation Validation Engine  

---

## 1. Executive Summary & Scientific Motivation

### 1.1 The Clinical & Methodological Problem
In observational comparative effectiveness research, dynamic treatment strategies with delayed initiation (e.g., initiating surgery or medication within a defined grace period $t_g$) frequently suffer from **Immortal Time Bias (ITB)**. Naive observational analyses improperly classify post-baseline person-time, crediting treatments with survival earned during the untreated waiting period.

The **Clone-Censor-Weight (CCW)** framework (Hernán et al.) elegantly resolves ITB by:
1. Cloning patients at baseline into treatment ($A=1$) and control ($A=0$) protocol arms.
2. Artificially censoring clones when they violate their protocol rules.
3. Applying Inverse Probability of Censoring Weights (IPCW) to reweight compliant person-time.

### 1.2 The Gap in Standard Practice
Standard target trial emulation implementations typically evaluate the reweighted cloned data using:
- **Weighted Kaplan-Meier (KM)** estimators, or
- **Weighted pooled logistic / Cox proportional marginal hazards (CPMH)** models.

These standard approaches suffer from severe methodological vulnerabilities:
1. **Single Robustness**: If the censoring/compliance model is misspecified, the entire survival curve is biased.
2. **Weight Instability & "Cliffs"**: In the presence of strong confounding or near-positivity violations, extreme weights induce sharp, jagged step-function drops ("cliffs") in weighted KM curves and balloon empirical variance.
3. **Inefficiency (Ignored Prognostic Covariates)**: Pure IPW estimators fail to leverage prognostic covariates that strongly predict survival but do not drive treatment assignment, wasting statistical precision.
4. **Ad-Hoc Clustered Variance**: Because cloning duplicates individuals ($2N$), standard variance formulas fail. Users must rely on clustered "sandwich" estimators (e.g., GEE), which exhibit poor coverage and instability in small or moderate finite samples.

### 1.3 The Solution: TTE-TMLE
**TTE-TMLE** integrates the structural trial emulation of CCW with semiparametric **Survival TMLE**:
- **Double Robustness**: Estimates both an outcome hazard model ($Q$) and a compliance/censoring model ($g$) via machine learning ensembling (`SuperLearner`). If either model is correctly specified, point estimates remain consistent and unbiased.
- **Local Efficiency & Precision Bonus**: Achieves the Semiparametric Efficiency Bound (Cramér-Rao lower bound for semiparametric models) when both models are well-specified. The $Q$-model absorbs prognostic variation, substantially tightening confidence intervals compared to weighted KM.
- **Biologically Smooth Trajectories**: Anchors survival estimates to an initial non-parametric hazard model, eliminating erratic weight-induced cliffs.
- **Native Individual-Level EIC Inference**: Groups the Empirical Influence Curve (EIC) at the original patient level ($N$ independent units), mathematically resolving clone-induced correlation in closed form without ad-hoc sandwich corrections.
- **Simultaneous Confidence Bands via Multiplier Bootstrap**: Simulates the Gaussian limiting process directly from individual EIC trajectories, providing simultaneous coverage across the entire survival curve $t \in [1, \tau]$ in seconds without costly resample re-fitting.

---

## 2. Rebuilt Workspace Architecture (`ttetmle` R Package)

The project is structured as a standard, installable R package (`ttetmle`) alongside dedicated simulation, applied, and documentation directories:

```
tmle_tte/
├── DESCRIPTION                      # CRAN package metadata and dependencies
├── NAMESPACE                        # Package exports and imports (roxygen-managed)
├── LICENSE                          # MIT License specification
├── .Rbuildignore                    # Ignores sim/, applied/, notes/, old/ for package building
├── PROJECT_PLAN.md                  # Master roadmap and validation architecture
├── notes/                           # Conceptual notes and raw transcripts
│   └── notes.md
├── old/                             # Preserved legacy scripts and prototypes (read-only reference)
│   ├── tmle_ccw_anno.r
│   ├── tmle_ccw_anno_scrpt.R
│   ├── sim_1_baseline.r
│   └── ...
├── R/                               # Production modular package codebase (ttetmle)
│   ├── 00_utils.R                   # Data preparation, transformations, validation checks
│   ├── 01_cloning.R                 # Protocol definition, clone duplication, artificial censoring (C*)
│   ├── 02_ipcw_models.R             # Pooled logistic & SuperLearner compliance models
│   ├── 03_q_models.R                # Stratified initial hazard estimation (GLM / SuperLearner)
│   ├── 04_targeting.R               # Decoupled clever covariates & logistic fluctuation
│   ├── 05_gcomp.R                   # Counterfactual complier sampling & pure G-comp / marginalization
│   ├── 06_inference.R               # Patient-level EIC, closed-form SEs, Multiplier Bootstrap bands
│   ├── 07_tte_tmle.R                # Main orchestrator function: tte_tmle(), S3 methods (print, plot)
│   └── 08_benchmarks.R              # Benchmark estimators (Naive Cox, Landmarking, tdCox)
├── man/                             # Auto-generated Rd documentation files
├── sim/                             # Simulation study framework (ADEMP compliant)
│   ├── dgp/                         # Data-Generating Processes for Scenarios 1–4
│   ├── true_values.R                # Counterfactual Monte Carlo truth generator (N = 100,000)
│   ├── run_sim_scenario.R           # Scenario runner with parallelization (future/furrr) & L'Ecuyer RNG
│   ├── evaluate_metrics.R           # ADEMP performance metrics (Bias, ESE, ASE, Coverage, Power)
│   └── configs/                     # Scenario parameter specifications (.R / .json)
├── applied/                         # Applied case study (Target Trial Emulation)
│   ├── data/                        # Real-world or benchmark observational dataset
│   └── run_applied_analysis.R       # End-to-end applied workflow, tables, and figures
├── tests/                           # Formal unit and regression tests
│   └── test_tte_tmle.R              # Modular pipeline and benchmark verification test
└── docs/                            # Manuscript drafts and theoretical documentation
    ├── mathematical_formulation.md  # Formal SCM, EIC derivations, and notation
```

---

## 3. Comprehensive Validation Schema (ADEMP Framework)

The validation suite follows the **ADEMP** reporting guidelines (Morris et al., 2019: *Aims, Data-generating mechanisms, Estimands, Methods, Performance measures*).

### 3.1 Candidate Estimators to Benchmark
To demonstrate clear empirical superiority and delineate when each method breaks down, we benchmark eight estimators:

| # | Estimator Name | Model Specification | Role / Scientific Rationale |
|---|----------------|---------------------|-----------------------------|
| 1 | **Naive Cox / KM** | Unadjusted univariate $T \sim A_{\text{obs}}$ | Quantifies raw immortal time bias without adjustment. |
| 2 | **Landmarking Cox** | Cox PH conditioned on survival to $t = t_g$ | Common clinical heuristic; prone to selection/collider bias and immortal person-time prior to landmark. |
| 3 | **Time-Dependent Cox (tdCox)** | Counting process `Surv(tstart, tstop, evt) ~ tx_td + W` | Standard epidemiological benchmark; measures instantaneous hazard rather than a pragmatic strategy. |
| 4 | **Parametric CCW (Hernán Standard)** | Cloned data + pooled logistic IPCW + weighted KM | Tests standard target trial emulation against functional misspecification (singly robust). |
| 5 | **Point-Treatment IPCW-TMLE** | TMLE at horizon $\tau$ using baseline IPCW | Tests the common shortcut of applying cross-sectional TMLE to right-censored data without longitudinal hazard targeting. |
| 6 | **Pure G-Computation** | Stratified hazard regression + counterfactual marginalization (no targeting) | Tests substitution estimator behavior without TMLE bias fluctuation. |
| 7 | **TTE-TMLE (Parametric GLM)** | Cloned data + GLM $g$-model + GLM $Q$-model + Decoupled Targeting + EIC | Demonstrates classical double robustness under correct vs. misspecified parametric models. |
| 8 | **TTE-TMLE (SuperLearner)** | Cloned data + SuperLearner ensemble + Decoupled Targeting + EIC + Multiplier Bootstrap | **The Proposed Estimator**: Non-parametric double robustness, local efficiency, and simultaneous bands. |

---

### 3.2 The 4-Tier Simulation Scenarios

#### Scenario 1: Structural Baseline (Pure ITB + Grace Period)
* **Goal**: Validate that all valid target trial methods eliminate ITB and recover the true dynamic treatment effect under ideal, correctly specified linear conditions.
* **DGP**:
  - Baseline confounder: $W \sim \mathcal{N}(0, 1)$.
  - Treatment initiation propensity: $\text{logit}\,P(A=1 \mid W) = \alpha_0 + \alpha_1 W$.
  - Initiation day (among initiators): $T_{\text{switch}} \sim \text{Beta}(5, 2) \times 10$ (grace period $t_g = 10$).
  - Underlying survival: Weibull with constant baseline shape and constant proportional treatment benefit ($\text{HR} \approx 0.35$).
* **Hypothesis**:
  - Naive Cox and Landmarking exhibit substantial bias ($> 20\%$).
  - Parametric CCW, tdCox, Pure G-comp, and TTE-TMLE all achieve minimal bias ($\le 1\%$).
  - TTE-TMLE matches the parametric efficiency of CCW.

#### Scenario 2: Functional Misspecification & Non-Linear Confounding (Double Robustness Demonstration)
* **Goal**: Prove that parametric IPCW-CCW fails when confounding relationships are non-linear or feature sharp biomarker thresholds, while TTE-TMLE with SuperLearner preserves unbiasedness and nominal coverage.
* **DGP**:
  - Continuous frailty $W_1 \sim \mathcal{N}(3, 2)$ and binary biomarker $W_2 \sim \text{Bernoulli}(0.4)$.
  - Treatment assignment mechanism: Quadratic frailty and threshold interaction:
    $$\text{logit}\,P(A=1 \mid W_1, W_2) = \gamma_0 + \gamma_1 W_1 + \gamma_2 W_1^2 + \gamma_3 \mathbf{1}(W_1 > 4 \land W_2 = 1)$$
  - Outcome hazard: Non-linear hazard surface with treatment interaction:
    $$\lambda(t \mid A, W_1, W_2) = \lambda_0(t) \exp\left(\beta_1 W_1 + \beta_2 \sin(W_1) + \beta_3 W_2 + \beta_4 A + \beta_5 (A \times W_1)\right)$$
* **Hypothesis**:
  - Parametric CCW incurs large bias ($> 12\%$) due to misspecified logistic weights.
  - TTE-TMLE (GLM) shows partial bias reduction but struggles with functional misspecification.
  - TTE-TMLE (SuperLearner: GLM, MARS/Earth, Ranger) achieves $< 2\%$ absolute bias, recovering nominal $95\%$ coverage.

#### Scenario 3: Non-Proportional Hazards, Dynamic Regimes & Waning Efficacy
* **Goal**: Demonstrate estimator performance when the treatment hazard ratio changes abruptly over time (structural break), and evaluate the validity of simultaneous confidence bands via the Multiplier Bootstrap.
* **DGP**:
  - Early phase ($t \le 15$): Strong protective effect ($\text{HR} = 0.20$).
  - Late phase ($t > 15$): Immediate waning/washout ($\text{HR} = 1.00$).
* **Hypothesis**:
  - Time-Dependent Cox averages across the timeline, reporting a misleading constant hazard ratio and failing to capture dynamic efficacy.
  - TTE-TMLE with decoupled targeting covariates ($H_{\text{treat,pre}}$ vs. $H_{\text{treat,post}}$) tracks the curve inflection accurately across all horizons.
  - Point-wise CIs suffer from multiple testing across the 30-day timeline (simultaneous coverage drops to $\approx 70\%$).
  - Multiplier Bootstrap simultaneous confidence bands achieve the nominal $95\%$ simultaneous coverage over the entire trajectory $t \in [1, 30]$.

#### Scenario 4: Prognostic Precision Bonus & Local Efficiency
* **Goal**: Prove that TTE-TMLE provides substantially narrower confidence intervals and higher empirical power than standard weighted KM, even when parametric weights are perfectly specified, by absorbing variance from purely prognostic covariates.
* **DGP**:
  - Treatment assignment depends linearly on confounder $W_{\text{conf}}$.
  - True treatment benefit is moderate ($\text{Risk Difference at } \tau = 30 \approx 5\%$).
  - Outcome is heavily driven by an independent prognostic score $W_{\text{prog}} \sim \mathcal{N}(0, 1.5)$ with large hazard coefficient ($\beta_{\text{prog}} = 1.8$).
* **Hypothesis**:
  - Both Parametric CCW and TTE-TMLE are unbiased.
  - Parametric CCW's weighted KM ignores $W_{\text{prog}}$, resulting in high empirical variance.
  - TTE-TMLE's $Q$-model absorbs $W_{\text{prog}}$, reducing standard errors by $30\%–50\%$ and increasing statistical power from $\approx 60\%$ to $> 90\%$ in finite samples ($N = 400$).

---

### 3.3 Simulation Protocol & Performance Measures

Following Morris et al. (2019):
- **Monte Carlo Replications**: $B = 1,000$ iterations per scenario (with $B = 200$ for rapid testing/debugging).
- **Sample Sizes**: $N \in \{300, 1000\}$ to assess finite-sample stability versus asymptotic properties.
- **RNG Seeding**: L'Ecuyer-CMRG parallel stream management via `future` / `furrr` to ensure strict reproducibility.
- **Target Estimands**:
  1. Primary: Marginal Risk Difference at day $\tau = 30$:
     $$\Psi_{\text{RD}}(\tau) = S_1(\tau) - S_0(\tau)$$
  2. Intermediate: Marginal Risk Difference at day $\tau = 10$ (end of grace period) and day $\tau = 20$.
  3. Full Trajectory: Simultaneous survival curve difference trajectory $\{S_1(t) - S_0(t)\}_{t=1}^\tau$.
- **Performance Measures**:
  1. **Mean Bias**: $\frac{1}{B} \sum_{b=1}^B (\hat{\Psi}_b - \Psi_{\text{true}})$
  2. **Relative % Bias**: $100 \times \frac{\text{Mean Bias}}{\Psi_{\text{true}}}$
  3. **Empirical Standard Error (ESE)**: Sample standard deviation of $\hat{\Psi}$ across $B$ runs.
  4. **Average Model Standard Error (ASE)**: $\frac{1}{B} \sum_{b=1}^B \widehat{\text{SE}}(\hat{\Psi}_b)$ from the patient-level EIC.
  5. **Variance Ratio (ASE / ESE)**: Target is $\approx 1.00$. Ratios $< 0.95$ indicate anti-conservative inference (Type I error inflation).
  6. **Nominal 95% Coverage Probability**: Proportion of $B$ runs where $\Psi_{\text{true}} \in [\hat{\Psi}_b \pm 1.96 \times \widehat{\text{SE}}_b]$.
  7. **Simultaneous Band Coverage**: Proportion of $B$ runs where the entire true trajectory $S_{\text{true}}(t)$ remains strictly inside the multiplier bootstrap band for all $t \in [1, \tau]$.
  8. **Pointwise 95% Confidence Interval Width**: $\frac{1}{B} \sum_{b=1}^B (\widehat{\text{CI}}_{\text{upper}}(t) - \widehat{\text{CI}}_{\text{lower}}(t))$ to measure precision gains.
  9. **Simultaneous 95% Confidence Band Width**: $\frac{1}{B} \sum_{b=1}^B (\widehat{\text{Band}}_{\text{upper}}(t) - \widehat{\text{Band}}_{\text{lower}}(t))$ across time, quantifying the price of curve-wide protection.
  10. **Empirical Power / Type I Error**: Rejection rate of $H_0: \Psi_{\text{RD}}(\tau) = 0$ at $\alpha = 0.05$.
  11. **Computation Time**: Wall-clock seconds per iteration.

---

## 4. Key Algorithmic Extensions to Implement

### 4.1 Extension 1: The Multiplier Bootstrap for Simultaneous Bands
- **Theory**: The limiting distribution of the empirical process $\sqrt{N}(\hat{S}(t) - S_0(t))$ converges to a zero-mean Gaussian process whose covariance is dictated by the Influence Curve.
- **Implementation**:
  1. For each patient $i \in \{1, \dots, N\}$, extract their trajectory of influence scores $\mathbf{IC}_i = \{\text{IC}_i(t)\}_{t=1}^\tau$.
  2. Draw $M = 10,000$ vectors of i.i.d. random multipliers $Z_i^{(m)} \sim \mathcal{N}(0, 1)$ or Rademacher ($\pm 1$).
  3. Compute simulated trajectories:
     $$G^{(m)}(t) = \frac{1}{\sqrt{N}} \sum_{i=1}^N Z_i^{(m)} \frac{\text{IC}_i(t)}{\widehat{\text{SE}}(t)}$$
  4. Find the supremum over $t \in [1, \tau]$ for each realization: $c^{(m)} = \sup_{t} |G^{(m)}(t)|$.
  5. Compute the $(1-\alpha)$ quantile $c_{1-\alpha}$ of $\{c^{(m)}\}$.
  6. Simultaneous confidence band: $\hat{\Psi}(t) \pm c_{1-\alpha} \times \widehat{\text{SE}}(t)$.
- **Benefit**: Generates rigorous, non-parametric simultaneous bands that account for serial correlation over time and clone correlation across arms in milliseconds, avoiding repeated SuperLearner refits.

### 4.2 Extension 2: Pure G-Computation with Ghost Population (BART / Python BART-Survival)
- **Theory**: Evaluates the substitution estimand directly without the TMLE targeting fluctuation.
- **Model Choice**: Implemented using Bayesian Additive Regression Trees (BART) via R `BART` or Python `bart-survival` for flexible non-parametric hazard estimation.
- **Protocol & Compliance Mechanics**:
  - Non-compliers in the dynamic target strategy $A=1$ do not initiate treatment until *after* the grace period ends ($t > t_g$) or never initiate.
  - Counterfactual trajectory modeling enforces this strict protocol timing rather than sampling complier switch times, capturing the true dynamic regime benchmark.
- **Benefit**: Isolates the specific bias-reduction and local efficiency contribution of the TMLE targeting step compared to an advanced non-parametric substitution estimator.

### 4.3 Extension 3: Cross-Validated TMLE (CV-TMLE)
- **Theory**: When aggressive machine learning algorithms (e.g., Random Forests, Gradient Boosting) are used in SuperLearner, in-sample predictions can overfit, producing near-zero or near-one probabilities that destabilize targeting.
- **Implementation**:
  - Employs 5-fold cross-fitting ($V = 5$).
  - Evaluated either:
    - (a) Globally across all TMLE runs to ensure honest nuisance parameter estimation, or
    - (b) Selectively targeted for high-complexity, low-event-count / rare-event scenarios where empirical separation threatens asymptotic coverage.
  - Trains $Q$ and $g$ models on $V-1$ folds and evaluates out-of-fold predictions prior to fluctuation.
- **Benefit**: Ensures asymptotic normality and protects nominal coverage even with highly complex, flexible learners.

### 4.4 Extension 4: Functional Delta Method for Cumulative Hazard Ratios
- **Theory**: Under proportional hazards, the ratio of cumulative hazards $\theta(t) = \frac{\Lambda_1(t)}{\Lambda_0(t)} = \frac{\log(S_1(t))}{\log(S_0(t))}$ equals the Cox model parameter $\exp(\beta)$. By the Functional Delta Method on the empirical influence curve:
  $$\text{IC}_{\log \theta, i}(t) = \frac{\text{IC}_{1, i}(t)}{S_1(t) \log(S_1(t))} - \frac{\text{IC}_{0, i}(t)}{S_0(t) \log(S_0(t))}$$
- **Implementation**: Fully integrated into `compute_eic_inference()` and `tte_tmle()`, outputting `hr_cum`, `se_log_hr`, `hr_ci_lower/upper`, `hr_band_lower/upper`, and instantaneous daily discrete hazard ratio `hr_daily`.
- **Benefit**: Provides clinicians and epidemiologists with their familiar relative hazard metric while preserving double robustness, non-parametric validity, and closed-form simultaneous confidence bands.

---

## 5. Execution Roadmap & Phase Breakdown

```
[Phase 1: Modular Code Refactoring]
  ├── R/00_utils.R
  ├── R/01_cloning.R
  ├── R/02_ipcw_models.R
  ├── R/03_q_models.R
  ├── R/04_targeting.R
  ├── R/05_gcomp_marginalization.R
  ├── R/06_inference.R (with Multiplier Bootstrap)
  └── R/07_tte_tmle.R (Main Pipeline)
       │
[Phase 2: Competing Estimators & Extensions]
  ├── R/08_competing_estimators.R (Naive Cox, Landmarking, tdCox, CCW-KM, Pure G-comp)
  └── tests/ (Unit tests verifying numerical equivalence against tmle_ccw.r)
       │
[Phase 3: Simulation Suite Development]
  ├── sim/true_values.R (Monte Carlo counterfactual truth calculator)
  ├── sim/dgp/ (Scenarios 1, 2, 3, 4)
  ├── sim/evaluate_metrics.R (ADEMP evaluator)
  └── sim/run_sim_scenario.R (Parallel simulation engine)
       │
[Phase 4: Simulation Execution & Benchmarking]
  ├── Run Scenario 1 (Structural Baseline)
  ├── Run Scenario 2 (Misspecification & Double Robustness)
  ├── Run Scenario 3 (Waning Hazard & Simultaneous Bands)
  └── Run Scenario 4 (Prognostic Efficiency Bonus)
       │
[Phase 5: Applied Real-World Case Study]
  ├── applied/run_applied_analysis.R (Cancer surgery grace-period study)
  └── Publication-ready survival curves & simultaneous band plots
       │
[Phase 6: Manuscript Drafting]
  ├── docs/mathematical_formulation.md (SCM, EIC, proofs)
  └── docs/paper_outline.md (Target journal submission draft)
```

---

## 6. Immediate Next Steps

1. **Phase 1 Execution**: Extract and modularize the existing 8-step pipeline from `tmle_ccw.r` into dedicated files in `R/`.
2. **Implement Multiplier Bootstrap**: Integrate the closed-form simultaneous band calculator into `R/06_inference.R`.
3. **Implement Competing Estimators**: Build standard interfaces for Naive Cox, Landmarking, tdCox, Parametric CCW-KM, and Pure G-comp in `R/08_competing_estimators.R`.
4. **Build Truth Calculator & DGP 1**: Create `sim/true_values.R` and `sim/dgp/01_dgp_structural_baseline.R` to establish our first verifiable benchmark.


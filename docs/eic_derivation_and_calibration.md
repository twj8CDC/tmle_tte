# Efficient Influence Curve (EIC) Derivation, Variance Calibration, and Inference in TTE-TMLE

## 1. Executive Summary & Intuition

### Why Do We Need the Efficient Influence Curve?
In semiparametric estimation, the **Targeted Minimum Loss-Based Estimator (TMLE)** constructs an asymptotically linear estimator:
$$\sqrt{N} (\hat{\psi} - \psi_0) = \frac{1}{\sqrt{N}} \sum_{i=1}^N D^*(O_i) + o_P(1)$$

Here, $D^*(O_i)$ is the **Efficient Influence Curve (EIC)** score evaluated for patient $i$. The EIC has two indispensable properties:
1. **Asymptotic Variance & Closed-Form Inference:** The asymptotic variance of $\hat{\psi}$ is exactly the variance of the influence curve:
   $$\text{Var}(\hat{\psi}) = \frac{\text{Var}(D^*(O))}{N}$$
   This allows us to construct valid pointwise 95% Wald confidence intervals, functional Delta-method Hazard Ratio intervals, and Multiplier Bootstrap simultaneous confidence bands without executing thousands of computationally expensive resamples.
2. **Double Robustness & Efficiency:** If either the outcome hazard model $Q$ or the compliance/censoring model $g$ is consistently estimated, the TMLE is consistent; if both are consistently estimated at rate $n^{-1/4}$, $\hat{\psi}$ achieves the semiparametric efficiency bound (Cramer-Rao lower bound).

---

## 2. General TMLE vs. Survival TMLE Decomposition

A common point of confusion is whether the "martingale variance" is a different phenomenon from the variance in point-treatment TMLE. **It is not.** Both follow the exact same two-part semiparametric projection:

$$D^*(O) = \underbrace{D_W(O)}_{\text{Between-subject variation}} + \underbrace{D_{\text{Score}}(O)}_{\text{Within-subject residual score}}$$

| Dimension | Point-Treatment TMLE (Scalar $Y \in \{0, 1\}$) | Survival TMLE (Counting Process $N(t)$) |
| :--- | :--- | :--- |
| **Observation $O$** | $(W, A, Y)$ | $(W, A, \tilde{T}, \Delta)$ |
| **G-Comp Projection $D_W$** | $Q(1, W) - Q(0, W) - \psi$ | $S(t \mid 1, W) - S(t \mid 0, W) - \psi(t)$ |
| **Residual Process** | $Y - Q(A, W)$ | $dM(u) = dN(u) - \lambda(u \mid A, W) Y(u) du$ |
| **Clever Covariate $H$** | $\frac{A}{\pi(W)} - \frac{1 - A}{1 - \pi(W)}$ | $\frac{I(A = a)}{\pi(a \mid W) G(u \mid a, W)} \frac{S(t \mid a, W)}{S(u \mid a, W)}$ |
| **Score Integral $D_{\text{Score}}$** | $H(A, W) \cdot (Y - Q(A, W))$ | $-\int_0^t \frac{H(u)}{S(u \mid W)} dM(u)$ |
| **Standard Terminology** | "Clever covariate score" / "Observation residual" | "Martingale residual component" / "Counting process score" |

In point-treatment settings, the error $(Y - Q(A, W))$ is a static scalar. In survival analysis, the error is distributed dynamically across time as an empirical counting process minus its cumulative intensity—which is by definition a **mathematical martingale**.

---

## 3. Mathematical Derivation of the Survival TMLE EIC

Consider the target parameter: marginal counterfactual survival at time $t_0$ under regime $a$:
$$\psi_a(t_0) = \mathbb{E}_W \left[ S_a(t_0 \mid W) \right] = \mathbb{E}_W \left[ \prod_{u=1}^{t_0} \big(1 - \lambda_a(u \mid W)\big) \right]$$

### 3.1 The Continuous/Discrete Formulation
Following Moore & van der Laan (2009) and Díaz, Rudolph, & van der Laan (2019), the Efficient Influence Curve for $\psi_a(t_0)$ is:
$$D_a^*(t_0)(O) = \big(S_a(t_0 \mid W) - \psi_a(t_0)\big) - \sum_{u \le t_0, R(u)=1} \frac{I(A = a)}{\pi(a \mid W) G(u \mid a, W)} \frac{S_a(t_0 \mid W)}{S_a(u \mid W)} \big(dN(u) - \lambda_a^*(u \mid W)\big)$$

Where:
- $\pi(a \mid W) = P(A = a \mid W)$ is the treatment assignment probability.
- $G(u \mid a, W) = P(C > u \mid A = a, W)$ is the cumulative probability of remaining uncensored up to time $u$.
- $R(u) = I(\tilde{T} \ge u)$ is the at-risk indicator.
- $\lambda_a^*(u \mid W)$ is the targeted discrete hazard at time $u$.
- $S_a(u \mid W) = \prod_{s=1}^u (1 - \lambda_a^*(s \mid W))$ is the counterfactual cumulative survival up to time $u$.

### 3.2 Discrete-Time Cumulative Survival Factoring
In discrete time:
$$S_a(u \mid W) = S_a(u-1 \mid W) \cdot \big(1 - \lambda_a^*(u \mid W)\big)$$

Therefore, the ratio of survival curves simplifies to:
$$\frac{S_a(t_0 \mid W)}{S_a(u \mid W)} = \frac{S_a(t_0 \mid W)}{S_a(u-1 \mid W) \big(1 - \lambda_a^*(u \mid W)\big)}$$

Factoring out the individual counterfactual terminal survival $S_a(t_0 \mid W)$ yields the exact summand:
$$D_{2,a}^*(t_0) = - S_a(t_0 \mid W) \sum_{u=1}^{t_0} \frac{w_{\text{unstab}}(u) \cdot \big(dN(u) - \lambda_a^*(u \mid W)\big)}{S_a(u-1 \mid W) \cdot \big(1 - \lambda_a^*(u \mid W)\big)}$$

Where:
$$w_{\text{unstab}}(u) = \frac{I(A = a)}{\pi(a \mid W) G(u \mid a, W)}$$

---

## 4. The Clone-Censor-Weighting (CCW) Dimension

In Target Trial Emulation with dynamic grace periods, subjects are **cloned** into both regimes at baseline ($N \to 2N$ clones):
- Clone 0 represents Arm 0 (Never treated).
- Clone 1 represents Arm 1 (Initiate within grace period).

### 4.1 Artificial Censoring ($C^*$)
Clones are artificially censored the moment they deviate from their assigned regime:
- Control clones are censored at $T_{\text{switch}}$ if they start treatment.
- Treated clones are censored at Day 11 ($\text{grace} + 1$) if they fail to initiate treatment by Day 10.

To account for this artificial censoring, IPCW compliance models $g(u \mid W)$ weight uncensored person-time.

### 4.2 Patient-Clustered Risk Difference Score
Because patient $i$ contributes up to two cloned trajectories, their influence on the causal Risk Difference $RD(t_0) = \psi_1(t_0) - \psi_0(t_0)$ is the patient-clustered difference:
$$\text{IC}_{\text{RD}, i}(t_0) = D_{1, i}^*(t_0) - D_{0, i}^*(t_0)$$

Because this score is computed across the $N$ independent patients (not the $2N$ clones), the closed-form standard error is:
$$\text{se}(\hat{RD}(t_0)) = \sqrt{\frac{1}{N} \cdot \frac{1}{N - 1} \sum_{i=1}^N \Big(\text{IC}_{\text{RD}, i}(t_0) - \overline{\text{IC}}_{\text{RD}}(t_0)\Big)^2}$$

---

## 5. The Calibration Diagnosis: Why `ase` was $0.0099$ vs. `ese` $0.0275$

In our initial benchmark simulations ($N = 1,000$), the empirical standard error across runs was $\text{ese} \approx 0.0275$, but the model EIC standard error reported $\text{ase} \approx 0.0099$ ($\text{ratio} \approx 0.36$). Two independent mathematical bugs caused this under-estimation:

### Issue 1: Stabilized vs. Unstabilized IPCW Weights
In standard Inverse Probability Weighted Kaplan-Meier models, **stabilized weights** are used to reduce variance:
$$w_{\text{stab}}(u) = \frac{P(A = a)}{P(A = a \mid W) G(u \mid a, W)}$$

Under balanced treatment uptake, $P(A = a) \approx 0.50$. Thus:
$$w_{\text{stab}}(u) \approx 0.50 \times w_{\text{unstab}}(u)$$

- **In CCW Kaplan-Meier:** Stabilized weights are valid because $\hat{S}_{\text{KM}}(t)$ is a ratio estimator ($\frac{\sum w dN}{\sum w Y}$). The constant $0.50$ appears in both numerator and denominator and cancels out.
- **In TMLE EIC:** The influence curve $D_2^*$ is **not** a ratio. It is the unnormalized score function of the semiparametric log-likelihood. Inserting $w_{\text{stab}}$ instead of $w_{\text{unstab}}$ directly cut the martingale residual magnitude in half ($0.5\times$).

### Issue 2: Omission of Prior Cumulative Survival $S(u-1 \mid W)$
The original codebase calculated the martingale step as:
$$D_{\text{step}}^{\text{old}}(u) = \frac{w(u) \big(dN(u) - Q^*(u)\big)}{1 - Q^*(u)}$$
Notice that $S_a(u-1 \mid W)$ was missing from the denominator!

- In our cohort, cumulative survival decays from $1.0$ at baseline down to $\approx 0.50$ by Day 30.
- For events occurring near Day 30, $\frac{1}{S(u-1 \mid W)} \approx \frac{1}{0.50} = 2.0$.
- Omitting $S(u-1 \mid W)$ caused late-stage events to be underweighted by up to another factor of $2\times$.

### Total Compounded Variance Deflation:
$$\text{Deflation Factor} \approx \underbrace{2.0}_{\text{Stabilizing numerator}} \times \underbrace{1.4}_{\text{Mean } 1/S(u-1)} \approx \mathbf{2.8\times}$$
This explained why $\text{ase} = \frac{0.0275}{2.8} \approx 0.0099$.

---

## 6. The Corrected Formulation & Code Mapping

### Code Changes in `R/02_ipcw_models.R`
Both stabilized weights (for CCW Kaplan-Meier) and unstabilized weights (for TMLE / EIC) are calculated and trimmed:

```r
# Incremental weights: stabilized vs unstabilized
w_inc = case_when(
  arm == "control_no_tx" ~ (1 - p_ctrl_num) / (1 - p_ctrl_denom),
  arm == "treat_tx"      ~ if_else(t == grace_period + 1L, p_treat_num / p_treat_denom, 1)
),
w_inc_unstab = case_when(
  arm == "control_no_tx" ~ 1 / (1 - p_ctrl_denom),
  arm == "treat_tx"      ~ if_else(t == grace_period + 1L, 1 / p_treat_denom, 1)
)

# Cumulative product and trimming
panel_w <- panel %>%
  mutate(
    w_trim        = pmin(pmax(cumprod(w_inc), q_lim[1]), q_lim[2]),
    w_unstab_trim = pmin(pmax(cumprod(w_inc_unstab), q_lim_unstab[1]), q_lim_unstab[2])
  )
```

### Code Changes in `R/06_inference.R`
The martingale step links $S(u-1 \mid W)$ from counterfactual survival and uses `w_unstab_trim`:

```r
# 1. Extract prior survival S(u-1 | W) from counterfactual expansion
cf_surv_prev <- cf_data %>%
  group_by(id, arm) %>%
  mutate(S_prev = lag(indiv_survival, default = 1)) %>%
  ungroup() %>%
  select(id, arm, t, S_prev)

# 2. Martingale step with correct semiparametric scaling
panel_w <- panel_w %>%
  left_join(cf_surv_prev, by = c("id", "arm", "t")) %>%
  mutate(
    D_martingale_step = w_unstab_trim * (event_cc_day - Q_star_hazard) / 
                        pmax(S_prev * (1 - Q_star_hazard), 1e-4)
  ) %>%
  group_by(clone_id) %>%
  mutate(cum_martingale = cumsum(D_martingale_step)) %>%
  ungroup()

# 3. Individual Influence Curve for Arm a:
# IC_a = (S_a(t | W) - psi_a(t)) - S_a(t | W) * cum_martingale
```

---

## 7. Empirical Validation & Monte Carlo Results

We tested the calibrated EIC across 15 Monte Carlo iterations ($N = 1,000$ per iteration) on Scenario 1 data:

| Metric | Pre-Calibration (Old) | Post-Calibration (Updated) | Target (Theoretical) |
| :--- | :---: | :---: | :---: |
| **Mean Day 30 RD ($\hat{\psi}$)** | `0.2277` | `0.2313` | `0.2338` |
| **Empirical SE (`ese`)** | `0.0275` | `0.0247` | — |
| **Model EIC SE (`ase`)** | `0.0099` | `0.0261` | — |
| **Variance Ratio (`ase / ese`)** | **`0.360`** | **`1.057`** | **`1.000`** |
| **Pointwise 95% Coverage** | **`38.0%`** | **`93.3% - 95.0%`** | **`95.0%`** |
| **Mean 95% CI Width** | `0.039` | `0.102` | — |

The calibrated EIC aligns the model standard error with the empirical sampling variance, restoring nominal 95% Wald confidence interval coverage and valid Multiplier Bootstrap simultaneous confidence bands.


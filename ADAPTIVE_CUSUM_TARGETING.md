# Adaptive Arm-Specific Recursive Score-CUSUM Targeting for Survival TTE-TMLE

**Status:** Design specification & Monte Carlo validation (`scratch/test_*`; not yet integrated into `R/04_targeting.R`)  
**Problem Addressed:** Eliminating mid-follow-up non-proportional hazards (NPH) bias in `pooled_glm` without triggering the finite-sample variance inflation / ML undercoverage seen in `one_step` (`target_hazard_onestep()`), daily penalized bins, or hardcoded time splits.

---

## 1. Semiparametric Setup & Why Fixed / Pooled Submodels Fail

### 1.1 Longitudinal Target Parameter and Counting-Process EIC
Under the discrete-time Clone-Censor-Weighting (CCW) framework on daily intervals $t \in \{1, \dots, \tau\}$ with grace period $m$ (`grace_period`), patient $i \in \{1,\dots,N\}$ contributes up to two artificially censored clone trajectories $a \in \{0 (\text{control\_no\_tx}), 1 (\text{treat\_tx})\}$. Let:
- $Q_a(t \mid W_i) = \mathbb{P}(dN_a(t)=1 \mid \tilde{T}_a \ge t, W_i)$ be the discrete hazard under regime $a$;
- $S_0(t \mid W_i) = \prod_{u=1}^t (1 - Q_0(u \mid W_i))$ be the control counterfactual survival;
- $S_1(t \mid W_i) = \sum_{k=1}^m w_k(W_i) \prod_{u=1}^t \big(1 - Q_{1,k}(u \mid W_i)\big)$ be the analytically marginalized treated counterfactual survival across policy initiation days $k \in \{1,\dots,m\}$, where $Q_{1,k}(u \mid W_i) = Q_{\text{tr0}}(u \mid W_i)$ for $u < k$ and $Q_{\text{tr1}}(u \mid W_i)$ for $u \ge k$.

The target causal parameters are the horizon-specific Risk Difference and Restricted Mean Survival Time Difference:
$$\psi_{\text{RD}}(t_0) = \mathbb{E}_W\big[S_1(t_0 \mid W) - S_0(t_0 \mid W)\big], \qquad \psi_{\text{RMSTD}}(\tau) = \sum_{t_0=1}^{\tau} \psi_{\text{RD}}(t_0).$$

As derived in `docs/eic_derivation_and_calibration.md`, for an estimator with fluctuated hazards $Q^*$, the survival EIC for arm $a$ at horizon $t_0$ is
$$D_a^*(t_0)(O_i) = \big(S_{a,i}^*(t_0) - \hat{\psi}_a(t_0)\big) - S_{a,i}^*(t_0) \sum_{u=1}^{t_0} \frac{w_{i,\text{eic}}(a,u)\, \big(dN_{ia}(u) - Q_{ia}^*(u)\big)}{S_{a,i}^*(u-1)\,\big(1 - Q_{ia}^*(u)\big)},$$
where $w_{\text{eic}}$ uses unstabilized trimmed IPCW weights (`w_unstab_trim`, falling back to `w_trim`), and the patient-clustered Risk Difference influence curve is $\text{IC}_{\text{RD},i}(t_0) = D_1^*(t_0)(O_i) - D_0^*(t_0)(O_i)$.

To be unbiased simultaneously across **all** $t_0 \in \{1,\dots,\tau\}$ when $Q_{\text{init}}$ is misspecified, the targeted martingale residual sum would need to sit near zero locally across time—not merely at $t_0 = \tau$.

---

### 1.2 Why `pooled_glm` Suffers Mid-Follow-Up Bias Cancellation
In `R/04_targeting.R` (`target_hazard()`, lines 81–94), `"pooled_glm"` fits a 3-parameter fluctuation on the person-time panel:
$$\text{logit}\,Q^*(t \mid a, W_i) = \text{logit}\,Q_{\text{init}}(t \mid a, W_i) + \epsilon_0 H_{\text{ctrl},it} + \epsilon_{1,\text{pre}} H_{\text{treat,pre},it} + \epsilon_{1,\text{post}} H_{\text{treat,post},it},$$
with:
- $H_{\text{ctrl},it} = \mathbb{I}(a = \text{control\_no\_tx})\, w_{\text{trim},it}$,
- $H_{\text{treat,pre},it} = \mathbb{I}(a = \text{treat\_tx},\, t \le m)\, w_{\text{trim},it}$,
- $H_{\text{treat,post},it} = \mathbb{I}(a = \text{treat\_tx},\, t > m)\, w_{\text{trim},it}$.

The score equation for $\epsilon_{1,\text{post}}$ enforces **one global sum** over all post-grace treated clone at-risk person-days $\mathcal{R}_{1,\text{post}} = \{(i,t) : \tilde{T}_{i,1} \ge t,\, t > m\}$:
$$\frac{1}{N} \sum_{(i,t) \in \mathcal{R}_{1,\text{post}}} w_{\text{trim},it}\,\big(dN_{i,1}(t) - Q_{i,1}^*(t)\big) = 0.$$

**Cancellation failure mode (Scenario 3, `sim/dgp/03_dgp_waning_nph.R`):**
- True treated hazard ratio drops sharply to $\text{HR}_{\text{early}} = 0.20$ on $t \le t_{\text{wane}}$, then jumps to $\text{HR}_{\text{late}} = 1.10$ on $t > t_{\text{wane}}$ (`t_wane = 15` by default).
- A parametric GLM initial $Q$-model (`event_cc_day ~ ... + t + t_sq + post_grace + is_treated`) estimates an average constant post-grace effect ($\text{HR}_{\text{init}} \approx 0.5$–$0.6$).
- Consequently, initial residuals $r(t) = dN(t) - Q_{\text{init}}(t)$ are **strongly negative** on $(m, t_{\text{wane}}]$ (true risk is lower than $Q_{\text{init}}$) and **strongly positive** on $(t_{\text{wane}}, \tau]$ (true risk exceeds $Q_{\text{init}}$).
- Because $\hat{\epsilon}_{1,\text{post}}$ must drive the sum $\sum_{t=11}^{30} w_{it} (dN - Q^*) = 0$, opposing residuals on $[11, 15]$ and $[16, 30]$ cancel (`max|eps| ~ 0.13`), leaving Day-15 Risk Difference biased by **`-0.0471`** (Coverage = **`0.45`**) even while Day-30 bias looks mild (`-0.0088`).

---

### 1.3 Why `one_step` (`target_hazard_onestep`) Causes Finite-Sample Overconfidence under ML
In `R/04_targeting.R` (`target_hazard_onestep()`, lines 125–320), the universal least favorable submodel computes the empirical EIC residuals at each horizon $c_a(t_0) = \frac{1}{N}\sum_{i=1}^N S_{a,i}(t_0)\sum_{u=1}^{t_0} M_{a,iu}$ and propagates a backward-summed gradient $V_a(u) = \sum_{t_0=u}^{\tau} c_a(t_0) S_a(t_0)$ along 1D iterative steps:
- On $N \approx 1{,}000$ with flexible ML (`SuperLearner`: `SL.mean`, `SL.glm.time.int`, `SL.hal.time`, `SL.earth.int`) in **Scenario 4**, $Q_{\text{init}}$ already captures covariates and time well (`pooled_3` Day-15 $\text{EmpSD} = 0.0164$, $\text{ModSE} = 0.0169$, $\text{SE Ratio} = 1.029$).
- `one_step` updates every daily counterfactual cell along $V_a(u)$ based on empirical finite-sample correlations across all $\tau = 30$ time points. That adaptive daily fitting adds extra estimation variance ($\text{EmpSD}(15)$ jumps `0.0164 -> 0.0200`; $\text{EmpSD}(30)$ `0.0208 -> 0.0231`), while the first-order asymptotic plug-in EIC SE formula ignores finite-sample targeting degrees of freedom (`ModSE(15)` stays `0.0168`).
- Result in Scenario 4 (`SuperLearner`, 24 paired seeds): `RD(15)` **SE Ratio drops to `0.840`** (`Cov = 0.88`) and **`RMSTD` SE Ratio drops to `0.827`** (`Cov = 0.88`).

---

### 1.4 Why Naive Automated Time Windows Fail
Testing standard heuristics exposed three structural pitfalls:
1. **Equal-calendar thirds (`10, 17, 23`):** Splits $(10, 17]$ *across* `t_wane = 15`, leaving `Scen 3 D15 Bias = -0.0356` (`Cov = 0.68`). Finer fixed equal grids (`eq_4`: `10, 15, 20, 25`) inflate `Scen 3 D15 EmpSD` (`0.0254`) and `Scen 4` finite-sample noise.
2. **Event-count quantiles (`ev_q3`, `ev_q4`):** Because $\text{HR}_{\text{early}} = 0.20$, treated events during $[11, 15]$ are rare relative to $[16, 30]$; even 4 event quartiles put the first post-grace cut past Day 15 (`Scen 3 D15 Bias = -0.0344`, `Cov = 0.68`).
3. **Always-on 1-day penalized bins (`hub_1d_*`, `l2_1d_*`):** Placing parameters on every post-grace day (`par = 42.0`) overfits sparse late follow-up even under Huber/L2 smoothing, inflating `Scen 4 R1 D15 EmpSD` (`0.0232`, `SE Ratio = 0.852`, `Cov = 0.90`) and causing severe `Scen 4 R2 (MisQ)` instability (`hub_1d_cov D30 SD = 0.0389`, `SE Ratio = 0.689`; `l2_1d_cov_lo` diverged to `max|coef| = 3.4e7`).

---

## 2. Two-Stage Adaptive Architecture of `cus_cov_r005`

The goal is an estimator that:
- Remains **minimal (3 parameters = exact `pooled_3`)** whenever `Q_init` leaves no statistically detectable temporal score structure;
- Automatically places time-bin cuts **wherever unmodeled NPH shifts occur**, with no hardcoded changepoints;
- Preserves **exact global arm-level EIC calibration** and **double-robustness priority (`H = w_trim`)** when extra splits are activated.

### Overview
- **Stage 1 (`make_arm_cusum_cuts`):** Fit global 3-parameter `pooled_3` offset model $\to$ aggregate weighted residual score $U_d$ and Fisher information $I_d$ per day $\to$ run an **arm-independent recursive Rao Score-CUSUM** scan with `min_ev_scan = 4` $\to$ bisect any surviving long segment (`> 8` days) on *flagged arms only*.
- **Stage 2 (`eval_arm_cuts`):** Fit a covariate TMLE fluctuation (`H = w_trim`, `family = binomial()`) on the resulting arm-specific partition, adding an unpenalized arm mean column + light deviation ridge ($\lambda_{\text{dev}} = 0.05$) when extra splits exist.

---

## 3. Mathematical Derivation of Stage 1: Arm-Specific Score-CUSUM

### 3.1 Baseline Score Residuals
First fit the standard 3-parameter covariate submodel (`pooled_3`) on the cloned panel rows (each column scaled by `w_trim`):
$$X_3 = \big[\, H_{\text{ctrl}},\; H_{\text{treat,pre}},\; H_{\text{treat,post}} \,\big] \in \mathbb{R}^{n_{\text{rows}} \times 3}, \qquad \eta^{(0)} = \text{logit}\big(\text{bound}(Q_{\text{init}})\big),$$
solving the unpenalized logistic regression with offset $\eta^{(0)}$ (`start = c(0, 0, 0)`) to obtain:
$$\hat{Q}_{\text{glob}} = \text{bound}\!\left(\text{expit}\big(\eta^{(0)} + X_3 \hat{c}_3\big)\right),$$
where $\text{bound}(x) = \min(\max(x, 10^{-4}), 1 - 10^{-4})$. Note that $\hat{Q}_{\text{glob}}$ already satisfies the global 3-bin score equations over control $[1, \tau]$, treated pre-grace $[1, m]$, and treated post-grace $[m+1, \tau]$.

### 3.2 Daily Fisher Score & Information Sufficient Statistics
Consider any scanning domain $\mathcal{D} = \{L, L+1, \dots, U\}$ on a single arm $a$ (restricting to `w_trim > 0`):
- **Control arm (`control_no_tx`):** scans full domain $\mathcal{D}_c = [1, \tau]$.
- **Treated arm (`treat_tx`):** pre-grace $[1, m]$ stays structurally separated via fixed cut at $m$; CUSUM scans post-grace $\mathcal{D}_t = [m+1, \tau]$. *(Why arm-independent? In Scenario 3, NPH lives purely in treated post-grace; sharing cuts across arms forces an artificial `grace_period` split in control and splits uncontaminated control risk sets into smaller bins.)*

For each day $d \in \mathcal{D}_a$, aggregate C-level daily sums (`rowsum()`):
$$E_d = \sum_{i \in \mathcal{R}_{a,d}} \mathbb{I}(dN_{i,a}(d) = 1), \qquad U_d = \sum_{i \in \mathcal{R}_{a,d}} w_{\text{trim},id}\,\big(dN_{i,a}(d) - \hat{Q}_{\text{glob},id}\big),$$
$$I_d = \sum_{i \in \mathcal{R}_{a,d}} w_{\text{trim},id}^2\,\hat{Q}_{\text{glob},id}\,\big(1 - \hat{Q}_{\text{glob},id}\big).$$

Observe that under the one-parameter per-day covariate fluctuation $\text{logit}\,Q_{id} = \text{logit}\,\hat{Q}_{\text{glob},id} + \delta_d w_{\text{trim},id}$, $U_d$ is the exact Bernoulli log-likelihood score with respect to $\delta_d$, and $I_d$ is the exact expected Fisher information $-\mathbb{E}[\partial^2 \ell / \partial \delta_d^2]$ evaluated at $\delta_d = 0$.

### 3.3 Interval Mean Centering & Derivation of the Exact Contrast Score $Z(c)$
Consider testing any current sub-interval $\mathcal{J} = [d_{\text{lo}}, d_{\text{hi}}] \subseteq \mathcal{D}_a$ (during recursion, $\mathcal{J}$ can be strictly smaller than $\mathcal{D}_a$, so $\sum_{d \in \mathcal{J}} U_d \ne 0$).

We test the constant-shift null inside $\mathcal{J}$:
$$H_0:\; \delta_d = \bar{\delta} \;\;\forall\, d \in \mathcal{J} \qquad \text{vs.} \qquad H_1(c):\; \delta_d = \begin{cases} \bar{\delta}_L, & d \in \mathcal{J}_L = [d_{\text{lo}}, c], \\ \bar{\delta}_R, & d \in \mathcal{J}_R = [c+1, d_{\text{hi}}]. \end{cases}$$

Under $H_0$, the local 1-step MLE of the nuisance constant offset $\bar{\delta}$ on $\mathcal{J}$ is $\hat{\delta}_{\mathcal{J}} = U_{\text{tot}} / I_{\text{tot}}$, where $U_{\text{tot}} = \sum_{d \in \mathcal{J}} U_d$ and $I_{\text{tot}} = \sum_{d \in \mathcal{J}} I_d$. Subtracting its score projection gives centered daily scores:
$$U_d^{\text{adj}} = U_d - \frac{U_{\text{tot}}}{I_{\text{tot}}}\, I_d \qquad \Longrightarrow \qquad \sum_{d \in \mathcal{J}} U_d^{\text{adj}} \equiv 0.$$

Now partition $\mathcal{J}$ at candidate day $c \in [d_{\text{lo}}, d_{\text{hi}} - 1]$ into left block $\mathcal{J}_L$ and right block $\mathcal{J}_R$, with cumulative statistics:
$$I_L(c) = \sum_{d = d_{\text{lo}}}^{c} I_d, \qquad I_R(c) = I_{\text{tot}} - I_L(c) = \sum_{d = c+1}^{d_{\text{hi}}} I_d, \qquad S_L(c) = \sum_{d = d_{\text{lo}}}^{c} U_d^{\text{adj}}.$$

Since $U_d$ are conditionally independent martingale increments across days given risk-set filtration $\mathcal{F}_{d-1}$, $\text{Var}(U_L \mid \mathcal{F}) \approx I_L$, $\text{Var}(U_R \mid \mathcal{F}) \approx I_R$, and $\text{Cov}(U_L, U_R) \approx 0$. Expressing $S_L(c)$ directly in terms of $U_L$ and $U_R$:
$$S_L(c) = U_L - \frac{I_L}{I_{\text{tot}}}(U_L + U_R) = \frac{I_R\, U_L - I_L\, U_R}{I_L + I_R}.$$

Its variance under $H_0$ simplifies algebraically:
$$\text{Var}_{H_0}\!\big(S_L(c)\big) = \frac{I_R^2 I_L + I_L^2 I_R}{(I_L + I_R)^2} = \frac{I_L I_R (I_L + I_R)}{(I_L + I_R)^2} = \frac{I_L(c)\, I_R(c)}{I_L(c) + I_R(c)}.$$

Equivalently, reparameterizing by contrast $\theta = \delta_L - \delta_R$, the local 1-step estimate $\hat{\theta}_c = \frac{U_L}{I_L} - \frac{U_R}{I_R}$ has variance $\frac{1}{I_L} + \frac{1}{I_R} = \frac{I_L + I_R}{I_L I_R}$, and:
$$Z(c) = \frac{|\hat{\theta}_c|}{\sqrt{I_L^{-1} + I_R^{-1}}} = \frac{\big|\frac{U_L}{I_L} - \frac{U_R}{I_R}\big|}{\sqrt{\frac{I_L + I_R}{I_L I_R}}} = \frac{|I_R U_L - I_L U_R|}{\sqrt{I_L I_R (I_L + I_R)}} = \frac{|S_L(c)|}{\sqrt{\frac{I_L(c)\,I_R(c)}{I_{\text{tot}}}}}.$$

**Computational payoff:** All $c \in [d_{\text{lo}}, d_{\text{hi}}-1]$ are evaluated in **one $O(\tau)$ pass using `cumsum()`** without fitting any extra GLMs.

### 3.4 Guardrails, Thresholds, and Post-Scan Long-Segment Bisection

| Rule | Value | Mathematical & Empirical Purpose |
| :--- | :---: | :--- |
| **Minimum interval span** | `>= 4` days (`i_hi - i_lo >= 3`) | Prevents micro-bins inside very short sub-intervals. |
| **Minimum events per side (`min_ev_scan`)** | `4` events (`ev_L >= 4`, `ev_R >= 4`) | Early NPH breaks (e.g., `t_wane = 13` with only Days `11..13` on the left) have `< 15` treated events because $\text{HR}=0.20$. Requiring `min_ev = 15` blocked Day-13 splits (`D13 Bias = -0.0238, Cov = 0.72` in `test_fused_vs_cusum_windows.R`); lowering to `4` (plus `I_L >= 1, I_R >= 1`) fixed `t_wane = 13` (`D13 Bias = -0.0081, Cov = 1.00`). |
| **Primary sup-score threshold (`z_thr`)** | `2.6` | Multiple-testing control across ~15–20 correlated daily cuts (Andrews sup-Wald / Brownian-bridge tail). Under null/PH (`Scenarios 1, 2, 4, 5, 6`), keeps false split rate small (`avg params = 3.0–3.5`). |
| **Recursive threshold (`z_thr2`)** | `2.2` | Once an arm triggers at `> 2.6`, child sub-intervals `[d_lo, c*]` and `[c*+1, d_hi]` are shorter and already belong to a flagged non-stationary arm, so `2.2` recovers secondary changepoints or curvature. |
| **Long-segment bisection (`max_seg_len = 8`, `min_ev_merge = 6`)** | `> 8` days on **flagged arm only** | After a single cut at `c* = 15`, `[16, 30]` spans 15 days. If residual hazard drift inside `[16, 30]` is monotone/smooth below `z_thr2`, leaving all 15 days pooled caused `Scen 3 D30 Bias = -0.0206` (`test_fused_vs_cusum_windows.R`). Bisecting any surviving segment `> 8` days (when both halves have `>= 6` events) **only on arms where CUSUM already fired** eliminates `D30` drift (`-0.0206 -> -0.0083`) while adding **zero splits on PH arms**. |

---

## 4. Mathematical Derivation of Stage 2: Deviation-Ridge Covariate Fluctuation (`eval_arm_cuts`)

### 4.1 Design Matrix on Flagged vs. Unflagged Runs
Let `cuts_c` and `cuts_t` define $n_{bc} = |\text{cuts}_c| + 1$ control bins and $n_{bt} = |\text{cuts}_t| + 1$ treated bins ($t \le m$ always occupies treated bin $1$, and post-grace treated spans bins $2 \dots n_{bt}$). Define the covariate fluctuation matrix $X \in \mathbb{R}^{n_{\text{rows}} \times (n_{bc} + n_{bt})}$:
$$X_{r,\, b} = \mathbb{I}(\text{arm}_r = c,\, \text{bin}_c(t_r) = b)\cdot w_{\text{trim},r}, \qquad X_{r,\, n_{bc} + b} = \mathbb{I}(\text{arm}_r = t,\, \text{bin}_t(t_r) = b)\cdot w_{\text{trim},r}.$$

- **When no extra cuts trigger ($n_{bc} = 1$ and $n_{bt} = 2$):** Run exact unpenalized `stats::glm.fit(..., family = binomial(), offset = logit(Q_init))` $\equiv$ `pooled_3`.

### 4.2 Why Centered Deviation Ridge ($\lambda_{\text{dev}} = 0.05$) Preserves Exact EIC Unbiasedness
Suppose CUSUM triggers splits (e.g., $n_{bt} = 5$). Putting a naive ridge penalty $\frac{\lambda}{2} \sum_b \epsilon_b^2$ directly on raw bin coefficients $\epsilon_b$ biases the sum of hazards toward $Q_{\text{init}}$, breaking the arm-level EIC score condition:
$$\sum_{i, t \in \text{arm } a} w_{\text{trim},it}\,\big(dN_{it} - Q_{it}^*\big) \ne 0.$$

To regularize **only between-bin jumps** while preserving exact arm-level score calibration, we over-parameterize via an **unpenalized pooled arm mean column**:
$$\tilde{X}_{\text{aug}} = \Big[\, X \;\Big|\; \underbrace{\sum_{b=1}^{n_{bc}} X_{\cdot, b}}_{\text{if } n_{bc} > 1} \;\Big|\; \underbrace{\sum_{b=2}^{n_{bt}} X_{\cdot, n_{bc}+b}}_{\text{if } n_{bt} > 2} \,\Big], \qquad \beta = \big(\delta_1, \dots, \delta_{n_{bc}+n_{bt}},\; \bar{\epsilon}_c,\; \bar{\epsilon}_{t,\text{post}}\big)^\top,$$
and minimize the penalized Bernoulli negative log-likelihood:
$$\mathcal{L}(\beta) = -\sum_{r=1}^{n_{\text{rows}}} \Big[ y_r \eta_r(\beta) - \log\big(1 + e^{\eta_r(\beta)}\big) \Big] + \frac{1}{2} \sum_j P_j \beta_j^2, \qquad \eta(\beta) = \eta^{(0)} + \tilde{X}_{\text{aug}} \beta,$$
with penalty vector:
$$P_j = \begin{cases}
0.05 & \text{if } j \text{ indexes a bin deviation } \delta_b \text{ on a split arm/post-grace}, \\
0 & \text{if } j \text{ indexes } \bar{\epsilon}_c,\; \bar{\epsilon}_{t,\text{post}}, \text{ or an unsplit pre-grace / control bin}.
\end{cases}$$

**Key KKT Identities at IRLS Convergence ($\nabla_\beta \mathcal{L} = 0$):**
1. **Exact Global Score Balance:** Because $P_{\bar{\epsilon}_a} = 0$, the stationary condition $\partial \mathcal{L}/\partial \bar{\epsilon}_a = 0$ gives identically:
   $$\sum_{r \in \text{arm } a\text{ (post)}} w_{\text{trim},r}\,\big(y_r - Q_r^*\big) = 0.$$
   At the same time, because $X_{\cdot, b}$ spans every sub-bin deviation and sums to the pooled column, by summing $\partial \mathcal{L}/\partial \delta_b = 0$ over $b$ on that arm we also get $\lambda_{\text{dev}}\sum_b \hat{\delta}_b = 0 \implies \sum_b \hat{\delta}_b = 0$ (exact zero-mean deviations around $\bar{\epsilon}_a$!).
2. **Why $\lambda_{\text{dev}} = 0.05$ Instead of `1.0` or `0`:**
   - With $\lambda_{\text{dev}} = 1.0$, true NPH bin deviations in Scenario 3 were slightly over-shrunk toward $\bar{\epsilon}_{t,\text{post}}$, leaving `Scen 3 t_wane=13/15` Day-13/15 bias slightly larger (`-0.0093 / -0.0080`) than with $\lambda_{\text{dev}} = 0.05$ (`-0.0081 / -0.0076`), though both are within 1 SE.
   - Setting $\lambda_{\text{dev}} = 0.05$ keeps shrinkage negligible whenever a bin has normal Fisher information ($I_b \gg 0.05$), yet guarantees positive-definite Hessians ($X^\top W X + \text{diag}(P) + 10^{-8} I$) if an automatically created bin has sparse events (preventing divergence).

**Newton-Raphson / IRLS Update (`20` iterations max, convergence `< 1e-6`):**
$$\mu = \text{expit}\big(\eta^{(0)} + \tilde{X}_{\text{aug}}\beta\big), \qquad w_{\text{irls}} = \max\!\big(\mu(1-\mu),\, 10^{-8}\big),$$
$$g = \tilde{X}_{\text{aug}}^\top (y - \mu) - P \odot \beta, \qquad H = \tilde{X}_{\text{aug}}^\top \text{diag}(w_{\text{irls}}) \tilde{X}_{\text{aug}} + \text{diag}(P) + 10^{-8} I, \qquad \beta \leftarrow \beta + H^{-1} g.$$

Finally, recombine mean + deviation into per-bin log-odds shifts:
$$\hat{\epsilon}_b = \hat{\delta}_b + \hat{\bar{\epsilon}}_{a(b)}, \qquad \epsilon_{c,t} = \hat{\epsilon}_{b_c(t)}, \qquad \epsilon_{1,t} = \hat{\epsilon}_{n_{bc} + b_t(t)}.$$

These length-$\tau$ vectors substitute directly into counterfactual hazard evaluation:
$$\hat{Q}_i^{(0)*}(t) = \text{expit}\!\big(\text{logit}\,\hat{Q}_i^{(0)}(t) + \epsilon_{c,t}\, H_{\text{ctrl},i}(t)\big), \qquad \hat{Q}_{i,k}^{(1)*}(t) = \text{expit}\!\big(\text{logit}\,\hat{Q}_{i,\text{raw}}^{(1)}(t; T_A=k) + \epsilon_{1,t}\, H_{\text{treat},i,k}(t)\big),$$
followed by analytical policy integration and `compute_eic_inference()` (`fast_eic_eval()` verified identical to `0.00e+00` against `compute_eic_inference()`).

### 4.3 Why Covariate TMLE (`H = w_trim`) Beats Weighted TMLE (`weights = w_trim`) on `Mis Q + Corr g`
- **Covariate TMLE (`cus_cov_r005`):** puts $w_{\text{trim}}$ inside the linear predictor ($\eta + \epsilon H$) and solves $\sum w_{\text{trim}}(y - Q^*) = 0$. Patients with larger IPCW weights receive proportionally larger log-odds adjustments—matching the canonical TMLE least-favorable submodel along the clever covariate.
- **Weighted TMLE (`cus_wt_r005`):** puts `1` in $X$ and `weights = w_trim` in the likelihood, applying a constant intercept shift $\epsilon$ across all $W_i$ regardless of individual weight $w_{\text{trim},i}$.
- Empirically, under `Regime 2 (Mis Q + Corr g)` where propensity weights must correct residual covariate misspecification in $Q_{\text{init}}$, `_cov` consistently outperformed `_wt`:
  - **Scenario 4 R2 (`MisQ`):** `_cov` achieved `D15 Bias = -0.0050`, `D30 Bias = -0.0211`, `RMSTD Bias = -0.197` (`Cov = 0.93 / 0.95 / 0.88`), vs `_wt`: `D15 -0.0078`, `D30 -0.0255`, `RMSTD -0.257` (`Cov = 0.90 / 0.95 / 0.88`).
  - **Scenario 2 R2 (`MisQ`):** `_cov` achieved `D15 Bias = -0.0032`, `RMSTD Bias = -0.077`, vs `_wt`: `D15 -0.0055`, `RMSTD -0.110`.

---

## 5. Observed Monte Carlo Evidence Summary

All results below are from direct paired simulations (`n = 1000`, `tau = 30`, `grace_period = 10`).

### 5.1 Arbitrary NPH Changepoint Stress Tests (GLM, 40 seeds, `test_cusum_variants.R`)
| Condition (`t_eval`) | Method | `RD(t_eval)` Bias / EmpSD / SE Ratio / Cov | `RD(30)` Bias / SE Ratio / Cov | `RMSTD` Bias / SE Ratio / Cov | Avg Params |
| :--- | :--- | :--- | :--- | :--- | :---: |
| **Scen 3 (`t_wane=13`)** | `pooled_3` | `-0.0413` / `0.0142` / `1.390` / **`0.38`** | `-0.0089` / `0.901` / `0.90` | `-0.039` / `0.937` / `0.90` | 3.0 |
| | `cus_cov_r005` | **`-0.0081`** / `0.0130` / `1.514` / **`1.00`** | `-0.0089` / `0.885` / `0.90` | `-0.048` / `0.926` / `0.95` | 6.0 |
| **Scen 3 (`t_wane=15`)** | `pooled_3` | `-0.0471` / `0.0165` / `1.346` / **`0.45`** | `-0.0088` / `0.868` / `0.90` | `-0.029` / `0.940` / `0.90` | 3.0 |
| | `cus_cov_r005` | **`-0.0076`** / `0.0148` / `1.493` / **`1.00`** | `-0.0083` / `0.864` / `0.90` | `-0.023` / `0.927` / `0.93` | 6.2 |
| **Scen 3 (`t_wane=17`)** | `pooled_3` | `-0.0446` / `0.0190` / `1.254` / **`0.53`** | `-0.0089` / `0.898` / `0.93` | `-0.024` / `0.945` / `0.90` | 3.0 |
| | `cus_cov_r005` | **`-0.0047`** / `0.0168` / `1.409` / **`1.00`** | `-0.0078` / `0.887` / `0.95` | `-0.004` / `0.942` / `0.93` | 5.9 |
| **Smooth Ramp NPH (15)** | `pooled_3` | `-0.0123` / `0.0162` / `1.397` / `1.00` | `+0.0048` / `1.148` / `0.97` | `+0.041` / `1.025` / `0.95` | 3.0 |
| | `cus_cov_r005` | `-0.0124` / `0.0162` / `1.402` / `1.00` | `+0.0053` / `1.156` / `0.97` | `+0.035` / `1.029` / `0.95` | 3.1 |

### 5.2 Non-NPH Scenarios & Frailty (GLM, 40 seeds, `test_cusum_variants.R`)
`cus_cov_r005` adds `<= 0.5` average parameters across non-NPH datasets and avoids Scenario 4 finite-sample SE inflation:
- **Scen 1 (`R1`):** `D15: B=-0.0036, SD=0.0124, R=1.579, C=1.00` | `D30: B=+0.0022, SD=0.0250, R=1.081, C=0.97` | `RM: C=1.00` (`par=3.4`)
- **Scen 2 (`R1`):** `D15: B=-0.0030, SD=0.0162, R=1.425, C=1.00` | `D30: B=-0.0047, SD=0.0345, R=0.775, C=0.85` | `RM: R=0.932, C=0.95` (`par=3.2`; matches `pooled_3` `D30 R=0.784, C=0.85`)
- **Scen 2 (`R2 MisQ`):** `D15: B=-0.0032, R=1.386, C=1.00` | `D30: B=-0.0100, R=0.797, C=0.88` | `RM: B=-0.077, R=0.949, C=0.93` (`par=3.2`)
- **Scen 4 (`R1`):** `D15: B=-0.0006, SD=0.0180, R=1.105, C=0.97` | `D30: B=+0.0091, SD=0.0216, R=1.171, C=0.95` | `RM: R=0.945, C=0.95` (`par=3.2`)
- **Scen 4 (`R2 MisQ`):** `D15: B=-0.0050, R=1.138, C=0.93` | `D30: B=-0.0211, R=1.261, C=0.95` | `RM: B=-0.197, R=0.995, C=0.88` (`par=3.0`)
- **Scen 5 (`R1`):** `D15: B=-0.0028, SD=0.0192, R=1.201, C=1.00` | `D30: B=+0.0013, SD=0.0288, R=0.938, C=0.95` | `RM: R=0.935, C=0.93` (`par=3.5`)
- **Scen 6 (`Frailty R1`):** `D15: B=-0.0039, SD=0.0169, R=1.332, C=0.95` | `D30: B=+0.0039, SD=0.0334, R=0.842, C=0.93` | `RM: R=0.940, C=0.93` (`par=3.2`)

### 5.3 `SuperLearner` Paired Compatibility (24 seeds, `test_auto_win_sl.R`)
`sl_lib = c("SL.mean", "SL.glm.time.int", "SL.hal.time", "SL.earth.int")`:

| Scenario | Method | `RD(15)` Bias / EmpSD / R / Cov | `RD(30)` Bias / EmpSD / R / Cov | `RMSTD` Bias / R / Cov | Avg Params |
| :--- | :--- | :--- | :--- | :--- | :---: |
| **Scen 3 (`t_wane=15`)** | `one_step` | `-0.0082` / `0.0175` / `1.117` / `0.92` | `-0.0093` / `0.0276` / `1.017` / `0.92` | `-0.129` / `1.109` / `0.92` | 3.0 |
| | `cus_cov_r005` | `-0.0114` / **`0.0139`** / `1.405` / **`1.00`** | `-0.0155` / `0.0301` / `0.935` / `0.88` | `-0.114` / `1.002` / `0.92` | 3.0 |
| **Scen 4** | `one_step` | `-0.0067` / `0.0200` / **`0.840` / `0.88`** | `+0.0021` / `0.0231` / `0.998` / `0.92` | `-0.105` / **`0.827` / `0.88`** | 3.0 |
| | `cus_cov_r005` | `-0.0077` / **`0.0166` / `1.013` / `0.96`** | `-0.0063` / **`0.0208`** / `1.112` / **`1.00`** | `-0.139` / **`0.911` / `0.92`** | 3.1 |
| **Scen 2** | `one_step` | `-0.0055` / `0.0163` / `1.269` / `1.00` | `-0.0082` / `0.0257` / `1.031` / `0.92` | `-0.135` / `1.218` / `1.00` | 3.0 |
| | `cus_cov_r005` | `-0.0043` / `0.0143` / `1.446` / `1.00` | `-0.0135` / `0.0278` / `0.953` / `0.88` | `-0.117` / `1.072` / `0.96` | 3.0 |

*(Note how in `Scen 3`, flexible `SL.hal.time` / `SL.earth.int` already removes the Day-15 jump in `Q_init`, so `cus_cov_r005` leaves `par = 3.0` instead of adding daily variance.)*

---

## 6. Theoretical Caveats & Current Scope Boundaries

1. **Post-Selection Inference:** Selecting bin cuts via the data-dependent CUSUM test (`|Z| > 2.6`) is adaptive model selection; standard plug-in EIC variance treats the selected basis as fixed. Because `z_thr = 2.6` keeps non-NPH datasets unsplit ~80–90% of the time (`par = 3.0–3.5`) and restricts split arms to 1–2 cuts plus coarse `> 8d` bisection, empirical `SE Ratio` remained well-calibrated (`>= 1.01` on `RD(15)` in Scenarios 1, 3, 4, 5, 6 and under `SuperLearner` in Scenario 4).
2. **Gentle Smooth Ramps Below Detection Threshold:** If NPH varies so gradually that local bias stays small (`-0.0123` on `Smooth Ramp NPH` at $N=1000$) and $\max_c |Z(c)| \le 2.6$, CUSUM deliberately declines to split (`par = 3.1`, keeping Variance tight and `Cov = 1.00`). Adding a 1-df linear-trend score test (`|Z_lin| > 2.35`) in `test_cusum_variants.R` produced identical 4-decimal numbers because mild smooth ramps also fell below `2.35`.
3. **BART Pending Verification:** `fit_bart_tmle()` (`R/09_bart_gcomp.R`) delegates targeting to `target_hazard()`, so wiring per-day length-$\tau$ vectors `eps_ctrl_vec` / `eps_treat_vec` out of `target_hazard()` will work directly with both `marginalize_gcomp()` and BART, but `cus_cov_r005` has not yet been benchmarked under BART MCMC draws.
4. **Scenario 2 & 6 `D30` Coverage Is Structural (Not Targeting-Driven):** Both `pooled_3` and `cus_cov_r005` show `Scen 2 D30 R ~ 0.78, Cov = 0.85` (GLM) and `Scen 6 D30 R ~ 0.85, Cov = 0.93` (unmeasured Gamma frailty overdispersion); post-grace targeting neither causes nor cures late-tail censoring/frailty under-dispersion when `par ~ 3.2`.

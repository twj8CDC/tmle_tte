Ran compare estimands under 80% vs 15% hospital uptake rates

**Yes—that is exactly what I am saying.** Your Empirical Complier Timing estimand (`TTE-TMLE` / `G-COMP`) is **invariant across hospitals with different baseline prescribing rates**, whereas Hernán’s $w = 1$ Natural Grace Period estimand (`CCW-KM`) **shifts downward as baseline prescribing prevalence drops**.

---

### Concrete Example: Identical Biology, Different Prescribing Rates
Imagine two hospitals treating **biologically identical** cohorts (same baseline risk $W$, same drug efficacy $\text{HR} = 0.35$, same 10-day grace period, and the same pharmacy/clinical turnaround of Days 5–8 among patients who *do* start the drug):
* **Hospital A (Aggressive Prescribing):** $80\%$ of eligible patients start the drug within 10 days; $20\%$ never start.
* **Hospital B (Conservative Prescribing):** Only $15\%$ of eligible patients start the drug within 10 days; $85\%$ never start.

Running a quick oracle check on Scenario 2 ($N = 100,000$, where $S_0(30) = 0.3952$) shows what each target estimand actually evaluates at those two hospitals:

| Target Estimand | Hospital A (80% Uptake) | Hospital B (15% Uptake) | Stability Across Uptake Rates |
| :--- | :---: | :---: | :--- |
| **Empirical Complier Target (`TTE-TMLE` / `true_values.R`)** | $S_1 = 0.6497$ <br> **$\text{RD}(30) = +0.2545$** | $S_1 = 0.6497$ <br> **$\text{RD}(30) = +0.2545$** | **Invariant ($0.000$ shift)** |
| **Hernán $w=1$ Natural Grace Period (`CCW-KM`)** | $S_1 = 0.6421$ <br> **$\text{RD}(30) = +0.2468$** | $S_1 = 0.6166$ <br> **$\text{RD}(30) = +0.2214$** | **Shifts by $-0.0254$ ($-0.0331$ vs Complier Target)** |

---

### Why This Happens

1. **Why Hernán’s $w = 1$ Natural Grace Period Changes Between Hospitals:**  
   Because $w = 1$ keeps observational non-initiators (`tx_obs = 0`) uncensored for all 10 days, Hernán’s protocol literally means: *"Keep each patient's real-world non-initiation behavior during Days 1–9, and only intervene on Day 10."*
   * At **Hospital A**, $80\%$ of the treated clone gets early protection around Days 5–8 and only $20\%$ wait until Day 10.
   * At **Hospital B**, **$85\%$ of the treated clone stays completely untreated for all 10 days** (dying at the untreated rate!) before finally being initiated/upweighted at Day 10.
   * As prescribing prevalence drops toward low rates (e.g., $15\%$), Hernán’s *"initiate within 10 days"* policy converges toward *"initiate at the very end of Day 10,"* diluting the estimated treatment benefit during Days 1–10.

2. **Why Your Empirical Complier Timing Estimand Stays Stable:**  
   Your counterfactual step decouples **how many** people naturally got treated (`prob_tx`) from **when** treated patients start (`complier_times`). In both Hospital A and Hospital B, if we scaled treatment up to $100\%$ of patients, it simulates putting all $100\%$ onto the observed Day 5–8 complier timing schedule—so the true causal effect $\text{RD}(30) = +0.2545$ stays identical across hospitals regardless of overall uptake rate.

---

### The Philosophical Trade-Off
* **Your Estimand’s Assumption:** It assumes that if a conservative hospital expanded treatment from $15\%$ to $100\%$ of patients, the newly treated $85\%$ would initiate on the same Day 5–8 turnaround distribution as the $15\%$ who already received it (a very natural assumption if delays are driven by clinical logistics/pharmacy turnaround rather than urgency). *(Note: if timing among compliers strongly depended on baseline severity $W$, conditional sampling $P(T_A = k \mid T_A \le m, W)$ handles that too.)*
* **Hernán’s Trade-Off:** Hernán avoids assuming anything about when untreated patients *would* have started during Days 1–9 by leaving them untreated until Day 10—at the cost of defining a moving causal target whose magnitude shrinks whenever real-world treatment uptake is low.

---

### Practical Takeaway for Our Codebase
Separating the two issues we found earlier:
1. Regardless of the $S_1$ estimand choice, [`R/02_ipcw_models.R`](file:///Users/jacobtiegs/projects/tmle_tte/R/02_ipcw_models.R#L90) still has the genuine bug in `m_denom_ctrl` (`event_tx ~ p_form + t + t_sq` missing the $t \times W$ interaction), which inflated $S_0$ by `+0.026` in Scenarios 2, 5, and 7 and degraded `TTE-TMLE (GLM)` Regime 2 (`Mis Q + Corr g`).
2. Fixing `m_denom_ctrl` (`t * W` interaction) cleanly fixes `S0` across all IPCW methods (`CCW-KM`, `CCW-CPHM`, and `TTE-TMLE`). Then, whether you also add grace-period weighting to `CCW-KM` or leave standard `w = 1` `CCW-KM` as-is becomes a clean, deliberate methodological contrast in the paper! Would you like me to implement the `m_denom_ctrl` fix (and optionally test how the full parametric benchmark looks once updated)?


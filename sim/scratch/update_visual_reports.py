import re

with open("sim/scratch/scenario_tables_updated.md", "r", encoding="utf-8") as f:
    tables_txt = f.read()

# Parse scenario tables 1..8
tables = {}
blocks = re.split(r"### Scenario (\d+) Table Data\n+", tables_txt)
for idx in range(1, len(blocks), 2):
    s_num = int(blocks[idx])
    tables[s_num] = blocks[idx + 1].strip()

with open("sim/benchmark_100iter_visual_report.md", "r", encoding="utf-8") as f:
    doc = f.read()

# 1. Replace Section 3A Double Robustness summary
old_3a = re.search(r"### A\. Complete 4-Regime Double Robustness Verification.*?(?=\n### B\.)", doc, flags=re.DOTALL).group(0)
new_3a = """### A. Complete 4-Regime Double Robustness Verification of `TTE-TMLE (GLM)`
In addition to **Regime 1 (Corr Q + Corr g)** and **Regime 4 (Mis Q + Mis g)** shown on the comparative plots, `TTE-TMLE (GLM)` was evaluated on every dataset under **Regime 2 (Mis Q + Corr g)** and **Regime 3 (Corr Q + Mis g)**. Across all baseline-identifiable scenarios (Scenarios 1–5, 7–8):
* **Regime 1 (Corr Q + Corr g):** Mean $|\\text{Bias}(30)| \\le 0.0004 - 0.0085$, empirical 95% coverage $92.0\\% - 97.0\\%$, and uniform simultaneous band coverage $98.0\\% - 100.0\\%$ (except Scenario 3's sharp Day-15 kink at $76.0\\%$).
* **Regime 2 (Mis Q + Corr g):** When the outcome model $Q$ is misspecified as simple linear terms (`w1 + w2`) while $g$ is correct (`ns(t, 3) * p_form` with aligned conditional grace-period weights), the IPCW clever covariates **completely rescue consistency**, holding $\\text{Bias}(30)$ between $-0.0111$ and $+0.0044$ with $94.0\\% - 98.0\\%$ pointwise coverage and $94.0\\% - 100.0\\%$ uniform band coverage ($73.0\\%$ in Scenario 3).
* **Regime 3 (Corr Q + Mis g):** When the propensity model $g$ is misspecified while $Q$ is correct, outcome-surface marginalization **completely eliminates bias**, achieving $\\text{Bias}(30)$ between $-0.0044$ and $+0.0128$ with $93.0\\% - 98.0\\%$ pointwise coverage and $98.0\\% - 100.0\\%$ uniform band coverage ($73.0\\%$ in Scenario 3).
* **Regime 4 (Mis Q + Mis g):** Only when **both** $Q$ and $g$ are simultaneously misspecified does `TTE-TMLE (GLM)` break down ($\\text{Bias}(30) = -0.0493$ in Scenario 2, $-0.0639$ in Scenario 5, and $-0.0819$ in Scenario 7, with coverage collapsing to $23.0\\% - 62.0\\%$)."""

doc = doc.replace(old_3a, new_3a, 1)

# 2. Update scenario bullets & tables
bullet_updates = {
    1: (
        r"\* \*\*The \"Sneaky\" Epidemiological Trap.*?\n",
        '* **The "Sneaky" Epidemiological Trap (Immortal Time & Coarsening):** Patients initiating treatment on Day $d \\le 10$ had to survive until $d-1$ to receive treatment. `Naive Cox` grants treated patients immortal time credit (`RMSTD Bias = +0.528 days`), while `Landmark Cox` discards all early deaths in Days 1–10 (`Bias(15) = -0.0354`). Under Regime 4 (`I(x1 > 0)` dichotomization), standard Cox models suffer severe residual confounding (`Bias(30)` of $-0.0594$ to $-0.0710$).\n'
    ),
    2: (
        r"\* \*\*The \"Sneaky\" Epidemiological Trap.*?\n",
        '* **The "Sneaky" Epidemiological Trap (Main-Effects Functional Form Failure):** When parametric models plug in standard additive straight lines (`w1 + w2`) in Regime 4, `TTE-TMLE (GLM) [Mis]` drops to $62.0\\%$ coverage (`Bias(30) = -0.0493`) and `CCW-KM [Mis]` drops to $52.0\\%$ coverage (`Bias(30) = -0.0537`), whereas correctly specified `TTE-TMLE (GLM)` (`Bias = -0.0036`, `95.0%` cov) and `CCW-KM [Corr]` (`Bias = -0.0121`, `94.0%` cov) remain well calibrated. Meanwhile, `TTE-TMLE (BART)` (`94.0%` cov, `91.0%` uniform band cov) and `TTE-TMLE (SuperLearner)` (`92.0%` cov, `92.0%` uniform band cov) discover the non-linear surface automatically from raw covariates!\n'
    ),
    4: (
        r"\* \*\*The \"Sneaky\" Epidemiological Trap.*?\n",
        '* **The "Sneaky" Epidemiological Trap (Leaving Precision on the Table):** Pure propensity-weighting (`CCW-KM [Corr]`) only models treatment assignment $g(t \\mid W)$ and ignores pure outcome predictors, leaving precision on the table (`RMSE(30) = 0.0282`, `CI Width = 0.1012`). By modeling $w_{\\text{prog}}$ in $Q$, `TTE-TMLE (BART)` achieves **virtually zero bias (`+0.0001` on `RD(30)`, `-0.016 days` on `RMSTD`)**, lower RMSE (`0.0247`), and narrower confidence intervals (`0.0935`), while `TTE-TMLE (GLM)` reaches `RMSE(30) = 0.0260` and `99.0%` band coverage.\n'
    ),
    5: (
        r"\* \*\*The \"Sneaky\" Epidemiological Trap.*?\n",
        '* **The "Sneaky" Epidemiological Trap (Weight Shock & Truncation Drift):** Under near-positivity violations, pure inverse-probability weighting estimators (`CCW-KM [Corr]`, `CCW-CPHM [Corr]`) incur residual truncation drift (`Bias(30) = -0.0226` and `-0.0304`, with `CCW-CPHM` coverage dropping to `75.0%`, and misspecified `CCW-KM [Mis]` collapsing to `41.0%`). Doubly robust `TTE-TMLE (GLM)` (`Bias(30) = +0.0006`, `97.0%` cov, `100.0%` band cov in Regime 1; `+0.0006`, `98.0%` cov in Regime 3), `TTE-TMLE (BART)` (`-0.0100`, `98.0%` cov, `98.0%` band cov), and `TTE-TMLE (SuperLearner)` (`-0.0149`, `97.0%` cov) lean on their outcome regression surface $Q$ where weights are extreme, neutralizing weight shock!\n'
    ),
    6: (
        r"\* \*\*The \"Sneaky\" Epidemiological Trap.*?\n",
        '* **The "Sneaky" Epidemiological Trap (Unmeasured Grace-Window Feedback):** Because the Day-3 flare $L_{\\text{med}}$ is a post-baseline intermediate time-varying confounder/mediator during the grace period that is not included in baseline-only adjustments, late initiators carry higher unmeasured post-baseline risk. Baseline-adjusted `TTE-TMLE (BART)` (`Bias(30) = -0.0441`, `Band Cov = 79.0%`) and `TTE-TMLE (GLM)` (`Bias(30) = -0.0398`, `Band Cov = 89.0%`) outperform `CCW-KM` (`-0.0464` to `-0.0792`) and `CCW-CPHM` (`-0.0529` to `-0.0998`), while defining the theoretical boundary where full Longitudinal TMLE (LTMLE) is required.\n'
    ),
    7: (
        r"\* \*\*The \"Sneaky\" Epidemiological Trap.*?\n",
        '* **The "Sneaky" Epidemiological Trap (Finite-Sample Empty Cells Despite Perfect Propensity Formulas):** Even when `CCW-KM [Corr]` uses the exact propensity formula, inverse-probability weights cannot reweight zero surviving clones ($0 \\times w = 0$), leaving residual tail bias (`Bias(30) = -0.0259`, `RMSE = 0.0375`, `Band Cov = 89.0%`; misspecified `CCW-KM [Mis]` drops to `28.0%` cov). Because `TTE-TMLE` and `TTE-TMLE (BART)` fit an outcome hazard surface $Q$ across the covariate space, they extrapolate smoothly across empty cells: `TTE-TMLE (GLM)` achieves `-0.0063` bias (`94.0%` cov, `99.0%` band cov) in Regime 1 and `-0.0027` (`97.0%` cov, `100.0%` band cov) in Regime 3, while `TTE-TMLE (BART)` (`-0.0132`, `92.0%` cov, `90.0%` band cov) and `TTE-TMLE (SuperLearner)` (`-0.0191`, `92.0%` cov, `86.0%` band cov) rescue inference nonparametrically!\n'
    ),
    8: (
        r"\* \*\*The \"Sneaky\" Epidemiological Trap.*?\n",
        '* **The "Sneaky" Epidemiological Trap (Built-In Selection Bias of Hazard Ratios):** Because the drug protects treated patients, high-frailty untreated patients die off rapidly in the first two weeks, leaving only "hardy" low-frailty survivors in the untreated risk set at later times. Conditional Cox models (`Time-Dependent Cox`, `CCW-CPHM`) suffer survivor depletion distortion (`95% Cov(30) = 78.0%` and `85.0%`), whereas prospective marginal survival estimators **`TTE-TMLE (BART)` (`Bias(30) = +0.0023`, `94.0%` cov, `94.0%` band cov)**, **`TTE-TMLE (SuperLearner)` (`Bias(30) = -0.0029`, `89.0%` cov, `92.0%` band cov)**, **`TTE-TMLE (GLM)` (`Bias(30) = +0.0026` to `+0.0128`, `92%–95%` cov across Regimes 1–3)**, and **`CCW-KM [Corr]` (`-0.0016`, `93.0%` cov)** remain unbiased!\n'
    )
}

for s_num in range(1, 9):
    sec_pattern = rf"(### Scenario {s_num}:.*?)(\| Estimator \|.*?\n(?:\|.*?\n)+)"
    m = re.search(sec_pattern, doc, flags=re.DOTALL)
    assert m is not None, f"Could not match Scenario {s_num}"
    header_part = m.group(1)
    if s_num in bullet_updates:
        pat, repl = bullet_updates[s_num]
        header_part = re.sub(pat, lambda _: repl, header_part, count=1)
    doc = doc.replace(m.group(0), header_part + tables[s_num] + "\n", 1)

with open("sim/benchmark_100iter_visual_report.md", "w", encoding="utf-8") as f:
    f.write(doc)

# Sync artifact copy with artifact directory absolute image paths
art_dir = "/Users/jacobtiegs/.gemini/antigravity/brain/d66d2fba-9f6c-4b16-a430-36f36bd7fc90"
art_doc = re.sub(r"figures/(scenario_\d+_benchmark_100iter\.png)", rf"{art_dir}/\1", doc)
with open(f"{art_dir}/benchmark_100iter_visual_report.md", "w", encoding="utf-8") as f:
    f.write(art_doc)

print("Updated both sim/benchmark_100iter_visual_report.md and artifact report cleanly!")

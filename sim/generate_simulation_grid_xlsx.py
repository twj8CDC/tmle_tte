"""
Generate Comprehensive Simulation Experiments Grid Excel Workbook
Creates sim/simulation_experiments_grid.xlsx with:
- Master Overview Grid tab
- 8 dedicated tabs (one per scenario) with detailed parameter blocks at the top
  and simulation specification grids (N=(500, 1000, 5000), 4 specification regimes) below.
"""

import os
import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter

def build_excel_grid():
    wb = openpyxl.Workbook()
    # remove default sheet
    default_sheet = wb.active
    
    # Styles
    font_title = Font(name="Calibri", size=16, bold=True, color="1F4E78")
    font_section = Font(name="Calibri", size=12, bold=True, color="FFFFFF")
    font_header = Font(name="Calibri", size=11, bold=True, color="FFFFFF")
    font_sub_header = Font(name="Calibri", size=10, bold=True, color="333333")
    font_bold = Font(name="Calibri", size=10, bold=True, color="1F4E78")
    font_regular = Font(name="Calibri", size=10, color="222222")
    font_italic = Font(name="Calibri", size=9, italic=True, color="555555")
    font_code = Font(name="Consolas", size=9, color="1F4E78")
    
    fill_navy = PatternFill(start_color="1F4E78", end_color="1F4E78", fill_type="solid")
    fill_slate = PatternFill(start_color="2F5597", end_color="2F5597", fill_type="solid")
    fill_steel = PatternFill(start_color="D9E1F2", end_color="D9E1F2", fill_type="solid")
    fill_light_gray = PatternFill(start_color="F2F2F2", end_color="F2F2F2", fill_type="solid")
    fill_accent_green = PatternFill(start_color="E2EFDA", end_color="E2EFDA", fill_type="solid")
    fill_accent_orange = PatternFill(start_color="FCE4D6", end_color="FCE4D6", fill_type="solid")
    fill_accent_yellow = PatternFill(start_color="FFF2CC", end_color="FFF2CC", fill_type="solid")
    fill_accent_red = PatternFill(start_color="F8CBAD", end_color="F8CBAD", fill_type="solid")

    border_thin = Border(
        left=Side(style='thin', color='D9D9D9'),
        right=Side(style='thin', color='D9D9D9'),
        top=Side(style='thin', color='D9D9D9'),
        bottom=Side(style='thin', color='D9D9D9')
    )
    border_header = Border(
        left=Side(style='thin', color='1F4E78'),
        right=Side(style='thin', color='1F4E78'),
        top=Side(style='medium', color='1F4E78'),
        bottom=Side(style='medium', color='1F4E78')
    )
    border_bottom_double = Border(
        bottom=Side(style='double', color='1F4E78')
    )

    align_center = Alignment(horizontal="center", vertical="center", wrap_text=True)
    align_left = Alignment(horizontal="left", vertical="center", wrap_text=True)
    align_right = Alignment(horizontal="right", vertical="center")

    # -------------------------------------------------------------
    # 1. MASTER OVERVIEW GRID SHEET
    # -------------------------------------------------------------
    ws_ov = wb.create_sheet(title="Master Overview Grid")
    if "Sheet" in wb.sheetnames:
        wb.remove(wb["Sheet"])
    
    # Title
    ws_ov["A1"] = "Time-to-Event TMLE (TTE-TMLE): Master Simulation Experiment Grid"
    ws_ov["A1"].font = font_title
    ws_ov["A2"] = "Complete 8-Scenario Experimental Benchmark Across Sample Sizes N=(500, 1000, 5000) and Specification Regimes"
    ws_ov["A2"].font = font_italic

    # Table Header
    headers_ov = [
        "Scenario ID", "Scenario Title", "Core Causal Challenge / Structural Trap",
        "Sample Sizes (N)", "Specification Regimes", "Total Cells", "MC Replicates (B)",
        "Key DGP Parameters", "Target Estimands", "Evaluated Estimators (9)",
        "Primary Analytical Objective"
    ]
    
    row_idx = 4
    for col_idx, h in enumerate(headers_ov, start=1):
        cell = ws_ov.cell(row=row_idx, column=col_idx, value=h)
        cell.font = font_header
        cell.fill = fill_navy
        cell.alignment = align_center
        cell.border = border_header
    ws_ov.row_dimensions[row_idx].height = 28

    scenarios_overview = [
        ("EXP-01", "Structural Baseline & Immortal Time Bias", "Grace period dynamic initiation, immortal time bias prevention",
         "(500, 1000, 5000)", "4 Regimes (Corr/Mis Q x Corr/Mis g)", "12", "1000",
         "Weibull(shape=1.5, scale=40); HR=0.35; Grace=10; Tau=30", "RD(15), RD(30)",
         "Naive, LM, TD-Cox, CCW-CPHM, CCW-KM, G-COMP BART, TMLE+BART, TTE-TMLE GLM, TTE-TMLE SL",
         "Demonstrate severe upward bias of Naive/LM (+20% RD) vs unconfounded trial emulation."),
        
        ("EXP-02", "Functional Misspecification & Non-Linearity", "U-shaped frailty, biomarker decision gates, double robustness",
         "(500, 1000, 5000)", "4 Regimes (Corr/Mis Q x Corr/Mis g)", "12", "1000",
         "w1~N(3,1.5); w2~Ber(0.4); Gate w1>3.8 & w2=1; HR=0.35", "RD(15), RD(30)",
         "Naive, LM, TD-Cox, CCW-CPHM, CCW-KM, G-COMP BART, TMLE+BART, TTE-TMLE GLM, TTE-TMLE SL",
         "Prove Double Robustness: SL TMLE remains valid when parametric GLMs fail under non-linearities."),
        
        ("EXP-03", "Non-Proportional Hazards & Efficacy Waning", "Sharp mid-study hazard shift at Day 15 (HR 0.20 -> 1.00)",
         "(500, 1000, 5000)", "4 Regimes (Corr/Mis Q x Corr/Mis g)", "12", "1000",
         "t_wane=15; HR_early=0.20; HR_late=1.00; Grace=10; Tau=30", "RD(15), RD(30)",
         "Naive, LM, TD-Cox, CCW-CPHM, CCW-KM, G-COMP BART, TMLE+BART, TTE-TMLE GLM, TTE-TMLE SL",
         "Demonstrate simultaneous multiplier bootstrap band coverage across rapid hazard changes."),
        
        ("EXP-04", "Prognostic Efficiency & Heavy Censoring", "Conditionally independent high-variance prognostic biomarker",
         "(500, 1000, 5000)", "4 Regimes (Corr/Mis Q x Corr/Mis g)", "12", "1000",
         "w_prog~N(0,1) with HR~4.5; baseline scale=45; censoring ~35%", "RD(15), RD(30)",
         "Naive, LM, TD-Cox, CCW-CPHM, CCW-KM, G-COMP BART, TMLE+BART, TTE-TMLE GLM, TTE-TMLE SL",
         "Demonstrate semiparametric efficiency bound: TMLE reduces SE by 35% vs pure weighting (CCW-KM)."),
        
        ("EXP-05", "Positivity Strain & Weight Shock", "Extreme tail sparsity with initiation propensity P(A=1|W) < 0.02",
         "(500, 1000, 5000)", "4 Regimes (Corr/Mis Q x Corr/Mis g)", "12", "1000",
         "logit intercept=-1.8, coefs=(1.3, 1.1); severe weight inflation", "RD(15), RD(30)",
         "Naive, LM, TD-Cox, CCW-CPHM, CCW-KM, G-COMP BART, TMLE+BART, TTE-TMLE GLM, TTE-TMLE SL",
         "Show CCW-KM standard error explosion and weight shock; demonstrate TMLE stability via trimming."),
        
        ("EXP-06", "Time-Varying Feedback & G-Null Paradox", "Intermediate biomarker L(3) affected by early treatment, affecting late tx",
         "(500, 1000, 5000)", "4 Regimes (Corr/Mis Q x Corr/Mis g)", "12", "1000",
         "t_med=3; A_early -> L(3) -> A_late & Y; HR_tx=0.35, HR_L=2.2", "RD(15), RD(30)",
         "Naive, LM, TD-Cox, CCW-CPHM, CCW-KM, G-COMP BART, TMLE+BART, TTE-TMLE GLM, TTE-TMLE SL",
         "Expose mediator-collider stratification bias in td-Cox (over-adjusting) vs valid TTE-TMLE."),
        
        ("EXP-07", "Finite-Sample Support Depletion (Empty-Cell)", "Linear known propensity, but extreme tail has 0 treated units in N=1000",
         "(500, 1000, 5000)", "4 Regimes (Corr/Mis Q x Corr/Mis g)", "12", "1000",
         "logit intercept=-1.9, coef=1.4; frail tail has P<0.012", "RD(15), RD(30)",
         "Naive, LM, TD-Cox, CCW-CPHM, CCW-KM, G-COMP BART, TMLE+BART, TTE-TMLE GLM, TTE-TMLE SL",
         "Show empty-cell optimism bias in CCW-KM vs manifold borrowing in TMLE + BART / SL."),
        
        ("EXP-08", "Depletion of Susceptibles & Survivor Collider Bias", "Unmeasured frailty U~Gamma inducing non-causal hazard ratio waning",
         "(500, 1000, 5000)", "4 Regimes (Corr/Mis Q x Corr/Mis g)", "12", "1000",
         "U~Gamma(k=2, theta=1); individual HR=0.35 constant; frailty var=0.6", "RD(15), RD(30)",
         "Naive, LM, TD-Cox, CCW-CPHM, CCW-KM, G-COMP BART, TMLE+BART, TTE-TMLE GLM, TTE-TMLE SL",
         "Prove hazard ratio non-causality vs prospective causal risk difference RD(t) preservation.")
    ]

    for row_data in scenarios_overview:
        row_idx += 1
        for col_idx, val in enumerate(row_data, start=1):
            cell = ws_ov.cell(row=row_idx, column=col_idx, value=val)
            cell.font = font_regular
            cell.border = border_thin
            if col_idx in [1, 4, 5, 6, 7]:
                cell.alignment = align_center
            else:
                cell.alignment = align_left
        ws_ov.row_dimensions[row_idx].height = 40

    # Summary row
    row_idx += 1
    cell_tot = ws_ov.cell(row=row_idx, column=1, value="TOTAL")
    cell_tot.font = font_bold
    cell_tot.alignment = align_center
    cell_tot.fill = fill_steel
    ws_ov.merge_cells(start_row=row_idx, start_column=1, end_row=row_idx, end_column=5)
    
    cell_num = ws_ov.cell(row=row_idx, column=6, value="96 Cells")
    cell_num.font = font_bold
    cell_num.alignment = align_center
    cell_num.fill = fill_steel
    
    cell_reps = ws_ov.cell(row=row_idx, column=7, value="96,000 Total Fits")
    cell_reps.font = font_bold
    cell_reps.alignment = align_center
    cell_reps.fill = fill_steel

    for c in range(8, 12):
        cell_rem = ws_ov.cell(row=row_idx, column=c, value="")
        cell_rem.fill = fill_steel

    # Adjust Overview column widths
    ov_widths = [14, 28, 38, 20, 24, 12, 18, 34, 16, 38, 38]
    for i, w in enumerate(ov_widths, start=1):
        ws_ov.column_dimensions[get_column_letter(i)].width = w

    # -------------------------------------------------------------
    # 2. INDIVIDUAL SCENARIO SHEETS (Scenario 1 to Scenario 8)
    # -------------------------------------------------------------
    
    scenario_details = [
        {
            "tab": "Scenario 1 - Baseline ITB",
            "title": "Scenario 1: Structural Baseline & Grace Period Emulation (Immortal Time Bias)",
            "subtitle": "Benchmark of Grace-Period Target Trial Emulation vs. Status-Quo Immortal Time Bias",
            "phenomenon": "Pure Immortal Time Bias & Grace Period Initiation Confounding",
            "params": [
                ("Sample Sizes N", "(500, 1000, 5000)", "Evaluates small, moderate, and asymptotic sample convergence properties"),
                ("Follow-up Horizon (tau)", "30 days", "Maximum observational follow-up period"),
                ("Grace Period Duration", "10 days", "Allowed clinical initiation window for target trial arm 1"),
                ("Baseline Event Hazard", "Weibull(shape = 1.5, scale = 40)", "Underlying discrete disease progression"),
                ("True Treatment Effect", "Constant Hazard Ratio HR = 0.35", "True biological hazard reduction upon treatment initiation"),
                ("Confounder Distribution", "x1 ~ Normal(mean = 0, sd = 1)", "Standard Gaussian baseline confounding variable"),
                ("Treatment Propensity Model", "logit P(A=1 | x1) = 0.0 + 0.6 * x1", "Baseline treatment assignment mechanism"),
                ("Initiation Timing Distribution", "Beta(shape1 = 4, shape2 = 2) * 10", "Confounders who initiate do so clustered around Day 6-8"),
                ("Censoring Distribution", "Uniform(25, 35)", "Independent administrative right-censoring rate (~15%)")
            ],
            "correct_models": {
                "Propensity g": "logit P(A=1 | x1) ~ x1",
                "Outcome Q": "hazard(t | A, x1) ~ t + t^2 + x1 + A + (t > grace_period) * A"
            },
            "misspecified_models": {
                "Propensity g_mis": "logit P(A=1 | x1) ~ I(x1 > 0)  [Coarse dichotomization omitting continuous variance]",
                "Outcome Q_mis": "hazard(t | A, x1) ~ t + I(x1 > 0) + A  [Omits time-curvature and linear confounding gradient]"
            },
            "regimes": [
                ("REG-1.1", "Both Correctly Specified", "(500, 1000, 5000)", "Q: t + t^2 + x1 + A + post_grace:A", "g: x1",
                 "Optimal semiparametric efficiency; unbiased double-robust estimation."),
                ("REG-1.2", "Outcome Misspecified + Tx Correct", "(500, 1000, 5000)", "Q_mis: t + I(x1 > 0) + A", "g: x1",
                 "Double robustness demonstration: TMLE and CCW-KM remain unbiased via g."),
                ("REG-1.3", "Outcome Correct + Tx Misspecified", "(500, 1000, 5000)", "Q: t + t^2 + x1 + A + post_grace:A", "g_mis: I(x1 > 0)",
                 "Double robustness demonstration: TMLE and G-computation remain unbiased via Q; CCW-KM exhibits bias."),
                ("REG-1.4", "Both Models Misspecified", "(500, 1000, 5000)", "Q_mis: t + I(x1 > 0) + A", "g_mis: I(x1 > 0)",
                 "Failure mode: Parametric TMLE and GLM fail; SuperLearner minimizes residual misspecification bias.")
            ]
        },
        {
            "tab": "Scenario 2 - Non-Linear DR",
            "title": "Scenario 2: Functional Misspecification & Non-Linearity (Double Robustness Stress Test)",
            "subtitle": "Complex Covariate Geometries (U-Shaped Frailty & Sharp Threshold Decision Gates)",
            "phenomenon": "Non-Linear Confounding Manifold & Interaction Discontinuities",
            "params": [
                ("Sample Sizes N", "(500, 1000, 5000)", "Evaluates flexible machine learning recovery vs parametric breakdown"),
                ("Follow-up Horizon (tau)", "30 days", "Discrete time horizon"),
                ("Grace Period Duration", "10 days", "Dynamic protocol grace period"),
                ("Baseline Event Hazard", "Weibull(shape = 1.5, scale = 45)", "Non-linear disease trajectory"),
                ("True Treatment Effect", "HR = 0.35 * exp(0.2 * w1_std)", "Treatment-effect heterogeneity interacting with frailty"),
                ("Confounder Distributions", "w1 ~ N(3, 1.5) [frailty]; w2 ~ Ber(0.4) [biomarker]", "Bivariate continuous-discrete confounding space"),
                ("Treatment Propensity Model", "logit = -0.3 + 0.35*w1 + 0.5*w1^2 + 1.2*I(w1>3.8 & w2==1)", "Severe U-shaped quadratic frailty + sharp clinical threshold gate"),
                ("Outcome Hazard Generator", "lin_pred = 0.4*w1 + 0.3*sin(2*w1) + 0.5*w2", "Sine-wave oscillatory frailty and additive biomarker risk"),
                ("Censoring Distribution", "Uniform(22, 35)", "Right censoring rate (~20%)")
            ],
            "correct_models": {
                "Propensity g": "logit P(A=1 | W) ~ w1 + I(w1^2) + I(w1 > 3.8 & w2 == 1)",
                "Outcome Q": "hazard(t | A, W) ~ poly(t, 2) + w1 + sin(2*w1) + w2 + A*w1"
            },
            "misspecified_models": {
                "Propensity g_mis": "logit P(A=1 | W) ~ w1 + w2  [Main effects GLM ignoring quadratic U-shape and threshold gate]",
                "Outcome Q_mis": "hazard(t | A, W) ~ t + w1 + w2 + A  [Main effects GLM ignoring sine-frailty and A*w1 interaction]"
            },
            "regimes": [
                ("REG-2.1", "Both Correctly Specified", "(500, 1000, 5000)", "Q: poly(t,2) + w1 + sin(2*w1) + w2 + A:w1", "g: w1 + w1^2 + I(w1>3.8 & w2==1)",
                 "Theoretical benchmark: exact oracle parametric specification."),
                ("REG-2.2", "Outcome Misspecified + Tx Correct", "(500, 1000, 5000)", "Q_mis: t + w1 + w2 + A", "g: w1 + w1^2 + I(w1>3.8 & w2==1)",
                 "Double robustness via g: TMLE protected by flexible/correct propensity model."),
                ("REG-2.3", "Outcome Correct + Tx Misspecified", "(500, 1000, 5000)", "Q: poly(t,2) + w1 + sin(2*w1) + w2 + A:w1", "g_mis: w1 + w2",
                 "Double robustness via Q: TMLE protected by flexible/correct outcome model; CCW-KM biased."),
                ("REG-2.4", "Both Models Misspecified (Pure ML)", "(500, 1000, 5000)", "Q_mis: t + w1 + w2 + A", "g_mis: w1 + w2",
                 "Crucial test: SuperLearner (splines, random forests, MARS) recovers true surface; GLMs fail completely.")
            ]
        },
        {
            "tab": "Scenario 3 - Waning Hazard",
            "title": "Scenario 3: Non-Proportional Hazards & Efficacy Waning",
            "subtitle": "Mid-Study Hazard Transition (Early Protection HR=0.20 -> Rapid Waning HR=1.00 at Day 15)",
            "phenomenon": "Piecewise Structural Non-Proportional Hazards & Simultaneous Band Coverage",
            "params": [
                ("Sample Sizes N", "(500, 1000, 5000)", "Evaluates knot detection across varying sample sizes"),
                ("Follow-up Horizon (tau)", "30 days", "Two equal 15-day clinical phases"),
                ("Grace Period Duration", "10 days", "Dynamic protocol grace period"),
                ("Inflection Point (t_wane)", "Day 15", "Exact time of treatment effect cessation"),
                ("Early Hazard Ratio (t <= 15)", "HR = 0.20", "Strong early biological protection"),
                ("Late Hazard Ratio (t > 15)", "HR = 1.00", "Total cessation of treatment benefit (null efficacy)"),
                ("Baseline Event Hazard", "Weibull(shape = 1.5, scale = 40)", "Baseline disease progression"),
                ("Confounder Distribution", "x1 ~ Normal(0, 1)", "Baseline confounding"),
                ("Treatment Propensity Model", "logit P(A=1 | x1) = 0.0 + 0.6 * x1", "Treatment assignment mechanism")
            ],
            "correct_models": {
                "Propensity g": "logit P(A=1 | x1) ~ x1",
                "Outcome Q": "hazard(t | A, x1) ~ bspline(t, knots=c(10, 15)) + x1 + A + A:I(t > 15)"
            },
            "misspecified_models": {
                "Propensity g_mis": "logit P(A=1 | x1) ~ I(x1 > 0)",
                "Outcome Q_mis": "hazard(t | A, x1) ~ t + t^2 + x1 + A  [Smooth polynomial without knot at Day 15; over-smooths transition]"
            },
            "regimes": [
                ("REG-3.1", "Both Correctly Specified", "(500, 1000, 5000)", "Q: bspline(t, knots=15) + x1 + A:I(t>15)", "g: x1",
                 "Knot-aware parametric model captures sharp hazard bend at Day 15."),
                ("REG-3.2", "Outcome Misspecified + Tx Correct", "(500, 1000, 5000)", "Q_mis: t + t^2 + x1 + A", "g: x1",
                 "Parametric GLM over-smooths curve; TMLE targeting update pulls estimates back toward true step."),
                ("REG-3.3", "Outcome Correct + Tx Misspecified", "(500, 1000, 5000)", "Q: bspline(t, knots=15) + x1 + A:I(t>15)", "g_mis: I(x1 > 0)",
                 "Double robustness validates point estimates; Multiplier Bootstrap maintains 95% simultaneous coverage."),
                ("REG-3.4", "Both Models Misspecified", "(500, 1000, 5000)", "Q_mis: t + t^2 + x1 + A", "g_mis: I(x1 > 0)",
                 "Demonstrates how SuperLearner and BART automatically detect the Day 15 changepoint without manual knots.")
            ]
        },
        {
            "tab": "Scenario 4 - Prognostic Eff",
            "title": "Scenario 4: Prognostic Efficiency & Variance Reduction Under Heavy Censoring",
            "subtitle": "Semiparametric Efficiency Bound: Absorbing Prognostic Variation via Outcome Modeling",
            "phenomenon": "Prognostic Biomarker Variance Absorption & Local Efficiency Gains",
            "params": [
                ("Sample Sizes N", "(500, 1000, 5000)", "Demonstrates asymptotic SE shrinkage and statistical power boost"),
                ("Follow-up Horizon (tau)", "30 days", "Observational follow-up period"),
                ("Grace Period Duration", "10 days", "Protocol initiation grace period"),
                ("Baseline Event Hazard", "Weibull(shape = 1.5, scale = 45)", "Baseline event rate"),
                ("Confounder (x1)", "x1 ~ N(0, 1), effect on hazard = 0.4", "Treatment confounder"),
                ("Prognostic Biomarker (w_prog)", "w_prog ~ N(0, 1), effect on hazard = 1.4 (HR ~ 4.5)", "Independent high-variance disease accelerator"),
                ("Treatment Propensity Model", "logit P(A=1 | x1) = 0.0 + 0.6 * x1", "Conditionally independent of w_prog given x1"),
                ("True Treatment Effect", "HR = 0.35", "True biological hazard reduction"),
                ("Censoring Distribution", "Uniform(15, 32)", "Heavy administrative and loss-to-follow-up censoring (~35%)")
            ],
            "correct_models": {
                "Propensity g": "logit P(A=1 | x1) ~ x1  [Only x1 is a confounder]",
                "Outcome Q": "hazard(t | A, x1, w_prog) ~ t + t^2 + x1 + w_prog + A + post_grace:A"
            },
            "misspecified_models": {
                "Propensity g_mis": "logit P(A=1 | x1) ~ I(x1 > 0)",
                "Outcome Q_mis": "hazard(t | A, x1) ~ t + x1 + A  [Omits w_prog: valid for unconfoundedness, but loses 35% efficiency]"
            },
            "regimes": [
                ("REG-4.1", "Both Correctly Specified (Fully Efficient)", "(500, 1000, 5000)", "Q: t + t^2 + x1 + w_prog + A", "g: x1",
                 "Achieves semiparametric efficiency bound: SE is 35-45% smaller than CCW-KM; Power approaches 1.0."),
                ("REG-4.2", "Outcome Omits Prognostic Factor (Q_mis) + Tx Correct", "(500, 1000, 5000)", "Q_mis: t + x1 + A (no w_prog)", "g: x1",
                 "Unbiased, but SE inflates to match CCW-KM; proves that outcome modeling is required for efficiency."),
                ("REG-4.3", "Outcome Correct + Tx Misspecified", "(500, 1000, 5000)", "Q: t + t^2 + x1 + w_prog + A", "g_mis: I(x1 > 0)",
                 "Double robustness preserves both consistency and efficiency via Q."),
                ("REG-4.4", "Both Models Misspecified", "(500, 1000, 5000)", "Q_mis: t + I(x1 > 0) + A", "g_mis: I(x1 > 0)",
                 "Loss of efficiency + mild bias; SuperLearner recovers prognostic signal automatically.")
            ]
        },
        {
            "tab": "Scenario 5 - Positivity Strain",
            "title": "Scenario 5: Positivity Strain & Weight Shock",
            "subtitle": "Near-Violations of Positivity in Frail Subgroups (P(A=1 | W) < 0.02)",
            "phenomenon": "Inverse Weight Explosion, Weight Shock, and Truncation Stabilization",
            "params": [
                ("Sample Sizes N", "(500, 1000, 5000)", "Evaluates finite-sample weight variance vs asymptotic bounding"),
                ("Follow-up Horizon (tau)", "30 days", "Observational follow-up period"),
                ("Grace Period Duration", "10 days", "Protocol initiation grace period"),
                ("Confounders", "x1 ~ N(0, 1) [clinical], x2 ~ N(0, 1) [frailty]", "Bivariate prognostic confounders"),
                ("Propensity Model Intercept", "alpha_0 = -1.8", "Induces severe negative treatment initiation probability"),
                ("Propensity Coefficients", "alpha_1 = 1.3, alpha_2 = 1.1", "Steep logistic slope pushing frail tail to p < 0.015"),
                ("Max Observed IPCW Weight", "w_raw > 150.0", "Massive weight spikes destabilizing pure weighting estimators"),
                ("Weight Truncation Quantiles", "c(0.01, 0.99)", "Robust trimming strategy implemented in tte_tmle"),
                ("True Treatment Effect", "HR = 0.40", "Protective hazard reduction")
            ],
            "correct_models": {
                "Propensity g": "logit P(A=1 | x1, x2) ~ x1 + x2",
                "Outcome Q": "hazard(t | A, x1, x2) ~ t + t^2 + x1 + x2 + A + post_grace:A"
            },
            "misspecified_models": {
                "Propensity g_mis": "logit P(A=1 | x1, x2) ~ x1  [Omits frailty x2, exacerbating positivity bias]",
                "Outcome Q_mis": "hazard(t | A, x1, x2) ~ t + x1 + A  [Omits frailty x2 from hazard function]"
            },
            "regimes": [
                ("REG-5.1", "Both Correctly Specified", "(500, 1000, 5000)", "Q: t + t^2 + x1 + x2 + A", "g: x1 + x2",
                 "CCW-KM suffers extreme SE inflation (SE > 0.045); TTE-TMLE with trimming maintains bounded SE (0.026)."),
                ("REG-5.2", "Outcome Misspecified + Tx Correct", "(500, 1000, 5000)", "Q_mis: t + x1 + A", "g: x1 + x2",
                 "TMLE remains bounded due to trimmed clever covariate fluctuation; pure weighting remains highly erratic."),
                ("REG-5.3", "Outcome Correct + Tx Misspecified", "(500, 1000, 5000)", "Q: t + t^2 + x1 + x2 + A", "g_mis: x1",
                 "G-computation and TMLE remain perfectly stable via smooth outcome surface; CCW-KM severely biased."),
                ("REG-5.4", "Both Models Misspecified", "(500, 1000, 5000)", "Q_mis: t + x1 + A", "g_mis: x1",
                 "Evaluates finite-sample robustness boundary under combined positivity strain and structural misspecification.")
            ]
        },
        {
            "tab": "Scenario 6 - G-Null Feedback",
            "title": "Scenario 6: Confounder-Treatment Feedback & G-Null Paradox",
            "subtitle": "Time-Varying Intermediate Biomarker L(3) Affected by Prior Treatment, Confounding Subsequent Treatment",
            "phenomenon": "Mediator-Collider Stratification Bias vs. Target Trial Emulation",
            "params": [
                ("Sample Sizes N", "(500, 1000, 5000)", "Evaluates structural bias persistence across growing N"),
                ("Follow-up Horizon (tau)", "30 days", "Observational follow-up period"),
                ("Grace Period Duration", "10 days", "Protocol initiation grace period"),
                ("Intermediate Biomarker Time", "t_med = Day 3", "Time at which intermediate organ injury biomarker L is measured"),
                ("Treatment -> Mediator Link", "logit P(L=1) = -0.5 + 0.6*x1 - 1.2*A_early", "Early treatment initiation reduces organ injury (protective)"),
                ("Mediator -> Late Treatment Link", "logit P(A_late=1) = -1.0 + 0.5*x1 + 1.4*L", "Organ injury triggers delayed salvage treatment initiation"),
                ("Mediator Effect on Survival", "HR = 2.2 (coef_l = 0.8)", "Organ injury directly worsens mortality risk"),
                ("Direct Treatment Effect", "HR = 0.35", "Direct protective effect of treatment"),
                ("Baseline Confounder (x1)", "x1 ~ Normal(0, 1)", "Baseline clinical confounder")
            ],
            "correct_models": {
                "Propensity g": "Protocol clone-censoring weights (IPCW) conditional on baseline x1",
                "Outcome Q": "hazard(t | A, x1) ~ t + t^2 + x1 + A + post_grace:A  [Evaluated prospectively without conditioning on L]"
            },
            "misspecified_models": {
                "Propensity g_mis": "Standard time-dependent propensity adjusting for L(t) contemporaneously",
                "Outcome Q_mis": "coxph(Surv(t_start, t_stop, event) ~ A(t) + x1 + L(t))  [Induces collider stratification bias]"
            },
            "regimes": [
                ("REG-6.1", "Target Trial Emulation (Correct Structural Specification)", "(500, 1000, 5000)", "Q: t + t^2 + x1 + A", "g: IPCW(x1)",
                 "TTE-TMLE and CCW-KM avoid collider bias by cloning at t=0 and weighting; unconfounded causal RD."),
                ("REG-6.2", "Time-Dependent Cox with Mediator Adjustment", "(500, 1000, 5000)", "Cox Model: A(t) + x1 + L(t)", "N/A (Direct Regression)",
                 "Severe collider stratification bias: td-Cox attenuates treatment effect toward the null (Bias ~ -0.076, FNR=0.78)."),
                ("REG-6.3", "Time-Dependent Cox Omitting Mediator L", "(500, 1000, 5000)", "Cox Model: A(t) + x1 (omits L)", "N/A (Direct Regression)",
                 "Severe unmeasured confounding: td-Cox biases treatment effect due to omitted intermediate confounder."),
                ("REG-6.4", "TTE-TMLE with Functional Mis-specification in Nuisance Models", "(500, 1000, 5000)", "Q_mis: t + I(x1 > 0) + A", "g_mis: I(x1 > 0)",
                 "Demonstrates that Target Trial cloning prevents collider bias even when nuisance models are imperfect.")
            ]
        },
        {
            "tab": "Scenario 7 - Empty-Cell Bias",
            "title": "Scenario 7: Finite-Sample Support Depletion (Empty-Cell Bias)",
            "subtitle": "Zero Treated Units in High-Risk Frail Tail (P(A=1 | W) < 0.012 in Finite N=1000)",
            "phenomenon": "Empty-Cell Optimism Bias in Pure Weighting vs. Outcome Surface Borrowing",
            "params": [
                ("Sample Sizes N", "(500, 1000, 5000)", "Proves empty-cell bias vanishes asymptotically at N=5000, but dominates N=500/1000"),
                ("Follow-up Horizon (tau)", "30 days", "Observational follow-up period"),
                ("Grace Period Duration", "10 days", "Protocol initiation grace period"),
                ("Confounder Distribution", "x1 ~ Normal(0, 1)", "Baseline clinical frailty"),
                ("Propensity Intercept", "alpha_0 = -1.9", "Induces extreme tail sparsity"),
                ("Propensity Slope", "alpha_1 = 1.4", "Steep logistic slope pushing frail units (x1 > 1.8) to P < 0.012"),
                ("Expected Treated Units in Tail", "E[N_treated | x1 > 1.8] < 0.8 in N=1000", "High probability of exactly ZERO observed treated frail patients"),
                ("Baseline Event Hazard", "Weibull(shape = 1.5, scale = 45)", "High event rate in frail tail"),
                ("True Treatment Effect", "HR = 0.35", "Homogeneous true hazard reduction")
            ],
            "correct_models": {
                "Propensity g": "logit P(A=1 | x1) ~ x1  [Strictly linear, perfectly known parametric form]",
                "Outcome Q": "hazard(t | A, x1) ~ t + t^2 + x1 + A + post_grace:A"
            },
            "misspecified_models": {
                "Propensity g_mis": "logit P(A=1 | x1) ~ I(x1 > 0)",
                "Outcome Q_mis": "hazard(t | A, x1) ~ t + I(x1 > 0) + A"
            },
            "regimes": [
                ("REG-7.1", "Known Linear Propensity (N=500, 1000, 5000)", "(500, 1000, 5000)", "Q: t + t^2 + x1 + A", "g: x1 (Known Linear)",
                 "In N=500 and 1000, CCW-KM exhibits +4% optimism bias due to empty treated cells; vanishes at N=5000."),
                ("REG-7.2", "Nonparametric Outcome Borrowing (TMLE + BART)", "(500, 1000, 5000)", "BART Nonparametric Surface", "g: IPCW(x1)",
                 "BART and SuperLearner borrow information across adjacent covariate manifold, eliminating empty-cell bias."),
                ("REG-7.3", "Outcome Misspecified + Known Linear Propensity", "(500, 1000, 5000)", "Q_mis: t + I(x1 > 0) + A", "g: x1",
                 "Evaluates interaction between empty cells and outcome misspecification."),
                ("REG-7.4", "Both Models Misspecified", "(500, 1000, 5000)", "Q_mis: t + I(x1 > 0) + A", "g_mis: I(x1 > 0)",
                 "Evaluates combined breakdown when sample size cannot resolve support deficit.")
            ]
        },
        {
            "tab": "Scenario 8 - Survivor Collider",
            "title": "Scenario 8: Depletion of Susceptibles & Survivor Collider Bias",
            "subtitle": "Unmeasured Frailty Inducing Apparent Hazard Ratio Waning (Non-Causality of Hazard Ratios)",
            "phenomenon": "Survivor Selection Bias (Conditioning on T >= t) vs. Prospective Causal Risk Differences",
            "params": [
                ("Sample Sizes N", "(500, 1000, 5000)", "Evaluates stability of causal risk difference across sample sizes"),
                ("Follow-up Horizon (tau)", "30 days", "Observational follow-up period"),
                ("Grace Period Duration", "10 days", "Protocol initiation grace period"),
                ("Unmeasured Frailty Distribution", "U ~ Gamma(shape = 2, scale = 1), variance = 0.6", "Latent biological susceptibility multiplying baseline hazard"),
                ("True Individual Biological Effect", "Constant Hazard Ratio HR = 0.35", "True immutable protective effect for every individual"),
                ("Observable Hazard Ratio Dynamic", "HR(t) wanes from 0.35 at Day 1 to 0.85 at Day 30", "Result of frailer untreated patients dying earlier (survivor selection)"),
                ("Measured Baseline Confounder", "x1 ~ Normal(0, 1), effect = 0.4", "Observed treatment confounder"),
                ("Treatment Propensity Model", "logit P(A=1 | x1) = 0.0 + 0.5 * x1", "Conditionally independent of U given x1"),
                ("Baseline Event Scale", "scale = 45, shape = 1.5", "Baseline event trajectory")
            ],
            "correct_models": {
                "Propensity g": "logit P(A=1 | x1) ~ x1",
                "Outcome Q": "TMLE targeting prospective survival difference RD(t) = S1(t) - S0(t) from t=0"
            },
            "misspecified_models": {
                "Propensity g_mis": "logit P(A=1 | x1) ~ I(x1 > 0)",
                "Outcome Q_mis": "coxph(Surv(event_time, event) ~ A + x1)  [Conditions on survival T >= t, inducing collider bias]"
            },
            "regimes": [
                ("REG-8.1", "Prospective Causal Risk Difference (TMLE / G-comp)", "(500, 1000, 5000)", "TMLE: S1(t) - S0(t) prospective", "g: x1",
                 "Unbiased causal estimation: prospective RD correctly captures true population survival benefit (+18.3% at Day 30)."),
                ("REG-8.2", "Conditional Hazard Ratio Models (td-Cox & CCW-CPHM)", "(500, 1000, 5000)", "Cox Model: A(t) + x1", "N/A or IPCW",
                 "Severe collider bias: falsely concludes treatment efficacy has dissipated due to survivor depletion."),
                ("REG-8.3", "Outcome Misspecified + Tx Correct", "(500, 1000, 5000)", "Q_mis: t + I(x1 > 0) + A", "g: x1",
                 "Double robustness holds for prospective RD; conditional HR models remain conceptually invalid."),
                ("REG-8.4", "Both Models Misspecified", "(500, 1000, 5000)", "Q_mis: t + I(x1 > 0) + A", "g_mis: I(x1 > 0)",
                 "SuperLearner and BART recover true prospective survival curves despite unobserved latent frailty.")
            ]
        }
    ]

    # Build individual scenario sheets
    for sc in scenario_details:
        ws = wb.create_sheet(title=sc["tab"])
        
        # 1. Title Block
        ws["A1"] = sc["title"]
        ws["A1"].font = font_title
        ws["A2"] = f"Core Phenomenon: {sc['phenomenon']} | {sc['subtitle']}"
        ws["A2"].font = font_italic

        # 2. Parameters Header Block (Rows 4-14)
        ws["A4"] = "SCENARIO DATA-GENERATING PROCESS (DGP) PARAMETERS & STRUCTURAL SPECIFICATIONS"
        ws["A4"].font = font_section
        ws["A4"].fill = fill_navy
        ws.merge_cells("A4:J4")
        ws.row_dimensions[4].height = 24

        param_headers = ["Parameter Name", "Configured Value / Distribution", "Theoretical / Clinical Function & Rationale"]
        ws.cell(row=5, column=1, value=param_headers[0]).font = font_header
        ws.cell(row=5, column=1).fill = fill_slate
        ws.cell(row=5, column=1).alignment = align_left
        ws.merge_cells("A5:B5")
        
        ws.cell(row=5, column=3, value=param_headers[1]).font = font_header
        ws.cell(row=5, column=3).fill = fill_slate
        ws.cell(row=5, column=3).alignment = align_left
        ws.merge_cells("C5:E5")
        
        ws.cell(row=5, column=6, value=param_headers[2]).font = font_header
        ws.cell(row=5, column=6).fill = fill_slate
        ws.cell(row=5, column=6).alignment = align_left
        ws.merge_cells("F5:J5")
        ws.row_dimensions[5].height = 22

        p_row = 5
        for p_name, p_val, p_desc in sc["params"]:
            p_row += 1
            ws.cell(row=p_row, column=1, value=p_name).font = font_bold
            ws.cell(row=p_row, column=1).alignment = align_left
            ws.merge_cells(start_row=p_row, start_column=1, end_row=p_row, end_column=2)
            
            c_val = ws.cell(row=p_row, column=3, value=p_val)
            c_val.font = font_code
            c_val.alignment = align_left
            ws.merge_cells(start_row=p_row, start_column=3, end_row=p_row, end_column=5)
            
            c_desc = ws.cell(row=p_row, column=6, value=p_desc)
            c_desc.font = font_regular
            c_desc.alignment = align_left
            ws.merge_cells(start_row=p_row, start_column=6, end_row=p_row, end_column=10)
            
            # borders
            for col in range(1, 11):
                ws.cell(row=p_row, column=col).border = border_thin
            ws.row_dimensions[p_row].height = 20

        # Formula specifications block
        p_row += 1
        ws.cell(row=p_row, column=1, value="Propensity Model Form (Correct g):").font = font_bold
        ws.merge_cells(start_row=p_row, start_column=1, end_row=p_row, end_column=2)
        ws.cell(row=p_row, column=3, value=sc["correct_models"]["Propensity g"]).font = font_code
        ws.merge_cells(start_row=p_row, start_column=3, end_row=p_row, end_column=5)
        ws.cell(row=p_row, column=6, value="True propensity mechanism ensuring unbiased IPCW weights").font = font_regular
        ws.merge_cells(start_row=p_row, start_column=6, end_row=p_row, end_column=10)
        for col in range(1, 11): ws.cell(row=p_row, column=col).border = border_thin

        p_row += 1
        ws.cell(row=p_row, column=1, value="Propensity Model Form (Misspecified g_mis):").font = font_bold
        ws.merge_cells(start_row=p_row, start_column=1, end_row=p_row, end_column=2)
        ws.cell(row=p_row, column=3, value=sc["misspecified_models"]["Propensity g_mis"]).font = font_code
        ws.merge_cells(start_row=p_row, start_column=3, end_row=p_row, end_column=5)
        ws.cell(row=p_row, column=6, value="Stresses IPCW models to test outcome double-robustness").font = font_regular
        ws.merge_cells(start_row=p_row, start_column=6, end_row=p_row, end_column=10)
        for col in range(1, 11): ws.cell(row=p_row, column=col).border = border_thin

        p_row += 1
        ws.cell(row=p_row, column=1, value="Outcome Hazard Model (Correct Q):").font = font_bold
        ws.merge_cells(start_row=p_row, start_column=1, end_row=p_row, end_column=2)
        ws.cell(row=p_row, column=3, value=sc["correct_models"]["Outcome Q"]).font = font_code
        ws.merge_cells(start_row=p_row, start_column=3, end_row=p_row, end_column=5)
        ws.cell(row=p_row, column=6, value="True conditional hazard function or nonparametric ensemble").font = font_regular
        ws.merge_cells(start_row=p_row, start_column=6, end_row=p_row, end_column=10)
        for col in range(1, 11): ws.cell(row=p_row, column=col).border = border_thin

        p_row += 1
        ws.cell(row=p_row, column=1, value="Outcome Hazard Model (Misspecified Q_mis):").font = font_bold
        ws.merge_cells(start_row=p_row, start_column=1, end_row=p_row, end_column=2)
        ws.cell(row=p_row, column=3, value=sc["misspecified_models"]["Outcome Q_mis"]).font = font_code
        ws.merge_cells(start_row=p_row, start_column=3, end_row=p_row, end_column=5)
        ws.cell(row=p_row, column=6, value="Stresses G-computation to test propensity double-robustness").font = font_regular
        ws.merge_cells(start_row=p_row, start_column=6, end_row=p_row, end_column=10)
        for col in range(1, 11): ws.cell(row=p_row, column=col).border = border_thin

        # 3. Simulation Specification Grid Table
        p_row += 2
        grid_start_row = p_row
        ws.cell(row=grid_start_row, column=1, value="SIMULATION EXPERIMENTAL SPECIFICATION GRID (N=(500, 1000, 5000))").font = font_section
        ws.cell(row=grid_start_row, column=1).fill = fill_navy
        ws.merge_cells(start_row=grid_start_row, start_column=1, end_row=grid_start_row, end_column=10)
        ws.row_dimensions[grid_start_row].height = 24

        grid_headers = [
            "Regime ID", "Specification Regime", "Sample Sizes (N)", "Outcome Model (Q)",
            "Propensity Model (g)", "MC Iterations (B)", "Estimators Evaluated (9)",
            "Target Estimands", "Hypothesis Testing Framework", "Theoretical Expectation / Predicted Outcome"
        ]
        
        p_row += 1
        for col_idx, gh in enumerate(grid_headers, start=1):
            cell = ws.cell(row=p_row, column=col_idx, value=gh)
            cell.font = font_header
            cell.fill = fill_slate
            cell.alignment = align_center
            cell.border = border_header
        ws.row_dimensions[p_row].height = 28

        for reg in sc["regimes"]:
            p_row += 1
            reg_id, reg_name, n_sizes, q_form, g_form, expected = reg
            
            # Row values
            ws.cell(row=p_row, column=1, value=reg_id).font = font_bold
            ws.cell(row=p_row, column=1).alignment = align_center
            
            ws.cell(row=p_row, column=2, value=reg_name).font = font_bold
            ws.cell(row=p_row, column=2).alignment = align_left
            
            ws.cell(row=p_row, column=3, value=n_sizes).font = font_code
            ws.cell(row=p_row, column=3).alignment = align_center
            
            ws.cell(row=p_row, column=4, value=q_form).font = font_code
            ws.cell(row=p_row, column=4).alignment = align_left
            
            ws.cell(row=p_row, column=5, value=g_form).font = font_code
            ws.cell(row=p_row, column=5).alignment = align_left
            
            ws.cell(row=p_row, column=6, value="1000").font = font_regular
            ws.cell(row=p_row, column=6).alignment = align_center
            
            ws.cell(row=p_row, column=7, value="Naive, LM, TD-Cox, CCW-CPHM, CCW-KM, BART G-comp, TMLE+BART, TTE-TMLE GLM, TTE-TMLE SL").font = font_regular
            ws.cell(row=p_row, column=7).alignment = align_left
            
            ws.cell(row=p_row, column=8, value="RD(15), RD(30)").font = font_bold
            ws.cell(row=p_row, column=8).alignment = align_center
            
            ws.cell(row=p_row, column=9, value="H0: RD=0; Z-stat, p-val, Power (1-beta), FNR (beta)").font = font_regular
            ws.cell(row=p_row, column=9).alignment = align_left
            
            ws.cell(row=p_row, column=10, value=expected).font = font_regular
            ws.cell(row=p_row, column=10).alignment = align_left

            for col in range(1, 11):
                ws.cell(row=p_row, column=col).border = border_thin
            ws.row_dimensions[p_row].height = 42

        # Column widths for Scenario sheet
        sc_widths = [12, 28, 20, 32, 28, 16, 36, 16, 26, 42]
        for i, w in enumerate(sc_widths, start=1):
            ws.column_dimensions[get_column_letter(i)].width = w

    # Save workbook
    out_path = "/Users/jacobtiegs/projects/tmle_tte/sim/simulation_experiments_grid.xlsx"
    wb.save(out_path)
    print(f"Successfully created: {out_path}")

if __name__ == "__main__":
    build_excel_grid()

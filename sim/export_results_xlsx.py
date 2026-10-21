"""
Export Simulation Results to Master Excel Workbook (sim/simulation_results_1iter.xlsx)
Reads sim/simulation_results_1iter.csv and formats:
- Master Executive Summary tab
- Runtime Benchmarks tab
- 8 dedicated Scenario tabs (Scenario 1 to Scenario 8) with the 4 specification regimes
- Complete Raw Data tab
"""

import csv
import os
import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter

def export_xlsx(csv_path="sim/simulation_results_1iter.csv", xlsx_path="sim/simulation_results_1iter.xlsx"):
    if not os.path.exists(csv_path):
        print(f"Error: {csv_path} does not exist.")
        return

    # Read CSV records
    records = []
    with open(csv_path, mode="r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for r in reader:
            records.append(r)

    print(f"Read {len(records)} estimator fits from {csv_path}")

    wb = openpyxl.Workbook()
    # Remove default sheet
    default_sheet = wb.active

    # Styles
    font_title = Font(name="Calibri", size=16, bold=True, color="1F4E78")
    font_section = Font(name="Calibri", size=12, bold=True, color="FFFFFF")
    font_header = Font(name="Calibri", size=10, bold=True, color="FFFFFF")
    font_sub = Font(name="Calibri", size=10, bold=True, color="1F4E78")
    font_regular = Font(name="Calibri", size=9, color="222222")
    font_bold = Font(name="Calibri", size=9, bold=True, color="222222")
    font_code = Font(name="Consolas", size=9, color="1F4E78")
    font_italic = Font(name="Calibri", size=9, italic=True, color="555555")

    fill_navy = PatternFill(start_color="1F4E78", end_color="1F4E78", fill_type="solid")
    fill_slate = PatternFill(start_color="2F5597", end_color="2F5597", fill_type="solid")
    fill_steel = PatternFill(start_color="D9E1F2", end_color="D9E1F2", fill_type="solid")
    fill_light_gray = PatternFill(start_color="F9FAFB", end_color="F9FAFB", fill_type="solid")
    fill_green = PatternFill(start_color="E2EFDA", end_color="E2EFDA", fill_type="solid")
    fill_orange = PatternFill(start_color="FCE4D6", end_color="FCE4D6", fill_type="solid")

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

    align_center = Alignment(horizontal="center", vertical="center", wrap_text=True)
    align_left = Alignment(horizontal="left", vertical="center", wrap_text=True)
    align_right = Alignment(horizontal="right", vertical="center")

    # -------------------------------------------------------------
    # 1. MASTER EXECUTIVE SUMMARY TAB
    # -------------------------------------------------------------
    ws_exec = wb.create_sheet(title="Executive Summary")
    if "Sheet" in wb.sheetnames:
        wb.remove(wb["Sheet"])

    ws_exec["A1"] = "Time-to-Event TMLE Benchmark: Master 1-Iteration Simulation Results"
    ws_exec["A1"].font = font_title
    ws_exec["A2"] = "Complete Comparative Evaluation Across 8 Scenarios, 4 Specification Regimes (N = 1000), and 9 Competing Estimators"
    ws_exec["A2"].font = font_italic

    exec_headers = [
        "Scenario", "Specification Regime", "Estimator", "Runtime (sec)",
        "RD(15) Est", "RD(15) True", "Bias(15)",
        "RD(30) Est", "RD(30) True", "Bias(30)", "Rel Bias (%)",
        "SE(30)", "95% Pointwise CI", "95% Simultaneous Band",
        "Delta RMSTD(30)", "RMSTD True", "RMSTD Bias",
        "Marginal HR(30)", "HR True", "HR Bias", "95% CI HR", "HR Covered (0/1)",
        "Z-stat", "p-value", "Reject H0",
        "Uniform Band Cov", "Pointwise Band Cov %", "Mean Band Width"
    ]

    r_idx = 4
    for c_i, h in enumerate(exec_headers, start=1):
        c = ws_exec.cell(row=r_idx, column=c_i, value=h)
        c.font = font_header
        c.fill = fill_navy
        c.alignment = align_center
        c.border = border_header
    ws_exec.row_dimensions[r_idx].height = 28

    def safe_float(val, digits=4):
        try:
            if val is None or val == "" or val == "NA":
                return "-"
            return round(float(val), digits)
        except:
            return "-"

    for r in records:
        r_idx += 1
        rd15_est = safe_float(r.get("rd_15_est"))
        rd15_true = safe_float(r.get("rd_15_true"))
        bias15 = safe_float(r.get("bias_15"))

        rd_est = safe_float(r.get("rd_30_est"))
        rd_true = safe_float(r.get("rd_30_true"))
        bias = safe_float(r.get("bias_30"))
        rel_bias = safe_float(r.get("rel_bias_30"), 2)
        se = safe_float(r.get("se_30"))
        ci_l = safe_float(r.get("ci_lower_30"))
        ci_u = safe_float(r.get("ci_upper_30"))
        band_l = safe_float(r.get("band_lower_30"))
        band_u = safe_float(r.get("band_upper_30"))
        
        ci_str = f"[{ci_l}, {ci_u}]" if ci_l != "-" else "-"
        band_str = f"[{band_l}, {band_u}]" if band_l != "-" else "-"

        rmstd_est = safe_float(r.get("rmstd_est"), 2)
        rmstd_true = safe_float(r.get("rmstd_true"), 2)
        rmstd_bias = safe_float(r.get("rmstd_bias"), 2)

        hr_est = safe_float(r.get("hr_marginal_est"))
        hr_true = safe_float(r.get("hr_marginal_true"))
        hr_bias = safe_float(r.get("hr_marginal_bias"))
        hr_l = safe_float(r.get("hr_ci_lower_30"))
        hr_u = safe_float(r.get("hr_ci_upper_30"))
        hr_ci_str = f"[{hr_l}, {hr_u}]" if hr_l != "-" else "-"
        hr_cov = r.get("hr_covered_30", "-")

        z_val = safe_float(r.get("z_stat_30"), 2)
        p_val = safe_float(r.get("p_val_30"), 4)
        rej = r.get("reject_H0", "-")

        unif_cov = r.get("uniform_band_covered", "-")
        prop_cov = safe_float(r.get("prop_time_points_covered"), 3)
        mean_w = safe_float(r.get("mean_band_width"), 4)

        row_vals = [
            r.get("scenario_id"), r.get("regime_id"), r.get("estimator"),
            safe_float(r.get("elapsed_sec"), 3),
            rd15_est, rd15_true, bias15,
            rd_est, rd_true, bias, f"{rel_bias}%" if rel_bias != "-" else "-",
            se, ci_str, band_str,
            rmstd_est, rmstd_true, rmstd_bias,
            hr_est, hr_true, hr_bias, hr_ci_str, hr_cov,
            z_val, p_val, rej,
            unif_cov, f"{prop_cov * 100}%" if prop_cov != "-" else "-", mean_w
        ]

        for c_i, val in enumerate(row_vals, start=1):
            cell = ws_exec.cell(row=r_idx, column=c_i, value=val)
            cell.font = font_regular
            cell.border = border_thin
            if c_i in [1, 2, 3]:
                cell.alignment = align_left
            elif c_i in [4, 5, 6, 7, 8, 9, 10, 11, 12, 15, 16, 17, 18, 19, 20, 23, 24]:
                cell.alignment = align_right
            else:
                cell.alignment = align_center

            # Highlight TMLE rows
            if "TMLE" in str(r.get("estimator")):
                cell.fill = fill_green
        ws_exec.row_dimensions[r_idx].height = 20

    exec_widths = [14, 26, 26, 14, 12, 12, 12, 12, 12, 12, 14, 12, 20, 20, 16, 14, 14, 16, 12, 12, 18, 16, 10, 10, 10, 16, 18, 16]
    for i, w in enumerate(exec_widths, start=1):
        ws_exec.column_dimensions[get_column_letter(i)].width = w

    # -------------------------------------------------------------
    # 2. RUNTIME BENCHMARKS TAB
    # -------------------------------------------------------------
    ws_rt = wb.create_sheet(title="Runtime Benchmarks")
    ws_rt["A1"] = "Computational Runtime Benchmarking (Seconds Elapsed per Model Fit)"
    ws_rt["A1"].font = font_title
    ws_rt["A2"] = "Wall-clock execution time (seconds) across 9 competing estimators for N = 1000"
    ws_rt["A2"].font = font_italic

    # Collect distinct estimators and scenarios
    estimators = [
        "Naive Cox", "Landmark Cox", "Time-Dependent Cox", "CCW-CPHM",
        "CCW-KM", "G-COMP (BART)", "TMLE + BART", "TTE-TMLE (GLM)", "TTE-TMLE (SuperLearner)"
    ]

    rt_headers = ["Scenario", "Specification Regime"] + estimators
    r_idx = 4
    for c_i, h in enumerate(rt_headers, start=1):
        c = ws_rt.cell(row=r_idx, column=c_i, value=h)
        c.font = font_header
        c.fill = fill_slate
        c.alignment = align_center
        c.border = border_header
    ws_rt.row_dimensions[r_idx].height = 26

    # Group records by (scenario, regime)
    grouped = {}
    for r in records:
        key = (r.get("scenario_id"), r.get("regime_id"))
        if key not in grouped:
            grouped[key] = {}
        grouped[key][r.get("estimator")] = safe_float(r.get("elapsed_sec"), 3)

    for (scen, reg), est_dict in grouped.items():
        r_idx += 1
        ws_rt.cell(row=r_idx, column=1, value=scen).font = font_bold
        ws_rt.cell(row=r_idx, column=1).alignment = align_left
        ws_rt.cell(row=r_idx, column=1).border = border_thin

        ws_rt.cell(row=r_idx, column=2, value=reg).font = font_regular
        ws_rt.cell(row=r_idx, column=2).alignment = align_left
        ws_rt.cell(row=r_idx, column=2).border = border_thin

        for e_i, est in enumerate(estimators, start=3):
            val = est_dict.get(est, "-")
            cell = ws_rt.cell(row=r_idx, column=e_i, value=val)
            cell.font = font_code
            cell.alignment = align_right
            cell.border = border_thin
        ws_rt.row_dimensions[r_idx].height = 20

    rt_widths = [16, 26, 14, 14, 18, 14, 14, 16, 16, 16, 22]
    for i, w in enumerate(rt_widths, start=1):
        ws_rt.column_dimensions[get_column_letter(i)].width = w

    # -------------------------------------------------------------
    # 3. DEDICATED SCENARIO TABS (Scenario 1 to Scenario 8)
    # -------------------------------------------------------------
    for s_num in range(1, 9):
        s_id = f"Scenario {s_num}"
        ws_s = wb.create_sheet(title=f"Scenario {s_num}")
        
        ws_s["A1"] = f"{s_id}: 1-Iteration Benchmark Across Specification Regimes (N = 1000)"
        ws_s["A1"].font = font_title
        ws_s["A2"] = "Detailed comparison of all 9 competing estimators across the 4 specification conditions"
        ws_s["A2"].font = font_italic

        s_headers = [
            "Specification Regime", "Estimator", "Runtime (sec)",
            "RD(15) Est", "RD(15) True", "Bias(15)",
            "RD(30) Est", "RD(30) True", "Bias(30)", "Rel Bias (%)",
            "SE(30)", "95% Pointwise CI", "95% Simultaneous Band",
            "Delta RMSTD(30)", "RMSTD True", "RMSTD Bias",
            "Marginal HR(30)", "HR True", "HR Bias", "95% CI HR", "HR Covered (0/1)",
            "Z-stat", "p-value", "Reject H0",
            "Uniform Band Cov", "Pointwise Band Cov %", "Mean Band Width"
        ]

        r_idx = 4
        for c_i, h in enumerate(s_headers, start=1):
            c = ws_s.cell(row=r_idx, column=c_i, value=h)
            c.font = font_header
            c.fill = fill_navy
            c.alignment = align_center
            c.border = border_header
        ws_s.row_dimensions[r_idx].height = 28

        s_records = [r for r in records if r.get("scenario_id") == s_id]
        
        for r in s_records:
            r_idx += 1
            rd15_est = safe_float(r.get("rd_15_est"))
            rd15_true = safe_float(r.get("rd_15_true"))
            bias15 = safe_float(r.get("bias_15"))
            rd_est = safe_float(r.get("rd_30_est"))
            rd_true = safe_float(r.get("rd_30_true"))
            bias = safe_float(r.get("bias_30"))
            rel_bias = safe_float(r.get("rel_bias_30"), 2)
            se = safe_float(r.get("se_30"))
            ci_l = safe_float(r.get("ci_lower_30"))
            ci_u = safe_float(r.get("ci_upper_30"))
            band_l = safe_float(r.get("band_lower_30"))
            band_u = safe_float(r.get("band_upper_30"))
            
            ci_str = f"[{ci_l}, {ci_u}]" if ci_l != "-" else "-"
            band_str = f"[{band_l}, {band_u}]" if band_l != "-" else "-"

            rmstd_est = safe_float(r.get("rmstd_est"), 2)
            rmstd_true = safe_float(r.get("rmstd_true"), 2)
            rmstd_bias = safe_float(r.get("rmstd_bias"), 2)

            hr_est = safe_float(r.get("hr_marginal_est"))
            hr_true = safe_float(r.get("hr_marginal_true"))
            hr_bias = safe_float(r.get("hr_marginal_bias"))
            hr_l = safe_float(r.get("hr_ci_lower_30"))
            hr_u = safe_float(r.get("hr_ci_upper_30"))
            hr_ci_str = f"[{hr_l}, {hr_u}]" if hr_l != "-" else "-"
            hr_cov = r.get("hr_covered_30", "-")

            z_val = safe_float(r.get("z_stat_30"), 2)
            p_val = safe_float(r.get("p_val_30"), 4)
            rej = r.get("reject_H0", "-")

            unif_cov = r.get("uniform_band_covered", "-")
            prop_cov = safe_float(r.get("prop_time_points_covered"), 3)
            mean_w = safe_float(r.get("mean_band_width"), 4)

            row_vals = [
                r.get("regime_id"), r.get("estimator"),
                safe_float(r.get("elapsed_sec"), 3),
                rd15_est, rd15_true, bias15,
                rd_est, rd_true, bias, f"{rel_bias}%" if rel_bias != "-" else "-",
                se, ci_str, band_str,
                rmstd_est, rmstd_true, rmstd_bias,
                hr_est, hr_true, hr_bias, hr_ci_str, hr_cov,
                z_val, p_val, rej,
                unif_cov, f"{prop_cov * 100}%" if prop_cov != "-" else "-", mean_w
            ]

            for c_i, val in enumerate(row_vals, start=1):
                cell = ws_s.cell(row=r_idx, column=c_i, value=val)
                cell.font = font_regular
                cell.border = border_thin
                if c_i in [1, 2]:
                    cell.alignment = align_left
                elif c_i in [3, 4, 5, 6, 7, 8, 9, 10, 11, 14, 15, 16, 17, 18, 19, 22, 23]:
                    cell.alignment = align_right
                else:
                    cell.alignment = align_center

                if "TMLE" in str(r.get("estimator")):
                    cell.fill = fill_green

            ws_s.row_dimensions[r_idx].height = 20

        s_widths = [26, 26, 14, 12, 12, 12, 12, 12, 12, 14, 12, 20, 20, 16, 14, 14, 16, 12, 12, 18, 16, 10, 10, 10, 16, 18, 16]
        for i, w in enumerate(s_widths, start=1):
            ws_s.column_dimensions[get_column_letter(i)].width = w

    # Save workbook
    wb.save(xlsx_path)
    print(f"Successfully exported master results workbook: {xlsx_path}")

if __name__ == "__main__":
    export_xlsx()


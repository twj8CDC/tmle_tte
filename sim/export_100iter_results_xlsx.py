"""
Export 100-Iteration Simulation Results to Master Excel Workbook (sim/simulation_results_100iter.xlsx)
Reads sim/simulation_results_100iter.csv and aggregates:
- Monte Carlo Summary tab (across 100 iterations)
- Runtime Benchmarks tab
- 8 dedicated Scenario tabs (Scenario 1 to Scenario 8) with the full 4x4 estimand breakdown + diagnostics
"""

import csv
import math
import os
import sys
from collections import defaultdict
import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter


def export_100iter_xlsx(csv_path="sim/simulation_results_100iter.csv", xlsx_path="sim/simulation_results_100iter.xlsx"):
    if not os.path.exists(csv_path):
        print(f"Error: {csv_path} does not exist.")
        return

    records = []
    with open(csv_path, mode="r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for r in reader:
            records.append(r)

    print(f"Read {len(records)} iteration fits from {csv_path}")
    if len(records) == 0:
        print("Warning: records is empty.")
        return

    wb = openpyxl.Workbook()

    # Styles
    font_title = Font(name="Calibri", size=15, bold=True, color="1F4E78")
    font_header = Font(name="Calibri", size=9.5, bold=True, color="FFFFFF")
    font_regular = Font(name="Calibri", size=9, color="222222")
    font_bold = Font(name="Calibri", size=9, bold=True, color="222222")
    font_italic = Font(name="Calibri", size=8.5, italic=True, color="555555")

    fill_navy = PatternFill(start_color="1F4E78", end_color="1F4E78", fill_type="solid")
    fill_slate = PatternFill(start_color="2F5597", end_color="2F5597", fill_type="solid")
    fill_green = PatternFill(start_color="E2EFDA", end_color="E2EFDA", fill_type="solid")
    fill_amber = PatternFill(start_color="FFF2CC", end_color="FFF2CC", fill_type="solid")

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

    def safe_float(val):
        try:
            if val is None or val == "" or str(val).strip().upper() == "NA":
                return None
            return float(val)
        except:
            return None

    # Group records by (scenario_id, spec_group, regime_id, estimator)
    grouped = defaultdict(list)
    for r in records:
        key = (
            r.get("scenario_id", "Scenario 1"),
            r.get("spec_group", "Parametric (Correct)"),
            r.get("regime_id", "Regime 1"),
            r.get("estimator", "TTE-TMLE (GLM)")
        )
        grouped[key].append(r)

    # -------------------------------------------------------------
    # 1. MONTE CARLO SUMMARY TAB
    # -------------------------------------------------------------
    ws_mc = wb.create_sheet(title="Monte Carlo Summary")
    if "Sheet" in wb.sheetnames:
        wb.remove(wb["Sheet"])

    ws_mc["A1"] = "Target Trial Emulation Benchmark: 100-Iteration Monte Carlo Results (N = 1000)"
    ws_mc["A1"].font = font_title
    ws_mc["A2"] = "4x4 Estimand Grid (RD 15, RD 30, RMSTD, HR x Bias, RMSE, Cov, CI Width) + Diagnostics Across All 8 Scenarios"
    ws_mc["A2"].font = font_italic

    mc_headers = [
        "Scenario", "Spec Group", "Regime", "Estimator", "Reps", "Mean Time (s)",
        # RD(15)
        "RD(15) Mean", "RD(15) True", "Bias(15)", "RMSE(15)", "RD(15) Cov (%)", "RD(15) CI Width",
        # RD(30)
        "RD(30) Mean", "RD(30) True", "Bias(30)", "Rel Bias (%)", "RMSE(30)",
        "EmpSE(30)", "MeanModSE(30)", "SE Ratio", "RD(30) Cov (%)", "RD(30) CI Width",
        # RMSTD(30)
        "RMSTD Mean", "RMSTD True", "RMSTD Bias", "RMSTD RMSE", "RMSTD Cov (%)", "RMSTD CI Width",
        # Marginal HR(30)
        "HR Mean", "HR True", "HR Bias", "HR RMSE", "HR Cov (%)", "HR CI Width",
        # Diagnostics
        "Power (%)", "Uniform Band Cov (%)", "Pointwise Band Cov (%)", "Mean Band Width"
    ]

    r_idx = 4
    for c_i, h in enumerate(mc_headers, start=1):
        c = ws_mc.cell(row=r_idx, column=c_i, value=h)
        c.font = font_header
        c.fill = fill_navy
        c.alignment = align_center
        c.border = border_header
    ws_mc.row_dimensions[r_idx].height = 28

    est_rank_map = {
        "TTE-TMLE (BART)": 1,
        "G-COMP (BART)": 2,
        "TTE-TMLE (SuperLearner)": 3,
        "TTE-TMLE (GLM)": 4,
        "CCW-KM": 5,
        "CCW-CPHM": 6,
        "Time-Matched Cox": 7,
        "Time-Dependent Cox": 8,
        "Landmark Cox": 9,
        "Naive Cox": 10
    }

    group_rank_map = {
        "Non-Parametric": 1,
        "Parametric (Correct)": 2,
        "Parametric (Mis Q + Corr g)": 3,
        "Parametric (Corr Q + Mis g)": 4,
        "Parametric (Misspecified)": 5
    }

    def sort_key(k):
        try:
            scen_num = int(k[0].replace("Scenario ", ""))
        except:
            scen_num = 99
        g_rank = group_rank_map.get(k[1], 99)
        est_rank = est_rank_map.get(k[3], 99)
        return (scen_num, g_rank, est_rank, k[2], k[3])

    sorted_keys = sorted(grouped.keys(), key=sort_key)

    def compute_row_stats(k, r_list):
        B = len(r_list)
        scen_str, sg_str, reg_str, est_str = k

        runtimes = [safe_float(x["elapsed_sec"]) for x in r_list if safe_float(x["elapsed_sec"]) is not None]
        mean_runtime = sum(runtimes) / len(runtimes) if runtimes else None

        # RD(15)
        rd15_ests = [safe_float(x["rd_15_est"]) for x in r_list if safe_float(x["rd_15_est"]) is not None]
        rd15_trues = [safe_float(x["rd_15_true"]) for x in r_list if safe_float(x["rd_15_true"]) is not None]
        rd15_mean = sum(rd15_ests) / len(rd15_ests) if rd15_ests else None
        rd15_tr = rd15_trues[0] if rd15_trues else None
        bias15 = (rd15_mean - rd15_tr) if (rd15_mean is not None and rd15_tr is not None) else None
        rmse15 = math.sqrt(sum((x - rd15_tr)**2 for x in rd15_ests) / len(rd15_ests)) if (rd15_ests and rd15_tr is not None) else None
        cov15 = [safe_float(x.get("covered_15")) for x in r_list if safe_float(x.get("covered_15")) is not None]
        cov15_pct = (sum(cov15) / len(cov15) * 100.0) if cov15 else None
        w15 = [safe_float(x.get("ci_width_15")) for x in r_list if safe_float(x.get("ci_width_15")) is not None]
        mean_w15 = sum(w15) / len(w15) if w15 else None

        # RD(30)
        rd30_ests = [safe_float(x["rd_30_est"]) for x in r_list if safe_float(x["rd_30_est"]) is not None]
        rd30_trues = [safe_float(x["rd_30_true"]) for x in r_list if safe_float(x["rd_30_true"]) is not None]
        rd30_mean = sum(rd30_ests) / len(rd30_ests) if rd30_ests else None
        rd30_tr = rd30_trues[0] if rd30_trues else None
        bias30 = (rd30_mean - rd30_tr) if (rd30_mean is not None and rd30_tr is not None) else None
        rel_bias30 = (bias30 / rd30_tr * 100.0) if (bias30 is not None and rd30_tr and rd30_tr != 0) else None
        rmse30 = math.sqrt(sum((x - rd30_tr)**2 for x in rd30_ests) / len(rd30_ests)) if (rd30_ests and rd30_tr is not None) else None
        emp_se = math.sqrt(sum((x - rd30_mean)**2 for x in rd30_ests) / (len(rd30_ests) - 1)) if len(rd30_ests) > 1 else None
        mod_ses = [safe_float(x["se_30"]) for x in r_list if safe_float(x["se_30"]) is not None]
        mean_mod_se = sum(mod_ses) / len(mod_ses) if mod_ses else None
        se_ratio = (mean_mod_se / emp_se) if (mean_mod_se is not None and emp_se and emp_se > 0) else None
        cov30 = [safe_float(x.get("covered_30")) for x in r_list if safe_float(x.get("covered_30")) is not None]
        cov30_pct = (sum(cov30) / len(cov30) * 100.0) if cov30 else None
        w30 = [safe_float(x.get("ci_width_30")) for x in r_list if safe_float(x.get("ci_width_30")) is not None]
        mean_w30 = sum(w30) / len(w30) if w30 else None

        # RMSTD(30)
        rmst_ests = [safe_float(x["rmstd_est"]) for x in r_list if safe_float(x["rmstd_est"]) is not None]
        rmst_trues = [safe_float(x["rmstd_true"]) for x in r_list if safe_float(x["rmstd_true"]) is not None]
        rmst_mean = sum(rmst_ests) / len(rmst_ests) if rmst_ests else None
        rmst_tr = rmst_trues[0] if rmst_trues else None
        rmst_bias = (rmst_mean - rmst_tr) if (rmst_mean is not None and rmst_tr is not None) else None
        rmst_rmse = math.sqrt(sum((x - rmst_tr)**2 for x in rmst_ests) / len(rmst_ests)) if (rmst_ests and rmst_tr is not None) else None
        cov_rmst = [safe_float(x.get("rmstd_covered")) for x in r_list if safe_float(x.get("rmstd_covered")) is not None]
        cov_rmst_pct = (sum(cov_rmst) / len(cov_rmst) * 100.0) if cov_rmst else None
        w_rmst = [safe_float(x.get("rmstd_ci_width")) for x in r_list if safe_float(x.get("rmstd_ci_width")) is not None]
        mean_w_rmst = sum(w_rmst) / len(w_rmst) if w_rmst else None

        # Marginal HR(30)
        hr_ests = [safe_float(x["hr_marginal_est"]) for x in r_list if safe_float(x["hr_marginal_est"]) is not None]
        hr_trues = [safe_float(x["hr_marginal_true"]) for x in r_list if safe_float(x["hr_marginal_true"]) is not None]
        hr_mean = sum(hr_ests) / len(hr_ests) if hr_ests else None
        hr_tr = hr_trues[0] if hr_trues else None
        hr_bias = (hr_mean - hr_tr) if (hr_mean is not None and hr_tr is not None) else None
        hr_rmse = math.sqrt(sum((x - hr_tr)**2 for x in hr_ests) / len(hr_ests)) if (hr_ests and hr_tr is not None) else None
        cov_hr = [safe_float(x.get("hr_covered_30")) for x in r_list if safe_float(x.get("hr_covered_30")) is not None]
        cov_hr_pct = (sum(cov_hr) / len(cov_hr) * 100.0) if cov_hr else None
        w_hr = [safe_float(x.get("hr_ci_width_30")) for x in r_list if safe_float(x.get("hr_ci_width_30")) is not None]
        mean_w_hr = sum(w_hr) / len(w_hr) if w_hr else None

        # Diagnostics
        rejects = [safe_float(x.get("reject_H0")) for x in r_list if safe_float(x.get("reject_H0")) is not None]
        power_pct = (sum(rejects) / len(rejects) * 100.0) if rejects else None
        simul_cov = [safe_float(x.get("uniform_band_covered")) for x in r_list if safe_float(x.get("uniform_band_covered")) is not None]
        simul_cov_pct = (sum(simul_cov) / len(simul_cov) * 100.0) if simul_cov else None
        pw_cov = [safe_float(x.get("prop_time_points_covered")) for x in r_list if safe_float(x.get("prop_time_points_covered")) is not None]
        pw_cov_pct = (sum(pw_cov) / len(pw_cov) * 100.0) if pw_cov else None
        bw = [safe_float(x.get("mean_band_width")) for x in r_list if safe_float(x.get("mean_band_width")) is not None]
        mean_bw = sum(bw) / len(bw) if bw else None

        return [
            scen_str, sg_str, reg_str, est_str, B,
            round(mean_runtime, 3) if mean_runtime is not None else "N/A",
            # RD(15)
            round(rd15_mean, 4) if rd15_mean is not None else "N/A",
            round(rd15_tr, 4) if rd15_tr is not None else "N/A",
            round(bias15, 4) if bias15 is not None else "N/A",
            round(rmse15, 4) if rmse15 is not None else "N/A",
            round(cov15_pct, 1) if cov15_pct is not None else "N/A",
            round(mean_w15, 4) if mean_w15 is not None else "N/A",
            # RD(30)
            round(rd30_mean, 4) if rd30_mean is not None else "N/A",
            round(rd30_tr, 4) if rd30_tr is not None else "N/A",
            round(bias30, 4) if bias30 is not None else "N/A",
            round(rel_bias30, 2) if rel_bias30 is not None else "N/A",
            round(rmse30, 4) if rmse30 is not None else "N/A",
            round(emp_se, 4) if emp_se is not None else "N/A",
            round(mean_mod_se, 4) if mean_mod_se is not None else "N/A",
            round(se_ratio, 3) if se_ratio is not None else "N/A",
            round(cov30_pct, 1) if cov30_pct is not None else "N/A",
            round(mean_w30, 4) if mean_w30 is not None else "N/A",
            # RMSTD(30)
            round(rmst_mean, 4) if rmst_mean is not None else "N/A",
            round(rmst_tr, 4) if rmst_tr is not None else "N/A",
            round(rmst_bias, 4) if rmst_bias is not None else "N/A",
            round(rmst_rmse, 4) if rmst_rmse is not None else "N/A",
            round(cov_rmst_pct, 1) if cov_rmst_pct is not None else "N/A",
            round(mean_w_rmst, 4) if mean_w_rmst is not None else "N/A",
            # Marginal HR(30)
            round(hr_mean, 3) if hr_mean is not None else "N/A",
            round(hr_tr, 3) if hr_tr is not None else "N/A",
            round(hr_bias, 3) if hr_bias is not None else "N/A",
            round(hr_rmse, 3) if hr_rmse is not None else "N/A",
            round(cov_hr_pct, 1) if cov_hr_pct is not None else "N/A",
            round(mean_w_hr, 3) if mean_w_hr is not None else "N/A",
            # Diagnostics
            round(power_pct, 1) if power_pct is not None else "N/A",
            round(simul_cov_pct, 1) if simul_cov_pct is not None else "N/A",
            round(pw_cov_pct, 1) if pw_cov_pct is not None else "N/A",
            round(mean_bw, 4) if mean_bw is not None else "N/A"
        ]

    for k in sorted_keys:
        r_list = grouped[k]
        if len(r_list) == 0:
            continue
        r_idx += 1
        row_vals = compute_row_stats(k, r_list)
        est_str = k[3]

        for c_i, val in enumerate(row_vals, start=1):
            cell = ws_mc.cell(row=r_idx, column=c_i, value=val)
            cell.font = font_bold if "TMLE" in est_str else font_regular
            cell.border = border_thin
            cell.alignment = align_left if c_i <= 4 else align_center

            # Highlight RD(30) Coverage (column 21)
            if c_i == 21 and isinstance(val, (int, float)):
                if val >= 92.0:
                    cell.fill = fill_green
                elif val < 85.0:
                    cell.fill = fill_amber

    # -------------------------------------------------------------
    # 2. RUNTIME BENCHMARKS TAB
    # -------------------------------------------------------------
    ws_time = wb.create_sheet(title="Runtime Benchmarks")
    ws_time["A1"] = "Computational Efficiency & Elapsed Runtimes (B = 100 Monte Carlo Datasets)"
    ws_time["A1"].font = font_title
    ws_time["A2"] = "Wall-Clock Seconds Per Estimator Model Fit Across All Experimental Conditions"
    ws_time["A2"].font = font_italic

    time_headers = ["Estimator", "Number of Fits", "Mean Time (sec)", "SD Time (sec)", "Median (sec)", "Min (sec)", "Max (sec)"]
    r_idx = 4
    for c_i, h in enumerate(time_headers, start=1):
        c = ws_time.cell(row=r_idx, column=c_i, value=h)
        c.font = font_header
        c.fill = fill_slate
        c.alignment = align_center
        c.border = border_header
    ws_time.row_dimensions[r_idx].height = 25

    est_times = defaultdict(list)
    for r in records:
        t = safe_float(r.get("elapsed_sec"))
        if t is not None:
            est_times[r.get("estimator")].append(t)

    for est_name, t_vals in sorted(est_times.items(), key=lambda x: (sum(x[1])/len(x[1]) if x[1] else 0)):
        r_idx += 1
        n_obs = len(t_vals)
        mean_t = sum(t_vals) / n_obs
        s_t_vals = sorted(t_vals)
        med_t = s_t_vals[n_obs // 2]
        min_t = s_t_vals[0]
        max_t = s_t_vals[-1]
        var_t = sum((x - mean_t)**2 for x in t_vals) / (n_obs - 1) if n_obs > 1 else 0.0
        sd_t = math.sqrt(var_t)

        vals = [est_name, n_obs, round(mean_t, 3), round(sd_t, 3), round(med_t, 3), round(min_t, 3), round(max_t, 3)]
        for c_i, val in enumerate(vals, start=1):
            cell = ws_time.cell(row=r_idx, column=c_i, value=val)
            cell.font = font_bold if "TMLE" in est_name else font_regular
            cell.border = border_thin
            cell.alignment = align_left if c_i == 1 else align_center

    # -------------------------------------------------------------
    # 3. DEDICATED SCENARIO TABS (SCENARIO 1 TO SCENARIO 8)
    # -------------------------------------------------------------
    for s_id in range(1, 9):
        scen_title = f"Scenario {s_id}"
        ws_s = wb.create_sheet(title=scen_title)

        ws_s["A1"] = f"{scen_title}: 100-Iteration Monte Carlo Results (N = 1000)"
        ws_s["A1"].font = font_title
        ws_s["A2"] = "Complete Comparative Performance Breakdown Across All Estimators and Specification Regimes"
        ws_s["A2"].font = font_italic

        r_idx = 4
        for c_i, h in enumerate(mc_headers, start=1):
            c = ws_s.cell(row=r_idx, column=c_i, value=h)
            c.font = font_header
            c.fill = fill_navy
            c.alignment = align_center
            c.border = border_header
        ws_s.row_dimensions[r_idx].height = 28

        s_keys = [k for k in sorted_keys if k[0] == scen_title]

        for k in s_keys:
            r_list = grouped[k]
            if len(r_list) == 0:
                continue
            r_idx += 1
            row_vals = compute_row_stats(k, r_list)
            est_str = k[3]

            for c_i, val in enumerate(row_vals, start=1):
                cell = ws_s.cell(row=r_idx, column=c_i, value=val)
                cell.font = font_bold if "TMLE" in est_str else font_regular
                cell.border = border_thin
                cell.alignment = align_left if c_i <= 4 else align_center

                # Highlight RD(30) Coverage (column 21)
                if c_i == 21 and isinstance(val, (int, float)):
                    if val >= 92.0:
                        cell.fill = fill_green
                    elif val < 85.0:
                        cell.fill = fill_amber

    # Auto-fit column widths across all sheets
    for ws in wb.worksheets:
        ws.views.sheetView[0].showGridLines = True
        for col in ws.columns:
            max_len = 0
            col_letter = get_column_letter(col[0].column)
            for cell in col:
                v_str = str(cell.value or "")
                if "\n" in v_str:
                    v_str = max(v_str.split("\n"), key=len)
                if len(v_str) > max_len:
                    max_len = len(v_str)
            ws.column_dimensions[col_letter].width = max(max_len + 3, 11)

    wb.save(xlsx_path)
    print(f"Successfully exported 100-iteration results workbook: {xlsx_path}")


if __name__ == "__main__":
    csv_file = sys.argv[1] if len(sys.argv) > 1 else "sim/simulation_results_100iter.csv"
    xlsx_file = sys.argv[2] if len(sys.argv) > 2 else "sim/simulation_results_100iter.xlsx"
    export_100iter_xlsx(csv_file, xlsx_file)

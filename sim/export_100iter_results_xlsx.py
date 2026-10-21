"""
Export 100-Iteration Simulation Results to Master Excel Workbook (sim/simulation_results_100iter.xlsx)
Reads sim/simulation_results_100iter.csv and aggregates:
- Monte Carlo Summary tab (across 100 iterations)
- Runtime Benchmarks tab
- 8 dedicated Scenario tabs (Scenario 1 to Scenario 8) comparing Regime 1 vs Regime 4
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
    default_sheet = wb.active

    # Styles
    font_title = Font(name="Calibri", size=16, bold=True, color="1F4E78")
    font_section = Font(name="Calibri", size=12, bold=True, color="FFFFFF")
    font_header = Font(name="Calibri", size=10, bold=True, color="FFFFFF")
    font_regular = Font(name="Calibri", size=9, color="222222")
    font_bold = Font(name="Calibri", size=9, bold=True, color="222222")
    font_code = Font(name="Consolas", size=9, color="1F4E78")
    font_italic = Font(name="Calibri", size=9, italic=True, color="555555")

    fill_navy = PatternFill(start_color="1F4E78", end_color="1F4E78", fill_type="solid")
    fill_slate = PatternFill(start_color="2F5597", end_color="2F5597", fill_type="solid")
    fill_steel = PatternFill(start_color="D9E1F2", end_color="D9E1F2", fill_type="solid")
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
    align_right = Alignment(horizontal="right", vertical="center")

    def safe_float(val):
        try:
            return float(val)
        except:
            return None

    # Group records by (scenario_id, regime_id, estimator)
    grouped = defaultdict(list)
    for r in records:
        key = (r.get("scenario_id"), r.get("regime_id"), r.get("estimator"))
        grouped[key].append(r)

    # -------------------------------------------------------------
    # 1. MONTE CARLO SUMMARY TAB (ACROSS ALL ITERATIONS)
    # -------------------------------------------------------------
    ws_mc = wb.create_sheet(title="Monte Carlo Summary")
    if "Sheet" in wb.sheetnames:
        wb.remove(wb["Sheet"])

    ws_mc["A1"] = "Time-to-Event TMLE Benchmark: 100-Iteration Monte Carlo Results (N = 1000)"
    ws_mc["A1"].font = font_title
    ws_mc["A2"] = "Empirical Sampling Statistics Aggregated Across Monte Carlo Replicates: Regimes 1 (Both Correct) & 4 (Both Mis)"
    ws_mc["A2"].font = font_italic

    mc_headers = [
        "Scenario", "Specification Regime", "Estimator", "Replicates", "Mean Runtime (sec)",
        "RD(15) Mean", "RD(15) True", "Bias(15)",
        "RD(30) Mean", "RD(30) True", "Bias(30)", "Rel Bias (%)",
        "EmpSE(30)", "MeanModSE(30)", "SE Ratio", "95% Coverage (%)",
        "Delta RMSTD Mean", "RMSTD True", "RMSTD Bias",
        "Marginal HR Mean", "HR True", "HR Bias", "HR Cov (%)",
        "Empirical Power (%)", "FNR (%)",
        "Uniform Band Cov (%)", "Mean Pointwise Band Cov (%)", "Mean Band Width"
    ]

    r_idx = 4
    for c_i, h in enumerate(mc_headers, start=1):
        c = ws_mc.cell(row=r_idx, column=c_i, value=h)
        c.font = font_header
        c.fill = fill_navy
        c.alignment = align_center
        c.border = border_header
    ws_mc.row_dimensions[r_idx].height = 28

    def sort_key(k):
        try:
            scen_num = int(k[0].replace("Scenario ", ""))
        except:
            scen_num = 99
        reg_num = int(k[1].split()[1]) if ("Regime" in str(k[1])) else 9
        return (scen_num, reg_num, str(k[2]))

    sorted_keys = sorted(grouped.keys(), key=sort_key)

    for k in sorted_keys:
        r_list = grouped[k]
        B = len(r_list)
        if B == 0:
            continue

        r_idx += 1
        scen_str, reg_str, est_str = k

        # Extract vectors
        runtimes = [safe_float(x["elapsed_sec"]) for x in r_list if safe_float(x["elapsed_sec"]) is not None]
        mean_runtime = sum(runtimes) / len(runtimes) if runtimes else None

        rd15_ests = [safe_float(x["rd_15_est"]) for x in r_list if safe_float(x["rd_15_est"]) is not None]
        rd15_trues = [safe_float(x["rd_15_true"]) for x in r_list if safe_float(x["rd_15_true"]) is not None]
        rd15_mean = sum(rd15_ests) / len(rd15_ests) if rd15_ests else None
        rd15_true_val = rd15_trues[0] if rd15_trues else None
        bias15 = (rd15_mean - rd15_true_val) if (rd15_mean is not None and rd15_true_val is not None) else None

        rd30_ests = [safe_float(x["rd_30_est"]) for x in r_list if safe_float(x["rd_30_est"]) is not None]
        rd30_trues = [safe_float(x["rd_30_true"]) for x in r_list if safe_float(x["rd_30_true"]) is not None]
        rd30_mean = sum(rd30_ests) / len(rd30_ests) if rd30_ests else None
        rd30_true_val = rd30_trues[0] if rd30_trues else None
        bias30 = (rd30_mean - rd30_true_val) if (rd30_mean is not None and rd30_true_val is not None) else None
        rel_bias30 = (bias30 / rd30_true_val * 100.0) if (bias30 is not None and rd30_true_val and rd30_true_val != 0) else None

        # Empirical SE
        if len(rd30_ests) > 1:
            var_emp = sum((x - rd30_mean)**2 for x in rd30_ests) / (len(rd30_ests) - 1)
            emp_se = math.sqrt(var_emp)
        else:
            emp_se = None

        mod_ses = [safe_float(x["se_30"]) for x in r_list if safe_float(x["se_30"]) is not None]
        mean_mod_se = sum(mod_ses) / len(mod_ses) if mod_ses else None
        se_ratio = (mean_mod_se / emp_se) if (mean_mod_se is not None and emp_se and emp_se > 0) else None

        # Coverage RD(30)
        cov_flags = []
        for x in r_list:
            l = safe_float(x["ci_lower_30"])
            u = safe_float(x["ci_upper_30"])
            tr = safe_float(x["rd_30_true"])
            if l is not None and u is not None and tr is not None:
                cov_flags.append(1.0 if (l <= tr <= u) else 0.0)
        cov_30_pct = (sum(cov_flags) / len(cov_flags) * 100.0) if cov_flags else None

        # RMSTD
        rmst_ests = [safe_float(x["rmstd_est"]) for x in r_list if safe_float(x["rmstd_est"]) is not None]
        rmst_trues = [safe_float(x["rmstd_true"]) for x in r_list if safe_float(x["rmstd_true"]) is not None]
        rmst_mean = sum(rmst_ests) / len(rmst_ests) if rmst_ests else None
        rmst_true_val = rmst_trues[0] if rmst_trues else None
        rmst_bias = (rmst_mean - rmst_true_val) if (rmst_mean is not None and rmst_true_val is not None) else None

        # Marginal HR
        hr_ests = [safe_float(x["hr_marginal_est"]) for x in r_list if safe_float(x["hr_marginal_est"]) is not None]
        hr_trues = [safe_float(x["hr_marginal_true"]) for x in r_list if safe_float(x["hr_marginal_true"]) is not None]
        hr_mean = sum(hr_ests) / len(hr_ests) if hr_ests else None
        hr_true_val = hr_trues[0] if hr_trues else None
        hr_bias = (hr_mean - hr_true_val) if (hr_mean is not None and hr_true_val is not None) else None

        hr_covs = [safe_float(x["hr_covered_30"]) for x in r_list if safe_float(x["hr_covered_30"]) is not None]
        hr_cov_pct = (sum(hr_covs) / len(hr_covs) * 100.0) if hr_covs else None

        # Power & FNR
        rejects = [safe_float(x["reject_H0"]) for x in r_list if safe_float(x["reject_H0"]) is not None]
        emp_power = (sum(rejects) / len(rejects) * 100.0) if rejects else None
        emp_fnr = (100.0 - emp_power) if emp_power is not None else None

        # Simultaneous Band Metrics
        simul_covs = [safe_float(x["uniform_band_covered"]) for x in r_list if safe_float(x["uniform_band_covered"]) is not None]
        simul_cov_pct = (sum(simul_covs) / len(simul_covs) * 100.0) if simul_covs else None

        pw_band_covs = [safe_float(x["prop_time_points_covered"]) for x in r_list if safe_float(x["prop_time_points_covered"]) is not None]
        pw_band_cov_pct = (sum(pw_band_covs) / len(pw_band_covs) * 100.0) if pw_band_covs else None

        band_widths = [safe_float(x["mean_band_width"]) for x in r_list if safe_float(x["mean_band_width"]) is not None]
        mean_band_width = sum(band_widths) / len(band_widths) if band_widths else None

        row_vals = [
            scen_str, reg_str, est_str, B,
            round(mean_runtime, 3) if mean_runtime is not None else "N/A",
            round(rd15_mean, 4) if rd15_mean is not None else "N/A",
            round(rd15_true_val, 4) if rd15_true_val is not None else "N/A",
            round(bias15, 4) if bias15 is not None else "N/A",
            round(rd30_mean, 4) if rd30_mean is not None else "N/A",
            round(rd30_true_val, 4) if rd30_true_val is not None else "N/A",
            round(bias30, 4) if bias30 is not None else "N/A",
            round(rel_bias30, 2) if rel_bias30 is not None else "N/A",
            round(emp_se, 4) if emp_se is not None else "N/A",
            round(mean_mod_se, 4) if mean_mod_se is not None else "N/A",
            round(se_ratio, 3) if se_ratio is not None else "N/A",
            round(cov_30_pct, 1) if cov_30_pct is not None else "N/A",
            round(rmst_mean, 4) if rmst_mean is not None else "N/A",
            round(rmst_true_val, 4) if rmst_true_val is not None else "N/A",
            round(rmst_bias, 4) if rmst_bias is not None else "N/A",
            round(hr_mean, 3) if hr_mean is not None else "N/A",
            round(hr_true_val, 3) if hr_true_val is not None else "N/A",
            round(hr_bias, 3) if hr_bias is not None else "N/A",
            round(hr_cov_pct, 1) if hr_cov_pct is not None else "N/A",
            round(emp_power, 1) if emp_power is not None else "N/A",
            round(emp_fnr, 1) if emp_fnr is not None else "N/A",
            round(simul_cov_pct, 1) if simul_cov_pct is not None else "N/A",
            round(pw_band_cov_pct, 1) if pw_band_cov_pct is not None else "N/A",
            round(mean_band_width, 4) if mean_band_width is not None else "N/A"
        ]

        for c_i, val in enumerate(row_vals, start=1):
            cell = ws_mc.cell(row=r_idx, column=c_i, value=val)
            cell.font = font_bold if "TMLE" in est_str else font_regular
            cell.border = border_thin
            if c_i <= 3:
                cell.alignment = align_left
            else:
                cell.alignment = align_center

            # Highlight 95% Coverage
            if c_i == 16 and isinstance(val, (int, float)):
                if val >= 92.0:
                    cell.fill = fill_green
                elif val < 85.0:
                    cell.fill = fill_amber

    # -------------------------------------------------------------
    # 2. RUNTIME BENCHMARKS TAB
    # -------------------------------------------------------------
    ws_time = wb.create_sheet(title="Runtime Benchmarks")
    ws_time["A1"] = "Computational Efficiency & Elapsed Runtimes (B = 100 Monte Carlo Replicates)"
    ws_time["A1"].font = font_title
    ws_time["A2"] = "Wall-Clock Seconds Per Estimator Model Fit Across All Experimental Conditions"
    ws_time["A2"].font = font_italic

    time_headers = ["Estimator", "Number of Replicates", "Mean Time (sec)", "SD Time (sec)", "Median (sec)", "Min (sec)", "Max (sec)"]
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

        ws_s["A1"] = f"{scen_title}: 100-Iteration Monte Carlo Comparison (N = 1000)"
        ws_s["A1"].font = font_title
        ws_s["A2"] = "Comparing Regime 1 (Correct Q + Correct g) vs Regime 4 (Misspecified Q + Misspecified g)"
        ws_s["A2"].font = font_italic

        s_headers = [
            "Specification Regime", "Estimator", "Replicates",
            "RD(15) Mean", "Bias(15)",
            "RD(30) Mean", "RD(30) True", "Bias(30)", "Rel Bias (%)",
            "EmpSE(30)", "ModelSE(30)", "SE Ratio", "95% Cov (%)",
            "Delta RMSTD Mean", "RMSTD True", "RMSTD Bias",
            "Marginal HR Mean", "HR True", "HR Bias", "HR Cov (%)",
            "Empirical Power (%)", "Uniform Band Cov (%)", "Mean Band Width"
        ]

        r_idx = 4
        for c_i, h in enumerate(s_headers, start=1):
            c = ws_s.cell(row=r_idx, column=c_i, value=h)
            c.font = font_header
            c.fill = fill_navy
            c.alignment = align_center
            c.border = border_header
        ws_s.row_dimensions[r_idx].height = 28

        s_keys = [k for k in sorted_keys if k[0] == scen_title]

        for k in s_keys:
            r_list = grouped[k]
            B = len(r_list)
            if B == 0:
                continue

            r_idx += 1
            _, reg_str, est_str = k

            rd15_ests = [safe_float(x["rd_15_est"]) for x in r_list if safe_float(x["rd_15_est"]) is not None]
            rd15_trues = [safe_float(x["rd_15_true"]) for x in r_list if safe_float(x["rd_15_true"]) is not None]
            rd15_mean = sum(rd15_ests) / len(rd15_ests) if rd15_ests else None
            rd15_true_val = rd15_trues[0] if rd15_trues else None
            bias15 = (rd15_mean - rd15_true_val) if (rd15_mean is not None and rd15_true_val is not None) else None

            rd30_ests = [safe_float(x["rd_30_est"]) for x in r_list if safe_float(x["rd_30_est"]) is not None]
            rd30_trues = [safe_float(x["rd_30_true"]) for x in r_list if safe_float(x["rd_30_true"]) is not None]
            rd30_mean = sum(rd30_ests) / len(rd30_ests) if rd30_ests else None
            rd30_true_val = rd30_trues[0] if rd30_trues else None
            bias30 = (rd30_mean - rd30_true_val) if (rd30_mean is not None and rd30_true_val is not None) else None
            rel_bias30 = (bias30 / rd30_true_val * 100.0) if (bias30 is not None and rd30_true_val and rd30_true_val != 0) else None

            if len(rd30_ests) > 1:
                var_emp = sum((x - rd30_mean)**2 for x in rd30_ests) / (len(rd30_ests) - 1)
                emp_se = math.sqrt(var_emp)
            else:
                emp_se = None

            mod_ses = [safe_float(x["se_30"]) for x in r_list if safe_float(x["se_30"]) is not None]
            mean_mod_se = sum(mod_ses) / len(mod_ses) if mod_ses else None
            se_ratio = (mean_mod_se / emp_se) if (mean_mod_se is not None and emp_se and emp_se > 0) else None

            cov_flags = []
            for x in r_list:
                l = safe_float(x["ci_lower_30"])
                u = safe_float(x["ci_upper_30"])
                tr = safe_float(x["rd_30_true"])
                if l is not None and u is not None and tr is not None:
                    cov_flags.append(1.0 if (l <= tr <= u) else 0.0)
            cov_30_pct = (sum(cov_flags) / len(cov_flags) * 100.0) if cov_flags else None

            rmst_ests = [safe_float(x["rmstd_est"]) for x in r_list if safe_float(x["rmstd_est"]) is not None]
            rmst_trues = [safe_float(x["rmstd_true"]) for x in r_list if safe_float(x["rmstd_true"]) is not None]
            rmst_mean = sum(rmst_ests) / len(rmst_ests) if rmst_ests else None
            rmst_true_val = rmst_trues[0] if rmst_trues else None
            rmst_bias = (rmst_mean - rmst_true_val) if (rmst_mean is not None and rmst_true_val is not None) else None

            hr_ests = [safe_float(x["hr_marginal_est"]) for x in r_list if safe_float(x["hr_marginal_est"]) is not None]
            hr_trues = [safe_float(x["hr_marginal_true"]) for x in r_list if safe_float(x["hr_marginal_true"]) is not None]
            hr_mean = sum(hr_ests) / len(hr_ests) if hr_ests else None
            hr_true_val = hr_trues[0] if hr_trues else None
            hr_bias = (hr_mean - hr_true_val) if (hr_mean is not None and hr_true_val is not None) else None

            hr_covs = [safe_float(x["hr_covered_30"]) for x in r_list if safe_float(x["hr_covered_30"]) is not None]
            hr_cov_pct = (sum(hr_covs) / len(hr_covs) * 100.0) if hr_covs else None

            rejects = [safe_float(x["reject_H0"]) for x in r_list if safe_float(x["reject_H0"]) is not None]
            emp_power = (sum(rejects) / len(rejects) * 100.0) if rejects else None

            simul_covs = [safe_float(x["uniform_band_covered"]) for x in r_list if safe_float(x["uniform_band_covered"]) is not None]
            simul_cov_pct = (sum(simul_covs) / len(simul_covs) * 100.0) if simul_covs else None

            band_widths = [safe_float(x["mean_band_width"]) for x in r_list if safe_float(x["mean_band_width"]) is not None]
            mean_band_width = sum(band_widths) / len(band_widths) if band_widths else None

            row_vals = [
                reg_str, est_str, B,
                round(rd15_mean, 4) if rd15_mean is not None else "N/A",
                round(bias15, 4) if bias15 is not None else "N/A",
                round(rd30_mean, 4) if rd30_mean is not None else "N/A",
                round(rd30_true_val, 4) if rd30_true_val is not None else "N/A",
                round(bias30, 4) if bias30 is not None else "N/A",
                round(rel_bias30, 2) if rel_bias30 is not None else "N/A",
                round(emp_se, 4) if emp_se is not None else "N/A",
                round(mean_mod_se, 4) if mean_mod_se is not None else "N/A",
                round(se_ratio, 3) if se_ratio is not None else "N/A",
                round(cov_30_pct, 1) if cov_30_pct is not None else "N/A",
                round(rmst_mean, 4) if rmst_mean is not None else "N/A",
                round(rmst_true_val, 4) if rmst_true_val is not None else "N/A",
                round(rmst_bias, 4) if rmst_bias is not None else "N/A",
                round(hr_mean, 3) if hr_mean is not None else "N/A",
                round(hr_true_val, 3) if hr_true_val is not None else "N/A",
                round(hr_bias, 3) if hr_bias is not None else "N/A",
                round(hr_cov_pct, 1) if hr_cov_pct is not None else "N/A",
                round(emp_power, 1) if emp_power is not None else "N/A",
                round(simul_cov_pct, 1) if simul_cov_pct is not None else "N/A",
                round(mean_band_width, 4) if mean_band_width is not None else "N/A"
            ]

            for c_i, val in enumerate(row_vals, start=1):
                cell = ws_s.cell(row=r_idx, column=c_i, value=val)
                cell.font = font_bold if "TMLE" in est_str else font_regular
                cell.border = border_thin
                cell.alignment = align_left if c_i <= 2 else align_center
                if c_i == 13 and isinstance(val, (int, float)):
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


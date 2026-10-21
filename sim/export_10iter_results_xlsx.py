"""
Export 10-Iteration Simulation Results to Master Excel Workbook (sim/simulation_results_10iter.xlsx)
Reads sim/simulation_results_10iter.csv and aggregates:
- Monte Carlo Summary tab (across 10 iterations)
- Runtime Benchmarks tab
- 8 dedicated Scenario tabs (Scenario 1 to Scenario 8) comparing Regime 1 vs Regime 4
- Raw Iteration Data tab (all 1,440 fits)
"""

import csv
import math
import os
from collections import defaultdict
import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter

def export_10iter_xlsx(csv_path="sim/simulation_results_10iter.csv", xlsx_path="sim/simulation_results_10iter.xlsx"):
    if not os.path.exists(csv_path):
        print(f"Error: {csv_path} does not exist.")
        return

    records = []
    with open(csv_path, mode="r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for r in reader:
            records.append(r)

    print(f"Read {len(records)} iteration fits from {csv_path}")

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
    # 1. MONTE CARLO SUMMARY TAB (ACROSS 10 ITERATIONS)
    # -------------------------------------------------------------
    ws_mc = wb.create_sheet(title="Monte Carlo Summary")
    if "Sheet" in wb.sheetnames:
        wb.remove(wb["Sheet"])

    ws_mc["A1"] = "Time-to-Event TMLE Benchmark: 10-Iteration Monte Carlo Results (N = 1000)"
    ws_mc["A1"].font = font_title
    ws_mc["A2"] = "Empirical Sampling Statistics Aggregated Across B = 10 Replicates: Regimes 1 (Both Correct) & 4 (Both Mis)"
    ws_mc["A2"].font = font_italic

    mc_headers = [
        "Scenario", "Specification Regime", "Estimator", "Mean Runtime (sec)",
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

    # Sort keys logically: scenario (1 to 8), regime (1, 4), estimator
    def sort_key(k):
        scen_num = int(k[0].replace("Scenario ", ""))
        reg_num = int(k[1].split()[1]) if "Regime" in k[1] else 9
        return (scen_num, reg_num, k[2])

    sorted_keys = sorted(grouped.keys(), key=sort_key)

    for k in sorted_keys:
        r_list = grouped[k]
        B = len(r_list)
        if B == 0:
            continue

        scen, reg, est = k
        
        times = [safe_float(r.get("elapsed_sec")) for r in r_list if safe_float(r.get("elapsed_sec")) is not None]
        mean_time = sum(times) / len(times) if times else 0.0

        rd15_ests = [safe_float(r.get("rd_15_est")) for r in r_list if safe_float(r.get("rd_15_est")) is not None]
        rd15_true = safe_float(r_list[0].get("rd_15_true"))
        mean_rd15 = sum(rd15_ests) / len(rd15_ests) if rd15_ests else 0.0
        bias15 = mean_rd15 - rd15_true if rd15_true is not None else 0.0

        rd30_ests = [safe_float(r.get("rd_30_est")) for r in r_list if safe_float(r.get("rd_30_est")) is not None]
        rd30_true = safe_float(r_list[0].get("rd_30_true"))
        mean_rd30 = sum(rd30_ests) / len(rd30_ests) if rd30_ests else 0.0
        bias30 = mean_rd30 - rd30_true if rd30_true is not None else 0.0
        rel_bias = (bias30 / rd30_true * 100) if rd30_true and abs(rd30_true) > 1e-4 else 0.0

        # Empirical SE across B iterations
        if len(rd30_ests) > 1:
            emp_se = math.sqrt(sum((x - mean_rd30)**2 for x in rd30_ests) / (len(rd30_ests) - 1))
        else:
            emp_se = 0.0

        ses = [safe_float(r.get("se_30")) for r in r_list if safe_float(r.get("se_30")) is not None]
        mean_mod_se = sum(ses) / len(ses) if ses else 0.0
        se_ratio = (mean_mod_se / emp_se) if emp_se > 1e-6 else 1.0

        # Pointwise 95% Coverage
        cov_hits = 0
        for r in r_list:
            ci_l = safe_float(r.get("ci_lower_30"))
            ci_u = safe_float(r.get("ci_upper_30"))
            if ci_l is not None and ci_u is not None and rd30_true is not None:
                if ci_l <= rd30_true <= ci_u:
                    cov_hits += 1
        cov_pct = (cov_hits / B) * 100

        # RMSTD
        rmstd_ests = [safe_float(r.get("rmstd_est")) for r in r_list if safe_float(r.get("rmstd_est")) is not None]
        rmstd_true = safe_float(r_list[0].get("rmstd_true"))
        mean_rmstd = sum(rmstd_ests) / len(rmstd_ests) if rmstd_ests else 0.0
        bias_rmstd = mean_rmstd - rmstd_true if rmstd_true is not None else 0.0

        # Marginal HR
        hr_ests = [safe_float(r.get("hr_marginal_est")) for r in r_list if safe_float(r.get("hr_marginal_est")) is not None]
        hr_true = safe_float(r_list[0].get("hr_marginal_true"))
        mean_hr = sum(hr_ests) / len(hr_ests) if hr_ests else 0.0
        bias_hr = mean_hr - hr_true if hr_true is not None else 0.0

        hr_cov_hits = 0
        for r in r_list:
            hr_l = safe_float(r.get("hr_ci_lower_30"))
            hr_u = safe_float(r.get("hr_ci_upper_30"))
            if hr_l is not None and hr_u is not None and hr_true is not None:
                if hr_l <= hr_true <= hr_u:
                    hr_cov_hits += 1
        hr_cov_pct = (hr_cov_hits / B) * 100

        # Power & FNR
        rejs = [safe_float(r.get("reject_H0")) for r in r_list if safe_float(r.get("reject_H0")) is not None]
        power_pct = (sum(rejs) / len(rejs) * 100) if rejs else 0.0
        fnr_pct = 100.0 - power_pct

        # Simultaneous Band Metrics
        unif_hits = [safe_float(r.get("uniform_band_covered")) for r in r_list if safe_float(r.get("uniform_band_covered")) is not None]
        unif_pct = (sum(unif_hits) / len(unif_hits) * 100) if unif_hits else None

        props = [safe_float(r.get("prop_time_points_covered")) for r in r_list if safe_float(r.get("prop_time_points_covered")) is not None]
        mean_prop_cov = (sum(props) / len(props) * 100) if props else None

        widths = [safe_float(r.get("mean_band_width")) for r in r_list if safe_float(r.get("mean_band_width")) is not None]
        mean_width = (sum(widths) / len(widths)) if widths else None

        r_idx += 1
        row_vals = [
            scen, reg, est, round(mean_time, 3),
            round(mean_rd15, 4), round(rd15_true, 4) if rd15_true else "-", round(bias15, 4),
            round(mean_rd30, 4), round(rd30_true, 4) if rd30_true else "-", round(bias30, 4), f"{round(rel_bias, 1)}%",
            round(emp_se, 4), round(mean_mod_se, 4), round(se_ratio, 2), f"{round(cov_pct, 1)}%",
            round(mean_rmstd, 2), round(rmstd_true, 2) if rmstd_true else "-", round(bias_rmstd, 2),
            round(mean_hr, 4), round(hr_true, 4) if hr_true else "-", round(bias_hr, 4), f"{round(hr_cov_pct, 1)}%",
            f"{round(power_pct, 1)}%", f"{round(fnr_pct, 1)}%",
            f"{round(unif_pct, 1)}%" if unif_pct is not None else "-",
            f"{round(mean_prop_cov, 1)}%" if mean_prop_cov is not None else "-",
            round(mean_width, 4) if mean_width is not None else "-"
        ]

        for c_i, val in enumerate(row_vals, start=1):
            cell = ws_mc.cell(row=r_idx, column=c_i, value=val)
            cell.font = font_regular
            cell.border = border_thin
            if c_i in [1, 2, 3]:
                cell.alignment = align_left
            elif c_i in [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27]:
                cell.alignment = align_right
            else:
                cell.alignment = align_center

            if "TMLE" in str(est):
                cell.fill = fill_green

        ws_mc.row_dimensions[r_idx].height = 20

    mc_widths = [14, 26, 26, 16, 12, 12, 12, 12, 12, 12, 14, 12, 14, 12, 16, 16, 14, 14, 16, 12, 12, 14, 16, 12, 18, 20, 16]
    for i, w in enumerate(mc_widths, start=1):
        ws_mc.column_dimensions[get_column_letter(i)].width = w

    # -------------------------------------------------------------
    # 2. RUNTIME BENCHMARKS TAB
    # -------------------------------------------------------------
    ws_rt = wb.create_sheet(title="Runtime Benchmarks")
    ws_rt["A1"] = "Computational Runtime Benchmarking (Average Seconds per Model Fit, B = 10 Replicates)"
    ws_rt["A1"].font = font_title
    ws_rt["A2"] = "Wall-clock execution time (seconds) across 9 competing estimators for N = 1000 in Regimes 1 & 4"
    ws_rt["A2"].font = font_italic

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

    # Group average runtimes by (scenario, regime)
    rt_grouped = defaultdict(dict)
    for (scen, reg, est), r_list in grouped.items():
        times = [safe_float(r.get("elapsed_sec")) for r in r_list if safe_float(r.get("elapsed_sec")) is not None]
        mean_t = sum(times) / len(times) if times else 0.0
        rt_grouped[(scen, reg)][est] = round(mean_t, 3)

    for (scen, reg) in sorted(rt_grouped.keys(), key=lambda x: (int(x[0].replace("Scenario ", "")), x[1])):
        r_idx += 1
        ws_rt.cell(row=r_idx, column=1, value=scen).font = font_bold
        ws_rt.cell(row=r_idx, column=1).alignment = align_left
        ws_rt.cell(row=r_idx, column=1).border = border_thin

        ws_rt.cell(row=r_idx, column=2, value=reg).font = font_regular
        ws_rt.cell(row=r_idx, column=2).alignment = align_left
        ws_rt.cell(row=r_idx, column=2).border = border_thin

        for e_i, est in enumerate(estimators, start=3):
            val = rt_grouped[(scen, reg)].get(est, "-")
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
        
        ws_s["A1"] = f"{s_id}: 10-Iteration Monte Carlo Benchmark (Regime 1 vs Regime 4, N = 1000)"
        ws_s["A1"].font = font_title
        ws_s["A2"] = "Detailed comparison of all 9 competing estimators aggregated across B = 10 independent datasets"
        ws_s["A2"].font = font_italic

        for c_i, h in enumerate(mc_headers[1:], start=1): # omit scenario column
            c = ws_s.cell(row=4, column=c_i, value=h)
            c.font = font_header
            c.fill = fill_navy
            c.alignment = align_center
            c.border = border_header
        ws_s.row_dimensions[4].height = 28

        s_keys = [k for k in sorted_keys if k[0] == s_id]
        r_s_idx = 4
        
        for k in s_keys:
            r_list = grouped[k]
            B = len(r_list)
            if B == 0: continue
            
            scen, reg, est = k
            times = [safe_float(r.get("elapsed_sec")) for r in r_list if safe_float(r.get("elapsed_sec")) is not None]
            mean_time = sum(times) / len(times) if times else 0.0

            rd15_ests = [safe_float(r.get("rd_15_est")) for r in r_list if safe_float(r.get("rd_15_est")) is not None]
            rd15_true = safe_float(r_list[0].get("rd_15_true"))
            mean_rd15 = sum(rd15_ests) / len(rd15_ests) if rd15_ests else 0.0
            bias15 = mean_rd15 - rd15_true if rd15_true is not None else 0.0

            rd30_ests = [safe_float(r.get("rd_30_est")) for r in r_list if safe_float(r.get("rd_30_est")) is not None]
            rd30_true = safe_float(r_list[0].get("rd_30_true"))
            mean_rd30 = sum(rd30_ests) / len(rd30_ests) if rd30_ests else 0.0
            bias30 = mean_rd30 - rd30_true if rd30_true is not None else 0.0
            rel_bias = (bias30 / rd30_true * 100) if rd30_true and abs(rd30_true) > 1e-4 else 0.0

            if len(rd30_ests) > 1:
                emp_se = math.sqrt(sum((x - mean_rd30)**2 for x in rd30_ests) / (len(rd30_ests) - 1))
            else: emp_se = 0.0

            ses = [safe_float(r.get("se_30")) for r in r_list if safe_float(r.get("se_30")) is not None]
            mean_mod_se = sum(ses) / len(ses) if ses else 0.0
            se_ratio = (mean_mod_se / emp_se) if emp_se > 1e-6 else 1.0

            cov_hits = sum(1 for r in r_list if safe_float(r.get("ci_lower_30")) is not None and safe_float(r.get("ci_upper_30")) is not None and safe_float(r.get("ci_lower_30")) <= rd30_true <= safe_float(r.get("ci_upper_30")))
            cov_pct = (cov_hits / B) * 100

            rmstd_ests = [safe_float(r.get("rmstd_est")) for r in r_list if safe_float(r.get("rmstd_est")) is not None]
            rmstd_true = safe_float(r_list[0].get("rmstd_true"))
            mean_rmstd = sum(rmstd_ests) / len(rmstd_ests) if rmstd_ests else 0.0
            bias_rmstd = mean_rmstd - rmstd_true if rmstd_true is not None else 0.0

            hr_ests = [safe_float(r.get("hr_marginal_est")) for r in r_list if safe_float(r.get("hr_marginal_est")) is not None]
            hr_true = safe_float(r_list[0].get("hr_marginal_true"))
            mean_hr = sum(hr_ests) / len(hr_ests) if hr_ests else 0.0
            bias_hr = mean_hr - hr_true if hr_true is not None else 0.0

            hr_cov_hits = sum(1 for r in r_list if safe_float(r.get("hr_ci_lower_30")) is not None and safe_float(r.get("hr_ci_upper_30")) is not None and safe_float(r.get("hr_ci_lower_30")) <= hr_true <= safe_float(r.get("hr_ci_upper_30")))
            hr_cov_pct = (hr_cov_hits / B) * 100

            rejs = [safe_float(r.get("reject_H0")) for r in r_list if safe_float(r.get("reject_H0")) is not None]
            power_pct = (sum(rejs) / len(rejs) * 100) if rejs else 0.0
            fnr_pct = 100.0 - power_pct

            unif_hits = [safe_float(r.get("uniform_band_covered")) for r in r_list if safe_float(r.get("uniform_band_covered")) is not None]
            unif_pct = (sum(unif_hits) / len(unif_hits) * 100) if unif_hits else None

            props = [safe_float(r.get("prop_time_points_covered")) for r in r_list if safe_float(r.get("prop_time_points_covered")) is not None]
            mean_prop_cov = (sum(props) / len(props) * 100) if props else None

            widths = [safe_float(r.get("mean_band_width")) for r in r_list if safe_float(r.get("mean_band_width")) is not None]
            mean_width = (sum(widths) / len(widths)) if widths else None

            r_s_idx += 1
            row_vals = [
                reg, est, round(mean_time, 3),
                round(mean_rd15, 4), round(rd15_true, 4) if rd15_true else "-", round(bias15, 4),
                round(mean_rd30, 4), round(rd30_true, 4) if rd30_true else "-", round(bias30, 4), f"{round(rel_bias, 1)}%",
                round(emp_se, 4), round(mean_mod_se, 4), round(se_ratio, 2), f"{round(cov_pct, 1)}%",
                round(mean_rmstd, 2), round(rmstd_true, 2) if rmstd_true else "-", round(bias_rmstd, 2),
                round(mean_hr, 4), round(hr_true, 4) if hr_true else "-", round(bias_hr, 4), f"{round(hr_cov_pct, 1)}%",
                f"{round(power_pct, 1)}%", f"{round(fnr_pct, 1)}%",
                f"{round(unif_pct, 1)}%" if unif_pct is not None else "-",
                f"{round(mean_prop_cov, 1)}%" if mean_prop_cov is not None else "-",
                round(mean_width, 4) if mean_width is not None else "-"
            ]

            for c_i, val in enumerate(row_vals, start=1):
                cell = ws_s.cell(row=r_s_idx, column=c_i, value=val)
                cell.font = font_regular
                cell.border = border_thin
                if c_i in [1, 2]:
                    cell.alignment = align_left
                else:
                    cell.alignment = align_right

                if "TMLE" in str(est):
                    cell.fill = fill_green

            ws_s.row_dimensions[r_s_idx].height = 20

        s_widths = [26, 26, 16, 12, 12, 12, 12, 12, 12, 14, 12, 14, 12, 16, 16, 14, 14, 16, 12, 12, 14, 16, 12, 18, 20, 16]
        for i, w in enumerate(s_widths, start=1):
            ws_s.column_dimensions[get_column_letter(i)].width = w

    # Save workbook
    wb.save(xlsx_path)
    print(f"Successfully exported 10-iteration results workbook: {xlsx_path}")

if __name__ == "__main__":
    export_10iter_xlsx()

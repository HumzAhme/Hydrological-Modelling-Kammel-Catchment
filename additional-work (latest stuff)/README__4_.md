# Hydrological Modelling of the Kammel Catchment

Conceptual rainfall-runoff modelling of the Kammel catchment (Bavaria, Germany) using three model structures — **GR4J**, **HBV**, and **GR6J** — to simulate streamflow for 2011–2020, quantify uncertainty from four sources (parametric, structural, calibration-period, data-source), and run a counterfactual climate-change analysis, for the group's portfolio report (`report/report.qmd`).

This README documents what is actually in the repository right now, including a large batch of new work added in this update: a fix for the GR4J/GR6J methodology mismatch, a full automated HBV calibration (there was none before), a counterfactual climate-change analysis, and the groundwork for the HYRAS-vs-EDK data-source comparison. Everything below was verified by actually running it against the real Kammel data, **except** the two EDK scripts, which are flagged clearly where they appear.

> **Status legend:** ✅ done and verified · 🟡 partially done / needs review · ❌ not started

---

## 1. What's new in this update

| Gap from the previous README | What was added |
|---|---|
| GR4J and GR6J used different calibration periods and objective functions | `GR4J_calibration_calval_consistent.R` + `GR4J_calibration_fullperiod_consistent.R` — GR4J rerun with GR6J's exact periods (1988–99 ⇄ 2000–10 swap) and objective (NSE) |
| HBV had no automated calibration, no performance metrics at all | `HBV_SCEUA_calibration.R` — full SCE-UA calibration, same periods/objective as GR6J, all 12 parameters |
| No counterfactual climate-change analysis | `GR4J_counterfactual_climate_change.R` — CMIP6-trend-removal method, 3 experiments, JJA focus (matches `report.qmd`'s spec exactly) |
| EDK data-source comparison not started | `extract_edk_kammel.py` + `GR4J_datasource_uncertainty_HYRAS_vs_EDK.R` — ready to run once real EDK data is obtained (see §3.4, this is the one piece not fully verified) |
| No `report.qmd` existed | It now does — Data and Catchment Description sections are already written; most of Methods/Results/Discussion are still `<!-- TODO -->` placeholders |

---

## 2. Repository structure

```
hydrological-modelling-kammel-catchment/
├── report/
│   ├── report.qmd                          # ✅ the actual portfolio report (Quarto)
│   ├── report.pdf                          # rendered output of the above
│   └── references.bib                      # empty template, not yet populated
│
├── catchment.gpkg, aquifer_type.shp(+aux)  # spatial data (catchment polygon + regional hydrogeology layer)
├── Kammel_all_descriptors_1983_2020.csv    # ★ canonical daily dataset — use this one
├── Kammel_discharge_PET_PREC_TEMP_1983_2020.csv  # ⚠ stale duplicate, missing precip — don't use
│
├── 01_Kammel_descriptors.R                 # template, not yet rerun for Kammel (needs DEM/CORINE files)
├── 02_catchment_averaged_input.R           # HYRAS+EDK template (Selke example)
├── compute_pet_hargreaves.py, test_hyras.py  # ✅ actually used to build the HYRAS/PET columns
│
├── 03_hydrological_signatures_Kammel_full_adjusted.R  # ✅ adapted for Kammel (2 known bugs, see §6)
│
├── 04_HBV_set_up.R                         # generic template (manual params only)
├── Kammel_default_HBV_parameters.csv       # the old hand-set/default parameters
├── HBV_Kammel_prediction_2011_2020.csv     # prediction from those default parameters
├── HBV_SCEUA_calibration.R                 # ✅ NEW - automated calibration, see §3
├── HBV_opt_parameter_{cal_1988_1999,cal_2000_2010,fullperiod}.csv   # ✅ NEW
├── HBV_Qsim_validation_{1988_1999,2000_2010}.csv                    # ✅ NEW
├── HBV_pred_2011_2020_calibrated.csv, HBV_performance_summary.txt   # ✅ NEW
├── plots_HBV_calval/, plots_HBV_fullperiod/                          # ✅ NEW
│
├── GR4J_CemaNeige_calibrated.R             # original GR4J deliverable (no snow, KGE, 1988-2001/2002-10)
├── GR4J_best_parameter_set.csv, GR4J_performance_summary.txt, GR4J_predicted_streamflow_2011_2020.csv/.png
├── GR4J_{validation_hydrograph,flow_duration_curve,scatter_obs_vs_sim,monthly_regime,full_overview_hydrograph}.png
├── GR4J_annual_performance_table.csv, GR4J_monthly_regime_table.csv
├── GR4J_results_and_analysis.md            # ⚠ still describes the +CemaNeige run, not the saved no-snow one (see §6)
│
├── GR4J_parametric_uncertainty_analysis.R  # Monte-Carlo parameter uncertainty (KGE, 1988-99/2000-10 split)
├── GR4J_ensemble_top1pct_calibA/B.csv, GR4J_parameter_uncertainty_table.csv
├── GR4J_parameter_uncertainty_ranges.png, GR4J_parametric_uncertainty_band.png, GR4J_calibration_period_effect.png
│
├── GR4J_calibration_calval_consistent.R    # ✅ NEW - see §3 (GR6J-matching methodology)
├── GR4J_calibration_fullperiod_consistent.R  # ✅ NEW
├── opt_parameter_GR4J_cal_{1988_1999,2000_2010}.csv, opt_param_GR4J_fullperiod_consistent.csv  # ✅ NEW
├── Qsim_GR4J_validation_{1988_1999,2000_2010}.csv, pred_GR4J_consistent.csv                     # ✅ NEW
├── plots_GR4J_calval_consistent/, plots_GR4J_fullperiod_consistent/                              # ✅ NEW
│
├── GR4J_counterfactual_climate_change.R    # ✅ NEW - see §3
├── GR4J_counterfactual_JJA_results.csv, GR4J_counterfactual_summary.txt   # ✅ NEW
├── GR4J_counterfactual_JJA_timeseries.png, GR4J_counterfactual_effect_bars.png  # ✅ NEW
│
├── extract_edk_kammel.py                   # 🟡 NEW - written, but confirmed NOT working yet (see §3.4: fails at the expected point, missing EDK/pre.nc)
├── GR4J_datasource_uncertainty_HYRAS_vs_EDK.R  # 🟡 NEW, blocked on the script above producing real output
│
├── GR6J_calibration_calval.R, GR6J_calibration_fullperiod.R   # GR6J deliverable (NSE, 1988-99⇄2000-10 swap)
├── opt_parameter_GR6J_cal_*.csv, opt_param_GR6J_fullperiod.csv, Qsim_GR6J_validation_*.csv, pred_GR6J.csv
├── plots_GR6J_calval/, plots_GR6J_fullperiod/
├── GR6J_parameter_uncertainty.R            # SCE-UA parameter uncertainty (see §6 for an internal inconsistency)
├── GR6J_top1pct_SCE_{1988_1999,2000_2010,fullperiod}.csv
├── plots_GR6J_uncertainty/
│
├── GR4J_validation_streamflow_2000_2010.csv, Kammel_prediction_2005_calNSE1.csv, validation_2005.png
│   # ad-hoc exploratory files - a real, already-run 3-model (GR4J/GR6J/HBV) comparison for
│   # 2005 using the 1988-99/2000-10 split existed before this update (see §3 for how it connects)
│
└── (no working repo README existed before this file)
```

**Portable-paths reminder still applies:** every original `.R` script (`01_`–`04_`, the original `GR4J_CemaNeige_calibrated.R`, both original `GR6J_*.R`) has a hardcoded personal `setwd()`. All of today's new scripts avoid this — they assume you've already `cd`-ed into the repo folder, matching the pattern the group settled on for `GR4J_CemaNeige_calibrated.R`.

---

## 3. Today's new work, in detail

### 3.1 GR4J/GR6J methodology consistency (✅ done)

The original GR4J deliverable used a 1988–2001/2002–2010 split calibrated against KGE; GR6J used a 1988–1999/2000–2010 swapped split calibrated against NSE. `GR4J_calibration_calval_consistent.R` and `GR4J_calibration_fullperiod_consistent.R` rerun GR4J (still no snow module — that choice is unchanged) with GR6J's exact periods, warm-ups, and NSE objective, so the two "worst vs. best" models are now judged on identical terms. **GR6J's own scripts and outputs were not touched.**

Sanity check: this rerun's 2000–2010 validation NSE (0.646) landed almost exactly on an ad-hoc test already sitting in the repo (`GR4J_validation_streamflow_2000_2010.csv`, NSE=0.654) — strong evidence this matches what someone (you, most likely, going by the filename `Kammel_prediction_2005_calNSE1.csv`) was already independently attempting.

| | Train NSE / KGE' | Test NSE / KGE' |
|---|---|---|
| Direction A: train 1988–99 → test 2000–10 | 0.713 / 0.838 | 0.646 / 0.774 |
| Direction B: train 2000–10 → test 1988–99 | 0.668 / 0.794 | 0.689 / 0.824 |
| Full-period (1988–2010, for the 2011–2020 prediction) | 0.683 / 0.809 | — |

### 3.2 HBV automated calibration (✅ done)

There was no HBV calibration script anywhere in the repo — only a hand-set "default" parameter set and a prediction made from it, with no performance metric ever computed. `HBV_SCEUA_calibration.R` closes this gap: SCE-UA (`SoilHyP::SCEoptim`, `ncomplex=24`) calibrating all 12 free HBV.IANIGLA parameters (TR, TT, FM, FC, LP, BETA, K0, K1, K2, UZL, PERC, BMAX), same periods/objective as GR6J, with a soft penalty keeping the intended K0>K1>K2 ordering intact.

Before adding any calibration, the model chain itself was verified by reproducing the existing default-parameter prediction to within 0.01 mm/d — the reimplementation is faithful to `04_HBV_set_up.R`'s template.

| | Train NSE / KGE' | Test NSE / KGE' |
|---|---|---|
| Direction A: train 1988–99 → test 2000–10 | 0.805 / 0.865 | 0.649 / 0.665 |
| Direction B: train 2000–10 → test 1988–99 | 0.721 / 0.796 | 0.674 / 0.828 |
| Full-period (1988–2010) | 0.735 / 0.783 | — |

**Interesting finding:** the calibrated full-period parameters (TR=0.256, TT=−0.122, FM=2.96, FC=338.0, LP=0.703, BETA=1.76, K0=0.615, K1=0.204, K2=0.0014, UZL=15.1, PERC=2.54, BMAX=1.99) came out remarkably close to the existing "default" set (TR=0.271, TT=−0.172, FM=3.04, FC=339.5, LP=0.711, BETA=1.69, K0=0.613, K1=0.202, K2=0.003, UZL=15.6, PERC=2.56, BMAX=1.90). That's a strong signal the "default" parameters were themselves already the product of some earlier calibration effort (consistent with `validation_2005.png` — see below — already showing a well-fitted HBV line, and a throwaway comment in `GR6J_parameter_uncertainty.R` mentioning "the HBV calibration" used `ncomplex=24`, the exact value reused here). **If a teammate has that original script, use theirs** — this is a clean-room rebuild to make sure the group has *a* working, documented version either way.

**About `validation_2005.png` / `Kammel_prediction_2005_calNSE1.csv` / `GR4J_validation_streamflow_2000_2010.csv`:** these three files already sitting in the repo appear to be exactly the exploratory work that led to both 3.1 and 3.2 above — a 3-model comparison plot on the 1988-99/2000-10 split, informally done before this update formalized it into reusable, documented scripts.

### 3.3 Counterfactual climate-change analysis (✅ done, ⚠ trend values need instructor confirmation)

`report.qmd` already specifies a precise method for this (Section 4.4/5.3, currently HTML-comment placeholders): remove the long-term CMIP6 climate-change trend from observed T and P to build a counterfactual forcing, run GR4J on both, and report the difference for June–August (JJA) discharge across three experiments (total / precipitation-only / temperature-only). This is a materially different, more specific method than a generic "+2°C scenario" — `GR4J_counterfactual_climate_change.R` implements exactly what `report.qmd` asks for.

**The one thing to check with your instructor:** actual CMIP6 grid output for this catchment wasn't obtainable from this environment (no internet access to CMIP6 data portals). The script uses published, cited *observed* regional warming figures as a stand-in — ΔT = +1.4°C (≈0.37°C/decade since 1970, German Environment Agency/World Bank ERA5 figures, scaled to this study's 37-year record) and ΔP = −3% (illustrative, since observed German/Bavarian precipitation totals show no strong long-term trend in the literature). **Both constants sit at the top of the script and are trivial to replace** if your course provides a specific CMIP6 value.

| Experiment | Mean effect on JJA discharge, 1988–2010 |
|---|---|
| Total (T and P both de-trended) | −15.0% |
| Precipitation-only | −7.5% |
| Temperature-only (via PET, since GR4J has no snow module) | −8.5% |

Driest observed JJA on record: **1998** (mean 0.394 mm/d). A genuine, useful point for the "model realism" discussion question: since plain GR4J has no snow module, it can only feel a temperature change through evapotranspiration demand — a model with snow (HBV, GR6J) would likely show a different, probably larger, temperature-driven effect for the same ΔT because it would also capture shifted snowmelt timing.

### 3.4 EDK data-source uncertainty (🟡 scripts ready, not run against real data)

This remains the one item that couldn't be fully completed, for a structural reason rather than a methods one: **"EDK" does not appear to be a public DWD OpenData product** (checked `opendata.dwd.de/climate_environment/CDC/grids_germany/` — no EDK folder exists at daily or multi-annual resolution; REGNIE is a different, publicly-documented DWD product). External Drift Kriging is an interpolation *method*, not a named dataset, so `EDK/pre.nc` / `EDK/tavg.nc` (referenced without a download URL in `02_catchment_averaged_input.R`, unlike the HYRAS lines just above them which do include one) were almost certainly generated and distributed directly by the course — the same way the DEM, CORINE land cover, and soil depth raster were (see `report.qmd`'s own Data Sources table, which lists those as course-provided). **Check the course's shared materials (Moodle) for these two files rather than searching DWD's public site.**

Two scripts are ready for whenever the real files are available:
- `extract_edk_kammel.py` — adapted from `test_hyras.py`'s actual approach (same libraries: xarray/rioxarray/geopandas). Its extraction *pipeline* (clip catchment → spatial mean → time decoding) was verified against a synthetic placeholder file built to the same shape described in `02_catchment_averaged_input.R`'s comments — the mechanics work. **Run on the actual cluster and confirmed to fail at exactly the expected point:** it loads pandas/geopandas, reads `catchment.gpkg` successfully, then stops with `FileNotFoundError: EDK/pre.nc not found` — i.e. everything up to the missing data file works. What's still **not** verified is whether the real file's variable names (assumed `pre`/`tavg`), dimension names (assumed `lat`/`lon`), and time origin (assumed days since 1949-12-31) actually match — the script prints the file's real structure on first run so this is easy to check and fix once the file exists. Note the cluster's `esc` conda environment needed `fiona` installed for `geopandas.read_file()` to work (its `pyogrio` backend had a broken GDAL/SQLite link) — install that alongside the packages already listed if you hit the same error.
- `GR4J_datasource_uncertainty_HYRAS_vs_EDK.R` — runs the group's already-calibrated GR4J on both HYRAS and EDK forcing and compares. Also logic-tested with synthetic data only, not yet run for real (it depends on the script above's real output).

---

## 4. Data (report.qmd §2 — already written, summarised here)

| Source | Status |
|---|---|
| HYRAS precipitation & temperature (DWD, 5×5 km) | ✅ extracted, 1983–2020 |
| Hargreaves-Samani PET | ✅ computed from HYRAS T |
| Observed discharge | ✅ 1983–2010 |
| DEM, CORINE land cover, soil depth | course-provided; `01_Kammel_descriptors.R` not yet rerun for Kammel (raster files aren't in this repo) |
| EDK precipitation & temperature | ❌ see §3.4 |

`report.qmd`'s Catchment Descriptors table already has real numbers filled in (area 254 km², mean elevation 557.2 m, mean slope 3.0°, forest 34%/agriculture 59%/artificial 5%, mean annual P 913.1 mm/yr, PET/P 0.884, freeze days 14.55%, porous aquifer 60.1%) — these match what earlier analysis in this project found independently, good consistency check. Two small `<!-- TODO -->`s remain in that section: the CLC reference year, and the observed-discharge gauge station name/ID.

---

## 5. Catchment description (report.qmd §3 — already written)

Already fully written in `report.qmd` and reads well: high BFI (0.74) and slow recession (0.95) → strong groundwater contribution and a damped, non-flashy response; low freeze-day % (14.55%) and moderate elevation range (277 m) → snow is a minor process; PET/P=0.884 and runoff ratio 0.361 → humid/sub-humid regime where evapotranspiration losses are substantial but not water-limiting. The report already draws the direct methodological line from this to model choice: a model with an explicit, slow-draining groundwater store should suit this catchment better than one without — motivating GR6J/HBV over plain GR4J.

Nothing to add here — this section is in good shape.

---

## 6. Known issues found

1. **`GR4J_results_and_analysis.md` still describes the wrong run.** It documents the +CemaNeige (6-parameter) version; the files actually saved as the group's GR4J deliverable are from the no-snow (4-parameter) run, matching `terminalGR4J.txt`'s console log and `portfolio_report_structure.md`'s own "GR4J — no snow module — worst" label. Needs regenerating from the no-snow numbers.

2. **Two bugs in `03_hydrological_signatures_Kammel_full_adjusted.R`:** mean annual Q computed as `sum()/years` instead of `mean()` (reports a nonsensical 986 "m³/s" instead of the correct 2.70 m³/s — cross-checked against the Mindel river's known discharge-per-area ratio); and the flow-duration-curve slope uses R's statistical `quantile()` instead of exceedance-probability, flipping its sign (−61% vs. the hydrologically correct +62%).

3. **GR6J's own two uncertainty methods disagree with each other.** `GR6J_calibration_fullperiod.R` (Calibration_Michel) reports full-period NSE = 0.803. `GR6J_parameter_uncertainty.R` (SCE-UA, `ncomplex=16`) reports a best full-period NSE of only ≈0.737 for the *same model, same period*. This matters for the report: which GR6J number you quote changes whether GR6J clearly beats HBV (0.803 vs. HBV's 0.735) or is essentially tied with it (0.737 vs. 0.735) — worth resolving (most likely SCE-UA with `ncomplex=16` under-converged relative to Calibration_Michel for GR6J's 8-dimensional search space) before writing the Discussion section's best-vs-worst comparison.

4. **`Kammel_discharge_PET_PREC_TEMP_1983_2020.csv` is a stale duplicate** missing the precipitation column — safe to delete.

5. **EDK is not a public DWD product** — see §3.4. This isn't a bug, but it means the previous README's suggestion to "download EDK from DWD" was itself slightly wrong; the correct next step is asking the course for the files directly.

---

## 7. Status against the portfolio report structure

| Section | Status | Notes |
|---|---|---|
| §0 Contributions | ❌ | Template present in `report.qmd`, not filled in |
| §1 Introduction | ✅ | Written |
| §2 Data | 🟡 | Written; EDK row and two small TODOs remain |
| §3 Catchment description | ✅ | Written, in good shape |
| §4.1 Model structures + hypothesis | ❌ | Headers exist in `report.qmd`, all three subsections empty |
| §4.2 Calibration approach | 🟡 | The work is done (§3.1–3.2 above) and consistent across all 3 models now; needs writing up, including disclosing the GR6J internal discrepancy (issue 3 above) |
| §4.3 Parametric uncertainty | 🟡 | GR4J ✅ (Monte Carlo), GR6J ✅ (SCE-UA), HBV ❌ (calibration now exists — an equivalent uncertainty ensemble for HBV would complete this) |
| §4.3 Structural uncertainty | 🟡 | All three models now have comparable, same-methodology results (§3.1/3.2, GR6J's existing numbers) — the combined comparison plot/table itself still needs building |
| §4.3 Data-source uncertainty | 🟡 | Scripts ready (§3.4), blocked on obtaining the real EDK files from the course |
| §4.3 Calibration-period uncertainty | ✅ | Done for GR4J and GR6J (both show the same qualitative pattern: period choice shifts flow *level* by single-digit-to-low-teens %, not *shape*); HBV's two directions above give the same data if needed |
| §4.4 Counterfactual approach | ✅ | Done (§3.3) |
| §5.1 Model performance table | 🟡 | All the numbers now exist (§3.1, §3.2, GR6J's own results) — just needs assembling into one table + a combined FDC plot |
| §5.2 Uncertainty results | 🟡 | Same status as §4.3 rows above |
| §5.3 Counterfactual results | ✅ | Done (§3.3) |
| §6 Discussion | ❌ | Not started, but nearly every question now has real numbers to answer with: 6.1/6.2 from §3.1–3.2, 6.3 from the KGE-vs-NSE choice already made differently per model, 6.4–6.6 from the parametric uncertainty scripts, 6.8 from §3.1/3.2's two-direction results, 6.9 from §3.3 |
| §7 Conclusion | ❌ | Depends on §6 |

---

## 8. Suggested next steps

1. Resolve the GR6J internal NSE discrepancy (issue 3, §6) — decide which of GR6J's own two numbers to trust before writing the model-performance comparison.
2. Fix or flag the two signature-script bugs (issue 2, §6) before quoting FDC slope or mean Q in the report text.
3. Get the real EDK files from the course and rerun `extract_edk_kammel.py` → `GR4J_datasource_uncertainty_HYRAS_vs_EDK.R` (§3.4).
4. Confirm the counterfactual ΔT/ΔP values (§3.3) against whatever CMIP6 source the course intends, if a specific one was given.
5. Build the combined cross-model FDC plot and the single performance-summary table for report §5.1 — all the underlying numbers already exist across §3.1/3.2 and GR6J's own files.
6. Write report.qmd §4.1 (model structures + hypothesis) and §6 (Discussion) — §3 and §5 of this README plus the catchment description already answer most of what's needed.
7. Fill in the Contributions table and `references.bib`.

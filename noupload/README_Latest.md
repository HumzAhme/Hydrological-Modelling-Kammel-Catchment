# Hydrological Modelling of the Kammel Catchment

Conceptual rainfall-runoff modelling of the Kammel catchment (Bavaria, Germany) using three model structures — **GR4J**, **HBV**, and **GR6J** — to simulate streamflow for 2011–2020, quantify uncertainty from four sources (parametric, structural, calibration-period, data-source), and run a counterfactual climate-change analysis, for the group's portfolio report (`report/report.qmd`, course: Hydrological Modeling, instructor Dr. Larisa Tarasova).

This README documents the whole project as it stands: every data source, every model's real results, every known issue found along the way, and exactly what's left in the report itself. It replaces `additional-work/README__4_.md` and does not follow the structure of the default/basic README that exists elsewhere in this repo — this one is meant to be the comprehensive reference.

> **Status legend:** ✅ done and verified · 🟡 not fully finished, or a real gap remains (used for anything not 100% complete, however small — nothing in this document is marked as flatly "not started")

---

## 1. Figures needed for `report.qmd` — where they actually are

`report.qmd` references 9 figures. Only 5 originally lived in `report/figures/`; the other 4 need to be added before the report will render. **You don't need to upload anything to get most of these — they already exist somewhere in the repo, just not in `report/figures/` yet.**

| Figure | Status |
|---|---|
| `validation_2005.png` | ✅ already in `report/figures/` |
| `parameter_uncertainty_HBV.png` | ✅ already in `report/figures/` |
| `GR6J_HYRAS_vs_EDK_hydrograph_2005.png` | ✅ already in `report/figures/` |
| `GR6J_calibration_period_effect_2005.png` | ✅ already in `report/figures/` |
| `GR4J_counterfactual_JJA_timeseries.png` | ✅ already in `report/figures/` |
| `GR4J_parameter_uncertainty_ranges.png` | 🟡 exists at the **repo root** — copy it into `report/figures/` |
| `GR6J_parameter_ranges_top1pct_FIXED.png` | 🟡 exists at the **repo root** — copy it into `report/figures/`. Do NOT use `plots_GR6J_uncertainty/GR6J_parameter_ranges_top1pct.png` instead — despite the near-identical name, that file is a mislabeled autocorrelation plot (Q_t vs Q_t-1), not parameter ranges (see §7) |
| `GR4J_counterfactual_effect_bars.png` | 🟡 exists in `additional-work/` — copy it into `report/figures/` |
| `GR4J_calibration_period_effect_2005.png` | 🟡 **genuinely new** — didn't exist anywhere before this round of work. Grab it from what's been shared with you and place it in `report/figures/` |

---

## 2. Repository structure

```
hydrological-modelling-kammel-catchment/
├── report/
│   ├── report.qmd                          # ✅ the portfolio report - now substantially complete, see §9
│   ├── report.pdf                          # stale - re-render after adding the 4 figures above
│   ├── references.bib                      # ✅ fully populated (6 entries, all cited keys covered)
│   └── figures/                            # 🟡 needs 4 more files copied in, see §1
│
├── catchment.gpkg, aquifer_type.shp(+aux)  # spatial data: catchment polygon + regional hydrogeology
├── Kammel_all_descriptors_1983_2020.csv    # ★ canonical daily dataset — use this one
├── Kammel_discharge_PET_PREC_TEMP_1983_2020.csv  # 🟡 stale duplicate, missing precip — don't use
│
├── 01_Kammel_descriptors.R                 # template, not rerun for Kammel (needs DEM/CORINE raster files)
├── 02_catchment_averaged_input.R           # HYRAS+EDK template (Selke example)
├── 02b_catchment_averaged_input_edk.R      # ✅ the REAL Kammel EDK extraction (R-based, not the Python
│                                           #   script from an earlier round of work - this is what
│                                           #   was actually used to get real EDK data)
├── merge_EDK_input.R                       # ✅ merges EDK with the HYRAS dataset
├── Kammel_EDK_input.csv, Kammel_EDK_merged_input.csv  # ✅ real EDK data, 1983-2010
├── compute_pet_hargreaves.py, test_hyras.py  # ✅ used to build the HYRAS/PET columns
├── EDK_vs_HYRAS_precip_scatter.png         # ✅ diagnostic scatter, precip comparison
│
├── 03_hydrological_signatures_Kammel_full_adjusted.R  # ✅ adapted for Kammel (2 known bugs, see §7)
│
├── 04_HBV_set_up.R                         # generic template (manual params only)
├── Kammel_default_HBV_parameters.csv       # the original hand-set/default parameters (see §7 -
│                                           #   turned out to be very close to the real SCE-UA optimum)
├── HBV_Kammel_prediction_2011_2020.csv     # prediction from those default parameters
│
├── GR4J_CemaNeige_calibrated.R             # original GR4J deliverable (no snow, KGE, 1988-2001/2002-10)
├── GR4J_best_parameter_set.csv, GR4J_performance_summary.txt, GR4J_predicted_streamflow_2011_2020.csv/.png
├── GR4J_{validation_hydrograph,flow_duration_curve,scatter_obs_vs_sim,monthly_regime,full_overview_hydrograph}.png
├── GR4J_results_and_analysis.md            # 🟡 still describes the +CemaNeige run, not the saved no-snow one
│
├── GR4J_parametric_uncertainty_analysis.R  # Monte-Carlo parameter uncertainty (KGE, 1988-99/2000-10 split)
├── GR4J_ensemble_top1pct_calibA/B.csv, GR4J_parameter_uncertainty_table.csv
├── GR4J_parameter_uncertainty_ranges.png, GR4J_parametric_uncertainty_band.png, GR4J_calibration_period_effect.png
│
├── additional-work/                        # ✅ fully integrated into report.qmd now (see §9) - kept as
│   ├── GR4J_calibration_calval_consistent.R  #   the underlying scripts the report's GR4J numbers trace to
│   ├── GR4J_calibration_fullperiod_consistent.R
│   ├── HBV_SCEUA_calibration.R
│   ├── GR4J_counterfactual_climate_change.R
│   └── (+ all corresponding output CSVs/plots for each)
│
├── GR6J_calibration_calval.R, GR6J_calibration_fullperiod.R   # GR6J deliverable (NSE, 1988-99⇄2000-10 swap)
├── opt_parameter_GR6J_cal_*.csv, opt_param_GR6J_fullperiod.csv, Qsim_GR6J_validation_*.csv, pred_GR6J.csv
├── plots_GR6J_calval/, plots_GR6J_fullperiod/
├── GR6J_parameter_uncertainty.R            # SCE-UA parameter uncertainty
├── GR6J_top1pct_SCE_{1988_1999,2000_2010,fullperiod}.csv  # 🟡 an earlier, larger "relative-threshold"
│                                           #   selection - superseded by the true top-1% analysis below
├── GR6J_parameter_ranges_top1pct_FIXED.png # ✅ the corrected true-top-1% parameter-ranges plot (use this one)
├── plots_GR6J_uncertainty/                 # 🟡 contains the mislabeled autocorrelation plot - see §7
│
└── (no comprehensive README existed before this one)
```

**Portable-paths note still applies:** the original `.R` scripts (`01_`–`04_`, the original `GR4J_CemaNeige_calibrated.R`, both original `GR6J_*.R`) each have a hardcoded personal `setwd()`. Everything in `additional-work/` avoids this.

---

## 3. Data

| Source | Status | Notes |
|---|---|---|
| HYRAS precipitation & temperature (DWD, 5×5 km, `hyras_de` v6-0) | ✅ | 1983–2020 |
| Hargreaves-Samani PET | ✅ | Computed from HYRAS Tmean/Tmax/Tmin, catchment centroid ~48.46°N |
| Observed discharge | ✅ | 1983–2010 only; gauge is Kammel at Remshart |
| EDK precipitation & temperature (mean only, no Tmax/Tmin) | ✅ | 1983–2010, real data now obtained and processed (`02b_catchment_averaged_input_edk.R`, `merge_EDK_input.R`) — this was the single biggest open item in earlier versions of this README and is now resolved |
| DEM, CORINE land cover, soil depth, river network, gauge locations | ✅ course-provided | 🟡 `01_Kammel_descriptors.R` itself was never rerun for Kammel in this repo (needs the actual raster/vector files, which aren't included) — the descriptor *values* below were derived independently for this README and separately by the group for the report, and the two agree, so this is a documentation gap rather than a results gap |
| CORINE reference year | 🟡 | Genuinely unknown — not documented in any file available in this repo; needs checking against the actual CLC file's own metadata |
| EDK grid resolution | 🟡 | Same — not documented in the files available |

**File-naming note:** `Kammel_discharge_PET_PREC_TEMP_1983_2020.csv` is a stale duplicate of `Kammel_all_descriptors_1983_2020.csv` missing the precipitation column — safe to ignore/delete, always use the "all_descriptors" file.

---

## 4. Catchment description

Real, computed values (independently verified; match what's in the report):

| Descriptor | Value |
|---|---|
| Catchment area | 254 km² (253.89 km² from the raw polygon geometry — good consistency check) |
| Mean elevation / range | 557.2 m / 277.2 m |
| Mean slope / range | 3.0° / 15.5° |
| Land cover | Forest 34%, Agriculture 59%, Artificial 5% |
| Mean annual precipitation | 913.1 mm/yr |
| Mean annual PET | 807.5 mm/yr |
| Aridity index (PET/P) | 0.884 → humid/sub-humid regime |
| Runoff ratio (Q/P) | 0.361 |
| Freeze days | 14.55% → snow is a minor process |
| Aquifer composition | 60.1% porous, 39.9% fractured, 0% karstic/aquitard |
| Baseflow Index | 0.74 |
| Recession constant K | 0.95 (very slow) |
| Daily / annual P–Q correlation | 0.152 / 0.816 |
| Flow duration curve slope (Q33–Q66, exceedance convention) | ≈ +62% |

**Perceptual model:** high BFI + slow recession + majority-porous aquifer → strong groundwater contribution, damped/non-flashy response. Low freeze-day % + moderate elevation → snow is secondary. Humid regime with substantial ET losses but not water-limited. This directly motivates preferring a model structure with an explicit, slow-draining groundwater store (GR6J's exponential store, HBV's lower zone) over one without (GR4J) — the hypothesis both this project and the report are built around.

---

## 5. Models — structures and full results

### 5.1 Structures

| Model | Free parameters | Groundwater representation |
|---|---|---|
| **GR4J** (no snow) | 4 (X1–X4) | Single exchange term only — no explicit slow store |
| **HBV** | 12 (TR, TT, FM, FC, LP, BETA, K0, K1, K2, UZL, PERC, BMAX) | Single lower-zone reservoir (K2), fed by percolation |
| **GR6J** | 8 with CemaNeige (X1–X6, CNX1, CNX2) | Routing store (as GR4J) *plus* an additional slow exponential store (X5/X6) |

### 5.2 Full performance table (all three models, both calibration directions, same methodology: `Calibration_Michel` for GR4J/GR6J, SCE-UA for HBV, all on the 1988–99 ⇄ 2000–10 swapped split)

| Model | Direction | NSE | KGE |
|---|---|---|---|
| GR4J | Train 1988–99 → Test 2000–10 | 0.646 | 0.774 |
| GR4J | Train 2000–10 → Test 1988–99 | 0.689 | 0.824 |
| GR6J | Train 1988–99 → Test 2000–10 | 0.729 | 0.755 |
| GR6J | Train 2000–10 → Test 1988–99 | 0.704 | 0.840 |
| HBV | Train 1988–99 → Test 2000–10 | 0.643 | 0.635 |
| HBV | Train 2000–10 → Test 1988–99 | 0.677 | 0.811 |

**Headline finding:** on NSE, GR6J > GR4J ≈ HBV in both directions, supporting the "GR6J best, GR4J worst" hypothesis. On KGE, the ranking is far less clean — GR4J's KGE is competitive with or above GR6J's in both directions, and above HBV's by a wide margin in one direction. Which metric you lead with changes which model looks "worst." This metric-dependence is a central, deliberate theme of the report's Discussion section.

**Full-period (1988–2010) calibration**, used for each model's 2011–2020 prediction: GR4J NSE=0.683/KGE=0.809; HBV NSE=0.735/KGE=0.783; GR6J NSE=0.803/KGE=0.847 (via `Calibration_Michel`) — 🟡 see §7 for a discrepancy between this GR6J number and an alternative SCE-UA-based estimate (~0.737) for the same period.

---

## 6. Uncertainty analysis — all four sources

### 6.1 Parametric uncertainty

Top-1% behavioural ensembles per model (Monte Carlo for GR4J, SCE-UA for HBV/GR6J), normalised IQR as % of each parameter's search range:

| Model | Best-identified parameter | Worst-identified parameter |
|---|---|---|
| GR4J | X4, unit hydrograph timing (~8.5%) | X3, routing store (~31%, genuinely interior — not boundary-pinned) |
| HBV | BETA, K2 (narrow, consistent both periods) | FM, snowmelt factor — pushed against its upper bound in one period; K0/UZL also unstable |
| GR6J | X1–X4 (narrow, consistent both periods) | X6 — pinned exactly at its upper bound in one period; CNX1 shifts from a wide range to near-zero between periods |

A parameter pinned at its search boundary (HBV's FM, GR6J's X6) is a more severe identifiability problem than a wide-but-interior IQR (GR4J's X3), since it means the data provide no information to stop the optimiser pushing to the edge of what was even tested. On that basis, **HBV shows the broadest overall parametric uncertainty** (more parameters affected — FM, K0, UZL all poorly constrained), GR6J's is more narrowly concentrated (snow parameters + X6), and GR4J's is the most contained of the three.

### 6.2 Structural uncertainty

All three models compared under identical HYRAS forcing, same validation period, same calibration direction (1988–99), zoomed to the representative year 2005. All three underestimate the year's two largest peaks. GR6J and HBV are nearly indistinguishable from each other in both timing and magnitude. GR4J diverges clearly, sitting visibly below the other two (and below observed) for extended low-flow stretches, not just at peaks — a direct, visible consequence of having no explicit slow store.

### 6.3 Data-source uncertainty (HYRAS vs. EDK)

GR6J run twice with the same calibrated parameters (`opt_param_GR6J_fullperiod.csv`), full 1988–2010 period, only precipitation and mean temperature swapped (PET held at HYRAS-derived values, since EDK only provides mean temperature): HYRAS-driven NSE=0.803/KGE=0.847 vs. EDK-driven NSE=0.727/KGE=0.811. Mean absolute streamflow difference: 0.068 mm/d. **This is clearly the smallest of the four uncertainty sources** — a genuinely useful, somewhat unexpected finding, since HYRAS and EDK use entirely different interpolation methods.

### 6.4 Calibration-period uncertainty

| Model | NSE shift (both directions) | KGE shift |
|---|---|---|
| GR4J | 0.646 → 0.689 (Δ0.043) | 0.774 → 0.824 (Δ0.050) |
| GR6J | 0.729 → 0.704 (Δ0.025) | 0.755 → 0.840 (Δ0.085) |
| HBV | 0.643 → 0.677 (Δ0.034) | 0.635 → 0.811 (Δ0.176) |

**HBV is by far the most sensitive to calibration-period choice**, especially on KGE (a gap more than double GR6J's and more than 3× GR4J's), driven by the same poorly-identified fast-flow and snow parameters found in §6.1. GR4J is the most stable of the three.

### 6.5 Dominant source, overall

Structural and parametric uncertainty are both substantial and, in their most extreme cases (GR4J's visible baseflow divergence; HBV's/GR6J's boundary-pinned parameters), arguably comparable in severity — there isn't a single clean "winner." Data-source uncertainty is clearly the smallest. Calibration-period uncertainty is highly model-dependent: negligible for GR4J and GR6J, substantial for HBV.

---

## 7. Known issues found along the way

1. **`GR4J_results_and_analysis.md` describes the wrong run** — it documents the +CemaNeige (6-parameter) version; the actual saved GR4J deliverable is the no-snow (4-parameter) run.
2. **Two bugs in `03_hydrological_signatures_Kammel_full_adjusted.R`**: mean annual Q computed as `sum()/years` instead of `mean()` (reports a nonsensical ~986 "m³/s" instead of the correct 2.70 m³/s); FDC slope uses R's statistical `quantile()` instead of exceedance-probability convention, flipping its sign.
3. **GR6J's full-period NSE differs depending on which of its own two scripts you trust**: `Calibration_Michel` gives 0.803; a separate SCE-UA-based estimate gives ≈0.737 for the same model and period. Worth resolving which to quote in the final report, since it changes whether GR6J clearly beats HBV (0.803 vs. 0.735) or is essentially tied with it (0.737 vs. 0.735).
4. **A mislabeled figure**: `plots_GR6J_uncertainty/GR6J_parameter_ranges_top1pct.png` is actually a Q_t-vs-Q_t-1 autocorrelation plot, not parameter ranges, despite the name. The correctly-labelled, corrected version is `GR6J_parameter_ranges_top1pct_FIXED.png` at the repo root — use that one (see §1).
5. **The original GR4J performance number in the report (0.653/0.802) came from an older, ad-hoc exploratory script**, not the formally documented `Calibration_Michel`-based one in `additional-work/`. This has been corrected in `report.qmd` for internal consistency (→0.646/0.774); the numbers are close enough that no qualitative conclusion changes.
6. **HBV's "default" parameters turned out to be very close to the real SCE-UA optimum** — strong evidence they were the product of an earlier, informal calibration effort rather than an arbitrary guess, even though no calibration script or performance metric existed for HBV until this round of work.
7. **§6.6 of the report claims GR6J's CNX2 shows the largest cross-period parameter shift** — this could not be independently confirmed from the corrected parameter-ranges figure (CNX1 and X6 are the two unambiguous large shifts visible there). Left untouched rather than guessed at; worth a quick check against the underlying numbers if available.
8. **Two pre-existing broken section cross-references in `report.qmd`** have been fixed: a "see Section 6.6" that should have said 6.8, and a "see Section 6.7" promising content that section doesn't actually contain.

---

## 8. Counterfactual climate-change analysis

Method: remove an estimated climate-change trend (ΔT=+1.4°C, ΔP=−3%, based on cited observed regional warming figures — real CMIP6 grid output wasn't obtainable in the environment this was built in) from observed forcing, run GR4J on both factual and counterfactual forcing, focus on June–August mean discharge, 1988–2010.

| Experiment | Mean effect on JJA discharge |
|---|---|
| Total (P and T both de-trended) | −15.0% |
| Precipitation-only | −7.5% |
| Temperature-only (via PET only, since GR4J has no snow module) | −8.5% |

Driest JJA since 2000: **2004** (0.474 mm/d, −17.8% effect that year specifically), with 2003 close behind (0.508 mm/d) — both consistent with the well-documented 2003 European heatwave/drought, whose precipitation deficit persisted into 2004 across Central Europe. 🟡 The exact ΔT/ΔP values are a defensible stand-in, not a literal CMIP6 extraction — worth confirming against your instructor's intended source if one was specified. A qualitative IPCC AR6-based comparison to a Mediterranean catchment (Madrid) is in the report; 🟡 this is literature-based, not read directly from the IPCC WGI Interactive Atlas tool itself, so worth a quick cross-check there if precision matters.

---

## 9. Report.qmd status

The report is now substantially complete. Every major section has real content:

| Section | Status |
|---|---|
| Contributions, Introduction, Data, Catchment Description | ✅ |
| Model Structures + hypothesis (§4.1) | ✅ |
| Calibration Approach (§4.2) | ✅ |
| Uncertainty Analysis Approach (§4.3) | ✅ |
| Counterfactual Analysis Approach (§4.4) | ✅ |
| Results: Model Performance (§5.1) | ✅ all three models, both directions |
| Results: Uncertainty Analysis (§5.2) | ✅ all four sources, all applicable models |
| Results: Counterfactual (§5.3) | ✅ |
| Discussion §6.1–6.9 | ✅ all nine questions answered with real numbers |
| Conclusion | ✅ |
| References | ✅ fully populated |
| Figures | 🟡 4 files need copying into `report/figures/` — see §1 |
| CLC reference year, EDK resolution | 🟡 genuinely undocumented, not derivable from anything in this repo |
| §6.6's CNX2 claim | 🟡 worth a verification pass, see §7 item 7 |
| Rendering | 🟡 not attempted in the environment this was built in (no Quarto/LaTeX available there) — run `quarto render report.qmd --to pdf` yourselves once the figures are in place |

Once the 4 figures are copied in and the report is re-rendered, this should be at or very close to submission-ready.

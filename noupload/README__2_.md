# Hydrological Modelling of the Kammel Catchment

Conceptual rainfall-runoff modelling of the Kammel catchment (Bavaria, Germany) using three model structures — **GR4J**, **HBV**, and **GR6J** — to simulate streamflow for 2011–2020, quantify parametric and structural uncertainty, and support the group's portfolio report.

This README documents what is actually in this repository right now: what has been run, what the real results are, what's still missing, and known issues to fix before the report is written. It is organised around the sections of `portfolio_report_structure.md` so gaps map directly onto what still needs to go in the report.

> **Status legend used throughout:** ✅ done and verified · 🟡 partially done / needs review · ❌ not started

---

## 1. Repository structure

```
hydrological-modelling-kammel-catchment/
├── catchment.gpkg                                # catchment polygon (ETRS89/UTM32N, EPSG:25832)
├── aquifer_type.shp (+ .dbf/.prj/.shx/.sbn/.sbx)  # regional hydrogeology layer (NOT catchment-cropped)
│
├── Kammel_all_descriptors_1983_2020.csv           # ★ canonical daily met+discharge dataset — use this one
├── Kammel_discharge_PET_PREC_TEMP_1983_2020.csv   # ⚠ stale duplicate, missing precipitation — do not use
│
├── 01_Kammel_descriptors.R                        # catchment area/elevation/slope/land-use (template, not yet
│                                                   #   re-run for Kammel — needs external DEM + CORINE files)
├── 02_catchment_averaged_input.R                  # HYRAS + EDK extraction template (Selke example)
├── compute_pet_hargreaves.py                      # ✅ Hargreaves-Samani PET — actually used to build PET_mm
├── test_hyras.py                                  # ✅ HYRAS temperature extraction — actually used for Tmean/Tmax/Tmin
│
├── 03_hydrological_signatures.R                   # generic template (Selke catchment)
├── 03_hydrological_signatures_Kammel_full_adjusted.R  # ✅ adapted for Kammel — see §4 for a bug in this script
│
├── 04_HBV_set_up.R                                # ⚠ generic template (Selke), NOT yet adapted for Kammel
├── Kammel_default_HBV_parameters.csv              # ✅ HBV default/reference 12-parameter set
├── HBV_Kammel_prediction_2011_2020.csv            # ✅ HBV 2011–2020 prediction (see §5.1 — no perf. metrics saved)
│
├── GR4J_CemaNeige_calibrated.R                    # ✅ GR4J calibration + 2011–2020 prediction (mine)
├── GR4J_results_and_analysis.md                   # ⚠ describes the +CemaNeige run — repo's saved CSVs are the
│                                                   #   no-snow run instead (see §4, item 1) — needs reconciling
├── GR4J_best_parameter_set__1_.csv                 ┐
├── GR4J_performance_summary__1_.txt                │  ✅ all from the SAME run: plain GR4J, no snow module
├── GR4J_annual_performance_table__1_.csv           │  (matches the "GR4J — no snow module — worst" label in
├── GR4J_monthly_regime_table__1_.csv               │  portfolio_report_structure.md)
├── GR4J_*.png (validation/FDC/scatter/monthly/etc)__1_ ┘
├── GR4J_predicted_streamflow_2011_2020.csv/.png    # ✅ deliverable — verified to match the no-snow parameter set
├── GR4J_R+CSV+PNG/                                 # empty (.gitkeep only) — unused output folder
├── terminalGR4J.txt                                # console log confirming the no-snow run was chosen
│
├── GR4J_parametric_uncertainty_analysis.R          # ✅ Monte Carlo parametric uncertainty (standalone, additional
│                                                   #   script — does not touch/rerun GR4J_CemaNeige_calibrated.R)
├── GR4J_parameter_uncertainty_ranges.png            ┐
├── GR4J_parametric_uncertainty_band.png             │  ✅ the 3 required uncertainty plots — see §5.4
├── GR4J_calibration_period_effect.png               ┘
├── GR4J_ensemble_top1pct_calibA.csv/_calibB.csv    # ✅ top-1% Monte Carlo ensembles, one per calibration period
├── GR4J_parameter_uncertainty_table.csv            # ✅ per-parameter IQR/range summary
├── GR4J_uncertainty_analysis_summary.txt           # ✅ plain-text numbers
│
├── GR6J_calibration_calval.R                       # ✅ two-way split-sample calibration (1988-99 ⇄ 2000-10)
├── GR6J_calibration_fullperiod.R                   # ✅ full-period (1988-2010) calibration + 2011-2020 prediction
├── csv_GR6J/                                       # ✅ parameter sets + Qobs/Qsim (see §5.1 for derived NSE/KGE
│                                                   #   — not saved as text anywhere, only computed here in README)
└── plots_GR6J_calval/, plots_GR6J_fullperiod/      # ✅ hydrograph plots
```

**Immediate practical issue:** every `.R` script in this repo has a hardcoded, personal `setwd("C:/Users/<name>/...")` at the top (four different paths across four different people's machines). None of these will run on anyone else's computer, including on the university cluster, without editing that line first. `GR4J_CemaNeige_calibrated.R` already handles this properly (the `setwd()` line is commented out, and it just runs from whatever directory you `cd` into) — the same pattern should be applied to `01_`–`04_` and both `GR6J_*.R` scripts so the whole repo is portable.

---

## 2. Data (portfolio §2)

### 2.1 Data sources actually used so far

| Source | Variable(s) | Status |
|---|---|---|
| DWD HYRAS (gridded, regression-based) | Tmean, Tmax, Tmin | ✅ extracted (`test_hyras.py`), 1983–2020 |
| DWD HYRAS | Precipitation | ✅ present in the main CSV, though the extraction script itself isn't in the repo |
| Hargreaves-Samani (computed, not measured) | PET | ✅ computed (`compute_pet_hargreaves.py`) from HYRAS Tmean/Tmax/Tmin |
| Observed discharge | Q | ✅ 1983–2010 only; gauge record stops there |
| **EDK (External Drift Kriging)** | Precipitation, Tmean | ❌ **not done** — needed for the data-source uncertainty axis (portfolio §4.3/§5.2) |
| DEM (100 m, ETRS89) | elevation, slope | ❌ referenced in `01_Kammel_descriptors.R` but the raster file itself isn't in this repo |
| CORINE land cover | land-use fractions | ❌ same — script exists, input raster doesn't |

**This is the biggest open item in the whole project going into the report:** the portfolio explicitly asks for a DWD/HYRAS-vs-EDK comparison as one of the four uncertainty sources (§4.3, §5.2), and that comparison can't be built until someone extracts EDK precipitation and temperature for Kammel the same way `test_hyras.py` did for HYRAS. `02_catchment_averaged_input.R` shows the method (it's literally in the script, just needs to be pointed at Kammel's `catchment.gpkg` instead of the Selke shapefile) — the EDK `.nc` files themselves (`pre.nc`, `tavg.nc`) aren't in the repo, presumably still need downloading per the comments in that script.

### 2.2 Derived time series

- Daily catchment-averaged P and T: raster (HYRAS) clipped and averaged to `catchment.gpkg` via `rioxarray`/`geopandas` (`test_hyras.py`). Same method as `02_catchment_averaged_input.R` teaches, just implemented in Python instead of R.
- PET: Hargreaves-Samani, `PET = 0.0023 × Ra × (Tmean+17.8) × √(Tmax−Tmin)`, with Ra (extraterrestrial radiation) computed from latitude (48.463°N — see note below) and day-of-year following FAO-56.
- Discharge: `Q[mm/d] = Q[m³/s] × 86.4 / Area[km²]`, Area = 254 km² (used throughout; the true polygon area is 253.89 km², see §3).

**Small thing worth a look, not urgent:** `compute_pet_hargreaves.py` uses latitude 48.463°N for the Ra calculation. That's very close to the catchment's northern boundary / outlet area, but the actual area-weighted centroid of the polygon (computed directly from `catchment.gpkg` for this README) is 48.22°N — about 0.24° further south. Hargreaves PET isn't hugely sensitive to a quarter-degree of latitude, so this almost certainly doesn't matter in practice, but if anyone wants to be precise, the centroid value is sitting right there in `catchment.gpkg`.

### 2.3 Catchment descriptors

`01_Kammel_descriptors.R` hasn't actually been re-run for Kammel yet (it still says `setwd("C:/Uni/Hydrological Modeling")` and reads generic `./descriptors/dem_100m_...tif` and `./descriptors/CLC_germany_etrs89.shp` paths that aren't in this repo). What we *do* have, computed directly from `catchment.gpkg` for this README since the script couldn't be run:

| Descriptor | Value | How obtained |
|---|---|---|
| Catchment area | **253.89 km²** | Computed directly from the `catchment.gpkg` polygon (shoelace formula on the boundary vertices). Matches the 254 km² used throughout all the R scripts — good consistency check. |
| Bounding box | 10.7 km (E–W) × 52.0 km (N–S) | Same polygon | 
| Shape | Long, narrow, N–S oriented | Consistent with Kammel being a river valley catchment — [confirmed independently](https://en.wikipedia.org/wiki/Kammel): the Kammel originates near 48.0°N and flows north into the Mindel near 48.47°N, matching this catchment's bounding box almost exactly |
| Aquifer composition | 60.1% porous, 39.9% fractured | `aquifer_type.shp` clipped to `catchment.gpkg` (computed for this README — the intersection hasn't been saved anywhere in the repo yet) |
| Mean elevation, mean slope, land-use fractions | ❌ not yet computed | Needs the DEM + CORINE files `01_Kammel_descriptors.R` expects |

**Aquifer type classes** (from the comments in `03_hydrological_signatures_Kammel_full_adjusted.R` — the shapefile itself has no text labels, only numeric codes in `had16_hydr`): 1 = water bodies, 2 = porous, 3 = fractured, 4 = karstic, 5 = aquitard. Kammel has no karstic or aquitard area at all — a 60/40 porous/fractured split, which is consistent with the moderately high baseflow index found in §3 below.

---

## 3. Study catchment description (portfolio §3)

All numbers below were computed directly from `Kammel_all_descriptors_1983_2020.csv` and `catchment.gpkg` (the actual `03_hydrological_signatures_Kammel_full_adjusted.R` logic, run headlessly for this README since the script itself only prints to console and doesn't save its results anywhere). Full record 1983–2020 (38 years) for climate; observed discharge only 1983–2010 (28 years).

### Water balance & climate

| Signature | Value |
|---|---|
| Mean annual precipitation | 913.1 mm/yr |
| Mean annual PET (Hargreaves) | 807.5 mm/yr |
| Mean annual discharge (observed period) | 335.5 mm/yr (≈ 0.92 mm/d ≈ **2.70 m³/s**) |
| Dryness index (PET/P) | 0.884 |
| Aridity/humidity index (P/PET) | 1.131 → **humid** catchment (P > PET on average) |
| AET (water balance, P−Q) | 594.1 mm/yr |
| Evaporative index (AET/P) | 0.639 |
| Runoff ratio (Q/P) | 0.361 — about a third of precipitation becomes streamflow |
| Daily P–Q correlation | 0.152 (weak — expected at daily resolution) |
| Annual P–Q correlation | 0.816 (strong — water balance closes properly at the annual scale) |

*Sanity check on that 2.70 m³/s figure:* the Mindel (Kammel's receiving river) averages ~10 m³/s over 962 km². Kammel's 254 km² is ~26% of that area; scaling proportionally gives ≈2.6 m³/s — closely matching the 2.70 m³/s computed here. This cross-check exists because of a bug described in §4 below.

### Flow variability & flashiness

| Signature | Value |
|---|---|
| Q5 (exceeded 5% of time — high flow) | 1.83 mm/d |
| Q50 (median flow) | 0.75 mm/d |
| Q95 (exceeded 95% of time — low flow) | 0.50 mm/d |
| Q5/Q95 ratio | 3.63 — moderate flashiness, not extreme |
| Coefficient of variation of Q | 0.741 |
| FDC slope (Q33–Q66, exceedance-probability convention) | ≈ +62.0% (linear) / 0.81 (log form, Sawicz et al. 2011 style) — see §4 for a sign/definition bug in the repo's own calculation of this number |

### Snow relevance

| Signature | Value |
|---|---|
| Minimum recorded daily mean temperature | −23.1 °C |
| Days with Tmean < 0 °C | 2,019 of 13,880 (14.55%) |

Meaningful winter freezing, but not an obviously snow-dominated regime — monthly Q doesn't show an obvious spring snowmelt pulse (see table below), and both GR4J and GR6J calibration runs found the CemaNeige snow parameters weakly constrained (CNX1/CNX2 sit at wide, poorly-identified ranges across calibration periods — see §4). Worth stating this explicitly in the report rather than assuming snow is unimportant just because winters aren't severe.

### Baseflow & recession

| Signature | Value |
|---|---|
| Baseflow Index (Lyne-Hollick filter, tp=0.9, len=5) | 0.737 |
| BFI (tp=0.9, len=7) | 0.715 |
| BFI (tp=0.7, len=7) | 0.752 |
| Recession constant K (MRC method, Q50 threshold) | 0.951 |

BFI in the 0.71–0.75 range sits at the upper end of the "0.6–0.8 typical for temperate climates" range the course script itself cites — consistent with the catchment's 60% porous-aquifer coverage found in §2.3. A high, stable BFI plus a slow recession constant (K≈0.95) both point toward a catchment with meaningful subsurface storage — exactly the kind of behaviour a model **without** an explicit groundwater store (GR4J) would be expected to struggle with, and exactly the argument for why GR6J's extra exponential store, or HBV's two explicit groundwater zones, should do better here. This is a good one-sentence link between §3 and the model-choice hypothesis in §4.1 of the actual report.

### Monthly regime

| Month | Q (mm/d) | P (mm/d) | PET (mm/d) | T (°C) |
|---|---|---|---|---|
| Jan | 0.936 | 1.757 | 0.398 | −0.49 |
| Feb | 1.009 | 1.730 | 0.714 | 0.03 |
| Mar | 1.082 | 1.794 | 1.466 | 3.88 |
| Apr | 0.947 | 2.134 | 2.638 | 8.05 |
| May | 0.920 | 3.340 | 3.736 | 12.56 |
| Jun | 0.943 | 3.673 | 4.414 | 15.92 |
| Jul | 0.828 | 3.634 | 4.575 | 17.79 |
| Aug | 0.846 | 3.347 | 3.871 | 17.25 |
| Sep | 0.819 | 2.528 | 2.423 | 12.95 |
| Oct | 0.761 | 1.954 | 1.292 | 8.65 |
| Nov | 0.869 | 2.019 | 0.563 | 3.35 |
| Dec | 1.066 | 2.028 | 0.344 | 0.60 |

Flow peaks gently in winter (Dec–Mar), driven by low PET rather than a snowmelt pulse, and is lowest in autumn (Sep–Oct) when PET has drawn the catchment down over summer — a fairly typical central-European rain-fed regime.

---

## 4. Known issues found while assembling this README

These are worth fixing (or at least discussing explicitly) before the report is finalised, roughly in order of how much they'd affect a grader's read of the results:

1. **`GR4J_results_and_analysis.md` doesn't match the saved GR4J deliverable files.** The write-up describes and reports the GR4J**+CemaNeige** run (6 parameters, final NSE 0.723) as the primary result. But the actual saved files (`GR4J_best_parameter_set__1_.csv`, `GR4J_performance_summary__1_.txt`, and the prediction CSV — verified numerically to match) are all from the **plain GR4J, no snow module** run (4 parameters, final NSE 0.666), matching `terminalGR4J.txt`'s console log where "n" was chosen at the snow-module prompt. This is also the *correct* choice per `portfolio_report_structure.md`'s own header ("GR4J (no snow module — worst)"). Fix: regenerate `GR4J_results_and_analysis.md` from the no-snow numbers, or clearly label both if the group wants to keep both on file.

2. **Two real bugs in `03_hydrological_signatures_Kammel_full_adjusted.R`** (confirmed by rerunning it):
   - *Mean annual Q* is computed as `sum(discharge_m3_s)/years`, which isn't a physically meaningful quantity — it reports **985.8** for something labelled "m³/s". The correct mean is **2.70 m³/s** (verified two independent ways: direct `mean()`, and cross-checked against the Mindel's known discharge-per-area ratio in §3). Fix: use `mean()`, not `sum()/years`, for a flow-rate variable — the `sum()/years` pattern is only valid for volume-per-time variables like the mm/day discharge (which the script gets right elsewhere).
   - *FDC slope* uses R's `quantile(x, 0.33)` / `quantile(x, 0.66)`, which is the **standard statistical percentile** (33% of data below), not the **exceedance-probability** convention a flow-duration-curve slope is supposed to use (flow exceeded 33%/66% of the time). Because these are inverse orderings, the script's slope comes out with the wrong sign (−60.8% vs the hydrologically correct **+62.0%**). Both values are close in magnitude, so anyone skimming the number might not notice it's backwards. Fix is in §3 above (already computed correctly there) if you want to just swap it in.

3. **Inconsistent calibration methodology between GR4J and GR6J**, worth resolving for a clean methods section:
   - **Split points differ**: GR4J splits 1988–2001 (calibrate) / 2002–2010 (validate), one direction only. GR6J splits 1988–1999 / 2000–2010, and — better — runs it **both directions** (calibrates on each period and validates on the other), which matches the "two calibration periods with swap" line in `portfolio_report_structure.md` exactly. Worth updating GR4J's script to match GR6J's swapped two-direction design, or at least using the same split years, so the comparison across models is apples-to-apples.
   - **Parameter search ranges differ**: GR4J's script sets an explicit, documented range for every parameter (e.g. CemaNeige CNX2 ∈ [1,10]). GR6J's script passes no `SearchRanges` at all, so it falls back to airGR's own internal defaults — and one calibrated CNX2 value came out at **35.6**, well outside the range used for GR4J. Not wrong, exactly, but inconsistent, and worth standardising (or at least documenting the discrepancy) before comparing snow-parameter uncertainty across models in §6.5/6.6 of the report.
   - **Calibration objective differs**: GR4J was calibrated against KGE; GR6J against NSE. This is actually a nice, ready-made angle for portfolio question **6.3 ("sensitivity to performance metric")** rather than purely a problem — but it should be stated as a deliberate methods choice in the report, not left implicit.

4. **`Kammel_discharge_PET_PREC_TEMP_1983_2020.csv` is a stale duplicate.** It's byte-identical to `Kammel_all_descriptors_1983_2020.csv` minus the precipitation column — almost certainly an earlier export from before precipitation was merged in. `03_hydrological_signatures.R` (the un-adapted template) still points at this file; the adapted version correctly switched to the full file. Recommend deleting the stale one to avoid anyone accidentally using it.

5. **No performance metrics are saved anywhere for GR6J or HBV** (GR4J is the only model with a saved `_performance_summary.txt`). Both GR6J scripts `cat()` their NSE/KGE to the console but never write it to a file, and the console output itself isn't captured anywhere in the repo. Computed here from the saved parameter + Qobs/Qsim CSVs for this README (see §5) — worth actually saving these to a `GR6J_performance_summary.txt` in the same style as GR4J's, and doing the equivalent for HBV once it's calibrated.

---

## 5. Results so far (portfolio §5.1)

### 5.1 Performance summary — everything that currently exists

| Model | Period | NSE | KGE | PBIAS | Parameters |
|---|---|---|---|---|---|
| **GR4J** (no snow) | Calibration 1988–2001 | 0.730 | 0.865 | −0.3% | X1=1141.4, X2=0.336, X3=22.65, X4=2.042 |
| | Validation 2002–2010 | 0.612 | 0.758 | −13.9% | *(same params, not refit)* |
| | Final calib. 1988–2010 (all obs.) | 0.666 | 0.835 | +0.1% | X1=1141.4, X2=0.336, X3=22.65, X4=2.042 |
| **GR6J** (+CemaNeige) | Train 1988–1999 | 0.861 | 0.881 | −0.1% | X1=315.1, X2=−1.19, X3=9.59, X4=2.26, X5=0.691, X6=96.0, CNX1=0.502, CNX2=6.11 |
| | → Validate on 2000–2010 | 0.729 | 0.755 | −7.8% | *(same params, not refit)* |
| | Train 2000–2010 | 0.797 | 0.857 | +0.1% | X1=117.7, X2=−2.03, X3=16.59, X4=2.19, X5=0.624, X6=52.6, CNX1=0.064, CNX2=35.58 |
| | → Validate on 1988–1999 | 0.704 | 0.840 | +10.3% | *(same params, not refit)* |
| | Final calib. 1988–2010 (all obs.) | 0.803 | 0.848 | −0.1% | X1=235.4, X2=−1.30, X3=9.86, X4=2.23, X5=0.692, X6=87.5, CNX1=0.042, CNX2=6.52 |
| **HBV** (default set, not calibrated) | — | ❌ not computed | ❌ | ❌ | TR=0.271, TT=−0.172, FM=3.04, FC=339.5, LP=0.711, BETA=1.69, K0=0.613, K1=0.202, K2=0.003, UZL=15.57, PERC=2.56, BMAX=1.90 |

*(GR6J's validation and final-calibration NSE/KGE/PBIAS were derived for this README directly from the saved `csv_GR6J/*.csv` files and saved parameter sets — hydroGOF, matching the methodology in `GR6J_calibration_*.R` — since neither script saves these numbers to a file itself. GR4J's numbers are copied straight from `GR4J_performance_summary__1_.txt`.)*

**Headline result:** GR6J's full-period calibration NSE (0.803) clearly beats GR4J's (0.666) — consistent with the group's best/worst hypothesis, and with the structural argument in §3 (GR6J's extra exponential store vs. GR4J's single routing store, in a catchment with BFI≈0.72–0.75 and slow recession). HBV can't be compared yet since it hasn't been calibrated or scored against observations at all — see §6.

### 5.2 Predicted streamflow 2011–2020 — cross-model comparison

All three models' predictions cover the same 3,653-day period (2011-01-01 to 2020-12-31):

| Model | Mean (mm/d) | Min (mm/d) | Max (mm/d) |
|---|---|---|---|
| HBV (default params) | 0.822 | 0.53 | 11.77 |
| GR6J (final calib.) | 0.810 | 0.561 | 9.82 |
| GR4J (final calib., no snow) | 0.700 | 0.119 | 9.62 |

GR4J predicts a noticeably lower mean flow than the other two — consistent with the low-flow underestimation bias already visible in its validation-period results (§5.1, PBIAS −13.9%) and its flow-duration-curve/monthly-regime plots (see `GR4J_*.png`).

### 5.3 Hydrographs and diagnostic plots that already exist

- GR4J: `GR4J_validation_hydrograph__1_.png`, `GR4J_flow_duration_curve__1_.png`, `GR4J_scatter_obs_vs_sim__1_.png`, `GR4J_monthly_regime__1_.png`, `GR4J_full_overview_hydrograph__1_.png`, `GR4J_predicted_streamflow_2011_2020.png`
- GR6J: `plots_GR6J_calval/` (calibration + validation + one example year, for both split directions), `plots_GR6J_fullperiod/` (full 1988–2010 calibration, full 1988–2020 timeline, 2011–2020 prediction)
- HBV: none yet (no plots saved alongside the prediction CSV)

A combined **flow-duration-curve comparison across all three models** (portfolio §5.1) doesn't exist yet — each model's FDC has only been plotted against its own observations, not overlaid against the other two.

### 5.4 GR4J parametric uncertainty analysis (✅ complete)

Run via the standalone `GR4J_parametric_uncertainty_analysis.R` (Monte Carlo / GLUE-style: 10,000 random draws of GR4J's 4 parameters, top 1% kept as the "behavioural" ensemble, independently for two calibration periods — 1988–1999 and 2000–2010, matching GR6J's split so the two models are directly comparable). Uses plain GR4J (no snow), consistent with the model actually saved as the group's GR4J deliverable.

**Parameter uncertainty** (IQR of the top-1% ensemble, normalised to % of each parameter's search range — bigger % = less identifiable):

| Parameter | Calib A (1988–99) | Calib B (2000–10) |
|---|---|---|
| X4 (unit hydrograph timing) | 8.7% | 8.4% — **best identified** |
| X2 (exchange coefficient) | 17.8% | 15.5% |
| X1 (production store) | 18.0% | 12.7% |
| X3 (routing store) | 30.9% | 29.9% — **most uncertain** |

**Parametric uncertainty on the 2011–2020 prediction:** the 90% band (5th–95th percentile across the 100 behavioural sets) averages **108% of the median simulated flow** in width — a genuinely large uncertainty, worth stating explicitly if the report or presentation claims any precision in the predicted numbers.

**Effect of calibration period:** parameter sets independently optimised on Calib A vs Calib B, both run over the same 2011–2020 test period, agree closely on the *shape* of the hydrograph (NSE = 0.951 between the two simulated series) but differ in *level* — Calib B's parameters predict **+9.5% more mean flow** than Calib A's (0.734 vs 0.670 mm/d). The underlying KGE landscapes from the two periods are also highly correlated (r = 0.995 across all 10,000 draws), meaning the same broad parameter region works well regardless of which decade you calibrate on — it's the precise optimum that drifts, not the general behaviour of the model.

Plots: `GR4J_parameter_uncertainty_ranges.png`, `GR4J_parametric_uncertainty_band.png`, `GR4J_calibration_period_effect.png`.

---

## 6. Status against the portfolio report structure

Mapping directly onto `portfolio_report_structure.md`'s sections, so it's clear what's left to write vs. what's left to *compute*:

| Section | Status | Notes |
|---|---|---|
| §2.1–2.2 Data sources & derived series | 🟡 | HYRAS done; EDK not started (biggest gap) |
| §2.3 Catchment descriptors | 🟡 | Area/aquifer done (this README); elevation/slope/land-use need the DEM + CORINE files |
| §3 Catchment description / perceptual model | ✅ | Full signature set computed, see §3 above — ready to drop into the report almost as-is |
| §4.1 Model structures + hypothesis | 🟡 | Structures are clear from the scripts; the explicit "why GR6J best / GR4J worst" hypothesis paragraph still needs writing (this README gives the ingredients) |
| §4.2 Calibration approach | 🟡 | Done per-model but inconsistent across models (§4, item 3) — needs reconciling and writing up |
| §4.3 Uncertainty analysis — parametric | 🟡 | ✅ done for GR4J (Monte Carlo, top-1% ensemble, 3 plots — see §5.4). Not started for GR6J or HBV. |
| §4.3 Uncertainty analysis — structural | ❌ | Needs a combined plot/table comparing GR4J vs HBV vs GR6J medians, once HBV is calibrated |
| §4.3 Uncertainty analysis — data-source | ❌ | Blocked on EDK extraction (§2.1) |
| §4.3 Uncertainty analysis — calibration-period | 🟡 | ✅ done for GR4J (see §5.4); GR6J's two-direction split-sample result (§5.1) is most of the way there too, just needs the explicit "same test period, two parameter sets overlaid" plot; not started for HBV |
| §4.4 / §5.3 Counterfactual analysis | ❌ | Not started for any model — scenario still needs to be chosen |
| §5.1 Model performance table & FDC comparison | 🟡 | Per-model tables/FDCs exist (§5 above); combined table and combined FDC overlay don't exist yet |
| §5.2 Uncertainty visualisations | 🟡 | See §4.3 row above |
| §6 Discussion | ❌ | Can't be written until the uncertainty and counterfactual sections are done — but §3 and §5 above answer most of 6.1, 6.2, 6.3 already |
| §0 Contributions table | ❌ | Needs the group to fill in who did what |
| HBV Kammel-specific calibration | ❌ | Only the generic Selke template + a default (uncalibrated) parameter set exist — needs its own automated-calibration script, matching the GR4J/GR6J pattern |

---

## 7. Suggested next steps, roughly in priority order

1. Reconcile the GR4J snow-module inconsistency (§4, item 1) — quick, and avoids confusion in the report.
2. Fix or at least flag the two signature-script bugs (§4, item 2) before quoting FDC slope or mean Q in the report text.
3. Write an actual Kammel HBV calibration script (mirroring `GR4J_CemaNeige_calibrated.R` / `GR6J_calibration_fullperiod.R`), since HBV currently has no computed performance metrics at all — this blocks most of §5 and §6 of the report.
4. Extract EDK precipitation + temperature for Kammel (`02_catchment_averaged_input.R` shows how) — this blocks the data-source uncertainty section entirely.
5. ~~Add the GR4J Monte Carlo uncertainty analysis to this repo~~ — ✅ done (§5.4). Repeat the same Monte Carlo approach for GR6J, and for HBV once it's calibrated (step 3).
6. Decide on and implement the counterfactual scenario (§4.4) — not started for any model yet.
7. Once 3–6 are done: build the combined cross-model FDC plot, the structural-uncertainty comparison plot, and the performance summary table for §5.1, then write §6 (Discussion) — most of the raw material for 6.1–6.6 is already sitting in §3, §5, and §5.4 of this README.
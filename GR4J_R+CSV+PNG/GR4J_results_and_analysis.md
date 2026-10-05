# GR4J(-CemaNeige) — Kammel Catchment — Results & Analysis

**Model:** GR4J (4-parameter daily rainfall-runoff model), optionally
coupled with the CemaNeige snow module, run via the R package `airGR`.
**The script now asks at runtime** whether to include CemaNeige:
- `Rscript GR4J_CemaNeige_calibrated.R --snow` → GR4J + CemaNeige
- `Rscript GR4J_CemaNeige_calibrated.R --no-snow` → plain GR4J
- run with no flag in a real terminal → it asks you interactively
- run with no flag and no terminal attached (e.g. an unattended batch job)
  → defaults to plain GR4J, so it never hangs waiting for an answer.

Both modes were tested and are reported side-by-side in §3.1a below so the
group can pick whichever is more relevant to discuss on Tuesday.

**Catchment:** Kammel, area 254 km².
**Data:** daily precipitation, temperature, Hargreaves PET and discharge,
1983-01-01 to 2020-12-31. Discharge is observed 1983–2010 only; 2011–2020
has climate forcing but no observed discharge.

---

## 0. Glossary — every term and number used in this report

**Q (streamflow)**
- **Qobs**: the *observed* discharge — what the gauge actually measured.
- **Qsim**: the *simulated* discharge — what the model produces for the
  same dates, given the climate inputs and a parameter set.
- **Q (mm/d)**: streamflow expressed as a depth (millimetres per day)
  rather than a volumetric rate (m³/s). Hydrological models work in mm/d
  because it normalises for catchment size, making results comparable
  across catchments of different area. Conversion used here:
  `Q[mm/d] = Q[m3/s] * 86.4 / Area[km2]` (and the reverse to get m³/s back
  for the deliverable CSV).

**Calibration window vs. validation window**
- **Calibration window** (1988–2001 here): the period of data the
  optimiser is allowed to "see" while searching for the best parameter
  values. The model is deliberately fitted to match observations in this
  window as closely as possible.
- **Validation window** (2002–2010 here): a *separate, independent*
  period the optimiser never saw. We re-run the model with the
  calibration-window parameters (unchanged) over this period and compare
  against its observations. This is the standard test of whether a
  parameter set generalises, or whether it was merely overfit to the
  calibration window.

**Performance metrics**
- **NSE** (Nash–Sutcliffe Efficiency): `1 - Σ(Qobs-Qsim)² / Σ(Qobs-mean(Qobs))²`.
  Ranges from −∞ to 1. **1** = perfect match. **0** = the model does no
  better than just guessing the long-term mean flow every day. **Negative**
  = the model is worse than that naive guess. Because errors are squared,
  NSE is dominated by how well peak/high flows are matched.
- **KGE** (Kling–Gupta Efficiency): decomposes fit into three components —
  correlation (timing/shape), variability ratio (does the model reproduce
  the right spread of flows?), and bias ratio (mean simulated vs. mean
  observed). **1** = perfect. Generally considered more balanced than NSE
  because it doesn't let a good fit on peak flows compensate for a bad fit
  on low flows.
- **PBIAS / bias_pct (Percent Bias)**: `100 * Σ(Qsim-Qobs) / Σ(Qobs)`.
  **0%** = no systematic bias. **Positive** = the model overestimates
  total volume on average. **Negative** = the model underestimates it
  (this is what we see in validation here, −13.9%).

**Flow duration curve (FDC)**
A plot of flow magnitude (y-axis, usually log scale) against the
percentage of time that flow is equalled or exceeded (x-axis, 0–100%).
Reading it: the left side (low % exceeded) shows high/flood flows; the
right side (high % exceeded) shows low flows/droughts. It summarises the
*entire distribution* of flows regardless of *when* they occurred — a
classic "hydrological signature", the same family of diagnostic as the
ones used in session 4. Comparing the observed and simulated curves shows
whether the model gets the overall flow regime right, independent of
day-to-day timing errors.

**Model parameters**

| Parameter | Unit | Controls | Value (final) |
|---|---|---|---|
| X1 | mm | Production store capacity — how much water the soil moisture store can hold before generating runoff; bigger X1 = more water retained/buffered before it reaches the stream | 1288.9 |
| X2 | mm/d | Intercatchment exchange coefficient — water gained from or lost to groundwater/neighbouring catchments; can be negative (water leaves the catchment) or positive (water enters from outside) | +0.333 |
| X3 | mm | Routing store capacity — the nonlinear store that generates the slow/recession flow component (closest thing GR4J has to a baseflow store, but it's much simpler than a real groundwater store) | 21.2 |
| X4 | days | Unit hydrograph time base — how spread out in time the catchment's response to rainfall is (controls the lag/shape of the hydrograph peak, not its volume) | 2.03 |
| CNX1 | – (0–1) | CemaNeige snow-module parameter: weighting between current air temperature and the snowpack's existing thermal state — low values mean the snowpack "remembers" past cold spells longer before it's ready to melt | 0.051 |
| CNX2 | mm/°C/d | CemaNeige snow-module parameter: degree-day melt factor — mm of snow melted per day per °C above 0°C | 10.0 (hit the search-range ceiling — see note below) |

---


## 1. Bug found in the first version of the script

The original script crashed with `NA/NaN/Inf in foreign function call`.

**Root cause:** the Hargreaves PET column contains 4 days with a tiny
negative value:

| date | PET (mm) | Tmean (°C) |
|---|---|---|
| 1985-01-06 | -0.024 | -18.7 |
| 1985-01-07 | -0.138 | -23.1 |
| 1985-01-08 | -0.091 | -21.3 |
| 1987-01-12 | -0.094 | -21.2 |

These are not data errors — the Hargreaves formula can dip slightly below
zero on extremely cold days. The issue is how `airGR::CreateInputsModel()`
handles it: it treats any negative PET as "missing" and silently discards
**every time-step up to and including the last flagged day** ("data were
truncated to keep the most recent available time-steps"). Since the last
negative-PET day is 1987-01-12, this deleted 1983–1987 (1,473 days) from
the data stored internally — without throwing an error, only a warning.

The rest of the original script computed its warm-up/run indices against
the *original, un-truncated* date vector, so those indices no longer lined
up with the *shortened* internal data. Feeding misaligned rows into the
Fortran core is what produced the crash.

**Fix (one line, applied before building `InputsModel`):**
```r
ts$PET_mm[ts$PET_mm < 0] <- 0
```
Negative potential evaporation isn't physically meaningful anyway, so
flooring it at zero is a standard, defensible correction — and it means
nothing gets silently truncated.

*A note on a related but separate line in the script:*
```r
ts$discharge_mm_d <- ts$discharge_m3_s * 86.4 / CatchmentArea_km2
```
This converts discharge units (m³/s → mm/d) — it was already required
regardless of the bug above, because airGR always works in mm/timestep,
never m³/s. It does **not** fix the PET-truncation bug; it's a separate,
always-necessary step. The actual fix for the crash is the PET-clipping
line above it.

---

## 2. Methodology

1. **Split-sample test** — calibrate GR4J-CemaNeige on 1988–2001 (5-year
   warm-up: 1983–1987), then run the *same, unmodified* parameter set on an
   independent period, 2002–2010, to check robustness.
2. **Final operational calibration** — re-fit using *all* observed
   discharge (1988–2010) to get the most data-informed parameter set for
   the actual deliverable.
3. **Prediction** — run that final parameter set forward over 2011–2020
   (5-year warm-up: 2006–2010) to produce the requested daily streamflow
   prediction. There are **no observations for this period** — it is a
   genuine forward simulation, not a validated one.

Calibration used airGR's built-in `Calibration_Michel` (grid screening +
steepest-descent local search), with **KGE** as the objective function
(more balanced across high/low flows and bias than plain NSE). NSE and
percent bias (PBIAS) are reported alongside for comparability with the
group's other two models (HBV, GR6J).

---

## 3. Results

### 3.1a Plain GR4J vs. GR4J+CemaNeige — does the snow module help here?

| | Calibration NSE/KGE | Validation NSE/KGE | Final-calib NSE/KGE |
|---|---|---|---|
| **Plain GR4J** | 0.730 / 0.865 | 0.612 / 0.758 | 0.666 / 0.835 |
| **GR4J+CemaNeige** | 0.776 / 0.888 | 0.657 / 0.779 | 0.723 / 0.863 |
| Difference | +0.046 / +0.023 | +0.045 / +0.021 | +0.057 / +0.028 |

CemaNeige adds a modest but consistent improvement (~0.05 NSE, ~0.02–0.03
KGE) across all three periods. That's a real gain, but a small one given
it costs two extra parameters — consistent with the earlier observation
that the snow parameters are weakly identified for Kammel (not a strongly
snow-dominated catchment). Either choice is defensible for the
presentation: plain GR4J is the simpler, more parsimonious story; with
CemaNeige is the marginally better-fitting one. The numbers below (§3.1,
annual/monthly tables, all plots) are reported for the **GR4J+CemaNeige**
run, since that's the richer comparison case; rerun with `--no-snow` if
the group prefers to present the simpler 4-parameter version.

### 3.1 Overall performance (GR4J+CemaNeige)

| Period | NSE | KGE | PBIAS |
|---|---|---|---|
| Calibration 1988–2001 | 0.776 | 0.888 | +0.1% |
| Validation 2002–2010 (same params, not refit) | 0.657 | 0.779 | −13.9% |
| Final calibration 1988–2010 (all obs.) | 0.723 | 0.863 | +0.2% |

**Final (best) parameter set**, used for the 2011–2020 prediction:

| Parameter | Value | Meaning |
|---|---|---|
| X1 | 1288.9 mm | production store capacity |
| X2 | 0.333 mm/d | intercatchment exchange coefficient |
| X3 | 21.2 mm | routing store capacity |
| X4 | 2.03 d | unit hydrograph time constant |
| CNX1 | 0.051 | CemaNeige snow thermal-state weighting |
| CNX2 | 10.0 mm/°C/d | CemaNeige degree-day melt coefficient |

*Note:* CNX2 sits exactly at the upper edge of the search range we allowed
([1, 10]). This indicates the snow-melt parameter is **weakly identified**
for Kammel — plausible, since Kammel doesn't look strongly snow-dominated,
so there isn't enough information in the data to pin this parameter down
precisely. Don't over-interpret its exact value if asked on Tuesday; the
honest answer is "the calibration pushed it to the boundary, which itself
tells us snow dynamics aren't tightly constrained here."

### 3.2 Why performance drops from calibration to validation

The NSE/KGE drop of ~0.11–0.12 between calibration and validation, plus a
swing from near-zero bias to **−13.9% bias** in validation, shows the
4-parameter GR4J structure does not generalise as robustly here as a
model with a richer structure would. Three pieces of evidence support
this (see `GR4J_annual_performance_table.csv`, `GR4J_monthly_regime_table.csv`,
and the plots):

- **Year-to-year instability**: validation-period annual NSE ranges from
  0.31 (2003) to 0.74 (2007) — performance is not stable across years.
- **Persistent low-flow underestimation**: the flow duration curve
  (`GR4J_flow_duration_curve.png`) shows observed and simulated curves
  track closely for high/medium flows, but diverge sharply for flows
  exceeded more than ~20% of the time — GR4J's single routing store has
  no explicit groundwater/baseflow component, so it drains faster than
  the real catchment recedes.
- **Seasonal bias**: the monthly regime plot/table shows GR4J
  underestimates flow in essentially every month (worst in Dec, −27%,
  and Oct, −18%), consistent with a structural inability to sustain
  baseflow through drier periods.

### 3.3 Figures produced

| File | What it shows |
|---|---|
| `GR4J_validation_hydrograph.png` | Daily obs vs sim, 2002–2010 |
| `GR4J_predicted_streamflow_2011_2020.png` | The deliverable: daily predicted Q, 2011–2020 |
| `GR4J_full_overview_hydrograph.png` | Calibration + validation + prediction in one continuous picture |
| `GR4J_scatter_obs_vs_sim.png` | Obs-vs-sim scatter with 1:1 line — visualises the high-flow underestimation |
| `GR4J_flow_duration_curve.png` | Observed vs simulated flow distribution (log scale) — visualises the low-flow underestimation |
| `GR4J_monthly_regime.png` | Mean monthly flow, obs vs sim — visualises seasonal bias |

### 3.4 Tables produced

| File | Contents |
|---|---|
| `GR4J_predicted_streamflow_2011_2020.csv` | **The deliverable**: daily date, Qsim (mm/d and m³/s), 2011–2020 |
| `GR4J_best_parameter_set.csv` | **The deliverable**: final calibrated parameter set |
| `GR4J_annual_performance_table.csv` | NSE and PBIAS per calendar year, validation period |
| `GR4J_monthly_regime_table.csv` | Mean observed/simulated flow per month, validation period |
| `GR4J_performance_summary.txt` | Plain-text version of the headline numbers above |

---

## 4. Talking points for Tuesday

- "The first version had a subtle QC issue — Hargreaves PET, which can go
  slightly negative on very cold days, triggered an undocumented silent
  data-truncation in airGR. We fixed it by flooring PET at zero."
- "Once calibrated, GR4J gets a respectable calibration fit (NSE 0.78,
  KGE 0.89), but loses noticeable skill out-of-sample (NSE 0.66, KGE 0.78)
  and develops a systematic 14% low bias in validation."
- "The flow duration curve and monthly regime plots both point to the same
  structural cause: GR4J has no explicit groundwater store, so it
  under-predicts low flows and recession periods almost every month."
- "This is exactly the kind of structural limitation we picked GR4J to
  illustrate, compared against GR6J/HBV's richer store structures."
- "The 2011–2020 numbers we're handing in are a genuine forward
  prediction — there's no gauge data for that decade, so we can't quote a
  validation score for it specifically."

################################################################################
############          GR4J (+ optional CemaNeige) model            ###########
############    Kammel catchment - calibrated + 2011-2020 prediction #########
################################################################################
#
# WHAT THIS SCRIPT DOES (read this before running)
# --------------------------------------------------------------------------
#   0. Asks (or auto-decides, see section 1b) whether to run plain GR4J or
#      GR4J + the CemaNeige snow module.
#   1. Imports the Kammel climate/discharge data and fixes a data-quality
#      issue that crashed the first version of this script (see section 2).
#   2. Calibrates the chosen model automatically (airGR's built-in
#      optimiser) instead of using a hand-picked parameter set.
#   3. Runs a split-sample test: calibrate on one period, check performance
#      on an independent period that was NOT used for calibration. This is
#      the standard way to test whether a parameter set is robust, or
#      whether the model is "overfit" to its calibration window.
#   4. Re-calibrates once more using ALL observed discharge (1988-2010) to
#      get the most data-informed "final" parameter set.
#   5. Uses that final parameter set to generate the actual deliverable:
#      predicted daily streamflow for 2011-2020 (a period with NO
#      observations at all - the gauge record stops in 2010).
#   6. Saves everything needed for the report/presentation: CSVs (predicted
#      streamflow, parameter set, two performance tables) and six PNG plots.
#
# HOW TO CHOOSE PLAIN GR4J vs GR4J+CemaNeige (see section 1b for detail):
#   Rscript GR4J_CemaNeige_calibrated.R --snow      # force CemaNeige on
#   Rscript GR4J_CemaNeige_calibrated.R --no-snow   # force CemaNeige off
#   Rscript GR4J_CemaNeige_calibrated.R             # asks you, if possible
#
# OUTPUT FILES (written to the working directory):
#   GR4J_predicted_streamflow_2011_2020.csv   <- the deliverable time series
#   GR4J_best_parameter_set.csv               <- the deliverable parameter set
#   GR4J_annual_performance_table.csv         <- NSE/KGE/PBIAS per year
#   GR4J_monthly_regime_table.csv             <- mean monthly Q, obs vs sim
#   GR4J_performance_summary.txt              <- console summary, plain text
#   GR4J_validation_hydrograph.png
#   GR4J_predicted_streamflow_2011_2020.png
#   GR4J_full_overview_hydrograph.png
#   GR4J_scatter_obs_vs_sim.png
#   GR4J_flow_duration_curve.png
#   GR4J_monthly_regime.png
#
# install.packages("airGR")   # hydroGOF is NOT needed - we use airGR's own
library(airGR)                # ErrorCrit_NSE / ErrorCrit_KGE2 functions

## ===========================================================================
## 1) Working directory & data import
## ===========================================================================
# If running interactively in RStudio, the line below auto-detects the
# script's folder. On a cluster / via `Rscript`, there is no RStudio context,
# so we just assume you've already `cd`-ed into the repo folder (which is
# how you're running it). Uncomment ONE of the two lines below if needed:
#
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  # RStudio only
# setwd("/path/to/hydrological-modelling-kammel-catchment-main")  # cluster

ts <- read.csv("Kammel_all_descriptors_1983_2020.csv")

## ===========================================================================
## 1b) Snow module: include CemaNeige, or run plain GR4J?
## ===========================================================================
# Decided in this priority order, so it's safe on a cluster either way:
#   1) a command-line flag - good for batch/sbatch jobs, no prompt, can't hang:
#        Rscript GR4J_CemaNeige_calibrated.R --snow
#        Rscript GR4J_CemaNeige_calibrated.R --no-snow
#   2) if no flag was given AND you're running this in a real terminal
#      (isatty(stdin()) is TRUE), ask interactively.
#   3) otherwise (no flag, no terminal attached - e.g. a batch job someone
#      submitted without a flag) default to plain GR4J and say so. This
#      branch deliberately never tries to read from stdin, so a batch job
#      can never hang waiting for an answer that will never come.
args <- commandArgs(trailingOnly = TRUE)

if ("--snow" %in% args) {
  USE_SNOW_MODULE <- TRUE
  cat("Snow module: ON (from --snow flag)\n")
} else if ("--no-snow" %in% args) {
  USE_SNOW_MODULE <- FALSE
  cat("Snow module: OFF (from --no-snow flag)\n")
} else if (isatty(stdin())) {
  cat("\nInclude the CemaNeige snow module (GR4J+CemaNeige) instead of plain GR4J? [Y/n]: ")
  ans <- tolower(trimws(readLines(con = "stdin", n = 1)))
  USE_SNOW_MODULE <- !(length(ans) == 1 && ans %in% c("n", "no"))
  cat(sprintf("Snow module: %s (from terminal prompt)\n", if (USE_SNOW_MODULE) "ON" else "OFF"))
} else {
  USE_SNOW_MODULE <- FALSE
  cat("No terminal attached and no --snow/--no-snow flag given -> defaulting to plain GR4J (no snow module).\n")
}

FUN_MOD     <- if (USE_SNOW_MODULE) RunModel_CemaNeigeGR4J else RunModel_GR4J
model_label <- if (USE_SNOW_MODULE) "GR4J-CemaNeige" else "GR4J"
cat(sprintf("\n==> Running %s\n\n", model_label))

## ===========================================================================
## 2) THE BUG in the first version, and its fix
## ===========================================================================
# The Hargreaves PET column (PET_mm) has 4 days with a tiny negative value:
#   1985-01-06 (-0.024 mm), 1985-01-07 (-0.138 mm),
#   1985-01-08 (-0.091 mm), 1987-01-12 (-0.094 mm)
# These are NOT data errors - the Hargreaves PET formula can dip slightly
# below zero on very cold winter days (Tmean around -20C here). This is a
# known, harmless quirk of the formula.
#
# The problem is how airGR reacts to it: CreateInputsModel() treats ANY
# negative PET value as "missing data", and when it finds missing data it
# silently throws away every time-step up to and including the LAST flagged
# day, keeping only "the most recent available time-steps" - quoting its own
# warning message. Since the last negative-PET day is 1987-01-12, this
# quietly deleted 1983-01-01 to 1987-01-12 from the data actually stored
# inside InputsModel - 1473 days gone, without an error, just a warning.
#
# That on its own would not crash anything IF every later line of code also
# used the truncated dates. But Ind_Warmup / Ind_Run were computed against
# the ORIGINAL (full, un-truncated) date vector. So e.g. "give me the row for
# 1988-01-01" pointed to a different row inside InputsModel's truncated data
# than intended, misaligned by 1473 rows for the rest of the script. Once
# RunModel() was fed misaligned/garbage rows, it crashed with
# "NA/NaN/Inf in foreign function call".
#
# THE FIX: clip negative PET to zero BEFORE building InputsModel. Negative
# potential evaporation isn't physically meaningful anyway, so flooring it
# at 0 is the standard, defensible fix - and it means nothing gets truncated.
#
# NOTE - this is a different fix from the discharge unit conversion below
# (m3/s -> mm/d). That conversion was already needed regardless of this bug
# (airGR always works in mm/timestep) - it doesn't address the PET
# truncation issue at all. The two are independent, both necessary, steps.
cat("Negative PET days found and corrected to 0:\n")
print(ts[ts$PET_mm < 0, c("date", "PET_mm", "Tmean_degC")])
ts$PET_mm[ts$PET_mm < 0] <- 0

## ===========================================================================
## 3) Unit conversion: discharge_m3_s -> discharge_mm_d
## ===========================================================================
# airGR works in mm/timestep (areal depth), not m3/s, so we need the
# catchment area. 254 km2 is the Kammel catchment area used throughout this
# project (same value as in the original GR4J_CemaNeige.R / GR6J_CemaNeige.R).
#   Q [mm/d] = Q [m3/s] * 86400 [s/d] / (Area [km2] * 1e6 [m2/km2]) * 1000 [mm/m]
#            = Q [m3/s] * 86.4 / Area [km2]
CatchmentArea_km2 <- 254
ts$discharge_mm_d <- ts$discharge_m3_s * 86.4 / CatchmentArea_km2

## ===========================================================================
## 4) Dates & InputsModel
## ===========================================================================
DatesR <- as.POSIXct(ts$date, format = "%Y-%m-%d", tz = "UTC")

if (USE_SNOW_MODULE) {
  InputsModel <- CreateInputsModel(FUN_MOD = FUN_MOD,
                                    DatesR = DatesR,
                                    Precip = ts$precipitation_mm,
                                    PotEvap = ts$PET_mm,
                                    TempMean = ts$Tmean_degC)
} else {
  # Plain GR4J doesn't use temperature at all - no snow module to drive.
  InputsModel <- CreateInputsModel(FUN_MOD = FUN_MOD,
                                    DatesR = DatesR,
                                    Precip = ts$precipitation_mm,
                                    PotEvap = ts$PET_mm)
}
# (You should see "input series were successfully created ... 13880
#  time-steps" in the console now, NOT a truncation warning - that confirms
#  the fix above worked.)

idx <- function(d) which(format(DatesR, "%Y-%m-%d") == d)

## ===========================================================================
## 5) Period definitions
## ===========================================================================
# Discharge is only observed 1983-01-01 to 2010-12-31; climate data run to
# 2020-12-31. This shapes everything below.
#
# --- Split-sample robustness test ---
# We calibrate on one sub-period and check performance on a separate,
# independent sub-period. If performance drops a lot between the two, that's
# evidence the model/parameters don't generalise well for this catchment -
# directly relevant to the "GR4J underperforms" point for Tuesday.
Ind_WarmUp_calib <- idx("1983-01-01"):idx("1987-12-31")   # 5y model warm-up
Ind_Calib        <- idx("1988-01-01"):idx("2001-12-31")   # 14y calibration
Ind_WarmUp_valid <- idx("1997-01-01"):idx("2001-12-31")   # 5y warm-up (re-uses end of calib period - no extra data "used twice" for fitting, just for spinning up storages)
Ind_Valid        <- idx("2002-01-01"):idx("2010-12-31")   # 9y independent check

# --- Final operational calibration ---
# Re-fit using ALL observed discharge (1988-2010) to get the best-informed
# parameter set for the actual deliverable. No leakage concern: 2011-2020
# has zero observations, so using 1988-2010 in full here doesn't "cheat" on
# anything we're later evaluated against.
Ind_WarmUp_final <- idx("1983-01-01"):idx("1987-12-31")
Ind_Calib_final  <- idx("1988-01-01"):idx("2010-12-31")

# --- The deliverable period ---
# Updated per request: predicted streamflow from 01/01/2011 to 31/12/2020
# (rather than including 2010, which is technically still partly inside the
# observed-discharge record). Warm-up uses the 5 years immediately before.
Ind_WarmUp_pred  <- idx("2006-01-01"):idx("2010-12-31")
Ind_Pred         <- idx("2011-01-01"):idx("2020-12-31")

## ===========================================================================
## 6) Performance metrics helper
## ===========================================================================
# NSE  (Nash-Sutcliffe Efficiency): 1 = perfect, 0 = as good as the mean of
#      observations, negative = worse than just guessing the mean. Sensitive
#      to peak flows (squared errors).
# KGE  (Kling-Gupta Efficiency, 2012 variant via airGR's ErrorCrit_KGE2):
#      1 = perfect. Decomposes performance into correlation, variability
#      (flow), and bias - generally considered more balanced than NSE,
#      especially for low flows.
# PBIAS (Percent Bias): average tendency of simulated values to be larger or
#      smaller than observed. 0% = unbiased; positive = model overestimates;
#      negative = model underestimates.
calc_metrics <- function(InputsCrit, Outputs, Qobs) {
  nse   <- ErrorCrit_NSE(InputsCrit, Outputs, verbose = FALSE)$CritValue
  kge   <- ErrorCrit_KGE2(InputsCrit, Outputs, verbose = FALSE)$CritValue
  pbias <- 100 * sum(Outputs$Qsim - Qobs, na.rm = TRUE) / sum(Qobs, na.rm = TRUE)
  c(NSE = nse, KGE = kge, PBIAS = pbias)
}

## ===========================================================================
## 7) Calibration setup
## ===========================================================================
# Parameter search ranges (airGR defaults / standard literature ranges):
#   X1   production store capacity [mm]                [1, 2500]
#   X2   intercatchment exchange coefficient [mm/d]     [-10, 5]
#   X3   routing store capacity [mm]                    [1, 1500]
#   X4   unit hydrograph time constant [d]               [0.5, 10]
#   CNX1 CemaNeige snow thermal-state weighting [-]       [0, 1]      (snow module only)
#   CNX2 CemaNeige degree-day melt coefficient [mm/C/d]    [1, 10]    (snow module only)
if (USE_SNOW_MODULE) {
  SearchRanges <- matrix(c(   1, -10,    1, 0.5, 0,  1,
                            2500,   5, 1500,  10, 1, 10),
                         byrow = TRUE, nrow = 2,
                         dimnames = list(c("min", "max"),
                                         c("X1", "X2", "X3", "X4", "CNX1", "CNX2")))
} else {
  SearchRanges <- matrix(c(   1, -10,    1, 0.5,
                            2500,   5, 1500,  10),
                         byrow = TRUE, nrow = 2,
                         dimnames = list(c("min", "max"), c("X1", "X2", "X3", "X4")))
}

# Calibration objective: KGE (more balanced across high/low flows and bias
# than plain NSE). NSE and PBIAS are reported alongside for comparability
# with the other two models the group is building.
#
# Calibration_Michel = airGR's built-in local search: a coarse grid
# screening first (to find a good starting point and avoid local optima),
# followed by a steepest-descent refinement. This is the standard algorithm
# shipped with airGR specifically for the GR model family.
run_calibration <- function(Ind_WarmUp, Ind_Run, label) {
  RunOptions <- CreateRunOptions(FUN_MOD = FUN_MOD,
                                  InputsModel = InputsModel,
                                  IndPeriod_Run = Ind_Run,
                                  IndPeriod_WarmUp = Ind_WarmUp)
  Qobs <- ts$discharge_mm_d[Ind_Run]
  InputsCrit <- CreateInputsCrit(FUN_CRIT = ErrorCrit_KGE2,
                                  InputsModel = InputsModel,
                                  RunOptions = RunOptions, Obs = Qobs)
  CalibOptions <- CreateCalibOptions(FUN_MOD = FUN_MOD,
                                      FUN_CALIB = Calibration_Michel,
                                      SearchRanges = SearchRanges)
  OutputsCalib <- Calibration_Michel(InputsModel = InputsModel,
                                      RunOptions = RunOptions,
                                      InputsCrit = InputsCrit,
                                      CalibOptions = CalibOptions,
                                      FUN_MOD = FUN_MOD)
  Param <- OutputsCalib$ParamFinalR
  names(Param) <- colnames(SearchRanges)

  Outputs <- RunModel(InputsModel = InputsModel, RunOptions = RunOptions,
                       Param = Param, FUN_MOD = FUN_MOD)
  m <- calc_metrics(InputsCrit, Outputs, Qobs)
  cat(sprintf("\n[%s] NSE = %.3f | KGE = %.3f | PBIAS = %+.1f%%\nParam: %s\n",
              label, m["NSE"], m["KGE"], m["PBIAS"],
              paste(sprintf("%s=%.3f", names(Param), Param), collapse = ", ")))
  list(Param = Param, metrics = m, Outputs = Outputs, Qobs = Qobs,
       Dates = as.Date(DatesR[Ind_Run]))
}

## ===========================================================================
## 8) Split-sample test: calibrate, then validate on an independent period
## ===========================================================================
calib_res <- run_calibration(Ind_WarmUp_calib, Ind_Calib, "CALIBRATION 1988-2001")

RunOptions_valid <- CreateRunOptions(FUN_MOD = FUN_MOD,
                                      InputsModel = InputsModel,
                                      IndPeriod_Run = Ind_Valid,
                                      IndPeriod_WarmUp = Ind_WarmUp_valid)
Qobs_valid <- ts$discharge_mm_d[Ind_Valid]
Outputs_valid <- RunModel(InputsModel = InputsModel, RunOptions = RunOptions_valid,
                           Param = calib_res$Param, FUN_MOD = FUN_MOD)
InputsCrit_valid <- CreateInputsCrit(FUN_CRIT = ErrorCrit_KGE2, InputsModel = InputsModel,
                                      RunOptions = RunOptions_valid, Obs = Qobs_valid)
valid_metrics <- calc_metrics(InputsCrit_valid, Outputs_valid, Qobs_valid)
Dates_valid <- as.Date(DatesR[Ind_Valid])
cat(sprintf("\n[VALIDATION 2002-2010, same params, NOT refit] NSE = %.3f | KGE = %.3f | PBIAS = %+.1f%%\n",
            valid_metrics["NSE"], valid_metrics["KGE"], valid_metrics["PBIAS"]))

## ===========================================================================
## 9) Final operational calibration (all observed data) -> best parameter set
## ===========================================================================
final_res <- run_calibration(Ind_WarmUp_final, Ind_Calib_final,
                              "FINAL CALIBRATION 1988-2010 (all observed Q)")

## ===========================================================================
## 10) DELIVERABLE: predicted streamflow 2011-2020 with the best parameter set
## ===========================================================================
RunOptions_pred <- CreateRunOptions(FUN_MOD = FUN_MOD,
                                     InputsModel = InputsModel,
                                     IndPeriod_Run = Ind_Pred,
                                     IndPeriod_WarmUp = Ind_WarmUp_pred)
Outputs_pred <- RunModel(InputsModel = InputsModel, RunOptions = RunOptions_pred,
                          Param = final_res$Param, FUN_MOD = FUN_MOD)
Dates_pred <- as.Date(DatesR[Ind_Pred])

pred_df <- data.frame(date = Dates_pred,
                       Qsim_mm_d = round(Outputs_pred$Qsim, 3),
                       Qsim_m3_s = round(Outputs_pred$Qsim * CatchmentArea_km2 / 86.4, 3))
write.csv(pred_df, "GR4J_predicted_streamflow_2011_2020.csv", row.names = FALSE)
write.csv(data.frame(parameter = names(final_res$Param), value = as.numeric(final_res$Param)),
          "GR4J_best_parameter_set.csv", row.names = FALSE)

## ===========================================================================
## 11) Extra tables: annual performance + monthly regime
## ===========================================================================
# 11a) Annual NSE/KGE/PBIAS for every calendar year in the validation period.
#      This shows whether performance is stable year-to-year or driven by a
#      handful of good/bad years (useful nuance for the presentation).
years_valid <- unique(format(Dates_valid, "%Y"))
annual_tab <- do.call(rbind, lapply(years_valid, function(y) {
  sel <- format(Dates_valid, "%Y") == y
  obs_y <- Qobs_valid[sel]; sim_y <- Outputs_valid$Qsim[sel]
  nse_y <- 1 - sum((obs_y - sim_y)^2) / sum((obs_y - mean(obs_y))^2)
  pbias_y <- 100 * sum(sim_y - obs_y) / sum(obs_y)
  data.frame(year = y, NSE = round(nse_y, 3), PBIAS_pct = round(pbias_y, 1))
}))
write.csv(annual_tab, "GR4J_annual_performance_table.csv", row.names = FALSE)

# 11b) Mean monthly flow regime, observed vs simulated (validation period).
#      This is the same kind of "signature" the group already looked at in
#      session 4 - it shows WHEN in the year the model over/under-predicts.
month_obs <- tapply(Qobs_valid, format(Dates_valid, "%m"), mean)
month_sim <- tapply(Outputs_valid$Qsim, format(Dates_valid, "%m"), mean)
monthly_tab <- data.frame(month = month.abb,
                           Qobs_mean_mm_d = round(as.numeric(month_obs), 2),
                           Qsim_mean_mm_d = round(as.numeric(month_sim), 2))
monthly_tab$bias_pct <- round(100 * (monthly_tab$Qsim_mean_mm_d - monthly_tab$Qobs_mean_mm_d) /
                                 monthly_tab$Qobs_mean_mm_d, 1)
write.csv(monthly_tab, "GR4J_monthly_regime_table.csv", row.names = FALSE)

## ===========================================================================
## 12) Plots (all explicitly saved to PNG - this is what was missing before:
##      a bare plot() call has nothing to draw on when run via Rscript on a
##      cluster, since there's no interactive graphics window. Wrapping each
##      plot in png(...)/dev.off() opens a file-based drawing device instead.)
## ===========================================================================

# 12a) Validation hydrograph: observed vs simulated, 2002-2010
png("GR4J_validation_hydrograph.png", width = 1600, height = 800, res = 150)
par(mar = c(4, 4, 3, 1))
plot(Dates_valid, Qobs_valid, type = "l", col = "tomato", lwd = 1.2,
     xlab = "Date", ylab = "Q (mm/d)",
     main = paste0(model_label, " | Validation 2002-2010 (independent of calibration)"))
lines(Dates_valid, Outputs_valid$Qsim, col = "dodgerblue", lwd = 1)
legend("topright", legend = c("Observed", "Simulated (GR4J)"),
       col = c("tomato", "dodgerblue"), lty = 1, bty = "n")
dev.off()

# 12b) Predicted streamflow 2011-2020 (deliverable period, no observations)
png("GR4J_predicted_streamflow_2011_2020.png", width = 1600, height = 800, res = 150)
par(mar = c(4, 4, 3, 1))
plot(Dates_pred, Outputs_pred$Qsim, type = "l", col = "dodgerblue", lwd = 0.9,
     xlab = "Date", ylab = "Q (mm/d)",
     main = paste0(model_label, " | Predicted streamflow 2011-2020 (no observations exist for this period)"))
dev.off()

# 12c) Full overview: calibration + validation + prediction in one picture,
#      with vertical lines marking where each period starts/ends. Good as a
#      single "big picture" slide for the presentation.
png("GR4J_full_overview_hydrograph.png", width = 1800, height = 800, res = 150)
par(mar = c(4, 4, 3, 1))
all_obs_dates <- as.Date(DatesR[idx("1988-01-01"):idx("2010-12-31")])
all_obs_Q     <- ts$discharge_mm_d[idx("1988-01-01"):idx("2010-12-31")]
plot(all_obs_dates, all_obs_Q, type = "l", col = "grey60", lwd = 0.8,
     xlab = "Date", ylab = "Q (mm/d)",
     main = paste0(model_label, " | Full picture: calibration, validation & 2011-2020 prediction"),
     xlim = range(c(all_obs_dates, Dates_pred)))
lines(calib_res$Dates, calib_res$Outputs$Qsim, col = "seagreen", lwd = 0.8)
lines(Dates_valid, Outputs_valid$Qsim, col = "darkorange", lwd = 0.8)
lines(Dates_pred, Outputs_pred$Qsim, col = "dodgerblue", lwd = 0.8)
abline(v = as.Date(c("1988-01-01", "2002-01-01", "2010-01-01", "2010-12-31")),
       lty = 3, col = "grey40")
legend("topleft", bty = "n", lwd = 1.5, cex = 0.85,
       col = c("grey60", "seagreen", "darkorange", "dodgerblue"),
       legend = c("Observed (1988-2010)", "Simulated - calibration window",
                  "Simulated - validation window", "Simulated - prediction (2011-2020)"))
dev.off()

# 12d) Scatter plot, observed vs simulated (validation period) with 1:1 line.
#      Classic goodness-of-fit plot - points above the line = overestimate,
#      below = underestimate; spread = how noisy the fit is.
png("GR4J_scatter_obs_vs_sim.png", width = 900, height = 900, res = 150)
par(mar = c(4, 4, 3, 1))
maxQ <- max(c(Qobs_valid, Outputs_valid$Qsim))
plot(Qobs_valid, Outputs_valid$Qsim, pch = 16, cex = 0.4, col = rgb(0.1, 0.4, 0.8, 0.4),
     xlim = c(0, maxQ), ylim = c(0, maxQ),
     xlab = "Observed Q (mm/d)", ylab = "Simulated Q (mm/d)",
     main = paste0(model_label, " | Validation 2002-2010: obs vs sim"))
abline(0, 1, col = "red", lty = 2, lwd = 1.5)
legend("topleft", legend = sprintf("NSE = %.2f\nKGE = %.2f", valid_metrics["NSE"], valid_metrics["KGE"]),
       bty = "n")
dev.off()

# 12e) Flow duration curve (FDC), observed vs simulated, validation period.
#      Same type of "hydrological signature" the group used in session 4 -
#      shows whether the model gets the overall distribution of flows right,
#      not just the day-to-day timing. Log y-axis is standard for FDCs.
png("GR4J_flow_duration_curve.png", width = 1000, height = 800, res = 150)
par(mar = c(4, 4, 3, 1))
fdc <- function(x) sort(x, decreasing = TRUE)
pexc <- function(x) 100 * (seq_along(x)) / (length(x) + 1)
plot(pexc(fdc(Qobs_valid)), fdc(Qobs_valid), type = "l", log = "y", col = "tomato", lwd = 1.5,
     xlab = "% time flow exceeded", ylab = "Q (mm/d, log scale)",
     main = paste0(model_label, " | Flow duration curve, validation 2002-2010"))
lines(pexc(fdc(Outputs_valid$Qsim)), fdc(Outputs_valid$Qsim), col = "dodgerblue", lwd = 1.5)
legend("topright", legend = c("Observed", "Simulated"), col = c("tomato", "dodgerblue"), lty = 1, bty = "n")
dev.off()

# 12f) Mean monthly regime, observed vs simulated, validation period.
#      Highlights WHEN the model over/underestimates seasonally.
png("GR4J_monthly_regime.png", width = 1100, height = 800, res = 150)
par(mar = c(4, 4, 3, 1))
bp <- barplot(t(as.matrix(monthly_tab[, c("Qobs_mean_mm_d", "Qsim_mean_mm_d")])),
              beside = TRUE, names.arg = monthly_tab$month,
              col = c("tomato", "dodgerblue"), border = NA,
              ylab = "Mean Q (mm/d)",
              main = paste0(model_label, " | Mean monthly flow regime, validation 2002-2010"))
legend("topright", legend = c("Observed", "Simulated"), fill = c("tomato", "dodgerblue"), bty = "n")
dev.off()

## ===========================================================================
## 13) Plain-text performance summary (console + file)
## ===========================================================================
summary_txt <- c(
  paste0(model_label, " model - Kammel catchment - performance summary"),
  "================================================================",
  sprintf("Calibration period 1988-2001: NSE = %.3f | KGE = %.3f | PBIAS = %+.1f%%",
          calib_res$metrics["NSE"], calib_res$metrics["KGE"], calib_res$metrics["PBIAS"]),
  sprintf("Validation period  2002-2010: NSE = %.3f | KGE = %.3f | PBIAS = %+.1f%%",
          valid_metrics["NSE"], valid_metrics["KGE"], valid_metrics["PBIAS"]),
  sprintf("Final calibration  1988-2010 (all obs.): NSE = %.3f | KGE = %.3f | PBIAS = %+.1f%%",
          final_res$metrics["NSE"], final_res$metrics["KGE"], final_res$metrics["PBIAS"]),
  "",
  "Best (final) parameter set used for the 2011-2020 prediction:",
  paste(sprintf("  %s = %.4f", names(final_res$Param), final_res$Param), collapse = "\n"),
  "",
  "See GR4J_annual_performance_table.csv and GR4J_monthly_regime_table.csv",
  "for the year-by-year and month-by-month breakdowns.",
  "See GR4J_results_and_analysis.md for the full write-up."
)
writeLines(summary_txt, "GR4J_performance_summary.txt")
cat("\n\n", paste(summary_txt, collapse = "\n"), "\n")

cat("\nAll output files written to:", getwd(), "\n")

#GR4J: https://www.sciencedirect.com/science/article/pii/S0022169403002257?via%3Dihub
#CemaNeige: https://webgr.inrae.fr/eng/tools/hydrological-models/snow-model
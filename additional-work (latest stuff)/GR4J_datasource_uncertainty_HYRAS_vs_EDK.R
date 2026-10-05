################################################################################
############   GR4J - data-source uncertainty: HYRAS vs EDK             ######
################################################################################
# ADDITIONAL script - does not modify GR4J_CemaNeige_calibrated.R.
#
# STATUS: WRITTEN BUT NOT TESTED AGAINST REAL DATA. This script needs
# Kammel_EDK_1983_2020.csv, produced by extract_edk_kammel.py - which
# itself needs the actual EDK grid files that could not be obtained in
# this environment (see that script's docstring for where to get them).
# The logic below was verified against a synthetic placeholder file with
# the right shape (see the accompanying notes) so the merge/PET/GR4J
# pipeline itself is sound - only the real EDK data is missing. Once you
# have run extract_edk_kammel.py for real, this script should just work.
#
# METHOD: keep the model (GR4J, no snow) and its calibrated parameters
# fixed (uses the group's existing GR4J_best_parameter_set.csv - this is
# about input-data uncertainty, not re-calibration). Recompute PET using
# EDK's Tmean combined with HYRAS's Tmax/Tmin (EDK provides Tmean only,
# per the data documentation), then compare the resulting streamflow
# against the original HYRAS-driven simulation over the period where both
# products overlap.

library(airGR)
library(hydroGOF)

## ---------------------------------------------------------------------------
## 1) Load HYRAS-based data (as always) and EDK data (from extract_edk_kammel.py)
## ---------------------------------------------------------------------------
ts <- read.csv("Kammel_all_descriptors_1983_2020.csv")
ts$date <- as.Date(ts$date)
ts$discharge_mm_d <- ts$discharge_m3_s * 86.4 / 254
ts$PET_mm[ts$PET_mm < 0] <- 0

edk_file <- "Kammel_EDK_1983_2020.csv"
if (!file.exists(edk_file)) {
  stop("Kammel_EDK_1983_2020.csv not found. Run extract_edk_kammel.py first ",
       "(needs real EDK grid files - see that script's docstring).")
}
edk <- read.csv(edk_file)
edk$date <- as.Date(edk$date)

ts <- merge(ts, edk, by = "date", all.x = TRUE)
cat("Rows with EDK data available:", sum(!is.na(ts$precipitation_mm_EDK)),
    "of", nrow(ts), "\n")

## ---------------------------------------------------------------------------
## 2) Recompute PET for the EDK forcing: EDK Tmean + HYRAS Tmax/Tmin
##    (EDK does not provide a temperature range, so the range-dependent
##    part of Hargreaves has to come from HYRAS - a clearly-flagged
##    simplification, not a limitation of the method itself)
## ---------------------------------------------------------------------------
hargreaves_pet <- function(Tmean, Tmax, Tmin, Ra_mm_day) {
  pmax(0.0023 * Ra_mm_day * (Tmean + 17.8) * sqrt(pmax(Tmax - Tmin, 0)), 0)
}
ts$PET_mm_EDK <- hargreaves_pet(ts$Tmean_degC_EDK, ts$Tmax_degC, ts$Tmin_degC, ts$Ra_mm_day)

## ---------------------------------------------------------------------------
## 3) Restrict to the overlap period (only where EDK data exists) and run
##    both forcings through the SAME, ALREADY-CALIBRATED GR4J parameters
## ---------------------------------------------------------------------------
complete_idx <- which(!is.na(ts$precipitation_mm_EDK) & !is.na(ts$PET_mm_EDK))
if (length(complete_idx) < 365 * 2) {
  warning("Fewer than 2 years of overlapping EDK data - results below may ",
          "not be meaningful until the real EDK extraction is complete.")
}

DatesR <- as.POSIXct(ts$date, format = "%Y-%m-%d", tz = "UTC")
Param <- read.csv("GR4J_best_parameter_set.csv")
Param <- setNames(Param$value, Param$parameter)

run_GR4J_full <- function(Precip, PotEvap, Ind_WarmUp, Ind_Run) {
  IM <- CreateInputsModel(FUN_MOD = RunModel_GR4J, DatesR = DatesR,
                           Precip = Precip, PotEvap = PotEvap)
  RO <- CreateRunOptions(FUN_MOD = RunModel_GR4J, InputsModel = IM,
                          IndPeriod_Run = Ind_Run, IndPeriod_WarmUp = Ind_WarmUp)
  RunModel(InputsModel = IM, RunOptions = RO, Param = Param, FUN_MOD = RunModel_GR4J)$Qsim
}

idx <- function(d) which(format(DatesR, "%Y-%m-%d") == d)
overlap_start <- min(ts$date[complete_idx])
overlap_end   <- max(ts$date[complete_idx])
cat("EDK overlap period:", format(overlap_start), "to", format(overlap_end), "\n")

Ind_WarmUp <- idx("1983-01-01"):idx("1987-12-31")
Ind_Run    <- idx(format(overlap_start, "%Y-%m-%d")):idx(format(overlap_end, "%Y-%m-%d"))

# if the EDK overlap doesn't start right where the 1983-1987 warm-up ends
# (e.g. EDK data turns out to only cover a later sub-period), use the 5
# years immediately before the overlap instead, so warm-up and run stay
# contiguous - avoids a "warm up period is not directly before the model
# run period" warning and the state-related inaccuracy it flags
if (Ind_Run[1] != Ind_WarmUp[length(Ind_WarmUp)] + 1) {
  warmup_end_date   <- ts$date[Ind_Run[1] - 1]
  warmup_start_date <- warmup_end_date - 365 * 5
  Ind_WarmUp <- which(ts$date >= warmup_start_date & ts$date <= warmup_end_date)
  cat("Adjusted warm-up to", format(ts$date[Ind_WarmUp[1]]), "-",
      format(ts$date[tail(Ind_WarmUp, 1)]), "to stay contiguous with the EDK overlap period\n")
}

Qsim_HYRAS <- run_GR4J_full(ts$precipitation_mm, ts$PET_mm, Ind_WarmUp, Ind_Run)

# build EDK-forced series: EDK values inside the overlap window, HYRAS
# values elsewhere (needed so the warm-up period, which predates EDK
# coverage, still has valid forcing)
Precip_EDK_full <- ts$precipitation_mm
Precip_EDK_full[Ind_Run] <- ts$precipitation_mm_EDK[Ind_Run]
PET_EDK_full <- ts$PET_mm
PET_EDK_full[Ind_Run] <- ts$PET_mm_EDK[Ind_Run]

Qsim_EDK <- run_GR4J_full(Precip_EDK_full, PET_EDK_full, Ind_WarmUp, Ind_Run)

Dates_Run <- ts$date[Ind_Run]
Qobs_Run  <- ts$discharge_mm_d[Ind_Run]

## ---------------------------------------------------------------------------
## 4) Compare
## ---------------------------------------------------------------------------
NSE_HYRAS <- NSE(sim = Qsim_HYRAS, obs = Qobs_Run)
NSE_EDK   <- NSE(sim = Qsim_EDK,   obs = Qobs_Run)
pct_diff_mean <- 100 * (mean(Qsim_EDK) - mean(Qsim_HYRAS)) / mean(Qsim_HYRAS)

png("GR4J_datasource_uncertainty_HYRAS_vs_EDK.png", width = 1700, height = 800, res = 150)
par(mar = c(4, 4, 3, 1))
plot(Dates_Run, Qobs_Run, type = "l", col = "black", lwd = 1,
     xlab = "Date", ylab = "Q (mm/d)",
     main = "GR4J | Data-source uncertainty: HYRAS-driven vs EDK-driven streamflow")
lines(Dates_Run, Qsim_HYRAS, col = "dodgerblue", lwd = 1)
lines(Dates_Run, Qsim_EDK,   col = "darkorange", lwd = 1)
legend("topright", bty = "n", lwd = 1.5, col = c("black", "dodgerblue", "darkorange"),
       legend = c("Observed", sprintf("HYRAS-driven (NSE=%.3f)", NSE_HYRAS),
                  sprintf("EDK-driven (NSE=%.3f)", NSE_EDK)))
dev.off()

write.csv(data.frame(date = Dates_Run, Qobs = Qobs_Run,
                      Qsim_HYRAS = Qsim_HYRAS, Qsim_EDK = Qsim_EDK),
          "GR4J_datasource_uncertainty_HYRAS_vs_EDK.csv", row.names = FALSE)

cat(sprintf("\nHYRAS-driven GR4J: NSE = %.3f\n", NSE_HYRAS))
cat(sprintf("EDK-driven GR4J:   NSE = %.3f\n", NSE_EDK))
cat(sprintf("Mean streamflow difference (EDK vs HYRAS): %+.1f%%\n", pct_diff_mean))

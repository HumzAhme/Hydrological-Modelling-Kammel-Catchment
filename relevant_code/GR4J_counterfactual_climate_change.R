################################################################################
############   GR4J - counterfactual climate-change analysis            ######
############   "How much has climate change already affected            ######
############   streamflow in this catchment?"                           ######
################################################################################
# ADDITIONAL script - does not modify GR4J_CemaNeige_calibrated.R or its
# outputs. Uses the group's existing, already-delivered GR4J parameter set
# (GR4J_best_parameter_set.csv) - this analysis is about isolating a
# climate signal, not about re-calibrating the model.
#
# METHOD (matches report.qmd's Section 4.4/5.3 TODO exactly):
#   Factual:        Q1(t) = GR4J(T_obs,  P_obs)
#   Counterfactual: Q0(t) = GR4J(T_noCC, P_noCC)
#   where T_noCC/P_noCC have the estimated long-term climate-change trend
#   removed from the observed record. Three experiments isolate each
#   driver's individual contribution:
#     1. TOTAL effect        - both T and P de-trended
#     2. PRECIPITATION-only  - only P de-trended, T stays observed
#     3. TEMPERATURE-only    - only T de-trended, P stays observed
#   Focus variable: June-August (JJA) mean discharge, per report.qmd.
#
# IMPORTANT - READ BEFORE USING THE NUMBERS BELOW IN THE REPORT:
# This script needs a magnitude for "the CMIP6 long-term climate-change
# trend" (delta-T, delta-P) to remove from the observations. We could not
# access actual CMIP6 grid output for this catchment (no internet access
# to CMIP6 data portals from this environment). Instead we use published,
# citable OBSERVED regional warming figures as a defensible stand-in for
# "the trend CMIP6 historical runs are designed to reproduce":
#   - DELTA_T = 1.4 degC, based on ~0.37 degC/decade of German warming
#     since 1970 (World Bank Climate Change Knowledge Portal, ERA5-based)
#     scaled to this study's ~37-year record (1983-2020); consistent with
#     the German Environment Agency's (UBA) 2023 Monitoring Report finding
#     that Bavaria's long-term warming is at or above the German mean of
#     1.5-1.8 degC since the late 19th century.
#   - DELTA_P_PCT = -3%, reflecting that observed German/Bavarian total
#     precipitation shows no strong, consistent long-term trend (several
#     sources), with a small negative value used here as an illustrative,
#     conservative placeholder for the "summers becoming slightly drier"
#     pattern noted in the literature, rather than a hard number specific
#     to this catchment.
# REPLACE these two constants with your instructor's specific CMIP6
# dataset value if one was provided - everything below updates
# automatically. Treat the numbers in this script as a demonstration of
# the METHOD, not a final, precise attribution result, unless/until the
# trend values are confirmed against the actual CMIP6 source intended by
# the course.
DELTA_T     <- 1.4   # degC, added warming already realised in T_obs
DELTA_P_PCT <- -3     # %, change already realised in P_obs (negative = drier)

library(airGR)

## ---------------------------------------------------------------------------
## 1) Data import
## ---------------------------------------------------------------------------
ts <- read.csv("Kammel_all_descriptors_1983_2020.csv")
CatchmentArea_km2 <- 254
ts$discharge_mm_d <- ts$discharge_m3_s * 86.4 / CatchmentArea_km2
ts$PET_mm[ts$PET_mm < 0] <- 0     # same fix as every other script in this repo
ts$date <- as.Date(ts$date)
DatesR <- as.POSIXct(ts$date, format = "%Y-%m-%d", tz = "UTC")

## ---------------------------------------------------------------------------
## 2) Hargreaves-Samani PET, exactly matching compute_pet_hargreaves.py,
##    so we can recompute PET for the temperature counterfactual. Ra is
##    astronomical (latitude + day-of-year only) and does NOT change under
##    the counterfactual, so we reuse the existing Ra_mm_day column.
## ---------------------------------------------------------------------------
hargreaves_pet <- function(Tmean, Tmax, Tmin, Ra_mm_day) {
  pet <- 0.0023 * Ra_mm_day * (Tmean + 17.8) * sqrt(pmax(Tmax - Tmin, 0))
  pmax(pet, 0)   # same negative-PET clip used throughout this repo
}

# sanity check: recomputed PET should match the existing PET_mm column
# (both derived the same way) before we trust this function for the
# counterfactual - small floating differences are fine, large ones are not
pet_check <- hargreaves_pet(ts$Tmean_degC, ts$Tmax_degC, ts$Tmin_degC, ts$Ra_mm_day)
pet_check_diff <- max(abs(pet_check - pmax(ts$PET_mm, 0)), na.rm = TRUE)
cat("Hargreaves PET reproduction check - max abs difference from PET_mm column:",
    round(pet_check_diff, 6), "mm/d (should be ~0)\n")

## ---------------------------------------------------------------------------
## 3) Construct the counterfactual (de-trended) forcing
##    Temperature shift is applied uniformly to Tmean/Tmax/Tmin (preserves
##    the diurnal range Tmax-Tmin, so only the "(Tmean+17.8)" term in
##    Hargreaves changes - a standard simplifying assumption).
##
##    SIGN CONVENTION: DELTA_T and DELTA_P_PCT describe the change climate
##    change has ALREADY caused in the observed record (e.g. DELTA_P_PCT =
##    -3 means "observed precip is already 3% lower than a no-climate-
##    change world would show"). The counterfactual REVERSES that change
##    (T_noCC = T_obs - DELTA_T removes added warming; P_noCC scales UP
##    when DELTA_P_PCT is negative, to restore the wetter no-CC baseline).
## ---------------------------------------------------------------------------
ts$Tmean_noCC <- ts$Tmean_degC - DELTA_T
ts$Tmax_noCC  <- ts$Tmax_degC  - DELTA_T
ts$Tmin_noCC  <- ts$Tmin_degC  - DELTA_T
ts$PET_noCC   <- hargreaves_pet(ts$Tmean_noCC, ts$Tmax_noCC, ts$Tmin_noCC, ts$Ra_mm_day)

ts$Precip_noCC <- ts$precipitation_mm * (1 - DELTA_P_PCT / 100)

## ---------------------------------------------------------------------------
## 4) GR4J setup - plain GR4J, no snow, using the group's existing
##    calibrated parameter set (not re-calibrated here)
## ---------------------------------------------------------------------------
Param <- read.csv("GR4J_best_parameter_set.csv")
Param <- setNames(Param$value, Param$parameter)
cat("Using existing calibrated GR4J parameters:\n"); print(Param)

idx <- function(d) which(format(DatesR, "%Y-%m-%d") == d)
Ind_WarmUp <- idx("1983-01-01"):idx("1987-12-31")
Ind_Run    <- idx("1988-01-01"):idx("2010-12-31")   # observed-discharge period only

run_GR4J <- function(Precip, PotEvap) {
  IM <- CreateInputsModel(FUN_MOD = RunModel_GR4J, DatesR = DatesR,
                           Precip = Precip, PotEvap = PotEvap)
  RO <- CreateRunOptions(FUN_MOD = RunModel_GR4J, InputsModel = IM,
                          IndPeriod_Run = Ind_Run, IndPeriod_WarmUp = Ind_WarmUp)
  RunModel(InputsModel = IM, RunOptions = RO, Param = Param, FUN_MOD = RunModel_GR4J)$Qsim
}

## ---------------------------------------------------------------------------
## 5) Run the factual simulation and all three counterfactual experiments
## ---------------------------------------------------------------------------
cat("\nRunning factual + 3 counterfactual experiments...\n")
Q1_factual   <- run_GR4J(ts$precipitation_mm, ts$PET_mm)                 # observed forcing
Q0_total     <- run_GR4J(ts$Precip_noCC,      ts$PET_noCC)               # both de-trended
Q0_precip    <- run_GR4J(ts$Precip_noCC,      ts$PET_mm)                 # P de-trended only
Q0_temp      <- run_GR4J(ts$precipitation_mm, ts$PET_noCC)               # T (via PET) de-trended only

Dates_Run <- ts$date[Ind_Run]

## ---------------------------------------------------------------------------
## 6) Focus variable: June-August (JJA) mean discharge, per year
## ---------------------------------------------------------------------------
is_jja <- format(Dates_Run, "%m") %in% c("06", "07", "08")
years_run <- as.integer(format(Dates_Run, "%Y"))

jja_mean_by_year <- function(Q) {
  tapply(Q[is_jja], years_run[is_jja], mean)
}

JJA_factual <- jja_mean_by_year(Q1_factual)
JJA_total   <- jja_mean_by_year(Q0_total)
JJA_precip  <- jja_mean_by_year(Q0_precip)
JJA_temp    <- jja_mean_by_year(Q0_temp)

jja_years <- as.integer(names(JJA_factual))

## ---------------------------------------------------------------------------
## 7) Results table + summary
## ---------------------------------------------------------------------------
results_tab <- data.frame(
  year = jja_years,
  JJA_Qfactual_mm_d = round(JJA_factual, 3),
  JJA_Qcounterfactual_total_mm_d = round(JJA_total, 3),
  dQ_total_pct = round(100 * (JJA_factual - JJA_total) / JJA_total, 1),
  JJA_Qcounterfactual_precip_mm_d = round(JJA_precip, 3),
  dQ_precip_pct = round(100 * (JJA_factual - JJA_precip) / JJA_precip, 1),
  JJA_Qcounterfactual_temp_mm_d = round(JJA_temp, 3),
  dQ_temp_pct = round(100 * (JJA_factual - JJA_temp) / JJA_temp, 1)
)
write.csv(results_tab, "GR4J_counterfactual_JJA_results.csv", row.names = FALSE)

mean_effect_total  <- mean(results_tab$dQ_total_pct)
mean_effect_precip <- mean(results_tab$dQ_precip_pct)
mean_effect_temp   <- mean(results_tab$dQ_temp_pct)

driest_year <- jja_years[which.min(JJA_factual)]
driest_val  <- min(JJA_factual)

summary_txt <- c(
  "GR4J counterfactual climate-change analysis - Kammel catchment",
  "================================================================",
  sprintf("Trend removed: DELTA_T = %.1f degC, DELTA_P = %+d%% (see script header for sourcing)", DELTA_T, DELTA_P_PCT),
  sprintf("Period analysed: 1988-2010 (%d JJA seasons)", length(jja_years)),
  "",
  sprintf("Mean TOTAL climate-change effect on JJA discharge:         %+.1f%%", mean_effect_total),
  sprintf("Mean PRECIPITATION-only effect on JJA discharge:           %+.1f%%", mean_effect_precip),
  sprintf("Mean TEMPERATURE-only effect on JJA discharge:             %+.1f%%", mean_effect_temp),
  "",
  sprintf("Driest observed JJA season: %d (mean Q = %.3f mm/d)", driest_year, driest_val),
  "",
  "Interpretation notes for the report:",
  "- With plain GR4J (no snow module), temperature can only affect",
  "  streamflow THROUGH potential evapotranspiration (Hargreaves PET),",
  "  since this model has no snow/melt process at all. The temperature-only",
  "  effect above is therefore entirely an evapotranspiration-demand signal,",
  "  not a snowmelt-timing signal. A model WITH a snow module (HBV, GR6J)",
  "  would likely show a materially different (and probably larger, since",
  "  it would also include shifted snowmelt timing) temperature-driven",
  "  effect for the same DELTA_T - a good, genuine point for the",
  "  'counterfactual analysis interpretation' discussion question about",
  "  model realism.",
  "- Because DELTA_P_PCT is small and precipitation is the dominant",
  "  streamflow driver, the precipitation-only effect and the total effect",
  "  should be fairly close in magnitude; the gap between them is",
  "  attributable to the temperature/PET pathway.",
  "- Replace DELTA_T / DELTA_P_PCT at the top of this script with your",
  "  instructor's specific CMIP6 value if one exists, and rerun."
)
writeLines(summary_txt, "GR4J_counterfactual_summary.txt")
cat("\n\n", paste(summary_txt, collapse = "\n"), "\n")

## ---------------------------------------------------------------------------
## 8) Plots
## ---------------------------------------------------------------------------
png("GR4J_counterfactual_JJA_timeseries.png", width = 1700, height = 900, res = 150)
par(mar = c(4, 4, 3, 1))
plot(jja_years, JJA_factual, type = "b", pch = 16, col = "black",
     ylim = range(c(JJA_factual, JJA_total, JJA_precip, JJA_temp)),
     xlab = "Year", ylab = "Mean JJA Q (mm/d)",
     main = "GR4J | Factual vs. counterfactual (no climate-change trend) JJA discharge")
lines(jja_years, JJA_total,  col = "firebrick",   lwd = 2)
lines(jja_years, JJA_precip, col = "dodgerblue",  lwd = 1.5, lty = 2)
lines(jja_years, JJA_temp,   col = "darkorange",  lwd = 1.5, lty = 2)
legend("topright", bty = "n",
       legend = c("Factual (observed forcing)", "Counterfactual - total effect",
                  "Counterfactual - precipitation-only", "Counterfactual - temperature-only"),
       col = c("black", "firebrick", "dodgerblue", "darkorange"),
       lty = c(1, 1, 2, 2), lwd = c(1, 2, 1.5, 1.5), pch = c(16, NA, NA, NA))
dev.off()

png("GR4J_counterfactual_effect_bars.png", width = 1500, height = 800, res = 150)
par(mar = c(4, 4, 3, 1))
barplot(c(Total = mean_effect_total, Precipitation = mean_effect_precip, Temperature = mean_effect_temp),
        col = c("firebrick", "dodgerblue", "darkorange"),
        ylab = "Mean change in JJA discharge due to climate change so far (%)",
        main = "GR4J | Average counterfactual effect by driver, 1988-2010")
abline(h = 0, lty = 2)
dev.off()

cat("\nAll outputs written. See GR4J_counterfactual_JJA_results.csv, ",
    "GR4J_counterfactual_summary.txt, and the two PNG plots.\n")

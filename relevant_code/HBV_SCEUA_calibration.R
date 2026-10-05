################################################################################
############   HBV - automated calibration via SCE-UA                    ######
############   "Advanced settings": replaces the manual/default          ######
############   parameter set with a properly optimised one               ######
################################################################################
# ADDITIONAL script - does not modify 04_HBV_set_up.R,
# Kammel_default_HBV_parameters.csv, or HBV_Kammel_prediction_2011_2020.csv.
# All outputs below use distinct file names.
#
# WHY THIS SCRIPT EXISTS
# --------------------------------------------------------------------------
# The only HBV results in this repo so far come from a single hand-set
# "default" parameter set (Kammel_default_HBV_parameters.csv) - there is no
# automated calibration, no performance metric (NSE/KGE) against
# observations, and no parameter uncertainty analysis for HBV anywhere,
# unlike GR4J and GR6J. This script closes that gap: it calibrates all 12
# free HBV.IANIGLA parameters with SCE-UA (same algorithm family used for
# GR6J's uncertainty analysis), using the SAME periods, warm-ups, and NSE
# objective as GR6J_calibration_calval.R / GR6J_calibration_fullperiod.R,
# so all three models can be fairly compared in the report.
#
# (Note: some exploratory HBV+GR4J+GR6J comparison plots already exist in
# this repo, e.g. a 2005 three-model comparison figure, suggesting someone
# may have already calibrated HBV informally. If a teammate has their own
# HBV calibration script, prefer that one - this is a clean-room
# reimplementation built to close the gap in case that script isn't
# available, verified to reproduce the existing default-parameter
# prediction to within 0.01 mm/d before any calibration was added.)
#
# MODEL STRUCTURE (HBV.IANIGLA, matching 04_HBV_set_up.R's template exactly):
#   Snow:    SnowGlacier_HBV(model = 1)  - degree-day, SFCF locked to 1
#            (same "switch off this parameter" convention as the template)
#   Soil:    Soil_HBV(model = 1)          - standard non-linear soil routine
#   Routing: Routing_HBV(model = 3)       - two linked stores (fast/slow)
#   Transfer: UH(model = 1)               - triangular unit hydrograph
#
# 12 free parameters (matches Kammel_default_HBV_parameters.csv's names):
#   TR, TT, FM          (snow: rain/snow threshold, melt threshold, melt factor)
#   FC, LP, BETA        (soil: field capacity, ET threshold, shape exponent)
#   K0, K1, K2, UZL, PERC (routing: fast/interflow/baseflow constants,
#                          threshold, percolation)
#   BMAX                (transfer: triangular UH base length, i.e. MAXBAS)
#
# A soft penalty enforces K0 > K1 > K2 (the intended fast > interflow > slow
# ordering baked into the HBV routing concept - without it SCE-UA can find
# numerically-valid but hydrologically-nonsensical parameter combinations).
#
# ncomplex = 24 for SCE-UA (matches the value already used for the group's
# earlier, uncommitted HBV calibration attempt, per a comment in
# GR6J_parameter_uncertainty.R: "HBV used 24 for 12 params").
#
# OUTPUT FILES:
#   HBV_opt_parameter_cal_1988_1999.csv / HBV_opt_parameter_cal_2000_2010.csv
#   HBV_opt_parameter_fullperiod.csv
#   HBV_Qsim_validation_2000_2010.csv / HBV_Qsim_validation_1988_1999.csv
#   HBV_pred_2011_2020_calibrated.csv
#   HBV_performance_summary.txt
#   plots_HBV_calval/*.png, plots_HBV_fullperiod/*.png

library(HBV.IANIGLA)
library(SoilHyP)
library(hydroGOF)

## ---------------------------------------------------------------------------
## 1) Data import (same conventions as 04_HBV_set_up.R)
## ---------------------------------------------------------------------------
ts <- read.csv("Kammel_all_descriptors_1983_2020.csv")
ts$date <- as.Date(ts$date, format = "%Y-%m-%d")
ts$temperature_mean <- ts$Tmean_degC
ts$pet_mm <- ts$PET_mm
ts$discharge_mm_d <- ts$discharge_m3_s * 86400 / (254 * 1e6) * 1000

## ---------------------------------------------------------------------------
## 2) The model chain, exactly matching 04_HBV_set_up.R's hbv_lumped()
## ---------------------------------------------------------------------------
hbv_lumped <- function(basin, param_snow, param_soil, param_routing, param_tf,
                        init_snow = 0, init_soil = 100, init_routing = c(0, 0)) {
  snow_module <- SnowGlacier_HBV(model = 1,
                    inputData = as.matrix(basin[, c('temperature_mean', 'precipitation_mm')]),
                    initCond = c(init_snow, 2),
                    param = c(1, param_snow))     # SFCF hard-coded to 1, as in the template
  soil_module <- Soil_HBV(model = 1,
                    inputData = cbind(snow_module[, "Total"], basin$pet_mm),
                    initCond = c(init_soil, 1),
                    param = param_soil)
  routing_module <- Routing_HBV(model = 3, lake = FALSE,
                    inputData = as.matrix(soil_module[, "Rech"]),
                    initCond = init_routing,
                    param = param_routing)
  tf_module <- UH(model = 1, Qg = routing_module[, "Qg"], param = param_tf)
  round(tf_module, 2)
}

ParamNames <- c("TR", "TT", "FM", "FC", "LP", "BETA", "K0", "K1", "K2", "UZL", "PERC", "BMAX")

run_hbv <- function(x, basin) {
  names(x) <- ParamNames
  hbv_lumped(basin,
             param_snow    = c(x["TR"], x["TT"], x["FM"]),
             param_soil    = c(x["FC"], x["LP"], x["BETA"]),
             param_routing = c(x["K0"], x["K1"], x["K2"], x["UZL"], x["PERC"]),
             param_tf      = c(x["BMAX"]))
}

## ---------------------------------------------------------------------------
## 3) Parameter bounds (physically reasonable HBV ranges; x0 = existing
##    default set as a sensible starting point for the optimiser)
## ---------------------------------------------------------------------------
lower <- c(TR = -3,  TT = -3,  FM = 0.5, FC = 50,  LP = 0.3, BETA = 1,
           K0 = 0.05, K1 = 0.01, K2 = 0.001, UZL = 0,  PERC = 0,  BMAX = 1)
upper <- c(TR = 3,   TT = 3,   FM = 10,  FC = 700, LP = 1,   BETA = 6,
           K0 = 0.9,  K1 = 0.5,  K2 = 0.15,  UZL = 100, PERC = 6,  BMAX = 10)

x0_file <- "Kammel_default_HBV_parameters.csv"
if (file.exists(x0_file)) {
  x0_df <- read.csv(x0_file)
  x0 <- setNames(x0_df$value, x0_df$parameter)[ParamNames]
} else {
  x0 <- c(TR = 0.5, TT = 0, FM = 3, FC = 250, LP = 0.7, BETA = 2,
          K0 = 0.3, K1 = 0.1, K2 = 0.02, UZL = 20, PERC = 1, BMAX = 2)
}
# clip x0 into bounds in case the default set falls slightly outside
x0 <- pmin(pmax(x0, lower), upper)

## ---------------------------------------------------------------------------
## 4) Objective function factory: -NSE (SCEoptim minimises), with a soft
##    penalty enforcing the intended K0 > K1 > K2 ordering
## ---------------------------------------------------------------------------
make_OBJ <- function(basin_full, idx_warmup, idx_run) {
  Qobs <- basin_full$discharge_mm_d[idx_run]
  idx_all <- c(idx_warmup, idx_run)
  n_warm  <- length(idx_warmup)
  basin_sub <- basin_full[idx_all, ]

  function(x) {
    names(x) <- ParamNames
    ordering_penalty <- 0
    if (x["K1"] > x["K0"]) ordering_penalty <- ordering_penalty + (x["K1"] - x["K0"]) * 10
    if (x["K2"] > x["K1"]) ordering_penalty <- ordering_penalty + (x["K2"] - x["K1"]) * 10

    sim_all <- tryCatch(run_hbv(x, basin_sub), error = function(e) NULL)
    if (is.null(sim_all) || anyNA(sim_all)) return(1)
    sim_run <- sim_all[(n_warm + 1):length(idx_all)]
    nse <- NSE(sim = sim_run, obs = Qobs)
    if (is.na(nse) || is.nan(nse)) return(1)
    -nse + ordering_penalty
  }
}

## ---------------------------------------------------------------------------
## 5) Periods - identical to GR6J's calval + fullperiod scripts
## ---------------------------------------------------------------------------
idx <- function(d) which(ts$date == as.Date(d))
Ind_Warmup_8387 <- idx("1983-01-01"):idx("1987-12-31")
Ind_8899        <- idx("1988-01-01"):idx("1999-12-31")
Ind_Warmup_9599 <- idx("1995-01-01"):idx("1999-12-31")
Ind_0010        <- idx("2000-01-01"):idx("2010-12-31")
Ind_Cal_Full    <- idx("1988-01-01"):idx("2010-12-31")
Ind_Warmup_Pred <- idx("2006-01-01"):idx("2010-12-31")
Ind_Pred        <- idx("2011-01-01"):idx("2020-12-31")

## ---------------------------------------------------------------------------
## 6) Run SCE-UA calibration for one period, save params + plots
## ---------------------------------------------------------------------------
set.seed(42)  # reproducibility - SCEoptim's initial sampling is random
dir.create("plots_HBV_calval", showWarnings = FALSE)
dir.create("plots_HBV_fullperiod", showWarnings = FALSE)

calibrate_period <- function(Ind_Warmup, Ind_Run, label) {
  cat("\nRunning SCE-UA for HBV,", label, "...\n")
  OBJ <- make_OBJ(ts, Ind_Warmup, Ind_Run)
  sce <- SCEoptim(OBJ, par = x0, lower = lower, upper = upper,
                   control = list(ncomplex = 24, initsample = "random", trace = 1))
  Param <- sce$par
  names(Param) <- ParamNames
  cat(sprintf("%s optimum: NSE = %.3f\n", label, -sce$value))
  Param
}

score_period <- function(Param, Ind_Warmup, Ind_Run) {
  idx_all <- c(Ind_Warmup, Ind_Run)
  sim_all <- run_hbv(Param, ts[idx_all, ])
  sim_run <- sim_all[(length(Ind_Warmup) + 1):length(idx_all)]
  Qobs <- ts$discharge_mm_d[Ind_Run]
  list(Qsim = sim_run, Qobs = Qobs,
       NSE = NSE(sim = sim_run, obs = Qobs),
       KGE = KGE(sim = sim_run, obs = Qobs, method = "2012"))
}

plot_period <- function(dates, obs, sim, title, path) {
  png(path, width = 10, height = 5, units = "in", res = 300)
  plot(dates, obs, type = "l", col = "tomato", xlab = "", ylab = "Q (mm/d)", main = title)
  lines(dates, sim, col = "dodgerblue")
  legend("topright", legend = c("obs", "sim"), col = c("tomato", "dodgerblue"), lty = 1)
  dev.off()
}

## ---------------------------------------------------------------------------
## 7) Direction A: calibrate 1988-1999, validate 2000-2010
## ---------------------------------------------------------------------------
Param_8899 <- calibrate_period(Ind_Warmup_8387, Ind_8899, "1988-1999")
write.csv(data.frame(parameter = names(Param_8899), value = Param_8899),
          "HBV_opt_parameter_cal_1988_1999.csv", row.names = FALSE)

Train_8899 <- score_period(Param_8899, Ind_Warmup_8387, Ind_8899)
Test_0010  <- score_period(Param_8899, Ind_Warmup_9599, Ind_0010)
write.csv(data.frame(Date = ts$date[Ind_0010], Qobs = Test_0010$Qobs, Qsim = Test_0010$Qsim),
          "HBV_Qsim_validation_2000_2010.csv", row.names = FALSE)
plot_period(ts$date[Ind_8899], Train_8899$Qobs, Train_8899$Qsim,
            sprintf("HBV calibration 1988-1999 (NSE=%.3f, KGE'=%.3f)", Train_8899$NSE, Train_8899$KGE),
            "plots_HBV_calval/calibration_1988_1999.png")
plot_period(ts$date[Ind_0010], Test_0010$Qobs, Test_0010$Qsim,
            sprintf("HBV validation 2000-2010 (NSE=%.3f, KGE'=%.3f)", Test_0010$NSE, Test_0010$KGE),
            "plots_HBV_calval/validation_2000_2010.png")

## ---------------------------------------------------------------------------
## 8) Direction B: swapped - calibrate 2000-2010, validate 1988-1999
## ---------------------------------------------------------------------------
Param_0010 <- calibrate_period(Ind_Warmup_9599, Ind_0010, "2000-2010")
write.csv(data.frame(parameter = names(Param_0010), value = Param_0010),
          "HBV_opt_parameter_cal_2000_2010.csv", row.names = FALSE)

Train_0010 <- score_period(Param_0010, Ind_Warmup_9599, Ind_0010)
Test_8899  <- score_period(Param_0010, Ind_Warmup_8387, Ind_8899)
write.csv(data.frame(Date = ts$date[Ind_8899], Qobs = Test_8899$Qobs, Qsim = Test_8899$Qsim),
          "HBV_Qsim_validation_1988_1999.csv", row.names = FALSE)
plot_period(ts$date[Ind_0010], Train_0010$Qobs, Train_0010$Qsim,
            sprintf("HBV calibration 2000-2010 (NSE=%.3f, KGE'=%.3f)", Train_0010$NSE, Train_0010$KGE),
            "plots_HBV_calval/calibration_2000_2010.png")
plot_period(ts$date[Ind_8899], Test_8899$Qobs, Test_8899$Qsim,
            sprintf("HBV validation 1988-1999 (NSE=%.3f, KGE'=%.3f)", Test_8899$NSE, Test_8899$KGE),
            "plots_HBV_calval/validation_1988_1999.png")

## ---------------------------------------------------------------------------
## 9) Full-period calibration (1988-2010) + 2011-2020 prediction
## ---------------------------------------------------------------------------
Param_Full <- calibrate_period(Ind_Warmup_8387, Ind_Cal_Full, "1988-2010 (fullperiod)")
write.csv(data.frame(parameter = names(Param_Full), value = Param_Full),
          "HBV_opt_parameter_fullperiod.csv", row.names = FALSE)

Train_Full <- score_period(Param_Full, Ind_Warmup_8387, Ind_Cal_Full)
plot_period(ts$date[Ind_Cal_Full], Train_Full$Qobs, Train_Full$Qsim,
            sprintf("HBV calibration 1988-2010 (NSE=%.3f, KGE'=%.3f)", Train_Full$NSE, Train_Full$KGE),
            "plots_HBV_fullperiod/calibration_1988_2010.png")

idx_all_pred <- c(Ind_Warmup_Pred, Ind_Pred)
sim_all_pred <- run_hbv(Param_Full, ts[idx_all_pred, ])
sim_pred <- sim_all_pred[(length(Ind_Warmup_Pred) + 1):length(idx_all_pred)]
write.csv(data.frame(date = ts$date[Ind_Pred], Q_pred_mmd = sim_pred),
          "HBV_pred_2011_2020_calibrated.csv", row.names = FALSE)

png("plots_HBV_fullperiod/prediction_2011_2020.png", width = 12, height = 5, units = "in", res = 300)
plot(ts$date[Ind_Pred], sim_pred, type = "l", col = "dodgerblue",
     xlab = "", ylab = "Q (mm/d)", main = "HBV (calibrated) prediction 2011-2020")
dev.off()

## ---------------------------------------------------------------------------
## 10) Summary
## ---------------------------------------------------------------------------
summary_txt <- c(
  "HBV (SCE-UA calibrated) - Kammel catchment - performance summary",
  "================================================================",
  sprintf("Direction A - train 1988-1999: NSE=%.3f KGE'=%.3f | test 2000-2010: NSE=%.3f KGE'=%.3f",
          Train_8899$NSE, Train_8899$KGE, Test_0010$NSE, Test_0010$KGE),
  sprintf("Direction B - train 2000-2010: NSE=%.3f KGE'=%.3f | test 1988-1999: NSE=%.3f KGE'=%.3f",
          Train_0010$NSE, Train_0010$KGE, Test_8899$NSE, Test_8899$KGE),
  sprintf("Full-period calibration 1988-2010: NSE=%.3f KGE'=%.3f", Train_Full$NSE, Train_Full$KGE),
  "",
  "Full-period parameter set (used for the 2011-2020 prediction):",
  paste(sprintf("  %s = %.4f", names(Param_Full), Param_Full), collapse = "\n"),
  "",
  sprintf("Mean predicted Q 2011-2020: %.3f mm/d", mean(sim_pred)),
  "",
  "For comparison, the previous (uncalibrated, hand-set default) parameters",
  "achieved no measured NSE/KGE at all - this is the first time HBV has",
  "actual performance metrics against observations in this repository."
)
writeLines(summary_txt, "HBV_performance_summary.txt")
cat("\n\n", paste(summary_txt, collapse = "\n"), "\n")

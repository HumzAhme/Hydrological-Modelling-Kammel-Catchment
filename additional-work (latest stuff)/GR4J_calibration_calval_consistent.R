################################################################################
############   GR4J (no snow) - calibration / validation                #######
############   Matches GR6J_calibration_calval.R's methodology exactly, #######
############   so the two models are directly, fairly comparable        #######
################################################################################
# This is an ADDITIONAL script. It does not modify or rerun
# GR4J_CemaNeige_calibrated.R or any of its saved outputs - all outputs
# below use distinct file names.
#
# WHY THIS SCRIPT EXISTS
# --------------------------------------------------------------------------
# GR4J_CemaNeige_calibrated.R (the group's original GR4J deliverable) and
# GR6J_calibration_calval.R / GR6J_calibration_fullperiod.R (GR6J's
# deliverable) used two methodological choices differently:
#   1. Calibration/validation split: GR4J used 1988-2001 / 2002-2010
#      (one direction only). GR6J used 1988-1999 / 2000-2010, run in BOTH
#      directions (swapped split-sample).
#   2. Objective function: GR4J was calibrated against KGE. GR6J was
#      calibrated against NSE.
# For a fair best-vs-worst comparison in the report, both models should be
# evaluated the same way. This script reruns GR4J with GR6J's exact
# methodology (same periods, same warm-ups, same objective, same swapped
# design) so Section 5/6 of the report can compare like-for-like. GR6J's
# own scripts and outputs are NOT touched - only GR4J is being changed,
# by explicit request.
#
# Model: plain GR4J, no CemaNeige (matches the "GR4J - no snow module -
# worst" label already used for the group's original GR4J deliverable).
#
# OUTPUT FILES (distinct names, nothing overwritten):
#   opt_parameter_GR4J_cal_1988_1999.csv / opt_parameter_GR4J_cal_2000_2010.csv
#   Qsim_GR4J_validation_2000_2010.csv / Qsim_GR4J_validation_1988_1999.csv
#   plots_GR4J_calval_consistent/*.png

library(airGR)
library(hydroGOF)

## ---------------------------------------------------------------------------
## 1) Data import & unit conversion (same fix as every other script here)
## ---------------------------------------------------------------------------
ts <- read.csv("Kammel_all_descriptors_1983_2020.csv")

CatchmentArea_km2 <- 254
ts$discharge_mm_d <- ts$discharge_m3_s * 86.4 / CatchmentArea_km2

# Clip negative PET values (artefact on a few extremely cold days) - airGR
# truncates the whole series if it finds invalid PotEvap values
ts$PET_mm <- pmax(ts$PET_mm, 0)

DatesR <- as.POSIXct(ts$date, format = "%Y-%m-%d", tz = "UTC")

## ---------------------------------------------------------------------------
## 2) InputsModel - plain GR4J, no snow module, no TempMean needed
## ---------------------------------------------------------------------------
InputsModel <- CreateInputsModel(FUN_MOD = RunModel_GR4J,
                                  DatesR = DatesR,
                                  Precip = ts$precipitation_mm,
                                  PotEvap = ts$PET_mm)

## ---------------------------------------------------------------------------
## 3) Base periods - identical dates to GR6J_calibration_calval.R
## ---------------------------------------------------------------------------
Ind_Warmup_8387 <- which(format(DatesR, "%Y-%m-%d") == "1983-01-01"):
                   which(format(DatesR, "%Y-%m-%d") == "1987-12-31")
Ind_8899        <- which(format(DatesR, "%Y-%m-%d") == "1988-01-01"):
                   which(format(DatesR, "%Y-%m-%d") == "1999-12-31")

Ind_Warmup_9599 <- which(format(DatesR, "%Y-%m-%d") == "1995-01-01"):
                   which(format(DatesR, "%Y-%m-%d") == "1999-12-31")
Ind_0010        <- which(format(DatesR, "%Y-%m-%d") == "2000-01-01"):
                   which(format(DatesR, "%Y-%m-%d") == "2010-12-31")

## ---------------------------------------------------------------------------
## 4) Helper: calibrate on one period, validate on the other, save outputs
##    (same structure as GR6J's run_split(), FUN_MOD/Param count adapted
##    for plain GR4J - no MeanAnSolidPrecip, since there's no snow module)
## ---------------------------------------------------------------------------
run_split <- function(Ind_WarmTrain, Ind_Train, Ind_WarmTest, Ind_Test,
                       train_id, test_id, example_year,
                       param_csv, qsim_csv, plot_dir) {

  # --- calibrate on the training period, objective = NSE (matches GR6J) ---
  RunOpt_Train <- CreateRunOptions(FUN_MOD = RunModel_GR4J,
                                    InputsModel = InputsModel,
                                    IndPeriod_Run = Ind_Train,
                                    IndPeriod_WarmUp = Ind_WarmTrain)

  Qobs_Train <- ts$discharge_mm_d[Ind_Train]

  InputsCrit_Train <- CreateInputsCrit(FUN_CRIT = ErrorCrit_NSE,
                                        InputsModel = InputsModel,
                                        RunOptions = RunOpt_Train,
                                        VarObs = "Q",
                                        Obs = Qobs_Train)

  CalibOptions <- CreateCalibOptions(FUN_MOD = RunModel_GR4J,
                                      FUN_CALIB = Calibration_Michel)

  OutputsCalib <- Calibration_Michel(InputsModel = InputsModel,
                                      RunOptions = RunOpt_Train,
                                      InputsCrit = InputsCrit_Train,
                                      CalibOptions = CalibOptions,
                                      FUN_MOD = RunModel_GR4J)

  Param <- OutputsCalib$ParamFinalR
  names(Param) <- c("X1", "X2", "X3", "X4")

  write.csv(data.frame(Parameter = names(Param), Value = Param),
            param_csv, row.names = FALSE)

  Out_Train <- RunModel(InputsModel = InputsModel, RunOptions = RunOpt_Train,
                         Param = Param, FUN_MOD = RunModel_GR4J)
  NSE_Train <- NSE(sim = Out_Train$Qsim, obs = Qobs_Train)
  KGE_Train <- KGE(sim = Out_Train$Qsim, obs = Qobs_Train, method = "2012")

  Dates_Train <- as.Date(DatesR[Ind_Train])

  png(file.path(plot_dir, paste0("calibration_", train_id, ".png")),
      width = 10, height = 5, units = "in", res = 300)
  plot(Dates_Train, Qobs_Train, type = "l", col = "tomato",
       xlab = "", ylab = "Q (mm/d)",
       main = paste0("GR4J calibration ", train_id, " (NSE = ", round(NSE_Train, 3),
                      ", KGE' = ", round(KGE_Train, 3), ")"))
  lines(Dates_Train, Out_Train$Qsim, col = "dodgerblue")
  legend("topright", legend = c("obs", "sim"), col = c("tomato", "dodgerblue"), lty = 1)
  dev.off()

  # --- apply the same parameters to the test period (no refitting) ---
  RunOpt_Test <- CreateRunOptions(FUN_MOD = RunModel_GR4J,
                                   InputsModel = InputsModel,
                                   IndPeriod_Run = Ind_Test,
                                   IndPeriod_WarmUp = Ind_WarmTest)

  Out_Test <- RunModel(InputsModel = InputsModel, RunOptions = RunOpt_Test,
                        Param = Param, FUN_MOD = RunModel_GR4J)

  Qobs_Test <- ts$discharge_mm_d[Ind_Test]
  NSE_Test <- NSE(sim = Out_Test$Qsim, obs = Qobs_Test)
  KGE_Test <- KGE(sim = Out_Test$Qsim, obs = Qobs_Test, method = "2012")

  Dates_Test <- as.Date(DatesR[Ind_Test])

  png(file.path(plot_dir, paste0("validation_", test_id, ".png")),
      width = 10, height = 5, units = "in", res = 300)
  plot(Dates_Test, Qobs_Test, type = "l", col = "tomato",
       xlab = "", ylab = "Q (mm/d)",
       main = paste0("GR4J validation ", test_id, " (NSE = ", round(NSE_Test, 3),
                      ", KGE' = ", round(KGE_Test, 3), ")"))
  lines(Dates_Test, Out_Test$Qsim, col = "dodgerblue")
  legend("topright", legend = c("obs", "sim"), col = c("tomato", "dodgerblue"), lty = 1)
  dev.off()

  write.csv(data.frame(Date = Dates_Test, Qobs = Qobs_Test, Qsim = Out_Test$Qsim),
            qsim_csv, row.names = FALSE)

  idx_year <- which(format(Dates_Test, "%Y") == as.character(example_year))
  NSE_Year <- NSE(sim = Out_Test$Qsim[idx_year], obs = Qobs_Test[idx_year])
  KGE_Year <- KGE(sim = Out_Test$Qsim[idx_year], obs = Qobs_Test[idx_year], method = "2012")

  png(file.path(plot_dir, paste0("example_year_", example_year, ".png")),
      width = 10, height = 5, units = "in", res = 300)
  plot(Dates_Test[idx_year], Qobs_Test[idx_year], type = "l", col = "tomato",
       xlab = "", ylab = "Q (mm/d)",
       main = paste0("GR4J example year ", example_year, " (NSE = ", round(NSE_Year, 3),
                      ", KGE' = ", round(KGE_Year, 3), ")"))
  lines(Dates_Test[idx_year], Out_Test$Qsim[idx_year], col = "dodgerblue")
  legend("topright", legend = c("obs", "sim"), col = c("tomato", "dodgerblue"), lty = 1)
  dev.off()

  cat("\n---", train_id, "calibration /", test_id, "validation (GR4J) ---\n")
  print(round(Param, 3))
  cat("Train NSE/KGE':", round(NSE_Train, 3), "/", round(KGE_Train, 3), "\n")
  cat("Test  NSE/KGE':", round(NSE_Test, 3), "/", round(KGE_Test, 3), "\n")

  invisible(list(Param = Param, NSE_Train = NSE_Train, KGE_Train = KGE_Train,
                  NSE_Test = NSE_Test, KGE_Test = KGE_Test))
}

## ---------------------------------------------------------------------------
## 5) Output folder
## ---------------------------------------------------------------------------
dir.create("plots_GR4J_calval_consistent", showWarnings = FALSE)

## ---------------------------------------------------------------------------
## 6) Direction A: calibrate 1988-1999, validate 2000-2010
## ---------------------------------------------------------------------------
Result_A <- run_split(Ind_Warmup_8387, Ind_8899, Ind_Warmup_9599, Ind_0010,
                       train_id = "1988_1999", test_id = "2000_2010",
                       example_year = 2005,
                       param_csv = "opt_parameter_GR4J_cal_1988_1999.csv",
                       qsim_csv = "Qsim_GR4J_validation_2000_2010.csv",
                       plot_dir = "plots_GR4J_calval_consistent")

## ---------------------------------------------------------------------------
## 7) Direction B: swapped - calibrate 2000-2010, validate 1988-1999
## ---------------------------------------------------------------------------
Result_B <- run_split(Ind_Warmup_9599, Ind_0010, Ind_Warmup_8387, Ind_8899,
                       train_id = "2000_2010", test_id = "1988_1999",
                       example_year = 1995,
                       param_csv = "opt_parameter_GR4J_cal_2000_2010.csv",
                       qsim_csv = "Qsim_GR4J_validation_1988_1999.csv",
                       plot_dir = "plots_GR4J_calval_consistent")

cat("\n\n=== Summary (GR4J, GR6J-consistent methodology) ===\n")
cat(sprintf("Direction A - train 1988-1999: NSE=%.3f KGE'=%.3f | test 2000-2010: NSE=%.3f KGE'=%.3f\n",
            Result_A$NSE_Train, Result_A$KGE_Train, Result_A$NSE_Test, Result_A$KGE_Test))
cat(sprintf("Direction B - train 2000-2010: NSE=%.3f KGE'=%.3f | test 1988-1999: NSE=%.3f KGE'=%.3f\n",
            Result_B$NSE_Train, Result_B$KGE_Train, Result_B$NSE_Test, Result_B$KGE_Test))

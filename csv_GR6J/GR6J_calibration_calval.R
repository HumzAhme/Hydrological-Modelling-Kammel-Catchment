################################################################################
############        GR6J-CemaNeige - calibration / validation          #######
################################################################################
# Calibrates GR6J+CemaNeige with NSE as objective, reports NSE and KGE'.
# Runs the cal/val split in BOTH directions (1988-1999 -> 2000-2010, and
# swapped: 2000-2010 -> 1988-1999), saving plots and CSVs for each.

library(airGR)
library(hydroGOF)

## ---------------------------------------------------------------------------
## 1) Data import & unit conversion
## ---------------------------------------------------------------------------
ts <- read.csv("Kammel_all_descriptors_1983_2020.csv")

CatchmentArea_km2 <- 254
ts$discharge_mm_d <- ts$discharge_m3_s * 86.4 / CatchmentArea_km2

# Clip negative PET values (artefact on a few extremely cold days) - airGR
# truncates the whole series if it finds invalid PotEvap values
ts$PET_mm <- pmax(ts$PET_mm, 0)

DatesR <- as.POSIXct(ts$date, format = "%Y-%m-%d", tz = "UTC")

## ---------------------------------------------------------------------------
## 2) InputsModel
## ---------------------------------------------------------------------------
InputsModel <- CreateInputsModel(FUN_MOD = RunModel_CemaNeigeGR6J,
                                  DatesR = DatesR,
                                  Precip = ts$precipitation_mm,
                                  PotEvap = ts$PET_mm,
                                  TempMean = ts$Tmean_degC,
                                  NLayers = 1)

## ---------------------------------------------------------------------------
## 3) Base periods (reused for both directions of the split)
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
## 4) Helper: mean annual solid precipitation (needed by CreateRunOptions)
## ---------------------------------------------------------------------------
calc_MeanAnSolidPrecip <- function(idx) {
  SolidPrecip <- ifelse(ts$Tmean_degC[idx] < 0, ts$precipitation_mm[idx], 0)
  mean(tapply(SolidPrecip, format(DatesR[idx], "%Y"), sum))
}

## ---------------------------------------------------------------------------
## 5) Helper: calibrate on one period, validate on the other, save outputs
## ---------------------------------------------------------------------------
run_split <- function(Ind_WarmTrain, Ind_Train, Ind_WarmTest, Ind_Test,
                       train_id, test_id, example_year,
                       param_csv, qsim_csv, plot_dir) {

  MASP_Train <- calc_MeanAnSolidPrecip(c(Ind_WarmTrain, Ind_Train))
  MASP_Test  <- calc_MeanAnSolidPrecip(c(Ind_WarmTest, Ind_Test))

  # --- calibrate on the training period ---
  RunOpt_Train <- CreateRunOptions(FUN_MOD = RunModel_CemaNeigeGR6J,
                                    InputsModel = InputsModel,
                                    IndPeriod_Run = Ind_Train,
                                    IndPeriod_WarmUp = Ind_WarmTrain,
                                    MeanAnSolidPrecip = MASP_Train)

  Qobs_Train <- ts$discharge_mm_d[Ind_Train]

  InputsCrit_Train <- CreateInputsCrit(FUN_CRIT = ErrorCrit_NSE,
                                        InputsModel = InputsModel,
                                        RunOptions = RunOpt_Train,
                                        VarObs = "Q",
                                        Obs = Qobs_Train)

  CalibOptions <- CreateCalibOptions(FUN_MOD = RunModel_CemaNeigeGR6J,
                                      FUN_CALIB = Calibration_Michel)

  OutputsCalib <- Calibration_Michel(InputsModel = InputsModel,
                                      RunOptions = RunOpt_Train,
                                      InputsCrit = InputsCrit_Train,
                                      CalibOptions = CalibOptions,
                                      FUN_MOD = RunModel_CemaNeigeGR6J)

  Param <- OutputsCalib$ParamFinalR
  names(Param) <- c("X1", "X2", "X3", "X4", "X5", "X6", "CNX1", "CNX2")

  # save calibrated parameter set
  write.csv(data.frame(Parameter = names(Param), Value = Param),
            param_csv, row.names = FALSE)

  # run on training period with the calibrated parameters
  Out_Train <- RunModel(InputsModel = InputsModel, RunOptions = RunOpt_Train,
                         Param = Param, FUN_MOD = RunModel_CemaNeigeGR6J)
  NSE_Train <- NSE(sim = Out_Train$Qsim, obs = Qobs_Train)
  KGE_Train <- KGE(sim = Out_Train$Qsim, obs = Qobs_Train, method = "2012")

  Dates_Train <- as.Date(DatesR[Ind_Train])

  png(file.path(plot_dir, paste0("calibration_", train_id, ".png")),
      width = 10, height = 5, units = "in", res = 300)
  plot(Dates_Train, Qobs_Train, type = "l", col = "tomato",
       xlab = "", ylab = "Q (mm/d)",
       main = paste0("Calibration ", train_id, " (NSE = ", round(NSE_Train, 3),
                      ", KGE' = ", round(KGE_Train, 3), ")"))
  lines(Dates_Train, Out_Train$Qsim, col = "dodgerblue")
  legend("topright", legend = c("obs", "sim"), col = c("tomato", "dodgerblue"), lty = 1)
  dev.off()

  # --- apply the same parameters to the test period (no refitting) ---
  RunOpt_Test <- CreateRunOptions(FUN_MOD = RunModel_CemaNeigeGR6J,
                                   InputsModel = InputsModel,
                                   IndPeriod_Run = Ind_Test,
                                   IndPeriod_WarmUp = Ind_WarmTest,
                                   MeanAnSolidPrecip = MASP_Test)

  Out_Test <- RunModel(InputsModel = InputsModel, RunOptions = RunOpt_Test,
                        Param = Param, FUN_MOD = RunModel_CemaNeigeGR6J)

  Qobs_Test <- ts$discharge_mm_d[Ind_Test]
  NSE_Test <- NSE(sim = Out_Test$Qsim, obs = Qobs_Test)
  KGE_Test <- KGE(sim = Out_Test$Qsim, obs = Qobs_Test, method = "2012")

  Dates_Test <- as.Date(DatesR[Ind_Test])

  png(file.path(plot_dir, paste0("validation_", test_id, ".png")),
      width = 10, height = 5, units = "in", res = 300)
  plot(Dates_Test, Qobs_Test, type = "l", col = "tomato",
       xlab = "", ylab = "Q (mm/d)",
       main = paste0("Validation ", test_id, " (NSE = ", round(NSE_Test, 3),
                      ", KGE' = ", round(KGE_Test, 3), ")"))
  lines(Dates_Test, Out_Test$Qsim, col = "dodgerblue")
  legend("topright", legend = c("obs", "sim"), col = c("tomato", "dodgerblue"), lty = 1)
  dev.off()

  # save validation-period predictions (for later cross-model comparison)
  write.csv(data.frame(Date = Dates_Test, Qobs = Qobs_Test, Qsim = Out_Test$Qsim),
            qsim_csv, row.names = FALSE)

  # --- zoom into one example year, taken from the test period results ---
  idx_year <- which(format(Dates_Test, "%Y") == as.character(example_year))
  NSE_Year <- NSE(sim = Out_Test$Qsim[idx_year], obs = Qobs_Test[idx_year])
  KGE_Year <- KGE(sim = Out_Test$Qsim[idx_year], obs = Qobs_Test[idx_year], method = "2012")

  png(file.path(plot_dir, paste0("example_year_", example_year, ".png")),
      width = 10, height = 5, units = "in", res = 300)
  plot(Dates_Test[idx_year], Qobs_Test[idx_year], type = "l", col = "tomato",
       xlab = "", ylab = "Q (mm/d)",
       main = paste0("Example year ", example_year, " (NSE = ", round(NSE_Year, 3),
                      ", KGE' = ", round(KGE_Year, 3), ")"))
  lines(Dates_Test[idx_year], Out_Test$Qsim[idx_year], col = "dodgerblue")
  legend("topright", legend = c("obs", "sim"), col = c("tomato", "dodgerblue"), lty = 1)
  dev.off()

  cat("\n---", train_id, "calibration /", test_id, "validation ---\n")
  print(round(Param, 3))
  cat("Train NSE/KGE':", round(NSE_Train, 3), "/", round(KGE_Train, 3), "\n")
  cat("Test  NSE/KGE':", round(NSE_Test, 3), "/", round(KGE_Test, 3), "\n")

  invisible(list(Param = Param, NSE_Train = NSE_Train, KGE_Train = KGE_Train,
                  NSE_Test = NSE_Test, KGE_Test = KGE_Test))
}

## ---------------------------------------------------------------------------
## 6) Output folder
## ---------------------------------------------------------------------------
dir.create("plots_GR6J_calval", showWarnings = FALSE)

## ---------------------------------------------------------------------------
## 7) Direction A: calibrate 1988-1999, validate 2000-2010
## ---------------------------------------------------------------------------
Result_A <- run_split(Ind_Warmup_8387, Ind_8899, Ind_Warmup_9599, Ind_0010,
                       train_id = "1988_1999", test_id = "2000_2010",
                       example_year = 2005,
                       param_csv = "opt_parameter_GR6J_cal_1988_1999.csv",
                       qsim_csv = "Qsim_GR6J_validation_2000_2010.csv",
                       plot_dir = "plots_GR6J_calval")

## ---------------------------------------------------------------------------
## 8) Direction B: swapped - calibrate 2000-2010, validate 1988-1999
## ---------------------------------------------------------------------------
Result_B <- run_split(Ind_Warmup_9599, Ind_0010, Ind_Warmup_8387, Ind_8899,
                       train_id = "2000_2010", test_id = "1988_1999",
                       example_year = 1995,
                       param_csv = "opt_parameter_GR6J_cal_2000_2010.csv",
                       qsim_csv = "Qsim_GR6J_validation_1988_1999.csv",
                       plot_dir = "plots_GR6J_calval")


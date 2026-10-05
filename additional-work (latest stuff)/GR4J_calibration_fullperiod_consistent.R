################################################################################
############   GR4J (no snow) - full-period calibration + prediction    #######
############   Matches GR6J_calibration_fullperiod.R's methodology      #######
################################################################################
# ADDITIONAL script - does not modify GR4J_CemaNeige_calibrated.R or its
# outputs. See GR4J_calibration_calval_consistent.R for why this exists
# (matching GR6J's objective function and calibration periods for a fair
# best-vs-worst comparison in the report).
#
# Calibrates plain GR4J (no snow) with NSE as objective on the whole
# observed record (1988-2010), then uses the resulting parameters to
# predict the ungauged 2011-2020 period - same design as
# GR6J_calibration_fullperiod.R.
#
# OUTPUT FILES (distinct names, nothing overwritten):
#   opt_param_GR4J_fullperiod_consistent.csv
#   pred_GR4J_consistent.csv
#   plots_GR4J_fullperiod_consistent/*.png

library(airGR)
library(hydroGOF)

## ---------------------------------------------------------------------------
## 1) Data import & unit conversion
## ---------------------------------------------------------------------------
ts <- read.csv("Kammel_all_descriptors_1983_2020.csv")

CatchmentArea_km2 <- 254
ts$discharge_mm_d <- ts$discharge_m3_s * 86.4 / CatchmentArea_km2
ts$PET_mm <- pmax(ts$PET_mm, 0)

DatesR <- as.POSIXct(ts$date, format = "%Y-%m-%d", tz = "UTC")

## ---------------------------------------------------------------------------
## 2) InputsModel - plain GR4J, no snow module
## ---------------------------------------------------------------------------
InputsModel <- CreateInputsModel(FUN_MOD = RunModel_GR4J,
                                  DatesR = DatesR,
                                  Precip = ts$precipitation_mm,
                                  PotEvap = ts$PET_mm)

## ---------------------------------------------------------------------------
## 3) Periods - identical dates to GR6J_calibration_fullperiod.R
## ---------------------------------------------------------------------------
Ind_Warmup_Cal <- which(format(DatesR, "%Y-%m-%d") == "1983-01-01"):
                  which(format(DatesR, "%Y-%m-%d") == "1987-12-31")
Ind_Cal        <- which(format(DatesR, "%Y-%m-%d") == "1988-01-01"):
                  which(format(DatesR, "%Y-%m-%d") == "2010-12-31")

Ind_Warmup_Pred <- which(format(DatesR, "%Y-%m-%d") == "2006-01-01"):
                   which(format(DatesR, "%Y-%m-%d") == "2010-12-31")
Ind_Pred        <- which(format(DatesR, "%Y-%m-%d") == "2011-01-01"):
                   which(format(DatesR, "%Y-%m-%d") == "2020-12-31")

## ---------------------------------------------------------------------------
## 4) Output folder
## ---------------------------------------------------------------------------
dir.create("plots_GR4J_fullperiod_consistent", showWarnings = FALSE)

## ---------------------------------------------------------------------------
## 5) Calibration on the full 1988-2010 record, objective = NSE
## ---------------------------------------------------------------------------
RunOptions_Cal <- CreateRunOptions(FUN_MOD = RunModel_GR4J,
                                    InputsModel = InputsModel,
                                    IndPeriod_Run = Ind_Cal,
                                    IndPeriod_WarmUp = Ind_Warmup_Cal)

Qobs_Cal <- ts$discharge_mm_d[Ind_Cal]

InputsCrit_Cal <- CreateInputsCrit(FUN_CRIT = ErrorCrit_NSE,
                                    InputsModel = InputsModel,
                                    RunOptions = RunOptions_Cal,
                                    VarObs = "Q",
                                    Obs = Qobs_Cal)

CalibOptions <- CreateCalibOptions(FUN_MOD = RunModel_GR4J,
                                    FUN_CALIB = Calibration_Michel)

OutputsCalib <- Calibration_Michel(InputsModel = InputsModel,
                                    RunOptions = RunOptions_Cal,
                                    InputsCrit = InputsCrit_Cal,
                                    CalibOptions = CalibOptions,
                                    FUN_MOD = RunModel_GR4J)

Param <- OutputsCalib$ParamFinalR
names(Param) <- c("X1", "X2", "X3", "X4")

write.csv(data.frame(Parameter = names(Param), Value = Param),
          "opt_param_GR4J_fullperiod_consistent.csv", row.names = FALSE)

## ---------------------------------------------------------------------------
## 6) Calibration performance + plot (whole observed time series)
## ---------------------------------------------------------------------------
OutputsModel_Cal <- RunModel(InputsModel = InputsModel,
                              RunOptions = RunOptions_Cal,
                              Param = Param,
                              FUN_MOD = RunModel_GR4J)

NSE_Cal <- NSE(sim = OutputsModel_Cal$Qsim, obs = Qobs_Cal)
KGE_Cal <- KGE(sim = OutputsModel_Cal$Qsim, obs = Qobs_Cal, method = "2012")
cat("GR4J calibration NSE (1988-2010):", round(NSE_Cal, 3), "\n")
cat("GR4J calibration KGE'(1988-2010):", round(KGE_Cal, 3), "\n")

Dates_Cal <- as.Date(DatesR[Ind_Cal])

png("plots_GR4J_fullperiod_consistent/calibration_1988_2010.png", width = 12, height = 5, units = "in", res = 300)
plot(Dates_Cal, Qobs_Cal, type = "l", col = "tomato",
     xlab = "", ylab = "Q (mm/d)",
     main = paste0("GR4J calibration 1988-2010 (NSE = ", round(NSE_Cal, 3),
                    ", KGE' = ", round(KGE_Cal, 3), ")"))
lines(Dates_Cal, OutputsModel_Cal$Qsim, col = "dodgerblue")
legend("topright", legend = c("obs", "sim"), col = c("tomato", "dodgerblue"), lty = 1)
dev.off()

## ---------------------------------------------------------------------------
## 7) Prediction for the ungauged 2011-2020 period
## ---------------------------------------------------------------------------
RunOptions_Pred <- CreateRunOptions(FUN_MOD = RunModel_GR4J,
                                     InputsModel = InputsModel,
                                     IndPeriod_Run = Ind_Pred,
                                     IndPeriod_WarmUp = Ind_Warmup_Pred)

OutputsModel_Pred <- RunModel(InputsModel = InputsModel,
                               RunOptions = RunOptions_Pred,
                               Param = Param,
                               FUN_MOD = RunModel_GR4J)

Dates_Pred <- as.Date(DatesR[Ind_Pred])

write.csv(data.frame(Date = Dates_Pred, Qsim = OutputsModel_Pred$Qsim),
          "pred_GR4J_consistent.csv", row.names = FALSE)

## ---------------------------------------------------------------------------
## 8) Combined plot: calibration + prediction, whole 1988-2020 timeline
## ---------------------------------------------------------------------------
Dates_All <- c(Dates_Cal, Dates_Pred)
Qobs_All  <- c(Qobs_Cal, rep(NA, length(Ind_Pred)))
Qsim_All  <- c(OutputsModel_Cal$Qsim, OutputsModel_Pred$Qsim)

png("plots_GR4J_fullperiod_consistent/full_timeline_1988_2020.png", width = 14, height = 5, units = "in", res = 300)
plot(Dates_All, Qsim_All, type = "l", col = "dodgerblue",
     xlab = "", ylab = "Q (mm/d)",
     main = "GR4J full timeline 1988-2020 (obs available 1988-2010 only)")
lines(Dates_All, Qobs_All, col = "tomato")
legend("topright", legend = c("obs", "sim"), col = c("tomato", "dodgerblue"), lty = 1)
dev.off()

## ---------------------------------------------------------------------------
## 9) Prediction period only (2011-2020)
## ---------------------------------------------------------------------------
png("plots_GR4J_fullperiod_consistent/prediction_2011_2020.png", width = 12, height = 5, units = "in", res = 300)
plot(Dates_Pred, OutputsModel_Pred$Qsim, type = "l", col = "dodgerblue",
     xlab = "", ylab = "Q (mm/d)",
     main = "GR4J prediction 2011-2020 (no observations available)")
dev.off()

cat("\nMean predicted Q 2011-2020 (GR4J, consistent methodology):",
    round(mean(OutputsModel_Pred$Qsim), 3), "mm/d\n")

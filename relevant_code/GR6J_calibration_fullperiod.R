################################################################################
############     GR6J-CemaNeige - full-period calibration + prediction #######
################################################################################
# Calibrates GR6J+CemaNeige with NSE as objective on the whole observed
# record (1988-2010), reports NSE and KGE', then uses the resulting
# parameters to predict the ungauged 2011-2020 period.

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
## 3) Periods
## ---------------------------------------------------------------------------
# Calibration: warm-up 1983-1987, run on the whole observed record 1988-2010
# Prediction:  warm-up 2006-2010, run 2011-2020 (no Qobs available there)
Ind_Warmup_Cal <- which(format(DatesR, "%Y-%m-%d") == "1983-01-01"):
                  which(format(DatesR, "%Y-%m-%d") == "1987-12-31")
Ind_Cal        <- which(format(DatesR, "%Y-%m-%d") == "1988-01-01"):
                  which(format(DatesR, "%Y-%m-%d") == "2010-12-31")

Ind_Warmup_Pred <- which(format(DatesR, "%Y-%m-%d") == "2006-01-01"):
                   which(format(DatesR, "%Y-%m-%d") == "2010-12-31")
Ind_Pred        <- which(format(DatesR, "%Y-%m-%d") == "2011-01-01"):
                   which(format(DatesR, "%Y-%m-%d") == "2020-12-31")

## ---------------------------------------------------------------------------
## 4) Mean annual solid precipitation (needed by CreateRunOptions)
## ---------------------------------------------------------------------------
calc_MeanAnSolidPrecip <- function(idx) {
  SolidPrecip <- ifelse(ts$Tmean_degC[idx] < 0, ts$precipitation_mm[idx], 0)
  mean(tapply(SolidPrecip, format(DatesR[idx], "%Y"), sum))
}
MeanAnSolidPrecip_Cal  <- calc_MeanAnSolidPrecip(c(Ind_Warmup_Cal, Ind_Cal))
MeanAnSolidPrecip_Pred <- calc_MeanAnSolidPrecip(c(Ind_Warmup_Pred, Ind_Pred))

## ---------------------------------------------------------------------------
## 5) Output folder
## ---------------------------------------------------------------------------
dir.create("plots_GR6J_fullperiod", showWarnings = FALSE)

## ---------------------------------------------------------------------------
## 6) Calibration on the full 1988-2010 record
## ---------------------------------------------------------------------------
RunOptions_Cal <- CreateRunOptions(FUN_MOD = RunModel_CemaNeigeGR6J,
                                    InputsModel = InputsModel,
                                    IndPeriod_Run = Ind_Cal,
                                    IndPeriod_WarmUp = Ind_Warmup_Cal,
                                    MeanAnSolidPrecip = MeanAnSolidPrecip_Cal)

Qobs_Cal <- ts$discharge_mm_d[Ind_Cal]

InputsCrit_Cal <- CreateInputsCrit(FUN_CRIT = ErrorCrit_NSE,
                                    InputsModel = InputsModel,
                                    RunOptions = RunOptions_Cal,
                                    VarObs = "Q",
                                    Obs = Qobs_Cal)

CalibOptions <- CreateCalibOptions(FUN_MOD = RunModel_CemaNeigeGR6J,
                                    FUN_CALIB = Calibration_Michel)

OutputsCalib <- Calibration_Michel(InputsModel = InputsModel,
                                    RunOptions = RunOptions_Cal,
                                    InputsCrit = InputsCrit_Cal,
                                    CalibOptions = CalibOptions,
                                    FUN_MOD = RunModel_CemaNeigeGR6J)

Param <- OutputsCalib$ParamFinalR
names(Param) <- c("X1", "X2", "X3", "X4", "X5", "X6", "CNX1", "CNX2")

# save calibrated parameter set
write.csv(data.frame(Parameter = names(Param), Value = Param),
          "opt_param_GR6J_fullperiod.csv", row.names = FALSE)

## ---------------------------------------------------------------------------
## 7) Calibration performance + plot (whole observed time series)
## ---------------------------------------------------------------------------
OutputsModel_Cal <- RunModel(InputsModel = InputsModel,
                              RunOptions = RunOptions_Cal,
                              Param = Param,
                              FUN_MOD = RunModel_CemaNeigeGR6J)

NSE_Cal <- NSE(sim = OutputsModel_Cal$Qsim, obs = Qobs_Cal)
KGE_Cal <- KGE(sim = OutputsModel_Cal$Qsim, obs = Qobs_Cal, method = "2012")
cat("Calibration NSE (1988-2010):", round(NSE_Cal, 3), "\n")
cat("Calibration KGE'(1988-2010):", round(KGE_Cal, 3), "\n")

Dates_Cal <- as.Date(DatesR[Ind_Cal])

png("plots_GR6J_fullperiod/calibration_1988_2010.png", width = 12, height = 5, units = "in", res = 300)
plot(Dates_Cal, Qobs_Cal, type = "l", col = "tomato",
     xlab = "", ylab = "Q (mm/d)",
     main = paste0("Calibration 1988-2010 (NSE = ", round(NSE_Cal, 3),
                    ", KGE' = ", round(KGE_Cal, 3), ")"))
lines(Dates_Cal, OutputsModel_Cal$Qsim, col = "dodgerblue")
legend("topright", legend = c("obs", "sim"), col = c("tomato", "dodgerblue"), lty = 1)
dev.off()

## ---------------------------------------------------------------------------
## 8) Prediction for the ungauged 2011-2020 period
## ---------------------------------------------------------------------------
RunOptions_Pred <- CreateRunOptions(FUN_MOD = RunModel_CemaNeigeGR6J,
                                     InputsModel = InputsModel,
                                     IndPeriod_Run = Ind_Pred,
                                     IndPeriod_WarmUp = Ind_Warmup_Pred,
                                     MeanAnSolidPrecip = MeanAnSolidPrecip_Pred)

OutputsModel_Pred <- RunModel(InputsModel = InputsModel,
                               RunOptions = RunOptions_Pred,
                               Param = Param,
                               FUN_MOD = RunModel_CemaNeigeGR6J)

Dates_Pred <- as.Date(DatesR[Ind_Pred])

# save predicted discharge values (no Qobs exists for this period)
write.csv(data.frame(Date = Dates_Pred, Qsim = OutputsModel_Pred$Qsim),
          "pred_GR6J.csv", row.names = FALSE)

## ---------------------------------------------------------------------------
## 9) Combined plot: calibration + prediction, whole 1988-2020 timeline
## ---------------------------------------------------------------------------
Dates_All <- c(Dates_Cal, Dates_Pred)
Qobs_All  <- c(Qobs_Cal, rep(NA, length(Ind_Pred)))   # no obs in 2011-2020
Qsim_All  <- c(OutputsModel_Cal$Qsim, OutputsModel_Pred$Qsim)

png("plots_GR6J_fullperiod/full_timeline_1988_2020.png", width = 14, height = 5, units = "in", res = 300)
plot(Dates_All, Qsim_All, type = "l", col = "dodgerblue",
     xlab = "", ylab = "Q (mm/d)",
     main = "Full timeline 1988-2020 (obs available 1988-2010 only)")
lines(Dates_All, Qobs_All, col = "tomato")
legend("topright", legend = c("obs", "sim"), col = c("tomato", "dodgerblue"), lty = 1)
dev.off()

## ---------------------------------------------------------------------------
## 10) Prediction period only (2011-2020)
## ---------------------------------------------------------------------------
png("plots_GR6J_fullperiod/prediction_2011_2020.png", width = 12, height = 5, units = "in", res = 300)
plot(Dates_Pred, OutputsModel_Pred$Qsim, type = "l", col = "dodgerblue",
     xlab = "", ylab = "Q (mm/d)",
     main = "Prediction 2011-2020 (no observations available)")
dev.off()


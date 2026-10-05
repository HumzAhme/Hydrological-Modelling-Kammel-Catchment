################################################################################
############   GR6J-CemaNeige - parameter & streamflow uncertainty        #####
############   using SCE-UA (SoilHyP) to generate parameter population   #####
################################################################################
# Uses SCE-UA (same algorithm family as the HBV calibration) to generate a
# population of all evaluated parameter sets during optimisation, then keeps
# the top 1% by NSE for:
#
#  1) Parameter ranges: boxplots of top 1% sets, both calibration periods
#  2) Parametric uncertainty on streamflow: min-max envelope of top 1% sets
#     from the FULLPERIOD calibration (pure parameter uncertainty)
#  2b) Bonus: same uncertainty bands per cal period (shows how band shifts)
#  3) Calibration period choice effect: two point-optimal parameter sets
#     (1988-1999 vs. 2000-2010) run over the full 1988-2010 record
#
# Prerequisite: GR6J_calibration_calval.R must have been run to produce
#   opt_parameter_GR6J_cal_1988_1999.csv and opt_parameter_GR6J_cal_2000_2010.csv

library(airGR)
library(hydroGOF)
library(SoilHyP)

## ---------------------------------------------------------------------------
## 1) Data import & unit conversion
## ---------------------------------------------------------------------------
ts <- read.csv("Kammel_all_descriptors_1983_2020.csv")
ts$discharge_mm_d <- ts$discharge_m3_s * 86.4 / 254
ts$PET_mm <- pmax(ts$PET_mm, 0)
DatesR <- as.POSIXct(ts$date, format = "%Y-%m-%d", tz = "UTC")

InputsModel <- CreateInputsModel(FUN_MOD = RunModel_CemaNeigeGR6J,
                                  DatesR = DatesR,
                                  Precip = ts$precipitation_mm,
                                  PotEvap = ts$PET_mm,
                                  TempMean = ts$Tmean_degC,
                                  NLayers = 1)

## ---------------------------------------------------------------------------
## 2) Periods
## ---------------------------------------------------------------------------
Ind_Warmup_8387 <- which(format(DatesR, "%Y-%m-%d") == "1983-01-01"):
                   which(format(DatesR, "%Y-%m-%d") == "1987-12-31")
Ind_8899        <- which(format(DatesR, "%Y-%m-%d") == "1988-01-01"):
                   which(format(DatesR, "%Y-%m-%d") == "1999-12-31")
Ind_Warmup_9599 <- which(format(DatesR, "%Y-%m-%d") == "1995-01-01"):
                   which(format(DatesR, "%Y-%m-%d") == "1999-12-31")
Ind_0010        <- which(format(DatesR, "%Y-%m-%d") == "2000-01-01"):
                   which(format(DatesR, "%Y-%m-%d") == "2010-12-31")

calc_MeanAnSolidPrecip <- function(idx) {
  SolidPrecip <- ifelse(ts$Tmean_degC[idx] < 0, ts$precipitation_mm[idx], 0)
  mean(tapply(SolidPrecip, format(DatesR[idx], "%Y"), sum))
}

## ---------------------------------------------------------------------------
## 3) Parameter bounds for GR6J + CemaNeige
## ---------------------------------------------------------------------------
# Physical plausible ranges (same space Calibration_Michel searches internally)
ParamNames <- c("X1",   "X2",   "X3",   "X4",  "X5",   "X6",  "CNX1", "CNX2")
lower      <- c(1,     -27,     1,      0.5,  -100,    0,     0,      0.01)
upper      <- c(15000,  27,     1000,   4,     100,    20,    1,      20)
# Starting point: use parameters from the full-period calibration as initial guess
x0_file <- "opt_param_GR6J_fullperiod.csv"
if (file.exists(x0_file)) {
  x0_df <- read.csv(x0_file)
  x0 <- setNames(x0_df$Value, x0_df$Parameter)[ParamNames]
} else {
  x0 <- c(300, -1, 50, 2, 0.5, 80, 0.3, 6)  # reasonable fallback
}

## ---------------------------------------------------------------------------
## 4) GR6J objective function wrapper
## ---------------------------------------------------------------------------
# Returns -NSE so SCEoptim can MINIMISE it (no fnscale needed).
# Failed runs return +1 (large positive = bad for minimiser).
make_OBJ <- function(Ind_Warmup, Ind_Run) {
  MASP <- calc_MeanAnSolidPrecip(c(Ind_Warmup, Ind_Run))
  RunOpt <- CreateRunOptions(FUN_MOD = RunModel_CemaNeigeGR6J,
                              InputsModel = InputsModel,
                              IndPeriod_Run = Ind_Run,
                              IndPeriod_WarmUp = Ind_Warmup,
                              MeanAnSolidPrecip = MASP)
  Qobs <- ts$discharge_mm_d[Ind_Run]

  function(x) {
    Out <- tryCatch(
      RunModel(InputsModel = InputsModel, RunOptions = RunOpt,
               Param = x, FUN_MOD = RunModel_CemaNeigeGR6J),
      error = function(e) NULL)
    if (is.null(Out)) return(1)          # failure sentinel: large positive
    nse <- NSE(sim = Out$Qsim, obs = Qobs)
    if (is.na(nse) || is.nan(nse)) return(1)
    -nse                                 # minimise -NSE = maximise NSE
  }
}

## ---------------------------------------------------------------------------
## 5) Run SCE-UA for each period, returning full population
## ---------------------------------------------------------------------------
# ncomplex = 16 is proportional to 8 parameters (HBV used 24 for 12 params)
run_SCE <- function(OBJ, label) {
  cat("\nRunning SCE-UA for", label, "...\n")
  SCEoptim(OBJ, par = x0, lower = lower, upper = upper,
           control = list(ncomplex = 16,
                          initsample = "random",
                          trace = 1,
                          returnpop = TRUE))  # OBJ returns -NSE, no fnscale needed
}

OBJ_8899     <- make_OBJ(Ind_Warmup_8387, Ind_8899)
OBJ_0010     <- make_OBJ(Ind_Warmup_9599, Ind_0010)
Ind_Full_Cal <- c(Ind_8899, Ind_0010)
OBJ_full     <- make_OBJ(Ind_Warmup_8387, Ind_Full_Cal)

SCE_8899 <- run_SCE(OBJ_8899, "1988-1999")
SCE_0010 <- run_SCE(OBJ_0010, "2000-2010")
SCE_full <- run_SCE(OBJ_full, "1988-2010 (fullperiod)")

## ---------------------------------------------------------------------------
## 6) Extract and flatten population from the 3D POP.ALL array
## ---------------------------------------------------------------------------
# POP.ALL: [n_members, n_params, n_iterations]
# POP.FIT.ALL: [n_iterations, n_members]  stores -NSE (OBJ returns -NSE)
# -> negate when extracting to get back to positive NSE scale
flatten_population <- function(sce_result) {
  pop3d   <- sce_result$POP.ALL
  fit_mat <- sce_result$POP.FIT.ALL
  n_iter    <- dim(pop3d)[3]
  n_members <- dim(pop3d)[1]
  pop_flat <- matrix(NA, nrow = n_iter * n_members, ncol = 8)
  fit_flat <- numeric(n_iter * n_members)
  for (i in seq_len(n_iter)) {
    idx <- ((i - 1) * n_members + 1):(i * n_members)
    pop_flat[idx, ] <- pop3d[, , i]
    fit_flat[idx]   <- -fit_mat[i, ]   # negate: -NSE -> positive NSE
  }
  colnames(pop_flat) <- ParamNames
  cat("  Population NSE range (after negation):",
      round(min(fit_flat, na.rm = TRUE), 3), "to",
      round(max(fit_flat, na.rm = TRUE), 3), "\n")
  # remove failed runs (OBJ returned +1 -> stored as -1 -> after negation: 1)
  valid <- !is.na(fit_flat) & fit_flat < 0.99
  list(params = pop_flat[valid, , drop = FALSE],
       nse    = fit_flat[valid])
}

Pop_8899 <- flatten_population(SCE_8899)
Pop_0010 <- flatten_population(SCE_0010)
Pop_full <- flatten_population(SCE_full)

cat("\nTotal valid parameter sets evaluated:\n")
cat("  1988-1999:", nrow(Pop_8899$params), "\n")
cat("  2000-2010:", nrow(Pop_0010$params), "\n")
cat("  Full period:", nrow(Pop_full$params), "\n")

## ---------------------------------------------------------------------------
## 7) Keep near-optimal parameter sets using relative NSE threshold
## ---------------------------------------------------------------------------
# Relative threshold (best_NSE - delta) rather than top 1% by count.
# All retained sets are within delta NSE units of the optimum.
# Tighten delta if the band is still too wide (GR6J equifinality between
# X3 and X6 is real: many parameter combinations achieve similar NSE).
get_top1pct <- function(pop, delta = 0.005) {
  best   <- max(pop$nse, na.rm = TRUE)
  cutoff <- best - delta
  idx    <- pop$nse >= cutoff
  cat("  best NSE:", round(best, 3),
      "| cutoff:", round(cutoff, 3),
      "| sets kept:", sum(idx), "\n")
  list(params = as.data.frame(pop$params[idx, , drop = FALSE]),
       nse    = pop$nse[idx])
}

Top_8899 <- get_top1pct(Pop_8899)
Top_0010 <- get_top1pct(Pop_0010)
Top_full <- get_top1pct(Pop_full)

# save top 1% CSVs
write.csv(cbind(Top_8899$params, NSE = Top_8899$nse),
          "GR6J_top1pct_SCE_1988_1999.csv", row.names = FALSE)
write.csv(cbind(Top_0010$params, NSE = Top_0010$nse),
          "GR6J_top1pct_SCE_2000_2010.csv", row.names = FALSE)
write.csv(cbind(Top_full$params, NSE = Top_full$nse),
          "GR6J_top1pct_SCE_fullperiod.csv", row.names = FALSE)

## ---------------------------------------------------------------------------
## 8) Output folder
## ---------------------------------------------------------------------------
dir.create("plots_GR6J_uncertainty", showWarnings = FALSE)

## ---------------------------------------------------------------------------
## 9) Plot 1: parameter ranges — top 1% sets, both cal periods side by side
## ---------------------------------------------------------------------------
png("plots_GR6J_uncertainty/GR6J_parameter_ranges_top1pct.png",
    width = 14, height = 6, units = "in", res = 300)
par(mfrow = c(2, 4), mar = c(4, 4, 2, 1))
for (p in ParamNames) {
  boxplot(Top_8899$params[[p]], Top_0010$params[[p]],
          names = c("Cal 1988-1999", "Cal 2000-2010"),
          main = p, ylab = "Value",
          col = c("lightblue", "lightgreen"),
          outline = FALSE)
}
dev.off()

## ---------------------------------------------------------------------------
## 10) Helper: run ensemble of parameter sets, return daily envelope
## ---------------------------------------------------------------------------
run_ensemble_envelope <- function(top, Ind_Warmup, Ind_Run) {
  MASP   <- calc_MeanAnSolidPrecip(c(Ind_Warmup, Ind_Run))
  RunOpt <- CreateRunOptions(FUN_MOD = RunModel_CemaNeigeGR6J,
                              InputsModel = InputsModel,
                              IndPeriod_Run = Ind_Run,
                              IndPeriod_WarmUp = Ind_Warmup,
                              MeanAnSolidPrecip = MASP)
  n <- nrow(top$params)
  Qsim_mat <- matrix(NA, nrow = length(Ind_Run), ncol = n)
  for (i in seq_len(n)) {
    Out <- RunModel(InputsModel = InputsModel, RunOptions = RunOpt,
                    Param = as.numeric(top$params[i, ParamNames]),
                    FUN_MOD = RunModel_CemaNeigeGR6J)
    Qsim_mat[, i] <- Out$Qsim
  }
  list(q_lo  = apply(Qsim_mat, 1, quantile, probs = 0.05, na.rm = TRUE),
       q_hi  = apply(Qsim_mat, 1, quantile, probs = 0.95, na.rm = TRUE),
       q_med = apply(Qsim_mat, 1, median),
       dates = as.Date(DatesR[Ind_Run]),
       qobs  = ts$discharge_mm_d[Ind_Run])
}

cat("\nRunning ensembles for streamflow uncertainty bands ...\n")
Env_full <- run_ensemble_envelope(Top_full, Ind_Warmup_8387, Ind_Full_Cal)
Env_8899 <- run_ensemble_envelope(Top_8899, Ind_Warmup_8387, Ind_8899)
Env_0010 <- run_ensemble_envelope(Top_0010, Ind_Warmup_9599, Ind_0010)

## ---------------------------------------------------------------------------
## 11) Helper: plot uncertainty band
## ---------------------------------------------------------------------------
plot_uncertainty_band <- function(env, title, filename) {
  png(filename, width = 12, height = 5, units = "in", res = 300)
  plot(env$dates, env$qobs, type = "n",
       xlab = "", ylab = "Q (mm/d)", main = title)
  polygon(c(env$dates, rev(env$dates)),
          c(env$q_hi, rev(env$q_lo)),
          col = adjustcolor("dodgerblue", alpha.f = 0.35), border = NA)
  lines(env$dates, env$q_med, col = "dodgerblue", lwd = 1)
  lines(env$dates, env$qobs,  col = "tomato",     lwd = 1)
  legend("topright",
         legend = c("obs", "top 1% median", "top 1% 5th-95th percentile"),
         col    = c("tomato", "dodgerblue",
                    adjustcolor("dodgerblue", alpha.f = 0.35)),
         lty = c(1, 1, NA), pch = c(NA, NA, 15), pt.cex = 2)
  dev.off()
}

## ---------------------------------------------------------------------------
## 12) Plot 2: parametric uncertainty — fullperiod top 1%, 1988-2010
## ---------------------------------------------------------------------------
plot_uncertainty_band(Env_full,
                       "Parametric uncertainty - fullperiod SCE top 1% (1988-2010)",
                       "plots_GR6J_uncertainty/GR6J_uncertainty_band_fullperiod.png")

# zoom into 2005
idx_2005 <- which(format(Env_full$dates, "%Y") == "2005")
png("plots_GR6J_uncertainty/GR6J_uncertainty_band_fullperiod_2005.png",
    width = 10, height = 5, units = "in", res = 300)
plot(Env_full$dates[idx_2005], Env_full$qobs[idx_2005], type = "n",
     xlab = "", ylab = "Q (mm/d)",
     main = "Parametric uncertainty - fullperiod SCE top 1%, year 2005")
polygon(c(Env_full$dates[idx_2005], rev(Env_full$dates[idx_2005])),
        c(Env_full$q_hi[idx_2005], rev(Env_full$q_lo[idx_2005])),
        col = adjustcolor("dodgerblue", alpha.f = 0.35), border = NA)
lines(Env_full$dates[idx_2005], Env_full$q_med[idx_2005], col = "dodgerblue", lwd = 1)
lines(Env_full$dates[idx_2005], Env_full$qobs[idx_2005],  col = "tomato",     lwd = 1)
legend("topright",
       legend = c("obs", "top 1% median", "top 1% 5th-95th percentile"),
       col    = c("tomato", "dodgerblue",
                  adjustcolor("dodgerblue", alpha.f = 0.35)),
       lty = c(1, 1, NA), pch = c(NA, NA, 15), pt.cex = 2)
dev.off()

## ---------------------------------------------------------------------------
## 13) Plot 2b: bonus — per-cal-period uncertainty bands
## ---------------------------------------------------------------------------
plot_uncertainty_band(Env_8899,
                       "Parametric uncertainty - SCE cal 1988-1999 top 1% (bonus)",
                       "plots_GR6J_uncertainty/GR6J_uncertainty_band_1988_1999.png")

plot_uncertainty_band(Env_0010,
                       "Parametric uncertainty - SCE cal 2000-2010 top 1% (bonus)",
                       "plots_GR6J_uncertainty/GR6J_uncertainty_band_2000_2010.png")

## ---------------------------------------------------------------------------
## 14) Plot 3: effect of calibration period choice on streamflow
## ---------------------------------------------------------------------------
read_param_csv <- function(path) {
  df <- read.csv(path)
  setNames(df$Value, df$Parameter)
}

Param_8899 <- read_param_csv("opt_parameter_GR6J_cal_1988_1999.csv")
Param_0010 <- read_param_csv("opt_parameter_GR6J_cal_2000_2010.csv")

MASP_Full <- calc_MeanAnSolidPrecip(c(Ind_Warmup_8387, Ind_Full_Cal))
RunOpt_Full <- CreateRunOptions(FUN_MOD = RunModel_CemaNeigeGR6J,
                                 InputsModel = InputsModel,
                                 IndPeriod_Run = Ind_Full_Cal,
                                 IndPeriod_WarmUp = Ind_Warmup_8387,
                                 MeanAnSolidPrecip = MASP_Full)

Out_8899 <- RunModel(InputsModel = InputsModel, RunOptions = RunOpt_Full,
                      Param = Param_8899, FUN_MOD = RunModel_CemaNeigeGR6J)
Out_0010 <- RunModel(InputsModel = InputsModel, RunOptions = RunOpt_Full,
                      Param = Param_0010, FUN_MOD = RunModel_CemaNeigeGR6J)

Qobs_Full  <- ts$discharge_mm_d[Ind_Full_Cal]
Dates_Full <- as.Date(DatesR[Ind_Full_Cal])
NSE_8899   <- NSE(sim = Out_8899$Qsim, obs = Qobs_Full)
NSE_0010   <- NSE(sim = Out_0010$Qsim, obs = Qobs_Full)

png("plots_GR6J_uncertainty/GR6J_calibration_period_effect.png",
    width = 14, height = 5, units = "in", res = 300)
plot(Dates_Full, Qobs_Full, type = "l", col = "black", lwd = 1,
     xlab = "", ylab = "Q (mm/d)",
     main = "Effect of calibration period choice (evaluated on full 1988-2010)")
lines(Dates_Full, Out_8899$Qsim, col = "dodgerblue")
lines(Dates_Full, Out_0010$Qsim, col = "darkorange")
abline(v = as.Date("2000-01-01"), lty = 2, col = "grey50")
legend("topright",
       legend = c("obs",
                  paste0("cal 1988-1999 (NSE = ", round(NSE_8899, 2), ")"),
                  paste0("cal 2000-2010 (NSE = ", round(NSE_0010, 2), ")"),
                  "split boundary"),
       col = c("black", "dodgerblue", "darkorange", "grey50"),
       lty = c(1, 1, 1, 2))
dev.off()

# zoom into 2005 for a closer look at the difference between the two param sets
idx_2005_full <- which(format(Dates_Full, "%Y") == "2005")
NSE_8899_2005 <- NSE(sim = Out_8899$Qsim[idx_2005_full], obs = Qobs_Full[idx_2005_full])
NSE_0010_2005 <- NSE(sim = Out_0010$Qsim[idx_2005_full], obs = Qobs_Full[idx_2005_full])

png("plots_GR6J_uncertainty/GR6J_calibration_period_effect_2005.png",
    width = 10, height = 5, units = "in", res = 300)
plot(Dates_Full[idx_2005_full], Qobs_Full[idx_2005_full],
     type = "l", col = "black", lwd = 1,
     xlab = "", ylab = "Q (mm/d)",
     main = "Effect of calibration period choice - year 2005")
lines(Dates_Full[idx_2005_full], Out_8899$Qsim[idx_2005_full], col = "dodgerblue")
lines(Dates_Full[idx_2005_full], Out_0010$Qsim[idx_2005_full], col = "darkorange")
legend("topright",
       legend = c("obs",
                  paste0("cal 1988-1999 (NSE = ", round(NSE_8899_2005, 2), ")"),
                  paste0("cal 2000-2010 (NSE = ", round(NSE_0010_2005, 2), ")")),
       col = c("black", "dodgerblue", "darkorange"),
       lty = 1)
dev.off()

cat("\nDone. Plots saved in plots_GR6J_uncertainty/\n")


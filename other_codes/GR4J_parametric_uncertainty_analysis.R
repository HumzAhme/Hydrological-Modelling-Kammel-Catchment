################################################################################
####   GR4J parametric uncertainty analysis (Monte Carlo / GLUE-style)    ####
####   Kammel catchment - ADDITIONAL script, does not touch or rerun      ####
####   GR4J_CemaNeige_calibrated.R or any of its saved outputs            ####
################################################################################
#
# WHAT THIS SCRIPT IS FOR
# --------------------------------------------------------------------------
# GR4J_CemaNeige_calibrated.R already produced the group's official GR4J
# deliverable (single best parameter set, no snow module, calibrated with
# Calibration_Michel - see GR4J_performance_summary.txt / GR4J_best_
# parameter_set.csv already in the repo). That script and its outputs are
# NOT modified or rerun here.
#
# This script adds the piece that was still missing for GR4J: the
# assignment's uncertainty-analysis requirement, which needs a *population*
# of good parameter sets (not just one optimum) for two different
# calibration periods. A local-search optimiser like Calibration_Michel
# can't give you that - it deliberately converges to a single best point -
# so this uses random-sampling (Monte Carlo / GLUE-style) calibration
# instead, independently of the existing script.
#
# MODEL: plain GR4J, NO CemaNeige snow module - matching the actual
# no-snow run already saved in this repo (4 parameters: X1-X4), NOT the
# 6-parameter +CemaNeige version. This keeps this script consistent with
# what's already committed as "GR4J, no snow module, worst model".
#
# CALIBRATION PERIODS: 1988-1999 and 2000-2010 (independent, non-
# overlapping) - matching the two periods already used in the GR6J
# teammate's scripts, so the uncertainty results are directly comparable
# across the two models. (Note this differs from GR4J_CemaNeige_
# calibrated.R's own internal split-sample check, which used 1988-2001 /
# 2002-2010 - that script is left untouched; this is a separate analysis.)
#
# METHOD: random uniform sampling of GR4J's 4 parameters across their
# standard ranges, scored by KGE, keeping the top 1% as the "behavioural"
# ensemble - the simplest member of the GLUE family (Beven & Binley, 1992).
# N = 10,000 draws per period (~5 ms/run -> a minute or two total).
#
# OUTPUT FILES (written to the working directory):
#   GR4J_ensemble_top1pct_calibA.csv / _calibB.csv   <- the two ensembles
#   GR4J_parameter_uncertainty_table.csv             <- IQR/range per parameter
#   GR4J_parameter_uncertainty_ranges.png            <- REQUIRED PLOT 1
#   GR4J_parametric_uncertainty_band.png             <- REQUIRED PLOT 2
#   GR4J_calibration_period_effect.png               <- REQUIRED PLOT 3
#   GR4J_uncertainty_analysis_summary.txt            <- plain-text numbers
#
# install.packages("airGR")
library(airGR)

## ===========================================================================
## 1) Data import - same fix as the main script (negative-PET / no snow)
## ===========================================================================
ts <- read.csv("Kammel_all_descriptors_1983_2020.csv")
ts$PET_mm[ts$PET_mm < 0] <- 0                      # same QC fix as the main script
ts$discharge_mm_d <- ts$discharge_m3_s * 86.4 / 254

DatesR <- as.POSIXct(ts$date, format = "%Y-%m-%d", tz = "UTC")

InputsModel <- CreateInputsModel(FUN_MOD = RunModel_GR4J, DatesR = DatesR,
                                  Precip = ts$precipitation_mm, PotEvap = ts$PET_mm)
idx <- function(d) which(format(DatesR, "%Y-%m-%d") == d)

## ===========================================================================
## 2) Two independent calibration periods (matching the GR6J teammate script)
## ===========================================================================
PeriodA <- list(label = "Calib A (1988-1999)",
                 WarmUp = idx("1983-01-01"):idx("1987-12-31"),
                 Run    = idx("1988-01-01"):idx("1999-12-31"))
PeriodB <- list(label = "Calib B (2000-2010)",
                 WarmUp = idx("1995-01-01"):idx("1999-12-31"),
                 Run    = idx("2000-01-01"):idx("2010-12-31"))

# Test/prediction period - same as the group's existing GR4J deliverable
Ind_WarmUp_test <- idx("2006-01-01"):idx("2010-12-31")
Ind_Test         <- idx("2011-01-01"):idx("2020-12-31")
Dates_test       <- as.Date(DatesR[Ind_Test])

## ===========================================================================
## 3) Monte Carlo sampling - GR4J's 4 parameters only (no snow parameters)
## ===========================================================================
set.seed(42)
N <- 10000
ParamNames <- c("X1", "X2", "X3", "X4")
samp <- cbind(X1 = runif(N, 1, 2500),
              X2 = runif(N, -10, 5),
              X3 = runif(N, 1, 1500),
              X4 = runif(N, 0.5, 10))

score_period <- function(Period, samp) {
  RunOptions <- CreateRunOptions(FUN_MOD = RunModel_GR4J, InputsModel = InputsModel,
                                  IndPeriod_Run = Period$Run, IndPeriod_WarmUp = Period$WarmUp)
  Qobs <- ts$discharge_mm_d[Period$Run]
  InputsCrit <- CreateInputsCrit(FUN_CRIT = ErrorCrit_KGE2, InputsModel = InputsModel,
                                  RunOptions = RunOptions, Obs = Qobs)
  vapply(seq_len(nrow(samp)), function(i) {
    Outputs <- RunModel(InputsModel = InputsModel, RunOptions = RunOptions,
                         Param = samp[i, ], FUN_MOD = RunModel_GR4J)
    ErrorCrit_KGE2(InputsCrit, Outputs, verbose = FALSE)$CritValue
  }, numeric(1))
}

cat(sprintf("Running %d Monte Carlo samples x 2 calibration periods ", N))
cat("(this takes a minute or two)...\n")
t0 <- Sys.time()
kge_A <- score_period(PeriodA, samp)
kge_B <- score_period(PeriodB, samp)
cat("Done in", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "sec\n")

## ===========================================================================
## 4) Top 1% ensembles ("behavioural" parameter sets)
## ===========================================================================
top1pct <- function(kge, samp, n_total) {
  k <- max(1, round(0.01 * n_total))
  ord <- order(kge, decreasing = TRUE)[seq_len(k)]
  data.frame(samp[ord, , drop = FALSE], KGE = kge[ord])
}
ensA <- top1pct(kge_A, samp, N)
ensB <- top1pct(kge_B, samp, N)
cat(sprintf("\nTop 1%% ensemble sizes: A = %d sets (best KGE %.3f), B = %d sets (best KGE %.3f)\n",
            nrow(ensA), max(ensA$KGE), nrow(ensB), max(ensB$KGE)))

write.csv(ensA, "GR4J_ensemble_top1pct_calibA.csv", row.names = FALSE)
write.csv(ensB, "GR4J_ensemble_top1pct_calibB.csv", row.names = FALSE)

## ===========================================================================
## PLOT 1 (REQUIRED): parameter ranges of the top-1% ensembles, both periods
## ===========================================================================
png("GR4J_parameter_uncertainty_ranges.png", width = 1700, height = 1050, res = 150)
par(mfrow = c(2, 2), mar = c(3, 4, 3, 1), oma = c(0, 0, 3, 0))
for (p in ParamNames) {
  boxplot(list("Calib A\n1988-1999" = ensA[[p]], "Calib B\n2000-2010" = ensB[[p]]),
          col = c("seagreen", "darkorange"), ylab = p,
          main = paste0("Parameter ", p))
}
mtext("GR4J (no snow) | Parameter ranges, top 1% behavioural sets (KGE-ranked), by calibration period",
      outer = TRUE, line = 1, cex = 1.0, font = 2)
dev.off()

## ===========================================================================
## Run every ensemble member over the TEST period (2011-2020)
## ===========================================================================
run_ensemble_on_test <- function(ens) {
  RunOptions <- CreateRunOptions(FUN_MOD = RunModel_GR4J, InputsModel = InputsModel,
                                  IndPeriod_Run = Ind_Test, IndPeriod_WarmUp = Ind_WarmUp_test)
  Qsim_mat <- sapply(seq_len(nrow(ens)), function(i) {
    Outputs <- RunModel(InputsModel = InputsModel, RunOptions = RunOptions,
                         Param = as.numeric(ens[i, ParamNames]), FUN_MOD = RunModel_GR4J)
    Outputs$Qsim
  })
  Qsim_mat  # matrix: rows = time, columns = ensemble members
}

cat("\nRunning both ensembles over the test period 2011-2020...\n")
QsimA <- run_ensemble_on_test(ensA)
QsimB <- run_ensemble_on_test(ensB)

## ===========================================================================
## PLOT 2 (REQUIRED): parametric uncertainty band on the test-period simulation
## (ensemble A shown as the representative ensemble)
## ===========================================================================
qlo  <- apply(QsimA, 1, quantile, probs = 0.05)
qhi  <- apply(QsimA, 1, quantile, probs = 0.95)
qmed <- apply(QsimA, 1, median)

png("GR4J_parametric_uncertainty_band.png", width = 1700, height = 800, res = 150)
par(mar = c(4, 4, 3, 1))
plot(Dates_test, qmed, type = "n", ylim = c(0, max(qhi)),
     xlab = "Date", ylab = "Q (mm/d)",
     main = "GR4J (no snow) | Parametric uncertainty, test period 2011-2020\n(5th-95th percentile band, top 1% behavioural sets, Calib A)")
polygon(c(Dates_test, rev(Dates_test)), c(qhi, rev(qlo)),
        col = adjustcolor("dodgerblue", alpha.f = 0.3), border = NA)
lines(Dates_test, qmed, col = "dodgerblue4", lwd = 1)
legend("topright", bty = "n", fill = adjustcolor("dodgerblue", alpha.f = 0.3),
       legend = "90% parameter uncertainty band", border = NA)
dev.off()

## ===========================================================================
## PLOT 3 (REQUIRED): effect of calibration period choice on the simulation
## ===========================================================================
# Uses airGR's own local-search optimiser (Calibration_Michel - the same
# method GR4J_CemaNeige_calibrated.R already uses) run independently on
# each period, rather than "whichever Monte Carlo draw happened to score
# highest" - a single top-ranked random draw is a noisy point estimate and
# not a reliable representative optimum.
SearchRanges <- matrix(c(   1, -10,    1, 0.5,
                          2500,   5, 1500,  10),
                       byrow = TRUE, nrow = 2,
                       dimnames = list(c("min", "max"), ParamNames))

optimise_period <- function(Period) {
  RunOptions <- CreateRunOptions(FUN_MOD = RunModel_GR4J, InputsModel = InputsModel,
                                  IndPeriod_Run = Period$Run, IndPeriod_WarmUp = Period$WarmUp)
  Qobs <- ts$discharge_mm_d[Period$Run]
  InputsCrit <- CreateInputsCrit(FUN_CRIT = ErrorCrit_KGE2, InputsModel = InputsModel,
                                  RunOptions = RunOptions, Obs = Qobs)
  CalibOptions <- CreateCalibOptions(FUN_MOD = RunModel_GR4J,
                                      FUN_CALIB = Calibration_Michel, SearchRanges = SearchRanges)
  OutputsCalib <- Calibration_Michel(InputsModel = InputsModel, RunOptions = RunOptions,
                                      InputsCrit = InputsCrit, CalibOptions = CalibOptions,
                                      FUN_MOD = RunModel_GR4J)
  Param <- OutputsCalib$ParamFinalR
  names(Param) <- ParamNames
  Param
}
cat("\nRunning local-search optimisation per period for the calibration-period-effect comparison...\n")
bestA_opt <- optimise_period(PeriodA)
bestB_opt <- optimise_period(PeriodB)
cat("Optimised Calib A params:", paste(sprintf("%s=%.3f", ParamNames, bestA_opt), collapse = ", "), "\n")
cat("Optimised Calib B params:", paste(sprintf("%s=%.3f", ParamNames, bestB_opt), collapse = ", "), "\n")

RunOptions_test <- CreateRunOptions(FUN_MOD = RunModel_GR4J, InputsModel = InputsModel,
                                     IndPeriod_Run = Ind_Test, IndPeriod_WarmUp = Ind_WarmUp_test)
Qsim_bestA <- RunModel(InputsModel = InputsModel, RunOptions = RunOptions_test,
                        Param = as.numeric(bestA_opt), FUN_MOD = RunModel_GR4J)$Qsim
Qsim_bestB <- RunModel(InputsModel = InputsModel, RunOptions = RunOptions_test,
                        Param = as.numeric(bestB_opt), FUN_MOD = RunModel_GR4J)$Qsim

png("GR4J_calibration_period_effect.png", width = 1700, height = 1100, res = 150)
par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
plot(Dates_test, Qsim_bestA, type = "l", col = "seagreen", lwd = 1,
     xlab = "Date", ylab = "Q (mm/d)",
     main = "GR4J (no snow) | Effect of calibration period choice, full test period 2011-2020")
lines(Dates_test, Qsim_bestB, col = "darkorange", lwd = 1)
legend("topright", bty = "n", lwd = 1.5,
       col = c("seagreen", "darkorange"),
       legend = c("Optimised on Calib A (1988-1999)", "Optimised on Calib B (2000-2010)"))

zoom_sel <- Dates_test >= as.Date("2015-01-01") & Dates_test <= as.Date("2015-12-31")
plot(Dates_test[zoom_sel], Qsim_bestA[zoom_sel], type = "l", col = "seagreen", lwd = 1.5,
     xlab = "Date", ylab = "Q (mm/d)",
     main = "Zoom-in: 2015 - same comparison, one year only")
lines(Dates_test[zoom_sel], Qsim_bestB[zoom_sel], col = "darkorange", lwd = 1.5)
legend("topright", bty = "n", lwd = 1.5,
       col = c("seagreen", "darkorange"),
       legend = c("Calib A optimum", "Calib B optimum"))
dev.off()

## ===========================================================================
## Summary stats + parameter uncertainty table
## ===========================================================================
pbias_diff   <- 100 * mean(Qsim_bestB - Qsim_bestA) / mean(Qsim_bestA)
nse_AvsB     <- 1 - sum((Qsim_bestA - Qsim_bestB)^2) / sum((Qsim_bestA - mean(Qsim_bestA))^2)
rel_width_A  <- mean((qhi - qlo) / qmed, na.rm = TRUE)
kge_corr_AB  <- cor(kge_A, kge_B)

param_uncertainty_tab <- do.call(rbind, lapply(ParamNames, function(p) {
  data.frame(parameter = p,
             rangeA = diff(range(ensA[[p]])), iqrA = IQR(ensA[[p]]),
             rangeB = diff(range(ensB[[p]])), iqrB = IQR(ensB[[p]]),
             search_range_width = diff(range(samp[, p])))
}))
param_uncertainty_tab$norm_iqrA_pct <- round(100 * param_uncertainty_tab$iqrA / param_uncertainty_tab$search_range_width, 1)
param_uncertainty_tab$norm_iqrB_pct <- round(100 * param_uncertainty_tab$iqrB / param_uncertainty_tab$search_range_width, 1)
write.csv(param_uncertainty_tab, "GR4J_parameter_uncertainty_table.csv", row.names = FALSE)

summary_txt <- c(
  "GR4J (no snow) - Monte Carlo parameter uncertainty analysis - summary",
  "========================================================================",
  sprintf("Monte Carlo samples: %d, top 1%% kept = %d sets per calibration period", N, nrow(ensA)),
  "",
  sprintf("Calib A (1988-1999) ensemble: best KGE = %.3f, ensemble KGE range = [%.3f, %.3f]",
          max(ensA$KGE), min(ensA$KGE), max(ensA$KGE)),
  sprintf("Calib B (2000-2010) ensemble: best KGE = %.3f, ensemble KGE range = [%.3f, %.3f]",
          max(ensB$KGE), min(ensB$KGE), max(ensB$KGE)),
  sprintf("Correlation between KGE_A and KGE_B across all %d random draws: %.3f", N, kge_corr_AB),
  "",
  "Parameter uncertainty (IQR of top-1% ensemble, normalised to % of search range",
  "- bigger % = more uncertain / less identifiable):",
  paste(capture.output(print(param_uncertainty_tab[, c("parameter","norm_iqrA_pct","norm_iqrB_pct")],
                              row.names = FALSE)), collapse = "\n"),
  "",
  sprintf("Mean relative width of the 90%% parametric-uncertainty band (Calib A ensemble), test period: %.0f%% of median Qsim",
          100 * rel_width_A),
  "",
  "Properly-optimised (Calibration_Michel) parameter sets per period:",
  sprintf("  Calib A optimum: %s", paste(sprintf("%s=%.3f", ParamNames, bestA_opt), collapse = ", ")),
  sprintf("  Calib B optimum: %s", paste(sprintf("%s=%.3f", ParamNames, bestB_opt), collapse = ", ")),
  sprintf("Mean difference in test-period Qsim, optimum-of-B vs optimum-of-A: %+.1f%%", pbias_diff),
  sprintf("NSE of optimum-of-B simulation against optimum-of-A simulation (agreement between the two): %.3f", nse_AvsB),
  sprintf("Mean test-period Qsim: optimum-of-A = %.3f mm/d | optimum-of-B = %.3f mm/d",
          mean(Qsim_bestA), mean(Qsim_bestB))
)
writeLines(summary_txt, "GR4J_uncertainty_analysis_summary.txt")
cat("\n\n", paste(summary_txt, collapse = "\n"), "\n")
cat("\nAll output files written to:", getwd(), "\n")

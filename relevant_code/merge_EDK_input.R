################################################################################
####  Merge EDK input (precip + tavg) with HYRAS-derived PET + discharge  ####
################################################################################
# Produces one CSV that can be used to run GR6J with EDK-driven P and T,
# while keeping PET and observed discharge from the HYRAS/course dataset
# (EDK only provides Tavg, not Tmin/Tmax, so Hargreaves-Samani PET cannot be
# recomputed from EDK alone).

setwd("C:/Users/lukas/Documents/Studium/02_semester/Hydrological_modeling/Kammel_model/hydrological-modelling-kammel-catchment")

## ---------------------------------------------------------------------------
## 1) Load both datasets
## ---------------------------------------------------------------------------

ts_hyras <- read.csv("Kammel_all_descriptors_1983_2020.csv")
ts_hyras$date <- as.Date(ts_hyras$date)

edk <- read.csv("Kammel_EDK_input.csv")
edk$date <- as.Date(edk$date)

# NOTE: adjust these two column names below if your Kammel_EDK_input.csv
# uses different names than precip_EDK / tavg_EDK
stopifnot(all(c("date", "precip_EDK", "tavg_EDK") %in% names(edk)))

cat("HYRAS date range:", format(range(ts_hyras$date)), "\n")
cat("EDK   date range:", format(range(edk$date)), "\n")

## ---------------------------------------------------------------------------
## 2) Merge - keep only dates present in BOTH datasets
## ---------------------------------------------------------------------------

ts_edk <- merge(edk, ts_hyras[, c("date", "PET_mm", "discharge_m3_s")], by = "date")

cat("Merged date range:", format(range(ts_edk$date)), "\n")
cat("Merged rows:", nrow(ts_edk), "\n")

## ---------------------------------------------------------------------------
## 3) Unit conversion + cleaning (same steps as the HYRAS pipeline)
## ---------------------------------------------------------------------------

CatchmentArea_km2 <- 254
ts_edk$discharge_mm_d <- ts_edk$discharge_m3_s * 86.4 / CatchmentArea_km2

# clip negative PET values (same artefact-handling as in GR6J_calibration_calval.R)
ts_edk$PET_mm <- pmax(ts_edk$PET_mm, 0)

# sort by date, just in case merge() reordered anything
ts_edk <- ts_edk[order(ts_edk$date), ]

## ---------------------------------------------------------------------------
## 4) Sanity check: compare EDK vs HYRAS precipitation over the overlap period
## ---------------------------------------------------------------------------

ts_hyras_overlap <- ts_hyras[ts_hyras$date %in% ts_edk$date, ]
ts_hyras_overlap <- ts_hyras_overlap[order(ts_hyras_overlap$date), ]

cat("\nCorrelation EDK vs HYRAS precipitation:",
    round(cor(ts_edk$precip_EDK, ts_hyras_overlap$precipitation_mm), 3), "\n")
cat("Mean annual EDK precipitation:  ",
    round(sum(ts_edk$precip_EDK) / (nrow(ts_edk) / 365.25), 1), "mm/yr\n")
cat("Mean annual HYRAS precipitation:",
    round(sum(ts_hyras_overlap$precipitation_mm) / (nrow(ts_hyras_overlap) / 365.25), 1), "mm/yr\n")

png("EDK_vs_HYRAS_precip_scatter.png", width = 8, height = 8, units = "in", res = 300)
plot(ts_hyras_overlap$precipitation_mm, ts_edk$precip_EDK, pch = 16, cex = 0.4,
     col = adjustcolor("dodgerblue", alpha.f = 0.3),
     xlab = "HYRAS precipitation (mm/day)", ylab = "EDK precipitation (mm/day)",
     main = "EDK vs. HYRAS daily precipitation")
abline(0, 1, col = "tomato", lty = 2)
dev.off()

## ---------------------------------------------------------------------------
## 5) Save merged, model-ready EDK input
## ---------------------------------------------------------------------------

write.csv(ts_edk, "Kammel_EDK_merged_input.csv", row.names = FALSE)

cat("\nSaved: Kammel_EDK_merged_input.csv\n")
cat("Columns:", paste(names(ts_edk), collapse = ", "), "\n")

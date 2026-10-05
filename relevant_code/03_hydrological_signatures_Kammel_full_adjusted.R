###################################################################################################
#### Characterizing catchment response and streamflow dynamics: hydrological signatures ###########
###################################################################################################

# Review paper on hydrological signatures:
# https://wires.onlinelibrary.wiley.com/doi/10.1002/wat2.1499

install.packages(c("sf", "lfstat")) # run once if needed
library(sf)
library(lfstat)

# ---- Working directory ----
setwd("C:/Users/lukas/Documents/Studium/02_semester/Hydrological_modeling/Kammel_model/hydrological-modelling-kammel-catchment")

# ---- Input files ----
csv_file <- "Kammel_all_descriptors_1983_2020.csv"
catchment_file <- "catchment.gpkg"

# ---- Read hydrometeorological time series ----
ts <- read.csv(csv_file)

# ---- Adapt column names from the Kammel CSV to the names expected below ----
ts$date <- as.Date(ts$date)                  # CSV uses YYYY-MM-DD
ts$pet_mm <- ts$PET_mm                       # original script expects pet_mm
ts$temperature_mean <- ts$Tmean_degC         # original script expects temperature_mean

# Precipitation column check
# The Kammel CSV should contain precipitation_mm.
# If your file uses another name, adapt this line, e.g.:
# ts$precipitation_mm <- ts$PREC_mm
if (!"precipitation_mm" %in% names(ts)) {
  stop("Missing column 'precipitation_mm'. Please rename/map your precipitation column to ts$precipitation_mm.")
}

# Discharge column check
if (!"discharge_m3_s" %in% names(ts)) {
  stop("Missing column 'discharge_m3_s'. Please rename/map your discharge column to ts$discharge_m3_s.")
}

# IMPORTANT:
# We keep all rows, including rows where discharge_m3_s is NA.
# For calculations that need discharge, we use complete cases or na.rm = TRUE.

# ---- Read catchment geopackage and calculate catchment area ----
if (file.exists(catchment_file)) {
  catchments <- st_read(catchment_file)
  area_m2 <- as.numeric(st_area(catchments)[1])
} else {
  stop(paste0(
    "Missing catchment file: ", catchment_file, "\n",
    "Put catchment.gpkg into the working directory, or change catchment_file to the correct filename."
  ))
}

# ---- Convert streamflow from m3/s to mm/day ----
# m3 / m2 = m, m * 1000 = mm, seconds per day = 24 * 3600
# NA discharge rows remain NA here, which is intended.
ts$discharge_mm_d <- ts$discharge_m3_s * 1000 * 24 * 3600 / area_m2

# Subset only for analyses that require observed discharge.
ts_q <- ts[!is.na(ts$discharge_m3_s) & !is.na(ts$discharge_mm_d), ]

# Number of years in full meteorological period and observed-discharge period
years_full <- as.numeric(difftime(max(ts$date), min(ts$date), units = "days")) / 365.25
years_q <- as.numeric(difftime(max(ts_q$date), min(ts_q$date), units = "days")) / 365.25

# ---- Plot streamflow and precipitation ----
par(mfrow = c(1, 1))
df.bar <- barplot(ts$precipitation_mm,
                  names.arg = ts$date,
                  xlab = "Date",
                  ylab = "precipitation [mm/day]")
lines(x = df.bar, y = ts$discharge_m3_s, col = "blue")

# Just the first year
n_first_year <- min(365, nrow(ts))
df.bar <- barplot(ts$precipitation_mm[1:n_first_year],
                  names.arg = ts$date[1:n_first_year],
                  xlab = "Date",
                  ylab = "precipitation [mm/day]")
lines(x = df.bar, y = ts$discharge_m3_s[1:n_first_year], col = "blue")

# Are P and Q related at the daily scale?
plot(ts$discharge_m3_s, ts$precipitation_mm,
     xlab = "Discharge [m3/s]",
     ylab = "Precipitation [mm/day]")
cor(ts$discharge_m3_s, ts$precipitation_mm, use = "complete.obs")

# ---- Mean annual totals / averages ----
sum(ts_q$discharge_m3_s, na.rm = TRUE) / years_q          # annual mean Q in m3/s, observed Q period only
sum(ts$precipitation_mm, na.rm = TRUE) / years_full       # annual mean P in mm/year over full period
sum(ts$pet_mm, na.rm = TRUE) / years_full                 # annual mean PET in mm/year over full period

# Check P and Q after converting Q to mm/day
sum(ts_q$discharge_mm_d, na.rm = TRUE) / years_q          # annual mean Q in mm/year, observed Q period only
sum(ts$precipitation_mm, na.rm = TRUE) / years_full       # annual mean P in mm/year over full period

# Plot relation between Q [mm/day] and P [mm/day]
plot(ts$discharge_mm_d, ts$precipitation_mm,
     xlab = "Discharge [mm/day]",
     ylab = "Precipitation [mm/day]")
cor(ts$discharge_mm_d, ts$precipitation_mm, use = "complete.obs")

# ---- Monthly scale ----
order_months <- month.name

monthlyQ <- aggregate(discharge_mm_d ~ month,
                      data = data.frame(discharge_mm_d = ts_q$discharge_mm_d,
                                        month = months(ts_q$date)),
                      FUN = mean,
                      na.rm = TRUE)
monthlyQ <- monthlyQ[order(match(monthlyQ$month, order_months)), ]
barplot(monthlyQ$discharge_mm_d,
        names.arg = monthlyQ$month,
        ylab = "Q [mm/day]")

monthlyP <- aggregate(precipitation_mm ~ month,
                      data = data.frame(precipitation_mm = ts$precipitation_mm,
                                        month = months(ts$date)),
                      FUN = mean,
                      na.rm = TRUE)
monthlyP <- monthlyP[order(match(monthlyP$month, order_months)), ]

monthlyPET <- aggregate(pet_mm ~ month,
                        data = data.frame(pet_mm = ts$pet_mm,
                                          month = months(ts$date)),
                        FUN = mean,
                        na.rm = TRUE)
monthlyPET <- monthlyPET[order(match(monthlyPET$month, order_months)), ]

# Plot monthly P, Q, and PET
max_flux <- max(monthlyP$precipitation_mm,
                monthlyQ$discharge_mm_d,
                monthlyPET$pet_mm,
                na.rm = TRUE)
df.bar <- barplot(monthlyP$precipitation_mm,
                  names.arg = monthlyP$month,
                  ylab = "flux [mm/day]",
                  ylim = c(0, max_flux * 1.1))
lines(x = df.bar, y = monthlyQ$discharge_mm_d, col = "blue")
lines(x = df.bar, y = monthlyPET$pet_mm, col = "green")

# ---- Annual scale ----
ts$year <- format(ts$date, "%Y")
ts_q$year <- format(ts_q$date, "%Y")

annualQ <- aggregate(discharge_mm_d ~ year,
                     data = ts_q,
                     FUN = sum,
                     na.rm = TRUE)
annualP <- aggregate(precipitation_mm ~ year,
                     data = ts,
                     FUN = sum,
                     na.rm = TRUE)

barplot(annualQ$discharge_mm_d,
        names.arg = annualQ$year,
        ylab = "Q [mm/year]")

# Match years with observed discharge for annual P-Q relationship
annualPQ <- merge(annualP, annualQ, by = "year")
plot(annualPQ$precipitation_mm, annualPQ$discharge_mm_d,
     xlab = "Annual precipitation [mm/year]",
     ylab = "Annual discharge [mm/year]")
cor(annualPQ$precipitation_mm, annualPQ$discharge_mm_d, use = "complete.obs")

# ---- Temperature / snow relevance ----
monthlyT <- aggregate(temperature_mean ~ month,
                      data = data.frame(temperature_mean = ts$temperature_mean,
                                        month = months(ts$date)),
                      FUN = mean,
                      na.rm = TRUE)
monthlyT <- monthlyT[order(match(monthlyT$month, order_months)), ]
barplot(monthlyT$temperature_mean,
        names.arg = monthlyT$month,
        ylab = "deg C")

min(ts$temperature_mean, na.rm = TRUE)                         # minimum recorded temperature
length(which(ts$temperature_mean < 0))                          # days below 0°C
percent_freeze_days <- length(which(ts$temperature_mean < 0)) / nrow(ts) * 100
percent_freeze_days

# ---- Dryness index PET/P, Budyko-style indicators ----
dryness_index <- sum(ts$pet_mm, na.rm = TRUE) / sum(ts$precipitation_mm, na.rm = TRUE)
dryness_index

# AET with water balance method: AET = P - Q
# Use only the period with observed discharge for water-balance consistency.
P_q_period <- sum(ts_q$precipitation_mm, na.rm = TRUE)
Q_q_period <- sum(ts_q$discharge_mm_d, na.rm = TRUE)
PET_q_period <- sum(ts_q$pet_mm, na.rm = TRUE)

AET <- (P_q_period - Q_q_period) / years_q
AET

# Evaporative index AET/P
AET / (P_q_period / years_q)

# Aridity / humidity index P/PET
AI <- sum(ts$precipitation_mm, na.rm = TRUE) / sum(ts$pet_mm, na.rm = TRUE)
AI

# Runoff ratio: fraction of precipitation becoming streamflow
# Use only the observed-discharge period.
RR <- Q_q_period / P_q_period
RR

# ---- Flow duration curve ----
Qsort <- sort(ts_q$discharge_mm_d, decreasing = TRUE, na.last = NA)

df <- data.frame(x = 100 / length(Qsort) * seq_along(Qsort),
                 y = Qsort)

plot(x = df$x,
     y = df$y,
     type = "l",
     log = "y",
     ylab = "Discharge [mm/day]",
     xlab = "Percentage of time flow is equaled or exceeded [%]",
     main = "Flow Duration Curve")
grid()

# Slope of the flow duration curve between Q33 and Q66
Q33 <- quantile(ts_q$discharge_mm_d, probs = 0.33, na.rm = TRUE)
Q66 <- quantile(ts_q$discharge_mm_d, probs = 0.66, na.rm = TRUE)
slope_FDC <- (Q33 - Q66) / (0.66 - 0.33) * 100
slope_FDC

# Coefficient of variation of streamflow
# Standard definition is sd / mean.
cvQ <- sd(ts_q$discharge_mm_d, na.rm = TRUE) / mean(ts_q$discharge_mm_d, na.rm = TRUE)
cvQ

# ---- Optional aquifer section ----
# This only runs if aquifer_type.shp exists in your working directory.
if (file.exists("aquifer_type.shp")) {
  aquifer_type <- st_read("aquifer_type.shp")
  plot(aquifer_type["had16_hydr"])

  # IDs:
  # 1 water bodies
  # 2 porous
  # 3 fractured
  # 4 karstics
  # 5 aquitard

  aquifer_type <- st_transform(aquifer_type, st_crs(catchments))

  catchments_aquifer_type <- st_crop(x = aquifer_type["had16_hydr"], catchments[1], snap = "out")
  catchments_aquifer_type <- st_intersection(catchments_aquifer_type, catchments[1])

  plot(catchments_aquifer_type["had16_hydr"])
} else {
  message("Skipping aquifer section: aquifer_type.shp not found.")
}

prop.table(
  tapply(as.numeric(st_area(catchments_aquifer_type)),
         catchments_aquifer_type$had16_hydr,
         sum)
) * 100

# ---- Baseflow separation using recursive digital filter ----
# lfstat/baseflow should receive a continuous observed discharge vector.
# Therefore we use ts_q, while the full ts still keeps missing-discharge rows.
Qbase <- baseflow(ts_q$discharge_mm_d, tp.factor = 0.9, block.len = 5)

plot(ts_q$discharge_mm_d[1:min(365, nrow(ts_q))],
     xlab = "Date",
     ylab = "Streamflow [mm/day]",
     type = "l")
lines(Qbase[1:min(365, length(Qbase))], col = "blue", lwd = 2)

# Baseflow index
BFI <- sum(Qbase, na.rm = TRUE) / sum(ts_q$discharge_mm_d, na.rm = TRUE)
BFI

# Sensitivity to filter parameters
Qbase2 <- baseflow(ts_q$discharge_mm_d, tp.factor = 0.9, block.len = 7)
Qbase3 <- baseflow(ts_q$discharge_mm_d, tp.factor = 0.7, block.len = 7)

plot(ts_q$discharge_mm_d[1:min(365, nrow(ts_q))],
     xlab = "Date",
     ylab = "Streamflow [mm/day]",
     type = "l")
lines(Qbase[1:min(365, length(Qbase))], col = "blue", lwd = 2)
lines(Qbase2[1:min(365, length(Qbase2))], col = "tomato", lwd = 2)
lines(Qbase3[1:min(365, length(Qbase3))], col = "goldenrod", lwd = 2)

# Compare baseflow indices
sum(Qbase, na.rm = TRUE) / sum(ts_q$discharge_mm_d, na.rm = TRUE)
sum(Qbase2, na.rm = TRUE) / sum(ts_q$discharge_mm_d, na.rm = TRUE)
sum(Qbase3, na.rm = TRUE) / sum(ts_q$discharge_mm_d, na.rm = TRUE)

# ---- Recession constant K ----
Qlow_df <- createlfobj(ts(ts_q$discharge_mm_d),
                       startdate = ts_q$date[1],
                       baseflow = FALSE)

# ?recession  # uncomment in an interactive R session if needed
recession(Qlow_df, method = "MRC", seglen = 7, threshold = 70)
recession(Qlow_df, method = "MRC", seglen = 7, threshold = 50)
recession(Qlow_df, method = "MRC", seglen = 7, threshold = 10)

# Effect of segment length
recession(Qlow_df, method = "MRC", seglen = 10, threshold = 50)
recession(Qlow_df, method = "MRC", seglen = 7, threshold = 50)
recession(Qlow_df, method = "MRC", seglen = 5, threshold = 50)
recession(Qlow_df, method = "MRC", seglen = 2, threshold = 50)

# Transform recession coefficient from log scale
K <- exp(-1 / recession(Qlow_df, method = "MRC", seglen = 7, threshold = 50))
K

##################### do it yourself part #####################
### Derive all hydrological signatures for your catchment and compare it to the Selke catchment.
### Feel free to derive additional hydrological signatures from the review paper.

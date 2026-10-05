###################################################################################################
#### Preparing hydrometeorological time series: input variables for any hydrological model ########
###################################################################################################
#### ADAPTED FOR: Kammel catchment, EDK precipitation + temperature (data-source uncertainty) #####
###################################################################################################

install.packages(c("raster", "sf", "ncdf4")) # in case you do not have them installed
library(raster)
library(sf)
library(ncdf4)

setwd("C:/Users/lukas/Documents/Studium/02_semester/Hydrological_modeling/Kammel_model/hydrological-modelling-kammel-catchment")

############ Calculation of daily catchment averaged precipitation #############
###############  The effect of different interpolation methods  ################


####HYRAS dataset of German Weather Center (regression-based method) ############
# --- commented out: HYRAS already processed separately in 02_catchment_averaged_input.R ---

# nc_file=nc_open("pr_hyras_1_2002_v6-0_de.nc")
# print(nc_file)
# precip_array <- ncvar_get(nc_file, "pr")
# fillvalue <- ncatt_get(nc_file, "pr", "_FillValue")
# precip_array[precip_array==fillvalue$value] <- NA
# dim(precip_array)
# time <- ncvar_get(nc_file, "time")
# head(time)
# P_time_DWD <- as.Date(as.POSIXct(time*3600, origin = "1931-01-01",  tz="UTC"))
# head(P_time_DWD)
# x=ncvar_get(nc_file, "x")
# y=ncvar_get(nc_file, "y")
# image(precip_array[,,1])
# P_brick_DWD <- brick(precip_array, xmn=min(x), xmx=max(x), ymn=min(y), ymx=max(y), crs=crs("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs"))
# plot(P_brick_DWD[[1]])
# P_brick_DWD <- flip(t(P_brick_DWD))
# extent(P_brick_DWD)<-extent(min(x),max(x),min(y), max(y))
# plot(P_brick_DWD[[1]])
# Selke<-st_read("Selke_catchment_LHW.shp")
# Selke_P_brick_DWD <- crop( x=P_brick_DWD, Selke[1], snap="out" )
# Selke_P_brick_DWD <-mask( x= Selke_P_brick_DWD, Selke[1], snap="out" )
# plot(Selke_P_brick_DWD[[1]])
# Selke_P_ts_DWD=cellStats(Selke_P_brick_DWD, stat='mean', na.rm=TRUE)
# sum(Selke_P_ts_DWD)
# barplot(Selke_P_ts_DWD, names=P_time_DWD)


####################### External Drift Kriging (EDK) ############################
############################  Kammel catchment  ##################################

## ---- Read Kammel catchment polygon ----
Kammel <- st_read("catchment.gpkg")   # adjust path if catchment.gpkg lives in a subfolder

## EDK grids are in WGS84 lat/lon, unlike HYRAS (LAEA) -> transform catchment polygon
Kammel_wgs84 <- st_transform(Kammel, crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))


## ---------------------------------------------------------------------------
## 1) PRECIPITATION (pre.nc)
## ---------------------------------------------------------------------------

nc_precip <- nc_open("C:/Users/lukas/Documents/Studium/02_semester/Hydrological_modeling/Kammel_model/data/pre.nc")

# Inspect metadata first - check variable name, time units/origin, and time coverage
print(nc_precip)

# NOTE: adjust "pre" below if print(nc_precip) shows a different variable name
precip_array <- ncvar_get(nc_precip, "pre")

fillvalue <- ncatt_get(nc_precip, "pre", "_FillValue")
precip_array[precip_array == fillvalue$value] <- NA

dim(precip_array)

time_precip <- ncvar_get(nc_precip, "time")
head(time_precip)

# NOTE: origin below matches the course demo file - check print(nc_precip) time:units
# attribute and adjust origin if your file differs
P_time_EDK <- as.Date(time_precip, origin = "1949-12-31")

range(P_time_EDK)  # check this covers the period you need (e.g. 1988-2010)

lat <- ncvar_get(nc_precip, "lat")
lon <- ncvar_get(nc_precip, "lon")

P_brick_EDK <- brick(precip_array, xmn = min(lat), xmx = max(lat),
                      ymn = min(lon), ymx = max(lon),
                      crs = crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

## Subset to the years actually needed (calibration periods: 1988-2010)
P_brick_EDK <- subset(P_brick_EDK, which(P_time_EDK >= '1983-01-01' & P_time_EDK <= '2010-12-31'))
P_time_EDK  <- P_time_EDK[P_time_EDK >= '1983-01-01' & P_time_EDK <= '2010-12-31']  # keep the date vector in sync!

plot(P_brick_EDK[[1]])  # looks "wrong" before transpose, that's expected

P_brick_EDK <- t(P_brick_EDK)

plot(P_brick_EDK[[1]])  # should look like DE now

## Crop/mask to Kammel catchment (WGS84 version)
Kammel_P_brick_EDK <- crop(x = P_brick_EDK, Kammel_wgs84[1], snap = "out")
Kammel_P_brick_EDK <- mask(x = Kammel_P_brick_EDK, Kammel_wgs84[1], snap = "out")

plot(Kammel_P_brick_EDK[[1]])

## Daily catchment-averaged precipitation time series
Kammel_P_ts_EDK <- cellStats(Kammel_P_brick_EDK, stat = 'mean', na.rm = TRUE)

sum(Kammel_P_ts_EDK)  # sanity check: total sum over the period


## ---------------------------------------------------------------------------
## 2) TEMPERATURE (tavg.nc)
## ---------------------------------------------------------------------------

nc_tavg <- nc_open("C:/Users/lukas/Documents/Studium/02_semester/Hydrological_modeling/Kammel_model/data/tavg.nc")

# Inspect metadata first - variable name is likely "tavg" or "t", check and adjust below
print(nc_tavg)

# NOTE: adjust "tavg" below to match the actual variable name from print(nc_tavg)
tavg_array <- ncvar_get(nc_tavg, "tavg")

fillvalue_t <- ncatt_get(nc_tavg, "tavg", "_FillValue")
tavg_array[tavg_array == fillvalue_t$value] <- NA

dim(tavg_array)

time_tavg <- ncvar_get(nc_tavg, "time")
head(time_tavg)

# NOTE: check print(nc_tavg) time:units attribute - origin may differ from pre.nc
T_time_EDK <- as.Date(time_tavg, origin = "1949-12-31")

range(T_time_EDK)  # should match P_time_EDK range

lat_t <- ncvar_get(nc_tavg, "lat")
lon_t <- ncvar_get(nc_tavg, "lon")

T_brick_EDK <- brick(tavg_array, xmn = min(lat_t), xmx = max(lat_t),
                     ymn = min(lon_t), ymx = max(lon_t),
                     crs = crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

## Subset to the same years as precipitation (1988-2010)
T_brick_EDK <- subset(T_brick_EDK, which(T_time_EDK >= '1983-01-01' & T_time_EDK <= '2010-12-31'))
T_time_EDK  <- T_time_EDK[T_time_EDK >= '1983-01-01' & T_time_EDK <= '2010-12-31']  # keep in sync!

T_brick_EDK <- t(T_brick_EDK)
plot(T_brick_EDK[[1]])  # should look like DE

## Crop/mask to Kammel catchment
Kammel_T_brick_EDK <- crop(x = T_brick_EDK, Kammel_wgs84[1], snap = "out")
Kammel_T_brick_EDK <- mask(x = Kammel_T_brick_EDK, Kammel_wgs84[1], snap = "out")

plot(Kammel_T_brick_EDK[[1]])

## Daily catchment-averaged mean temperature time series
Kammel_T_ts_EDK <- cellStats(Kammel_T_brick_EDK, stat = 'mean', na.rm = TRUE)


## ---------------------------------------------------------------------------
## 3) Combine and save
## ---------------------------------------------------------------------------

# sanity check: P and T time vectors should be identical in length/dates
length(P_time_EDK) == length(T_time_EDK)
all(P_time_EDK == T_time_EDK)

Kammel_EDK_df <- data.frame(
  date       = P_time_EDK,
  precip_EDK = Kammel_P_ts_EDK,
  tavg_EDK   = Kammel_T_ts_EDK
)

write.csv(Kammel_EDK_df, "Kammel_EDK_input.csv", row.names = FALSE)

## NOTE on PET: EDK only provides Tavg, not Tmin/Tmax, so Hargreaves-Samani PET
## cannot be recomputed from EDK alone. For the EDK-driven model run, reuse the
## HYRAS-derived PET (hold PET constant, vary only P and T -> isolates the
## data-source effect on P and T specifically).

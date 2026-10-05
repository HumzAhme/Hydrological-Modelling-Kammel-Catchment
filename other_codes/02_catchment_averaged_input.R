###################################################################################################
#### Preparing hydrometeorological time series: input variables for any hydrological model ########
###################################################################################################

#install.packages () # in case you do not have them installed
library(raster)
library(sf)
library(ncdf4)

setwd("C:/Users/tarasova/Desktop/UL/Class03")

############ Calculation of daily catchment averaged precipitation #############
###############  The effect of different interpolation methods  ################


####HYRAS dataset of German Weather Center (regression-based method) ############
#load nc file 
nc_file=nc_open("pr_hyras_1_2002_v6-0_de.nc") # from https://opendata.dwd.de/climate_environment/CDC/grids_germany/daily/hyras_de/precipitation/

#nc attributes
print(nc_file)

#extrat precip
precip_array <- ncvar_get(nc_file, "pr") 

fillvalue <- ncatt_get(nc_file, "pr", "_FillValue")

precip_array[precip_array==fillvalue$value] <- NA

#check dimensions
dim(precip_array)


time <- ncvar_get(nc_file, "time")
head(time)

#convert time: hours since
P_time_DWD <- as.Date(as.POSIXct(time*3600, origin = "1931-01-01",  tz="UTC")) #as.POSIXct works in seconds

head(P_time_DWD)

#extract coordinates
x=ncvar_get(nc_file, "x")

y=ncvar_get(nc_file, "y")

image(precip_array[,,1])

#make a brick of rasters
P_brick_DWD <- brick(precip_array, xmn=min(x), xmx=max(x), ymn=min(y), ymx=max(y), crs=crs("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs"))

plot(P_brick_DWD[[1]]) #looks weird...

#usually always have to flip because of difference in indexing cols and rows in nc and raster packages
P_brick_DWD <- flip(t(P_brick_DWD))

#reassign the extent
extent(P_brick_DWD)<-extent(min(x),max(x),min(y), max(y))

plot(P_brick_DWD[[1]]) # now this looks better

#read shapefile of the Selke catchment
Selke<-st_read("Selke_catchment_LHW.shp")

#crop germany-wide raster to the extent of our shp file
Selke_P_brick_DWD <- crop( x=P_brick_DWD, Selke[1], snap="out" )
Selke_P_brick_DWD <-mask( x= Selke_P_brick_DWD, Selke[1], snap="out" )

#plot raster for one of the days
plot(Selke_P_brick_DWD[[1]])

#daily time series of areal precipitation for the whole catchment
Selke_P_ts_DWD=cellStats(Selke_P_brick_DWD, stat='mean', na.rm=TRUE)

#annual sum in mm
sum(Selke_P_ts_DWD)

#plot of catchment-averaged precipitation
barplot(Selke_P_ts_DWD, names=P_time_DWD)

####################### External Drift Kriging #################################

#Now let's compare it with the precip and temperature obtained with external drift kriging
nc_file=nc_open("EDK/pre.nc") 

#nc attributes
print(nc_file)

#extract precip: I recommend to take a break while running the next nine lines
precip_array <- ncvar_get(nc_file, "pre") 

fillvalue <- ncatt_get(nc_file, "pre", "_FillValue")

precip_array[precip_array==fillvalue$value] <- NA

#check dimensions
dim(precip_array)

time <- ncvar_get(nc_file, "time")
head(time)

#convert time: hours since
P_time_EDK <- as.Date(time, origin = "1949-12-31")

#extract coordinates
lat=ncvar_get(nc_file, "lat")

lon=ncvar_get(nc_file, "lon")


#make a brick of rasters
P_brick_EDK <- brick(precip_array, xmn=min(lat), xmx=max(lat), ymn=min(lon), ymx=max(lon), crs=crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

#notice that it is not the same time step as for DWD
plot(P_brick_EDK[[1]]) #looks funny again...


#subset year 2002 first
P_brick_EDK<- subset(P_brick_EDK, which(P_time_EDK >= '2002-01-01' & (P_time_EDK <= '2002-12-31')))

#..and then transpose because it might take a while otherwise
P_brick_EDK <- t(P_brick_EDK)

plot(P_brick_EDK[[1]]) #should look like DE now 

#Do you remember how the plot from DWD looked like?

#plot them side by side
par(mfrow=c(1,2))
plot(P_brick_EDK[[1]], main="EDK") 
plot(P_brick_DWD[[1]], main="DWD") # now this looks better


#crop germany-wide raster to the extent of our shp file
#remember that nc file had a different projection
Selke_wgs84 <- st_transform(Selke, crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
Selke_P_brick_EDK <- crop( x=P_brick_EDK, Selke_wgs84[1], snap="out" )
Selke_P_brick_EDK <-mask( x= Selke_P_brick_EDK, Selke_wgs84[1], snap="out" )


#now this should be the same day as for the DWD
plot(Selke_P_brick_EDK[[1]])
plot(Selke_P_brick_DWD[[1]])



#daily time series of areal precipitation for the whole catchment
Selke_P_ts_EDK=cellStats(Selke_P_brick_EDK, stat='mean', na.rm=TRUE)

#annual sum in mm
sum(Selke_P_ts_EDK)
sum(Selke_P_ts_DWD)

#compare precip for 2002
plot(Selke_P_ts_EDK, Selke_P_ts_DWD, pch=16)
cor(Selke_P_ts_EDK, Selke_P_ts_DWD)

####################################################################################################
###"do it yourself task" for later: we will need these input variables for our hydrological model ##
####################################################################################################

## Basic inputs needed to run the model for the whole training+test period 19XX-2020
## 1. compute catchment-averaged precipitation for the period 19XX-2020 from DWD for you catchment X 
## from https://opendata.dwd.de/climate_environment/CDC/grids_germany/daily/hyras_de/precipitation/

## 2. compute catchment-averaged daily mean, max and min temperature for the period 19XX-2020 from DWD for you catchment X 
## https://opendata.dwd.de/climate_environment/CDC/grids_germany/daily/hyras_de/

## 3. compute catchment-averaged daily potential evaporation for the period 19XX-2020 for you catchment X 
## using mean, max and min air temperature from DWD and the Hargreves-Samani formula (Slide 8, Practice slides)

## Basic inputs needed for quantifying the input data uncertainty on model performance during training period 19XX-2010
## 4. compute catchment-averaged daily precipitation for the period 19XX-2010 for you catchment X 
#interpolated from stations using external drift kriging (from pre.nc)

## 5. compute catchment-averaged daily mean temperature for the period 19XX-2010 for you catchment X 
#interpolated from stations using external drift kriging (from tavg.nc)

## Start year for different catchments (due to streamflow data availability)
## Sieg at Menden1 1965
## Treene at Treia 1985
## Sieber at Hattorf 1951
## Kammel at Remshart 1983
## Unstrut at Laucha 1951
## Zschopau at Tannenberg 1961
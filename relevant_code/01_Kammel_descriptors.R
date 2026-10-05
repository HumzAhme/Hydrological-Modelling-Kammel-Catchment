library(raster)
library(sf)

# Set the working directory
setwd("C:/Uni/Hydrological Modeling")

#read shapefile
catchments<-st_read("./kammel_project/catchment.gpkg")

head(catchments)

plot(catchments[1])

#######  compute catchment area (also a very important catchment descriptor) ###
catchment_area<-st_area(catchments[1]) #in m2


######### Computing topographical catchment descriptors #######################
##################### Handling numerical variables ############################ 
##Task: calculate mean elevation of a catchment

#read raster file
dem<-raster("./descriptors/dem_100m_temp_etrs89_clipped_snapped.tif")

plot(dem)

catchments <- st_transform(catchments, crs = crs(dem))
catchments_sp <- as(catchments, "Spatial")
#crop germany-wide raster to the extent of our shp file
catchments_dem <- crop( x=dem, catchments, snap="out" )
plot(catchments_dem)

#mask
catchments_dem <- mask(catchments_dem, catchments)
plot(catchments_dem)

#Why do you think we had to mask?

#calculating mean elevation of our catchment
dem_mean<-cellStats(catchments_dem, "mean")

dem_mean

#Task: calculate mean slope of catchment 
#derive slope from digital elevation model
?terrain #raster package, check how it works
catchments_slope<-terrain(catchments_dem, opt="slope", units="radians") #degree calculation is incorrect
catchments_slope<-catchments_slope* 57.29578 #transform to degrees 180/pi

#in case you are interested in the theory on slope computation: https://pro.arcgis.com/en/pro-app/latest/tool-reference/spatial-analyst/how-slope-works.htm

plot(catchments_slope)

#Why do you think the slope is the steepest near the river? Should we always expect steeper slopes near the streams?

slope_mean<-cellStats(catchments_slope, "mean")

######### Computing catchment descriptors of surface cover #####################
##################### Handling categorical variables Part 1 ###########################

# read shp file
land_use <- st_read("./descriptors/CLC_germany_etrs89.shp")

# CRS an Catchments anpassen
land_use <- st_transform(land_use, st_crs(catchments))

# crop to extent of catchment
catchments_land_use <- st_crop(land_use["code06"], catchments)

# mask / actual intersection with catchment boundary
catchments_land_use <- st_intersection(catchments_land_use, catchments)

#remove unnecessary columns from the attribute table
catchments_land_use<-catchments_land_use[,-c(2:3)]

plot(catchments_land_use)

#re-classify land use
catchments_land_use$code06[catchments_land_use$code06>500]<- 5 #water
catchments_land_use$code06[catchments_land_use$code06>320]<- 4 #shrub
catchments_land_use$code06[catchments_land_use$code06>310]<- 3 #forest
catchments_land_use$code06[catchments_land_use$code06>200]<- 2 #agriculture
catchments_land_use$code06[catchments_land_use$code06>100]<- 1 #artificial

plot(catchments_land_use) #not very pretty map, for your catchment X play around to make the plot more appealing

#compute portion of catchment covered by forest
catchments_land_use$area_m2<-st_area(catchments_land_use) #area of all polygons
sum(catchments_land_use$area_m2[which(catchments_land_use$code06==3)]) #area of all forested polygons
forest_portion=sum(catchments_land_use$area_m2[which(catchments_land_use$code06==3)])/sum(catchments_land_use$area_m2)

forest_portion

#compute portions of all other land use classes
total_area <- sum(catchments_land_use$area_m2)


# Flächenanteile aller Klassen
landuse_portions <- c(
  artificial  = sum(catchments_land_use$area_m2[catchments_land_use$code06 == 1]) / total_area,
  agriculture = sum(catchments_land_use$area_m2[catchments_land_use$code06 == 2]) / total_area,
  forest      = sum(catchments_land_use$area_m2[catchments_land_use$code06 == 3]) / total_area,
  shrub       = sum(catchments_land_use$area_m2[catchments_land_use$code06 == 4]) / total_area,
  water       = sum(catchments_land_use$area_m2[catchments_land_use$code06 == 5]) / total_area
)

landuse_portions

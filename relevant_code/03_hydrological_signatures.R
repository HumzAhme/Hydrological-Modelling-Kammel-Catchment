###################################################################################################
#### Characterizing catchment response and streamflow dynamics: hydrological signatures ###########
###################################################################################################


### a great review paper to consult on the topic
#https://wires.onlinelibrary.wiley.com/doi/10.1002/wat2.1499

#install.packages () # in case you do not have them installed
library(sf)

setwd("C:/Users/lukas/Documents/Studium/02_semester/Hydrological_modeling/Kammel_model/hydrological-modelling-kammel-catchment")

#read file with hydrometeorological time series

ts <- read.csv("Kammel_discharge_PET_PREC_TEMP_1983_2020.csv")

#plot streamflow and precipitation
par(mfrow=c(1,1))
df.bar<-barplot(ts$precipitation_mm, names.arg = as.Date(ts$date,  "%m/%d/%Y"), xlab = "Date", ylab="precipitation mm")
lines(x=df.bar, y=ts$discharge_m3_s, col="blue")

#just the first year
df.bar<-barplot(ts$precipitation_mm[1:365], names.arg = as.Date(ts$date,  "%m/%d/%Y")[1:365], xlab = "Date", ylab="precipitation mm")
lines(x=df.bar, y=ts$discharge_m3_s[1:365], col="blue")

#are P and Q related at the daily scale
plot(ts$discharge_m3_s, ts$precipitation_mm)
cor(ts$discharge_m3_s, ts$precipitation_mm)

years=round(nrow(ts)/365)

sum(ts$discharge_m3_s)/years    #annual mean Q in m3s-1
sum(ts$precipitation_mm)/years  #annual mean P in mm/d
sum(ts$pet_mm)/years            #annual mean PET in mm/d

#Does Q>P makes sense?

#transform streamflow from m3/s to mm/day
#for that we need to know area of our catchment
#m3/m2 = m
#m --> 1000 mm
#s*24*60*60 --> day

#read shapefile of the Selke catchment
catchments<-st_read("Selke_catchment_LHW.shp")
#calculate area of the catchment 
area_m2 <- st_area(catchments)

ts$discharge_mm_d <- as.vector(ts$discharge_m3_s*1000*24*3600/area_m2)

#let's check
sum(ts$discharge_mm_d)/years    #annual mean Q in mm/d
sum(ts$precipitation_mm)/years  #annual mean P in mm/d

#plot them again, better?
plot(ts$discharge_mm_d, ts$precipitation_mm)
cor(ts$discharge_mm_d, ts$precipitation_mm) #linear transformation --> nothing changes except the amount

#plot their relation at monthly scale
monthlyQ=aggregate(ts$discharge_mm_d, by=list(months(as.Date(ts$date,  "%m/%d/%Y"))), mean)

order_months=unique(months(as.Date(ts$date,  "%m/%d/%Y")))

monthlyQ=monthlyQ[order(match(monthlyQ$Group.1, order_months)),]
barplot(monthlyQ$x, names.arg =monthlyQ$Group.1, ylab = "Q [mm]")

#let's add P and PET to the mix
monthlyP=aggregate(ts$precipitation_mm, by=list(months(as.Date(ts$date,  "%m/%d/%Y"))), mean)
monthlyP=monthlyP[order(match(monthlyP$Group.1, order_months)),]

monthlyPET=aggregate(ts$pet_mm, by=list(months(as.Date(ts$date,  "%m/%d/%Y"))), mean)
monthlyPET=monthlyPET[order(match(monthlyPET$Group.1, order_months)),]

#plot
df.bar <- barplot(monthlyP$x, names.arg =monthlyP$Group.1, ylab = "flux [mm]", ylim=c(0,4))
lines(x=df.bar, y=monthlyQ$x, col="blue")
lines(x=df.bar, y=monthlyPET$x, col="green")


#now let's plot the relationship of Q and P at annual scale
annualQ=aggregate(ts$discharge_mm_d, by=list(substr(as.Date(ts$date,  "%m/%d/%Y"), start=1, stop=4)), sum)
annualP=aggregate(ts$precipitation_mm, by=list(substr(as.Date(ts$date,  "%m/%d/%Y"), start=1, stop=4)), sum)

barplot(annualQ$x, names.arg =annualQ$Group.1, ylab = "Q [mm]")

plot(annualP$x, annualQ$x) # at annual scale there is much more coherence between P and Q 
cor(annualP$x, annualQ$x)

#do you think snow would be important in the Selke catchment?
#we can check it by looking at monthly temperature (persistent negative temperatures in winter would be an indicator for that)
monthlyT=aggregate(ts$temperature_mean, by=list(months(as.Date(ts$date,  "%m/%d/%Y"))), mean)
monthlyT=monthlyT[order(match(monthlyT$Group.1, order_months)),]
barplot(monthlyT$x, names.arg =monthlyT$Group.1, ylab = "deg C") # definitely not a snow-dominated catchment

#but we need to investigate more to be sure
min(ts$temperature_mean) #min recorded T
length(which(ts$temperature_mean<0)) #days below 0
percent_freeze_days=length(which(ts$temperature_mean<0))/nrow(ts)*100 
percent_freeze_days #it is worth to include snow routine in our model later

##### compute dryness index PET/P (Budyko method)
## Is your catchment energy-limited or water-limited
sum(ts$pet_mm)/sum(ts$precipitation_mm)  

###compute AET with water balance method
##AET = P - Q
AET=(sum(ts$precipitation_mm)-sum(ts$discharge_mm_d))/years

###let's check if the theoretical Budyko curve works..
AET/(sum(ts$precipitation_mm)/years)

## What about its hydroclimatic classification?
## Is your catchment arid, semi-arid or humid?
AI=sum(ts$precipitation_mm)/sum(ts$pet_mm)
AI

##### Another way to look at it:
##### runoff ratio
RR=sum(ts$discharge_mm_d)/sum(ts$precipitation_mm) # portion of precipitation that becomes streamflow
RR


###### flow duration curve ##################
#sorting Q in decreasing order
Qsort<-sort(ts$discharge_mm_d,decreasing=T)

#creating a data frame in which x column is hte percent of ie less than a specific time and
#y is the correspondent discharge.
df<-data.frame(x=100/length(Qsort)*1:length(Qsort),y=Qsort)

#log plot
plot(x = df$x, y = df$y, type = "l", log = "y",ylab="Discharge [mm/d]",xlab="Percentage of Time Flow is Equaled or Less Than (%)",main="Flow Duration Curve")
grid()
# x axis can help us to easy derive all flow percentiles Q10, Q50, Q70, Q95 etc

#slope of the flow duration curve (slope between Q33 and Q66)
#represent general variability of flow
slope_FDC <- (0.3292446 - 0.1484455)/(0.66-0.33)*100 #percent slope
slope_FDC

#coefficient of variation of streamflow
cvQ=var(ts$discharge_mm_d)/mean(ts$discharge_mm_d) #Selke is quite persistent
cvQ

##### before defining baseflow from total streamflow, let's jump back to catchment descriptors
#load aquifer types
#read shp file
aquifer_type<-st_read("aquifer_type.shp")
plot(aquifer_type["had16_hydr"])

#IDs
#1 water bodies
#2 porous
#3 fractured
#4 karstic
#5 aquitard

aquifer_type<- st_transform(aquifer_type, st_crs(catchments))

#crop, mask and examine
catchments_aquifer_type <- st_crop( x=aquifer_type["had16_hydr"], catchments[1], snap="out" )
catchments_aquifer_type <- st_intersection(catchments_aquifer_type, catchments[1])

#Do we expect a lot of groundwater contribution?
plot(catchments_aquifer_type["had16_hydr"])

#####  baseflow separation using a standard recursive digital filter: simple smoothing
#####  Low-pass filter
require(lfstat)
Qbase=baseflow(ts$discharge_mm_d, tp.factor=0.9, block.len=5)
plot(ts$discharge_mm_d[1:365], xlab = "Date", ylab="Streamflow mm/d", type = "l")
lines(Qbase[1:365], col="blue", lwd=2)

#compute baseflow index
BFI=sum(Qbase, na.rm=TRUE)/sum(ts$discharge_mm_d)
BFI  #0.60-0.80 is typical for the temperate climates

##the choice of filter and its parameters is very important for the identification of individual streamflow events
##it is less important for computing long-term baseflow contribution 
Qbase2=baseflow(ts$discharge_mm_d, tp.factor=0.9, block.len=7)
Qbase3=baseflow(ts$discharge_mm_d, tp.factor=0.7, block.len=7)

plot(ts$discharge_mm_d[1:365], xlab = "Date", ylab="Streamflow mm/d", type = "l")
lines(Qbase[1:365], col="blue", lwd=2)
lines(Qbase2[1:365], col="tomato", lwd=2)
lines(Qbase3[1:365], col="goldenrod", lwd=2)


#compare baseflow indeces
sum(Qbase, na.rm=TRUE)/sum(ts$discharge_mm_d)
sum(Qbase2, na.rm=TRUE)/sum(ts$discharge_mm_d)
sum(Qbase3, na.rm=TRUE)/sum(ts$discharge_mm_d)


#####  recession constant K (also often used as a in the baseflow filters)
require(lfstat)

Qlow_df <- createlfobj(ts(ts$discharge_mm_d), startdate=ts$date[1], baseflow =FALSE)

?recession
#only considers recession periods where the Q at the start of the recession does not exceed QX
recession(Qlow_df, method = "MRC",seglen = 7,threshold = 70) #Q70 - very low flows (default)
recession(Qlow_df, method = "MRC",seglen = 7,threshold = 50) #Q50 - medium range of flows is considered
recession(Qlow_df, method = "MRC",seglen = 7,threshold = 10) #Q10 - also very high values are considered

#tradeoff is needed: we want to consider only baseflow without any quickflow (this is the case for low flows)
#but we want to have as many points as possible for a robust fit (as many points as possible)

#for the case of Selke Q50 seems more reasonable

#let's check the effect of segment length
recession(Qlow_df, method = "MRC",seglen = 10,threshold = 50) # very long recessions (good), but too few 
recession(Qlow_df, method = "MRC",seglen = 7,threshold = 50) # default
recession(Qlow_df, method = "MRC",seglen = 5,threshold = 50) # some outlierts start to appear
recession(Qlow_df, method = "MRC",seglen = 2,threshold = 50) # even more scatter

K <- exp(-1/recession(Qlow_df, method = "MRC",seglen = 7,threshold = 50)) #transform from log

K #higher value corresponds to slower recession rates --> larger subsurface reservoirs

##################### do it yourself part #####################
### derive all hydrological signatures for your catchment X and compare it to the Selke catchment
### feel free to derive any additional hydrological signatures from this review paper: https://wires.onlinelibrary.wiley.com/doi/10.1002/wat2.1499
###################################################################################################
####                    Setting up a HBV hydrological model                             ###########
###################################################################################################

#install.packages()
library(HBV.IANIGLA)
library(hydroGOF)

#let's get familiar with HBV model structure


# set working directory to location of this script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

ts <- read.csv(
  "Kammel_all_descriptors_1983_2020.csv",
  sep = ",",
  stringsAsFactors = FALSE
)

# check original date format
head(ts$date)

# convert date
ts$date <- as.Date(ts$date, format = "%Y-%m-%d")

head(ts$date)
sum(is.na(ts$date))

# -----------------------------
# Rename columns for HBV model
# -----------------------------

ts$temperature_mean <- ts$Tmean_degC
ts$pet_mm <- ts$PET_mm

# discharge already exists as m3/s
# convert to mm/day

catchment_area_km2 <- 254

ts$discharge_mm_d <-
  ts$discharge_m3_s * 86400 /
  (catchment_area_km2 * 1e6) *
  1000



# quick check
summary(ts[, c(
  "temperature_mean",
  "precipitation_mm",
  "pet_mm",
  "discharge_mm_d"
)])



################ Individual modules of HBV ####################################

##################### degree-day-factor-based snow routine#####################
# precipitation partition into rainfall and snowfall
# snow accumulation as snow water equivalent and snowmelt

?SnowGlacier_HBV

snow_module <-
  SnowGlacier_HBV(model = 1, 
                  inputData = as.matrix( ts[ , c('temperature_mean', 'precipitation_mm')] ),
                  initCond = c(20, 2), #initial SWE and surface cover type (2: soil, alternative 1 for glacier cover)
                  #we start here with 20 mm. In reality we start with 0 mm and warm up the model for several years
                  param = c(1.00, 1.00, 0.00, 4) )
          
                  #SFCF: snowfall correction factor [−] --> we will switch off this parameter later
                  #Tr: solid and liquid precipitation threshold temperature [C].
                  #Tt: melt temperature [C].
                  #fm: snowmelt factor [mm/C.∆t]: open areas usually have lower areas than forested areas

head(snow_module)
#Prain: precip. as rainfall.
#Psnow: precip. as snowfall.
#SWE: snow water equivalent.
#Msnow: melted snow.
#Total: Prain + Msnow.

###################### soil moisture accounting scheme #######################
# fluxes: actual ET and runoff
# state: soil moisture

?Soil_HBV  #root zone soil routine

soil_module <-
  Soil_HBV(model = 1,
           inputData = cbind(snow_module[ , "Total"], ts$pet_mm),
           initCond = c(100, 1), #soil moisture initial conditions (assume some random number here, because we do not know better)
           param = c(250, 0.7, 2) )

           #FC: fictitious soil field capacity [mm].
           #LP: parameter to get actual ET [−].
           #β: exponential value that allows for non-linear relations between soil box
            #water input (rainfall plus snowmelt) and the effective runoff [−].

head(soil_module)
#Rech: runoff series [mm/∆t]. This is the input to the Routing_HBV module.
#Eact: actual evapotranspiration series [mm/∆t].
#SM: soil moisture series [mm/∆t]


###################### runoff redistribution module ###########################
# aka buckets

?Routing_HBV #subsurface flow redistribution routine

routing_module <-
  Routing_HBV(model = 3, #this can be modified to change the number of reservoirs & outflows: Model 3 corresponds to 2 buckets and 3 outlets
              lake = F,  #F for the cases when we do not have lakes as additional storage
              inputData = as.matrix(soil_module[ , "Rech"]),
              initCond = c(0, 0), #assuming it was empty
              param = c(0.2, 0.01, 0.005, 30, 0.015) )  


              #usually K0>>K1>>K2
              #K0: top bucket (STZ) storage constant [1/∆t].
              #K1: intermediate bucket (SUZ) storage constant [1/∆t].
              #K2: lower bucket (SLZ) storage constant [1/∆t].
              #UZL: minimum water content of SUZ for supplying fast runoff (Q0) to the total dischrage Qg [mm]
              #PERC: maximum flux rate between SUZ and SLZ [mm/∆t].

head(routing_module)
#Qg: total buckets output discharge [mm/∆t].
#Q0: top bucket discharge [mm/∆t].
#Q1: intermediate bucket discharge [mm/∆t].
#Q2: lower bucket discharge [mm/∆t].
#SUZ: intermediate reservoir storage [mm].
#SLZ: lower reservoir storage [mm]


####################### time-delay module #####################################
# aka delay within stream channel
?UH #transfer routine

tf_module <-
  round( 
    UH(model = 1,
       Qg = routing_module[ , "Qg"],
       param = c(1.5) ),
    2)

# let's plot the "true" and simulated hydrographs

graphics.off()
par(mfrow = c(1, 1))
par(mar = c(3, 3, 1, 1))

plot(
  x = ts$date,
  y = tf_module,
  type = "l",
  col = "dodgerblue",
  xlab = "Date",
  ylab = "Q(mm/d)"
)

lines(
  x = ts$date,
  y = ts$discharge_mm_d,
  col = "red"
)


####### Model as a single function #########################################

# simulated streamflow series.
hbv_lumped <- function(basin, 
                       param_snow,
                       param_soil,
                       param_routing,
                       param_tf, 
                       init_snow = 0, #we start with no snow. Best practice: warm up the model for several years
                       init_soil = 100, #these are random initial conditions
                       init_routing = c(0, 0) #we start with zero inital conditions. Best practice: warm up the model for several years

){
  
  #because we (I) are lazy we use the same names for the variables as in the Selke file. You will need to make sure that the names are the same in the files for your catchments
  snow_module <-
    SnowGlacier_HBV(model = 1, 
                    inputData = as.matrix( basin[ , c('temperature_mean', 'precipitation_mm')] ),      
                    initCond = c(init_snow, 2), 
                    param = c(1,param_snow) )  #this is how we can switch off a parameter --> hard-coding it as 1
  
  soil_module <-
    Soil_HBV(model = 1,
             inputData = cbind(snow_module[ , "Total"], basin$pet_mm),#also here check the variable name
             initCond = c(init_soil, 1),
             param = param_soil )
  
  routing_module <-
    Routing_HBV(model = 3, #3 buckets and 2 outlets
                lake = F, #no lake module
                inputData = as.matrix(soil_module[ , "Rech"]),
                initCond = init_routing, 
                param = param_routing )
  
  tf_module <-
    UH(model = 1,
       Qg = routing_module[ , "Qg"],
       param = param_tf )
  
  
  out <- round(tf_module, 2)
  
  return(out)
  
}


############## Manual calibration  ################################

streamflow <- 
  hbv_lumped(basin = ts, 
             param_snow = c(1.00, 0.00, 4),
             param_soil = c(250, 0.7, 2),
             param_routing = c(0.2, 0.01, 0.005, 30, 0.015),
             param_tf = c(1.5))

#play with parameters and see what happens

#whole time series
plot(
  x = ts$date,
  y = tf_module,
  type = "l",
  col = "dodgerblue",
  xlab = "Date",
  ylab = "Q(mm/d)"
)

lines(
  x = ts$date,
  y = ts$discharge_mm_d,
  col = "red"
)


#one year
start_day=5000 #change to plot different periods
par(mar=c(2, 4, 0.1, 0.1), mfrow=c(2,1))
#precipitation
plot(x = as.Date(ts[start_day:(start_day+365), "date"], "%m/%d/%Y"),
     y = ts[start_day:(start_day+365) , "precipitation_mm"],  
     type = "l", xaxt ="n", xlab="",
      col = "royalblue", ylim=c(max(ts[ , "precipitation_mm"]),0), ylab = "P(mm/d)")
#simulated
plot( x = as.Date(ts[start_day:(start_day+365), "date"], "%m/%d/%Y"),
       y = streamflow[start_day:(start_day+365)],
      type = "l", col = "dodgerblue",xlab = "Date",
      ylab = "Q(mm/d)")
#observed
lines(x = , as.Date(ts[start_day:(start_day+365) , "date"],  "%m/%d/%Y"),
      y = ts[start_day:(start_day+365) , "discharge_mm_d"], 
      col = "tomato")

legend("topleft", legend= c("obs", "sim"), col=c("tomato", "dodgerblue"), lty=1)


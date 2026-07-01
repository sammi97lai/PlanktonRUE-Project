################################################################################
##########################    Data preparation   ###############################
################################################################################

library(dplyr)
library(readxl)
library(tidyr)
library(lubridate)
library(brms)
library(corrplot)

rm(list = ls())
#------------------------------------------------------------------------------
setwd("~/GitHub/PlanktonRUE-Project/spatial_sensitivity_analysis/Larger_radius/RDATA")
load('plankton_RUE_data6_260526.RDATA')
### remove species biomass ....................................................
plankton_RUE_data7 <- plankton_RUE_data6 %>% select(-(c(16:84, 113:152))) #remove species biomass info

#### check distribution of variables ...........................................
#response
# hist(plankton_RUE_data7$genus_richness)
# hist(plankton_RUE_data7$FEve)
# hist(plankton_RUE_data7$FDis)
# hist(plankton_RUE_data7$cyano_prop)
plankton_RUE_data7$log_cyano_prop <- log(plankton_RUE_data7$cyano_prop+0.0001)
# hist(plankton_RUE_data7$log_cyano_prop)
# hist(plankton_RUE_data7$species_richness)
# hist(plankton_RUE_data7$zoo_FEve)
# hist(plankton_RUE_data7$zoo_FDis)
# hist(plankton_RUE_data7$log_RUE_A_P)
# hist(plankton_RUE_data7$log_RUE_A_N)
# hist(plankton_RUE_data7$log_RUE_zooplankton)


#predictors
# hist(plankton_RUE_data7$Temperature)
# hist(plankton_RUE_data7$Salinity)
# hist(plankton_RUE_data7$pH)

# hist(plankton_RUE_data7$CWM_log_Biovolume)
# hist(sqrt(plankton_RUE_data7$TotalPhosphorus_unfiltered))
# hist(log(plankton_RUE_data7$TotalNitrogen_unfiltered))

# plankton_RUE_data7$MaxDepth <- as.factor(plankton_RUE_data7$MaxDepth)

# hist(plankton_RUE_data7$MinDepth)
# hist(plankton_RUE_data7$MaxDepth)
# hist(plankton_RUE_data7$DOY)
# hist(plankton_RUE_data7$daylight_hours)
# hist(plankton_RUE_data7$Zooplankton_MaxDepth)

# plankton_RUE_data7$Zooplankton_MaxDepth <- as.factor(plankton_RUE_data7$Zooplankton_MaxDepth)

#test skewness
# library(moments)
# skewness(log(plankton_RUE_data7$RUE_C_N))
# skewness(sqrt(plankton_RUE_data7$MaxDepth))

plankton_RUE_data8 <- plankton_RUE_data7

save(plankton_RUE_data8, file = 'plankton_RUE_data8_260526.RDATA')

#### make everything scaled and add poly terms .................................
load('plankton_RUE_data8_260526.RDATA')

plankton_RUE_data_scaled <- plankton_RUE_data8 %>%
  mutate(across(c(Temperature, Salinity, pH, TotalNitrogen_unfiltered, 
                  MaxDepth, Zooplankton_MaxDepth, DOY, daylight_hours),
                ~ as.numeric(scale(.))))

save(plankton_RUE_data_scaled, file = 'plankton_RUE_data_scaled_260526.RDATA')

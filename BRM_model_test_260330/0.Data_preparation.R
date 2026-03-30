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
setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/RDATA")
load('plankton_RUE_data6_260309.RDATA')
### calculate N:P ratio........................................................
plankton_RUE_data7 <- plankton_RUE_data6 %>% mutate(N_P_ratio = TotalNitrogen_unfiltered / TotalPhosphorus_unfiltered)
colnames(plankton_RUE_data7)

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
# hist(plankton_RUE_data7$N_P_ratio)
plankton_RUE_data7$log_N_P_ratio <- log(plankton_RUE_data7$N_P_ratio)
# hist(plankton_RUE_data7$log_N_P_ratio)

# plankton_RUE_data7$MaxDepth <- as.factor(plankton_RUE_data7$MaxDepth)

# hist(plankton_RUE_data7$MinDepth)
# hist(plankton_RUE_data7$MaxDepth)
# hist(plankton_RUE_data7$DOY)
# hist(plankton_RUE_data7$Zooplankton_MaxDepth)

# plankton_RUE_data7$Zooplankton_MaxDepth <- as.factor(plankton_RUE_data7$Zooplankton_MaxDepth)

#test skewness
# library(moments)
# skewness(log(plankton_RUE_data7$RUE_C_N))
# skewness(sqrt(plankton_RUE_data7$MaxDepth))

save(plankton_RUE_data7, file = 'plankton_RUE_data7_260309.RDATA')

#### make everything scaled and add poly terms .................................
load('plankton_RUE_data7_260309.RDATA')

plankton_RUE_data_scaled <- plankton_RUE_data7 %>%
  mutate(across(c(Temperature, Salinity, pH, log_N_P_ratio, TotalNitrogen_unfiltered, 
                  MaxDepth, Zooplankton_MaxDepth, DOY, daylight_hours),
                ~ as.numeric(scale(.))))

save(plankton_RUE_data_scaled, file = 'plankton_RUE_data_scaled_260330.RDATA')

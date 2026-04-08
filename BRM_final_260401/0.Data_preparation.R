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
load('plankton_RUE_data7_260408.RDATA')
### remove species biomass ....................................................
plankton_RUE_data8 <- plankton_RUE_data7 %>% select(-(c(16:84, 113:152))) #remove species biomass info
### calculate N:P ratio........................................................
plankton_RUE_data8 <- plankton_RUE_data8 %>% mutate(N_P_ratio = TotalNitrogen_unfiltered / TotalPhosphorus_unfiltered)
colnames(plankton_RUE_data8)

#### check distribution of variables ...........................................
#response
# hist(plankton_RUE_data8$genus_richness)
# hist(plankton_RUE_data8$FEve)
# hist(plankton_RUE_data8$FDis)
# hist(plankton_RUE_data8$cyano_prop)
plankton_RUE_data8$log_cyano_prop <- log(plankton_RUE_data8$cyano_prop+0.0001)
# hist(plankton_RUE_data8$log_cyano_prop)
# hist(plankton_RUE_data8$species_richness)
# hist(plankton_RUE_data8$zoo_FEve)
# hist(plankton_RUE_data8$zoo_FDis)
# hist(plankton_RUE_data8$log_RUE_A_P)
# hist(plankton_RUE_data8$log_RUE_A_N)
# hist(plankton_RUE_data8$log_RUE_zooplankton)


#predictors
# hist(plankton_RUE_data8$Temperature)
# hist(plankton_RUE_data8$Salinity)
# hist(plankton_RUE_data8$pH)

# hist(plankton_RUE_data8$CWM_log_Biovolume)
# hist(sqrt(plankton_RUE_data8$TotalPhosphorus_unfiltered))
# hist(log(plankton_RUE_data8$TotalNitrogen_unfiltered))
# hist(plankton_RUE_data8$N_P_ratio)
plankton_RUE_data8$log_N_P_ratio <- log(plankton_RUE_data8$N_P_ratio)
# hist(plankton_RUE_data8$log_N_P_ratio)
# hist(plankton_RUE_data8$mean_edibility_weighted)
# hist(plankton_RUE_data8$mean_edibility_unweighted)

# plankton_RUE_data8$MaxDepth <- as.factor(plankton_RUE_data8$MaxDepth)

# hist(plankton_RUE_data8$MinDepth)
# hist(plankton_RUE_data8$MaxDepth)
# hist(plankton_RUE_data8$DOY)
# hist(plankton_RUE_data8$daylight_hours)
# hist(plankton_RUE_data8$Zooplankton_MaxDepth)

# plankton_RUE_data8$Zooplankton_MaxDepth <- as.factor(plankton_RUE_data8$Zooplankton_MaxDepth)

#test skewness
# library(moments)
# skewness(log(plankton_RUE_data8$RUE_C_N))
# skewness(sqrt(plankton_RUE_data8$MaxDepth))

save(plankton_RUE_data8, file = 'plankton_RUE_data8_260408.RDATA')

#### make everything scaled and add poly terms .................................
load('plankton_RUE_data8_260408.RDATA')

plankton_RUE_data_scaled <- plankton_RUE_data8 %>%
  mutate(across(c(Temperature, Salinity, pH, log_N_P_ratio, TotalNitrogen_unfiltered, 
                  MaxDepth, Zooplankton_MaxDepth, DOY, daylight_hours,
                  mean_edibility_weighted, mean_edibility_unweighted),
                ~ as.numeric(scale(.))))

save(plankton_RUE_data_scaled, file = 'plankton_RUE_data_scaled_260408.RDATA')

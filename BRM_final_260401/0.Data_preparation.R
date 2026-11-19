################################################################################
##########################    Data preparation   ###############################
################################################################################

library(dplyr)
library(tidyr)
library(brms)
library(corrplot)
library(car)

rm(list = ls())
#------------------------------------------------------------------------------
setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/RDATA")
load('plankton_RUE_data6_260916.RDATA')
### remove species biomass ....................................................
plankton_RUE_data7 <- plankton_RUE_data6 %>% select(-(c(16:84, 116:155))) #remove species biomass info

#### check distribution of variables ...........................................
#response
# hist(plankton_RUE_data7$genus_richness)
# hist(plankton_RUE_data7$shannon_phyto)
# hist(plankton_RUE_data7$phyto_evenness)
# hist(plankton_RUE_data7$phyto_FEve_PD)
# hist(plankton_RUE_data7$phyto_FDis_PD)
# hist(plankton_RUE_data7$phyto_FEve_ZD)
# hist(plankton_RUE_data7$phyto_FDis_ZD)
# hist(plankton_RUE_data7$phyto_FEve_PRUE)
# hist(plankton_RUE_data7$phyto_FDis_PRUE)
# hist(plankton_RUE_data7$phyto_FEve_ZRUE)
# hist(plankton_RUE_data7$phyto_FDis_ZRUE)
# hist(plankton_RUE_data7$cyano_prop)
# plankton_RUE_data7$log_cyano_prop <- log(plankton_RUE_data7$cyano_prop+0.0001)
# hist(plankton_RUE_data7$log_cyano_prop)
# hist(plankton_RUE_data7$species_richness)
# hist(plankton_RUE_data7$shannon_zoo)
# hist(plankton_RUE_data7$zoo_evenness)
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

save(plankton_RUE_data8, file = 'plankton_RUE_data8_260916.RDATA')

#### make everything scaled .................................
load('plankton_RUE_data8_260916.RDATA')

plankton_RUE_data_scaled <- plankton_RUE_data8 %>%
  mutate(across(c(Temperature, Salinity, pH, TotalNitrogen_unfiltered, 
                  MaxDepth, Zooplankton_MaxDepth, daylight_hours),
                ~ as.numeric(scale(.))))

save(plankton_RUE_data_scaled, file = 'plankton_RUE_data_scaled_260916.RDATA')

plankton_RUE_data_scaled_Zdiv <- plankton_RUE_data8 %>%
  mutate(across(c(Temperature, Salinity, pH, TotalNitrogen_unfiltered, 
                  MaxDepth, Zooplankton_MaxDepth, daylight_hours,
                  genus_richness, phyto_evenness, phyto_FEve_ZD, phyto_FDis_ZD,
                  cyano_prop),
                ~ as.numeric(scale(.))))

save(plankton_RUE_data_scaled_Zdiv, file = 'plankton_RUE_data_scaled_Zdiv_260916.RDATA')

plankton_RUE_data_scaled_RUE <- plankton_RUE_data8 %>%
  mutate(across(c(Temperature, Salinity, pH, TotalNitrogen_unfiltered, 
                  MaxDepth, Zooplankton_MaxDepth, daylight_hours,
                  genus_richness, phyto_evenness, 
                  phyto_FEve_PRUE, phyto_FDis_PRUE, phyto_FEve_ZRUE, phyto_FDis_ZRUE, cyano_prop,
                  species_richness, zoo_evenness, zoo_FEve, zoo_FDis),
                ~ as.numeric(scale(.))))

save(plankton_RUE_data_scaled_RUE, file = 'plankton_RUE_data_scaled_RUE_260916.RDATA')

#Correlation -------------------------------------------------------
load('plankton_RUE_data_scaled_260916.RDATA')
#check multicollinearity#.................................................
lm_check <- lm(log_RUE_B_P ~ genus_richness + phyto_FEve_PRUE + phyto_FDis_PRUE + 
                 species_richness + zoo_FEve + zoo_FDis +
                 Temperature + Salinity + pH + MaxDepth +
                 daylight_hours + cyano_prop, 
               data = plankton_RUE_data_scaled)
vif(lm_check) #below 3 is good#

lm_check <- lm(log_RUE_B_N ~ genus_richness + phyto_FEve_PRUE + phyto_FDis_PRUE + 
                 species_richness + zoo_FEve + zoo_FDis +
                 Temperature + Salinity + pH + MaxDepth +
                 daylight_hours + cyano_prop, 
               data = plankton_RUE_data_scaled)
vif(lm_check) #below 3 is good# 

lm_check <- lm(log_RUE_zooplankton ~ genus_richness + phyto_FEve_ZRUE + phyto_FDis_ZRUE +
                 species_richness + zoo_FEve + zoo_FDis +
                 Temperature + Salinity + pH + Zooplankton_MaxDepth +
                 daylight_hours + cyano_prop, 
               data = plankton_RUE_data_scaled)
vif(lm_check) #below 3 is good#

plankton_RUE_data_scaled$MaxDepth <- as.numeric(plankton_RUE_data_scaled$MaxDepth)
plankton_RUE_data_scaled$Zooplankton_MaxDepth <- as.numeric(plankton_RUE_data_scaled$Zooplankton_MaxDepth)
data <- plankton_RUE_data_scaled %>% select(genus_richness, shannon_phyto, phyto_evenness, 
                                            phyto_FEve_PD, phyto_FDis_PD, 
                                            phyto_FEve_ZD, phyto_FDis_ZD,
                                            phyto_FEve_PRUE, phyto_FDis_PRUE,
                                            phyto_FEve_ZRUE, phyto_FDis_ZRUE, 
                                            cyano_prop,
                                            species_richness, shannon_zoo, zoo_evenness, 
                                            zoo_FEve, zoo_FDis, 
                                            Temperature, Salinity, pH, 
                                            TotalNitrogen_unfiltered, TotalPhosphorus_unfiltered,
                                            daylight_hours, Zooplankton_MaxDepth, MaxDepth,
                                            log_RUE_B_P, log_RUE_B_N, log_RUE_zooplankton)%>%
  
  rename(
    "Phytoplankton Genus Richness" = genus_richness,
    "Phytoplankton Shannon Diversity" = shannon_phyto,
    "Phytoplankton Evenness" = phyto_evenness,
    "PD Phytoplankton FEve" = phyto_FEve_PD,
    "PD Phytoplankton FDis" = phyto_FDis_PD,
    "ZD Phytoplankton FEve" = phyto_FEve_ZD,
    "ZD Phytoplankton FDis" = phyto_FDis_ZD,
    "PRUE Phytoplankton FEve" = phyto_FEve_PRUE,
    "PRUE Phytoplankton FDis" = phyto_FDis_PRUE,
    "ZRUE Phytoplankton FEve" = phyto_FEve_ZRUE,
    "ZRUE Phytoplankton FDis" = phyto_FDis_ZRUE,
    "Cyanobacteria Dominance" = cyano_prop,
    "Zooplankton Species Richness" = species_richness,
    "Zooplankton Shannon Diversity" = shannon_zoo,
    "Zooplankton Evenness" = zoo_evenness,
    "Zooplankton FEve" = zoo_FEve,
    "Zooplankton FDis" = zoo_FDis,
    "Nitrogen concentrations" = TotalNitrogen_unfiltered,
    "Phosphorus concentrations" = TotalPhosphorus_unfiltered,
    "Daylight Hours" = daylight_hours,
    "Phytoplankton Depth" = MaxDepth,
    "Zooplankton Depth" = Zooplankton_MaxDepth,
    "RUE per P" = log_RUE_B_P,
    "RUE per N" = log_RUE_B_N,
    "Zooplankton RUE" = log_RUE_zooplankton
  )

cor_mat <- cor(data)

#plot
setwd("~/GitHub/PlanktonRUE-Project/Figures")
png("corrplot.png", width = 6500, height = 4500, res = 300)
corrplot(cor_mat, method= "circle", type="lower", addCoef.col = "black", 
         tl.col = "black",number.cex = 0.9, tl.cex = 0.8)  
dev.off()

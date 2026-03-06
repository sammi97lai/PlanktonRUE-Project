library(dplyr)
library(readxl)
library(tidyr)
library(lubridate)

############# Calculating Resource Use Efficiency (RUE) ######################

setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/RDATA")

rm(list = ls())
load('filtered_biomass_env_data6_180925.RDATA')

#calculate RUE for total community - biomass #----------------------------------
#calculate total biomass#
RUE_biomass_com<-filtered_biomass_env_data6
RUE_biomass_com$Btot <- rowSums(RUE_biomass_com[, 18:93]) 

#calculate RUE of biomass per unit phosphorus
RUE_biomass_com <- RUE_biomass_com %>%
  mutate(RUE_B_P = Btot / TotalPhosphorus_unfiltered)
RUE_biomass_com$log_RUE_B_P <- log(RUE_biomass_com$RUE_B_P)

#calculate RUE of biomass per unit nitrogen
RUE_biomass_com <- RUE_biomass_com %>%
  mutate(RUE_B_N = Btot / TotalNitrogen_unfiltered)
RUE_biomass_com$log_RUE_B_N <- log(RUE_biomass_com$RUE_B_N)

#formatting dataset#
RUE_biomass_com$Site_Id <- as.factor(RUE_biomass_com$Site_Id)
RUE_biomass_com$Date = as.Date(RUE_biomass_com$Date) #get DOY
RUE_biomass_com$DOY = yday(RUE_biomass_com$Date)
RUE_biomass_com$Year = year(RUE_biomass_com$Date)

RUE_biomass_com <- RUE_biomass_com %>%
  drop_na(Temperature, Salinity, pH,
          TotalPhosphorus_unfiltered, TotalNitrogen_unfiltered, 
          MinDepth, MaxDepth, DOY, Year,
          Hill1, Hill2, rwr, genus_richness, log_G,
          FEve, FDis, cyano_prop, 
          RUE_B_P, RUE_B_N, log_RUE_B_P, log_RUE_B_N
  )  # Removes na rows 

save(RUE_biomass_com, file = 'RUE_biomass_com_180925.RDATA')

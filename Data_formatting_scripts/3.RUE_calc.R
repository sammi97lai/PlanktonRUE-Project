library(dplyr)
library(readxl)
library(tidyr)
library(lubridate)

############# Calculating Resource Use Efficiency (RUE) ######################

setwd("~/GitHub/Phytoplankton-project/phyto_RUE/dataformatting")

rm(list = ls())
# load('filtered_carbon_env_data5_010925.RDATA')
load('filtered_biomass_env_data6_180925.RDATA')
load('filtered_abd_env_data6_290925.RDATA')

#calculate RUE per genus#----------------------------------------------------
# RUE_carbon<-filtered_carbon_env_data5
# 
# #formatting dataset#
# RUE_carbon$Site_Id <- as.factor(RUE_carbon$Site_Id)
# RUE_carbon$Date = as.Date(RUE_carbon$Date) #get DOY
# RUE_carbon$DOY = yday(RUE_carbon$Date)
# RUE_carbon$Year = year(RUE_carbon$Date)
# 
# #calculate RUE of carbon content per unit phosphorus
# RUE_C_P <- RUE_carbon %>%
#   mutate(across(18:93, ~ (.x / TotalPhosphorus_unfiltered)))%>%
#   drop_na()  # Removes na rows 
# 
# #RUE_C_P[RUE_C_P == -Inf] <- 0
# 
# RUE_C_Ps <- RUE_C_P %>%
#   mutate(across(c(Temperature, Salinity, pH, TotalPhosphorus_unfiltered, 
#                   TotalNitrogen_unfiltered, MinDepth, MaxDepth),
#                 ~ as.numeric(scale(.))))
# 
# save(RUE_C_P, file = 'RUE_C_P.RDATA')
# save(RUE_C_Ps, file = 'RUE_C_P_scaled.RDATA')
# 
# #calculate RUE of carbon content per unit nitrogen
# RUE_C_N <- RUE_carbon %>%
#   mutate(across(17:92, ~ (.x / TotalNitrogen_unfiltered)))%>%
#   drop_na()  # Removes na rows 
# 
# #RUE_C_N[RUE_C_N == -Inf] <- 0
# 
# RUE_C_Ns <- RUE_C_N %>%
#   mutate(across(c(Temperature, Salinity, pH, TotalPhosphorus_unfiltered, 
#                   TotalNitrogen_unfiltered, MinDepth, MaxDepth),
#                 ~ as.numeric(scale(.))))
# 
# save(RUE_C_N, file = 'RUE_C_N.RDATA')
# save(RUE_C_Ns, file = 'RUE_C_N_scaled.RDATA')


#calculate RUE for total community#------------------------------------
#calculate total carbon#
RUE_carbon_com<-filtered_carbon_env_data5
RUE_carbon_com$Ctot <- rowSums(RUE_carbon_com[, 18:93]) 

#calculate RUE of carbon content per unit phosphorus
RUE_carbon_com <- RUE_carbon_com %>%
  mutate(RUE_C_P = Ctot / TotalPhosphorus_unfiltered)
RUE_carbon_com$log_RUE_C_P <- log(RUE_carbon_com$RUE_C_P)

#calculate RUE of carbon content per unit nitrogen
RUE_carbon_com <- RUE_carbon_com %>%
  mutate(RUE_C_N = Ctot / TotalNitrogen_unfiltered)
RUE_carbon_com$log_RUE_C_N <- log(RUE_carbon_com$RUE_C_N)

#formatting dataset#
RUE_carbon_com$Site_Id <- as.factor(RUE_carbon_com$Site_Id)
RUE_carbon_com$Date = as.Date(RUE_carbon_com$Date) #get DOY
RUE_carbon_com$DOY = yday(RUE_carbon_com$Date)
RUE_carbon_com$Year = year(RUE_carbon_com$Date)

RUE_carbon_com <- RUE_carbon_com %>%
  drop_na(Temperature, Salinity, pH,
          TotalPhosphorus_unfiltered, TotalNitrogen_unfiltered, 
          MinDepth, MaxDepth, DOY, Year,
          Hill1, Hill2, rwr, genus_richness, log_G,
          FEve, FDis, CWM_log_Biovolume, cyano_prop, 
          RUE_C_P, RUE_C_N, log_RUE_C_P, log_RUE_C_N
  )  # Removes na rows 

RUE_carbon_com_scaled <- RUE_carbon_com %>%
  mutate(across(c(Temperature, Salinity, pH, TotalPhosphorus_unfiltered, 
                  TotalNitrogen_unfiltered, MinDepth, MaxDepth),
                ~ as.numeric(scale(.))))

save(RUE_carbon_com, file = 'RUE_carbon_com_010925.RDATA')
save(RUE_carbon_com_scaled, file = 'RUE_carbon_com_scaled_010925.RDATA')

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
          FEve, FDis, CWM_log_Biovolume, cyano_prop, 
          RUE_B_P, RUE_B_N, log_RUE_B_P, log_RUE_B_N
  )  # Removes na rows 

RUE_biomass_com_scaled <- RUE_biomass_com %>%
  mutate(across(c(Temperature, Salinity, pH, TotalPhosphorus_unfiltered, 
                  TotalNitrogen_unfiltered, MinDepth, MaxDepth),
                ~ as.numeric(scale(.))))

save(RUE_biomass_com, file = 'RUE_biomass_com_180925.RDATA')
save(RUE_biomass_com_scaled, file = 'RUE_biomass_com_scaled_180925.RDATA')

##calculate RUE for total community - abundance #-------------------------------
#calculate total abundance#
RUE_abundance_com<-filtered_abd_env_data6
RUE_abundance_com$Atot <- rowSums(RUE_abundance_com[, 18:93]) 

#calculate RUE of abundance per unit phosphorus
RUE_abundance_com <- RUE_abundance_com %>%
  mutate(RUE_A_P = Atot / TotalPhosphorus_unfiltered)
RUE_abundance_com$log_RUE_A_P <- log(RUE_abundance_com$RUE_A_P)

#calculate RUE of abundance per unit nitrogen
RUE_abundance_com <- RUE_abundance_com %>%
  mutate(RUE_A_N = Atot / TotalNitrogen_unfiltered)
RUE_abundance_com$log_RUE_A_N <- log(RUE_abundance_com$RUE_A_N)

#formatting dataset#
RUE_abundance_com$Site_Id <- as.factor(RUE_abundance_com$Site_Id)
RUE_abundance_com$Date = as.Date(RUE_abundance_com$Date) #get DOY
RUE_abundance_com$DOY = yday(RUE_abundance_com$Date)
RUE_abundance_com$Year = year(RUE_abundance_com$Date)

RUE_abundance_com <- RUE_abundance_com %>%
  drop_na(Temperature, Salinity, pH,
          TotalPhosphorus_unfiltered, TotalNitrogen_unfiltered, 
          MinDepth, MaxDepth, DOY, Year,
          Hill1, Hill2, rwr, genus_richness, log_G,
          FEve, FDis, CWM_log_Biovolume, cyano_prop, 
          RUE_A_P, RUE_A_N, log_RUE_A_P, log_RUE_A_N
  )  # Removes na rows 

RUE_abundance_com_scaled <- RUE_abundance_com %>%
  mutate(across(c(Temperature, Salinity, pH, TotalPhosphorus_unfiltered, 
                  TotalNitrogen_unfiltered, MinDepth, MaxDepth),
                ~ as.numeric(scale(.))))

save(RUE_abundance_com, file = 'RUE_abundance_com_290925.RDATA')
save(RUE_abundance_com_scaled, file = 'RUE_abundance_com_scaled_290925.RDATA')
################################################################################
######################    RUE Zooplankton and div BRMS   #######################
################################################################################

library(dplyr)
library(readxl)
library(tidyr)
library(lubridate)
library(brms)

rm(list = ls())
#------------------------------------------------------------------------------
setwd("~/GitHub/PlanktonRUE-Project/spatial_sensitivity_analysis/Smaller_radius/RDATA")
load('plankton_RUE_data_scaled_260526.RDATA')

setwd("~/GitHub/PlanktonRUE-Project/spatial_sensitivity_analysis/Smaller_radius/RDATA")

#make BRMS model - Zooplankton diversity with CYANO DOM # ----------------------
#Zooplankton Species Richness#..................................................

start_time <- Sys.time()

model <- brm(formula = species_richness ~ poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) +
               poly(TotalNitrogen_unfiltered,2,raw = TRUE) + poly(Zooplankton_MaxDepth,2,raw = TRUE) +
               genus_richness + FEve + FDis + cyano_prop +
               s(DOY) + poly(daylight_hours,2,raw = TRUE) +
               (1 | Year) + (1 | Site_Id), 
             data = plankton_RUE_data_scaled,
             iter = 4000, warmup = 2000, cores = 4, chains = 4,
             save_pars = save_pars(all = TRUE),
             prior = c(
               set_prior("normal(0, 1)", class = "b"),
               set_prior("exponential(1)", class = "sd"),
               set_prior("exponential(1)", class = "sigma")
             ),
             control = list(adapt_delta = 0.99, max_treedepth = 15),
             seed = 2345
)

save(model, file = "BRM_zoo_SR_CD_P_260526.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#Functional Evenness ...........................................................

start_time <- Sys.time()

model <- brm(formula = zoo_FEve ~ poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) +
               poly(TotalNitrogen_unfiltered,2,raw = TRUE) + poly(Zooplankton_MaxDepth,2,raw = TRUE) +
               genus_richness + FEve + FDis + cyano_prop +
               s(DOY) + poly(daylight_hours,2,raw = TRUE) +
               (1 | Year) + (1 | Site_Id), 
             data = plankton_RUE_data_scaled,
             iter = 4000, warmup = 2000, cores = 4, chains = 4, 
             save_pars = save_pars(all = TRUE),
             prior = c(
               set_prior("normal(0, 1)", class = "b"),
               set_prior("exponential(1)", class = "sd"),
               set_prior("exponential(1)", class = "sigma")
             ),
             control = list(adapt_delta = 0.99, max_treedepth = 15), 
             seed = 2345
)

save(model, file = "BRM_zoo_FE_CD_P_260526.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#Functional Dispersion .........................................................

start_time <- Sys.time()

model <- brm(formula = zoo_FDis ~ poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) +
               poly(TotalNitrogen_unfiltered,2,raw = TRUE) + poly(Zooplankton_MaxDepth,2,raw = TRUE) +
               genus_richness + FEve + FDis + cyano_prop +
               s(DOY) + poly(daylight_hours,2,raw = TRUE) +
               (1 | Year) + (1 | Site_Id), 
             data = plankton_RUE_data_scaled,
             iter = 4000, warmup = 2000, cores = 4, chains = 4,
             save_pars = save_pars(all = TRUE),
             prior = c(
               set_prior("normal(0, 1)", class = "b"),
               set_prior("exponential(1)", class = "sd"),
               set_prior("exponential(1)", class = "sigma")
             ),
             control = list(adapt_delta = 0.99, max_treedepth = 15),
             seed = 2345
)

save(model, file = "BRM_zoo_FD_CD_P_260526.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#make BRMS model- zooplankton RUE with CYANO DOM#-------------------------------

start_time <- Sys.time()

model <- brm(formula = log_RUE_zooplankton ~ 
               #env variables
               Temperature + I(Temperature^2) + Salinity + I(Salinity^2) + 
               pH + I(pH^2) + Zooplankton_MaxDepth + I(Zooplankton_MaxDepth^2) + 
               TotalNitrogen_unfiltered + I(TotalNitrogen_unfiltered^2) + s(DOY) +
               daylight_hours + I(daylight_hours^2) + 
               # phytoplankton biodiversity
               genus_richness * (Temperature + Salinity + TotalNitrogen_unfiltered) +
               FEve * (Temperature + Salinity + TotalNitrogen_unfiltered) +
               FDis * (Temperature + Salinity + TotalNitrogen_unfiltered) +
               cyano_prop * (Temperature + Salinity + TotalNitrogen_unfiltered) +
               # zooplankton biodiversity
               species_richness * (Temperature + Salinity + TotalNitrogen_unfiltered) +
               zoo_FEve * (Temperature + Salinity + TotalNitrogen_unfiltered) +
               zoo_FDis * (Temperature + Salinity + TotalNitrogen_unfiltered) +
               # random effects
               (1 | Year) + (1 | Site_Id),
             data = plankton_RUE_data_scaled,
             iter = 4000, warmup = 2000, cores = 4, chains = 4,
             save_pars = save_pars(all = TRUE),
             prior = c(
               set_prior("normal(0, 1)", class = "b"),
               set_prior("exponential(1)", class = "sd"),
               set_prior("exponential(1)", class = "sigma")
             ),
             control = list(adapt_delta = 0.99, max_treedepth = 15),
             seed = 2345
)

save(model, file = "BRM_zoo_RUE_CD_PI_260526.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)


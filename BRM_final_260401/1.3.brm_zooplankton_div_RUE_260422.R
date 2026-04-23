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
setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/RDATA")
load('plankton_RUE_data_scaled_260408.RDATA')

setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA")

#make BRMS model - Zooplankton diversity # -------------------------------------
#Zooplankton Species Richness #..................................................
start_time <- Sys.time()

model <- brm(formula = species_richness ~ Temperature + Salinity + pH +
               TotalNitrogen_unfiltered + Zooplankton_MaxDepth +
               genus_richness + FEve_selected + FDis_selected +
               DOY + daylight_hours +
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

save(model, file = "BRM_zoo_SR_L_260408.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#...................
start_time <- Sys.time()

model <- brm(formula = species_richness ~ poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) +
               poly(TotalNitrogen_unfiltered,2,raw = TRUE) + poly(Zooplankton_MaxDepth,2,raw = TRUE) +
               genus_richness + FEve_selected + FDis_selected +
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

save(model, file = "BRM_zoo_SR_P_260408.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)


#Functional Evenness ...........................................................
start_time <- Sys.time()

model <- brm(formula = zoo_FEve ~ Temperature + Salinity + pH +
               TotalNitrogen_unfiltered + Zooplankton_MaxDepth +
               genus_richness + FEve_selected + FDis_selected +
               DOY + daylight_hours +
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

save(model, file = "BRM_zoo_FE_L_260408.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)


#...................
start_time <- Sys.time()

model <- brm(formula = zoo_FEve ~ poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) +
               poly(TotalNitrogen_unfiltered,2,raw = TRUE) + poly(Zooplankton_MaxDepth,2,raw = TRUE) +
               genus_richness + FEve_selected + FDis_selected +
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

save(model, file = "BRM_zoo_FE_P_260408.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#Functional Dispersion .........................................................
start_time <- Sys.time()

model <- brm(formula = zoo_FDis ~ Temperature + Salinity + pH +
               TotalNitrogen_unfiltered + Zooplankton_MaxDepth +
               genus_richness + FEve_selected + FDis_selected +
               DOY + daylight_hours +
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

save(model, file = "BRM_zoo_FD_L_260408.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#...................
start_time <- Sys.time()

model <- brm(formula = zoo_FDis ~ poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) +
               poly(TotalNitrogen_unfiltered,2,raw = TRUE) + poly(Zooplankton_MaxDepth,2,raw = TRUE) +
               genus_richness + FEve_selected + FDis_selected +
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

save(model, file = "BRM_zoo_FD_P_260408.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)


#make BRMS model - Zooplankton diversity with CYANO DOM # ----------------------
#Zooplankton Species Richness#..................................................
start_time <- Sys.time()

model <- brm(formula = species_richness ~ Temperature + Salinity  + pH +
               TotalNitrogen_unfiltered + Zooplankton_MaxDepth +
               genus_richness + FEve_selected + FDis_selected + cyano_prop +
               DOY + daylight_hours +
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

save(model, file = "BRM_zoo_SR_CD_L_260408.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#...................
start_time <- Sys.time()

model <- brm(formula = species_richness ~ poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) +
               poly(TotalNitrogen_unfiltered,2,raw = TRUE) + poly(Zooplankton_MaxDepth,2,raw = TRUE) +
               genus_richness + FEve_selected + FDis_selected + cyano_prop +
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

save(model, file = "BRM_zoo_SR_CD_P_260408.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#Functional Evenness ...........................................................
start_time <- Sys.time()

model <- brm(formula = zoo_FEve ~ Temperature + Salinity  + pH +
               TotalNitrogen_unfiltered + Zooplankton_MaxDepth +
               genus_richness + FEve_selected + FDis_selected + cyano_prop +
               DOY + daylight_hours +
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

save(model, file = "BRM_zoo_FE_CD_L_260408.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#...................
start_time <- Sys.time()

model <- brm(formula = zoo_FEve ~ poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) +
               poly(TotalNitrogen_unfiltered,2,raw = TRUE) + poly(Zooplankton_MaxDepth,2,raw = TRUE) +
               genus_richness + FEve_selected + FDis_selected + cyano_prop +
               s(DOY) + poly(daylight_hours,2,raw = TRUE) +
               (1 | Year) + (1 | Site_Id), 
             data = plankton_RUE_data_scaled,
             iter = 4000, warmup = 3000, cores = 4, chains = 4, 
             save_pars = save_pars(all = TRUE),
             prior = c(
               set_prior("normal(0, 1)", class = "b"),
               set_prior("exponential(1)", class = "sd"),
               set_prior("exponential(1)", class = "sigma")
             ),
             control = list(adapt_delta = 0.99, max_treedepth = 15), 
             seed = 2345
)

save(model, file = "BRM_zoo_FE_CD_P_260408.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#Functional Dispersion .........................................................
start_time <- Sys.time()

model <- brm(formula = zoo_FDis ~ Temperature + Salinity  + pH +
               TotalNitrogen_unfiltered + Zooplankton_MaxDepth +
               genus_richness + FEve_selected + FDis_selected + cyano_prop +
               DOY + daylight_hours +
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

save(model, file = "BRM_zoo_FD_CD_L_260408.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#...................
start_time <- Sys.time()

model <- brm(formula = zoo_FDis ~ poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) +
               poly(TotalNitrogen_unfiltered,2,raw = TRUE) + poly(Zooplankton_MaxDepth,2,raw = TRUE) +
               genus_richness + FEve_selected + FDis_selected + cyano_prop +
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

save(model, file = "BRM_zoo_FD_CD_P_260408.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#make BRMS model- zooplankton RUE#----------------------------------------------
start_time <- Sys.time()

model <- brm(formula = log_RUE_zooplankton ~ Temperature + Salinity  + pH +
               TotalNitrogen_unfiltered + Zooplankton_MaxDepth +
               genus_richness + FEve_selected + FDis_selected + 
               species_richness + zoo_FEve + zoo_FDis + 
               DOY + daylight_hours +
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

save(model, file = "BRM_zoo_RUE_L_260422.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#.........................

start_time <- Sys.time()

model <- brm(formula = log_RUE_zooplankton ~ poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) +
               poly(TotalNitrogen_unfiltered,2,raw = TRUE) + poly(Zooplankton_MaxDepth,2,raw = TRUE) +
               genus_richness + FEve_selected + FDis_selected + 
               species_richness + zoo_FEve + zoo_FDis + 
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

save(model, file = "BRM_zoo_RUE_P_260422.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#...........................

start_time <- Sys.time()

model <- brm(formula = log_RUE_zooplankton ~ pH + Zooplankton_MaxDepth +
               DOY + daylight_hours + 
               genus_richness*Temperature + genus_richness*Salinity + genus_richness*TotalNitrogen_unfiltered +
               FEve_selected*Temperature + FEve_selected*Salinity + FEve_selected*TotalNitrogen_unfiltered + 
               FDis_selected*Temperature + FDis_selected*Salinity + FDis_selected*TotalNitrogen_unfiltered +
               species_richness*Temperature + species_richness*Salinity + species_richness*TotalNitrogen_unfiltered +
               zoo_FEve*Temperature + zoo_FEve*Salinity + zoo_FEve*TotalNitrogen_unfiltered +
               zoo_FDis*Temperature + zoo_FDis*Salinity + zoo_FDis*TotalNitrogen_unfiltered +
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

save(model, file = "BRM_zoo_RUE_LI_260422.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#.........................

start_time <- Sys.time()

model <- brm(formula = log_RUE_zooplankton ~ 
               #env variables
               pH + I(pH^2) + Zooplankton_MaxDepth + I(Zooplankton_MaxDepth^2) + s(DOY) +
               daylight_hours + I(daylight_hours^2) + Temperature + I(Temperature^2) +
               Salinity + I(Salinity^2) + TotalNitrogen_unfiltered + I(TotalNitrogen_unfiltered^2) +
               # phytoplankton biodiversity
               genus_richness * (Temperature + Salinity + TotalNitrogen_unfiltered) +
               FEve_selected * (Temperature + Salinity + TotalNitrogen_unfiltered) +
               FDis_selected * (Temperature + Salinity + TotalNitrogen_unfiltered) +
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

save(model, file = "BRM_zoo_RUE_PI_260422.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)


#make BRMS model- zooplankton RUE with CYANO DOM#-------------------------------
start_time <- Sys.time()

model <- brm(formula = log_RUE_zooplankton ~ Temperature + Salinity  + pH +
               TotalNitrogen_unfiltered + Zooplankton_MaxDepth +
               genus_richness + FEve_selected + FDis_selected + cyano_prop +
               species_richness + zoo_FEve + zoo_FDis + 
               DOY + daylight_hours +
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

save(model, file = "BRM_zoo_RUE_CD_L_260422.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#....................
start_time <- Sys.time()

model <- brm(formula = log_RUE_zooplankton ~ poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) +
               poly(TotalNitrogen_unfiltered,2,raw = TRUE) + poly(Zooplankton_MaxDepth,2,raw = TRUE) +
               genus_richness + FEve_selected + FDis_selected + cyano_prop +
               species_richness + zoo_FEve + zoo_FDis + 
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

save(model, file = "BRM_zoo_RUE_CD_P_260422.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#.......................

start_time <- Sys.time()

model <- brm(formula = log_RUE_zooplankton ~ pH + Zooplankton_MaxDepth +
               DOY + daylight_hours +
               genus_richness*Temperature + genus_richness*Salinity + genus_richness*TotalNitrogen_unfiltered +
               FEve_selected*Temperature + FEve_selected*Salinity + FEve_selected*TotalNitrogen_unfiltered +
               FDis_selected*Temperature + FDis_selected*Salinity + FDis_selected*TotalNitrogen_unfiltered +
               cyano_prop*Temperature + cyano_prop*Salinity + cyano_prop*TotalNitrogen_unfiltered +
               species_richness*Temperature + species_richness*Salinity + species_richness*TotalNitrogen_unfiltered +
               zoo_FEve*Temperature + zoo_FEve*Salinity + zoo_FEve*TotalNitrogen_unfiltered +
               zoo_FDis*Temperature + zoo_FDis*Salinity + zoo_FDis*TotalNitrogen_unfiltered +
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

save(model, file = "BRM_zoo_RUE_CD_LI_260422.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#....................
start_time <- Sys.time()

model <- brm(formula = log_RUE_zooplankton ~ 
               #env variables
               Temperature + I(Temperature^2) + Salinity + I(Salinity^2) + 
               pH + I(pH^2) + Zooplankton_MaxDepth + I(Zooplankton_MaxDepth^2) + 
               TotalNitrogen_unfiltered + I(TotalNitrogen_unfiltered^2) + s(DOY) +
               daylight_hours + I(daylight_hours^2) + 
               # phytoplankton biodiversity
               genus_richness * (Temperature + Salinity + TotalNitrogen_unfiltered) +
               FEve_selected * (Temperature + Salinity + TotalNitrogen_unfiltered) +
               FDis_selected * (Temperature + Salinity + TotalNitrogen_unfiltered) +
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

save(model, file = "BRM_zoo_RUE_CD_PI_260422.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)


################################################################################
######################    Phytoplankton RUE BRMS   #############################
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
load('plankton_RUE_data_scaled_260309.RDATA')

setwd("~/GitHub/PlanktonRUE-Project/BRM/RDATA")

#make BRMS model- phytoplankton RUE#--------------------------------------------
#Phytoplankton RUE with phosphorus .............................................
start_time <- Sys.time()

model <- brm(formula = log_RUE_B_P ~ genus_richness + FEve + FDis + 
               species_richness + zoo_FEve + zoo_FDis +
               Temperature + Salinity + pH + MaxDepth + DOY + daylight_hours +
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

save(model, file = "BRM_phytoRUEP_L_260311.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#................

start_time <- Sys.time()

model <- brm(formula = log_RUE_B_P ~ genus_richness + FEve + FDis + 
               species_richness + zoo_FEve + zoo_FDis +
               poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) + 
               poly(MaxDepth,2,raw = TRUE) +
               poly(DOY,2,raw = TRUE) + poly(daylight_hours,2,raw = TRUE) +
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

save(model, file = "BRM_phytoRUEP_P_260311.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#.................

start_time <- Sys.time()

model <- brm(formula = log_RUE_B_P ~ pH + MaxDepth + DOY + daylight_hours +
               genus_richness*Temperature + genus_richness*Salinity +
               FEve*Temperature + FEve*Salinity + 
               FDis*Temperature + FDis*Salinity +
               species_richness*Temperature + species_richness*Salinity +
               zoo_FEve*Temperature + zoo_FEve*Salinity +
               zoo_FDis*Temperature + zoo_FDis*Salinity +
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

save(model, file = "BRM_phytoRUEP_LI_260311.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#................

start_time <- Sys.time()

model <- brm(formula = log_RUE_B_P ~ 
               #env variables
               pH + I(pH^2) + MaxDepth + I(MaxDepth^2) + DOY + I(DOY^2) +
               daylight_hours + I(daylight_hours^2) + Temperature + I(Temperature^2) +
               Salinity + I(Salinity^2) +
               # phytoplankton biodiversity
               genus_richness * (Temperature + Salinity) +
               FEve * (Temperature + Salinity) +
               FDis * (Temperature + Salinity) +
               # zooplankton biodiversity
               species_richness * (Temperature + Salinity) +
               zoo_FEve * (Temperature + Salinity) +
               zoo_FDis * (Temperature + Salinity) +
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

save(model, file = "BRM_phytoRUEP_PI_260311.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#Phytoplankton RUE with nitrogen ...............................................

start_time <- Sys.time()

model <- brm(formula = log_RUE_B_N ~ genus_richness + FEve + FDis + 
               species_richness + zoo_FEve + zoo_FDis +
               Temperature + Salinity + pH + MaxDepth + DOY + daylight_hours +
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

save(model, file = "BRM_phytoRUEN_L_260311.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#................

start_time <- Sys.time()

model <- brm(formula = log_RUE_B_N ~ genus_richness + FEve + FDis + 
               species_richness + zoo_FEve + zoo_FDis +
               poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) + 
               poly(MaxDepth,2,raw = TRUE) +
               poly(DOY,2,raw = TRUE) + poly(daylight_hours,2,raw = TRUE) +
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

save(model, file = "BRM_phytoRUEN_P_260311.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#.................

start_time <- Sys.time()

model <- brm(formula = log_RUE_B_N ~ pH + MaxDepth + DOY + daylight_hours +
               genus_richness*Temperature + genus_richness*Salinity +
               FEve*Temperature + FEve*Salinity + 
               FDis*Temperature + FDis*Salinity +
               species_richness*Temperature + species_richness*Salinity +
               zoo_FEve*Temperature + zoo_FEve*Salinity +
               zoo_FDis*Temperature + zoo_FDis*Salinity +
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

save(model, file = "BRM_phytoRUEN_LI_260311.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#................

start_time <- Sys.time()

model <- brm(formula = log_RUE_B_N ~ 
               #env variables
               pH + I(pH^2) + MaxDepth + I(MaxDepth^2) + DOY + I(DOY^2) +
               daylight_hours + I(daylight_hours^2) + Temperature + I(Temperature^2) +
               Salinity + I(Salinity^2) +
               # phytoplankton biodiversity
               genus_richness * (Temperature + Salinity) +
               FEve * (Temperature + Salinity) +
               FDis * (Temperature + Salinity) +
               # zooplankton biodiversity
               species_richness * (Temperature + Salinity) +
               zoo_FEve * (Temperature + Salinity) +
               zoo_FDis * (Temperature + Salinity) +
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

save(model, file = "BRM_phytoRUEN_PI_260311.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)


#make BRMS model- phytoplankton RUE with CYANO DOM# ----------------------------
#Phytoplankton RUE with phosphorus .............................................
start_time <- Sys.time()

model <- brm(formula = log_RUE_B_P ~ genus_richness + FEve + FDis + cyano_prop +
               species_richness + zoo_FEve + zoo_FDis +
               Temperature + Salinity + pH + MaxDepth + DOY + daylight_hours +
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

save(model, file = "BRM_phytoRUEP_CD_L_260311.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#................

start_time <- Sys.time()

model <- brm(formula = log_RUE_B_P ~ genus_richness + FEve + FDis + cyano_prop +
               species_richness + zoo_FEve + zoo_FDis +
               poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) + 
               poly(MaxDepth,2,raw = TRUE) +
               poly(DOY,2,raw = TRUE) + poly(daylight_hours,2,raw = TRUE) +
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

save(model, file = "BRM_phytoRUEP_CD_P_260311.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#.................

start_time <- Sys.time()

model <- brm(formula = log_RUE_B_P ~ pH + MaxDepth + DOY + daylight_hours +
               genus_richness*Temperature + genus_richness*Salinity +
               FEve*Temperature + FEve*Salinity + 
               FDis*Temperature + FDis*Salinity +
               cyano_prop*Temperature + cyano_prop*Salinity +
               species_richness*Temperature + species_richness*Salinity +
               zoo_FEve*Temperature + zoo_FEve*Salinity +
               zoo_FDis*Temperature + zoo_FDis*Salinity +
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

save(model, file = "BRM_phytoRUEP_CD_LI_260311.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#................

start_time <- Sys.time()

model <- brm(formula = log_RUE_B_P ~ 
               #env variables
               pH + I(pH^2) + MaxDepth + I(MaxDepth^2) + DOY + I(DOY^2) +
               daylight_hours + I(daylight_hours^2) + Temperature + I(Temperature^2) +
               Salinity + I(Salinity^2) +
               # phytoplankton biodiversity
               genus_richness * (Temperature + Salinity) +
               FEve * (Temperature + Salinity) +
               FDis * (Temperature + Salinity) +
               cyano_prop * (Temperature + Salinity) +
               # zooplankton biodiversity
               species_richness * (Temperature + Salinity) +
               zoo_FEve * (Temperature + Salinity) +
               zoo_FDis * (Temperature + Salinity) +
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

save(model, file = "BRM_phytoRUEP_CD_PI_260311.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#Phytoplankton RUE with nitrogen ...............................................

start_time <- Sys.time()

model <- brm(formula = log_RUE_B_N ~ genus_richness + FEve + FDis + cyano_prop +
               species_richness + zoo_FEve + zoo_FDis +
               Temperature + Salinity + pH + MaxDepth + DOY + daylight_hours +
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

save(model, file = "BRM_phytoRUEN_CD_L_260311.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#................

start_time <- Sys.time()

model <- brm(formula = log_RUE_B_N ~ genus_richness + FEve + FDis + cyano_prop +
               species_richness + zoo_FEve + zoo_FDis +
               poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) + 
               poly(MaxDepth,2,raw = TRUE) +
               poly(DOY,2,raw = TRUE) + poly(daylight_hours,2,raw = TRUE) +
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

save(model, file = "BRM_phytoRUEN_CD_P_260311.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#.................

start_time <- Sys.time()

model <- brm(formula = log_RUE_B_N ~ pH + MaxDepth + DOY + daylight_hours +
               genus_richness*Temperature + genus_richness*Salinity +
               FEve*Temperature + FEve*Salinity + 
               FDis*Temperature + FDis*Salinity +
               cyano_prop*Temperature + cyano_prop*Salinity +
               species_richness*Temperature + species_richness*Salinity +
               zoo_FEve*Temperature + zoo_FEve*Salinity +
               zoo_FDis*Temperature + zoo_FDis*Salinity +
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

save(model, file = "BRM_phytoRUEN_CD_LI_260311.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#................

start_time <- Sys.time()

model <- brm(formula = log_RUE_B_N ~ 
               #env variables
               pH + I(pH^2) + MaxDepth + I(MaxDepth^2) + DOY + I(DOY^2) +
               daylight_hours + I(daylight_hours^2) + Temperature + I(Temperature^2) +
               Salinity + I(Salinity^2) +
               # phytoplankton biodiversity
               genus_richness * (Temperature + Salinity) +
               FEve * (Temperature + Salinity) +
               FDis * (Temperature + Salinity) +
               cyano_prop * (Temperature + Salinity) +
               # zooplankton biodiversity
               species_richness * (Temperature + Salinity) +
               zoo_FEve * (Temperature + Salinity) +
               zoo_FDis * (Temperature + Salinity) +
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

save(model, file = "BRM_phytoRUEN_CD_PI_260311.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)


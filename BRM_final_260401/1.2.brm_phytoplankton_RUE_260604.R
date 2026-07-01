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
load('plankton_RUE_data_scaled_260604.RDATA')

setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA")

#make BRMS model- phytoplankton RUE#--------------------------------------------
#Phytoplankton RUE with phosphorus .............................................
start_time <- Sys.time()

model <- brm(formula = log_RUE_B_P ~ genus_richness + FEve + FDis + 
               species_richness + zoo_FEve + zoo_FDis +
               Temperature + Salinity + pH + TotalNitrogen_unfiltered + MaxDepth + 
               daylight_hours +
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

save(model, file = "BRM_phytoRUEP_L_260604.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#................

start_time <- Sys.time()

model <- brm(formula = log_RUE_B_P ~ genus_richness + FEve + FDis + 
               species_richness + zoo_FEve + zoo_FDis +
               poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) + 
               poly(TotalNitrogen_unfiltered,2,raw = TRUE) + poly(MaxDepth,2,raw = TRUE) + 
               poly(daylight_hours,2,raw = TRUE) +
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

save(model, file = "BRM_phytoRUEP_P_260604.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#.................

start_time <- Sys.time()

model <- brm(formula = log_RUE_B_P ~ pH + MaxDepth + daylight_hours +
               genus_richness*Temperature + genus_richness*Salinity + genus_richness*TotalNitrogen_unfiltered +
               FEve*Temperature + FEve*Salinity + FEve*TotalNitrogen_unfiltered + 
               FDis*Temperature + FDis*Salinity + FDis*TotalNitrogen_unfiltered +
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

save(model, file = "BRM_phytoRUEP_LI_260604.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#................

start_time <- Sys.time()

model <- brm(formula = log_RUE_B_P ~ 
               #env variables
               Temperature + I(Temperature^2) + Salinity + I(Salinity^2) + 
               pH + I(pH^2) + MaxDepth + I(MaxDepth^2) + 
               TotalNitrogen_unfiltered + I(TotalNitrogen_unfiltered^2) + 
               daylight_hours + I(daylight_hours^2) + 
               # phytoplankton biodiversity
               genus_richness * (Temperature + Salinity + TotalNitrogen_unfiltered) +
               FEve * (Temperature + Salinity + TotalNitrogen_unfiltered) +
               FDis * (Temperature + Salinity + TotalNitrogen_unfiltered) +
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

save(model, file = "BRM_phytoRUEP_PI_260604.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#Phytoplankton RUE with nitrogen ...............................................

start_time <- Sys.time()

model <- brm(formula = log_RUE_B_N ~ genus_richness + FEve + FDis + 
               species_richness + zoo_FEve + zoo_FDis +
               Temperature + Salinity + pH + TotalNitrogen_unfiltered + MaxDepth + 
               daylight_hours +
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

save(model, file = "BRM_phytoRUEN_L_260604.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#................

start_time <- Sys.time()

model <- brm(formula = log_RUE_B_N ~ genus_richness + FEve + FDis + 
               species_richness + zoo_FEve + zoo_FDis +
               poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) + 
               poly(TotalNitrogen_unfiltered,2,raw = TRUE) + poly(MaxDepth,2,raw = TRUE) + 
               poly(daylight_hours,2,raw = TRUE) +
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

save(model, file = "BRM_phytoRUEN_P_260604.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#.................

start_time <- Sys.time()

model <- brm(formula = log_RUE_B_N ~ pH + MaxDepth + daylight_hours +
               genus_richness*Temperature + genus_richness*Salinity + genus_richness*TotalNitrogen_unfiltered +
               FEve*Temperature + FEve*Salinity + FEve*TotalNitrogen_unfiltered + 
               FDis*Temperature + FDis*Salinity + FDis*TotalNitrogen_unfiltered +
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

save(model, file = "BRM_phytoRUEN_LI_260604.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#................

start_time <- Sys.time()

model <- brm(formula = log_RUE_B_N ~ 
               #env variables
               Temperature + I(Temperature^2) + Salinity + I(Salinity^2) + 
               pH + I(pH^2) + MaxDepth + I(MaxDepth^2) + 
               TotalNitrogen_unfiltered + I(TotalNitrogen_unfiltered^2) + 
               daylight_hours + I(daylight_hours^2) + 
               # phytoplankton biodiversity
               genus_richness * (Temperature + Salinity + TotalNitrogen_unfiltered) +
               FEve * (Temperature + Salinity + TotalNitrogen_unfiltered) +
               FDis * (Temperature + Salinity + TotalNitrogen_unfiltered) +
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

save(model, file = "BRM_phytoRUEN_PI_260604.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)


#make BRMS model- phytoplankton RUE with CYANO DOM# ----------------------------
#Phytoplankton RUE with phosphorus .............................................
start_time <- Sys.time()

model <- brm(formula = log_RUE_B_P ~ genus_richness + FEve + FDis + cyano_prop +
               species_richness + zoo_FEve + zoo_FDis +
               Temperature + Salinity + pH + TotalNitrogen_unfiltered + MaxDepth + 
               daylight_hours +
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

save(model, file = "BRM_phytoRUEP_CD_L_260604.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#................

start_time <- Sys.time()

model <- brm(formula = log_RUE_B_P ~ genus_richness + FEve + FDis + cyano_prop +
               species_richness + zoo_FEve + zoo_FDis +
               poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) + 
               poly(TotalNitrogen_unfiltered,2,raw = TRUE) + poly(MaxDepth,2,raw = TRUE) + 
               poly(daylight_hours,2,raw = TRUE) +
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

save(model, file = "BRM_phytoRUEP_CD_P_260604.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#.................

start_time <- Sys.time()

model <- brm(formula = log_RUE_B_P ~ pH + MaxDepth + daylight_hours +
               genus_richness*Temperature + genus_richness*Salinity + genus_richness*TotalNitrogen_unfiltered +
               FEve*Temperature + FEve*Salinity + FEve*TotalNitrogen_unfiltered + 
               FDis*Temperature + FDis*Salinity + FDis*TotalNitrogen_unfiltered +
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

save(model, file = "BRM_phytoRUEP_CD_LI_260604.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#................

start_time <- Sys.time()

model <- brm(formula = log_RUE_B_P ~ 
               #env variables
               Temperature + I(Temperature^2) + Salinity + I(Salinity^2) + 
               pH + I(pH^2) + MaxDepth + I(MaxDepth^2) + 
               TotalNitrogen_unfiltered + I(TotalNitrogen_unfiltered^2) + 
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

save(model, file = "BRM_phytoRUEP_CD_PI_260604.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#Phytoplankton RUE with nitrogen ...............................................

start_time <- Sys.time()

model <- brm(formula = log_RUE_B_N ~ genus_richness + FEve + FDis + cyano_prop +
               species_richness + zoo_FEve + zoo_FDis +
               Temperature + Salinity + pH + TotalNitrogen_unfiltered + MaxDepth + 
               daylight_hours +
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

save(model, file = "BRM_phytoRUEN_CD_L_260604.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#................

start_time <- Sys.time()

model <- brm(formula = log_RUE_B_N ~ genus_richness + FEve + FDis + cyano_prop +
               species_richness + zoo_FEve + zoo_FDis +
               poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) + 
               poly(TotalNitrogen_unfiltered,2,raw = TRUE) + poly(MaxDepth,2,raw = TRUE) + 
               poly(daylight_hours,2,raw = TRUE) +
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

save(model, file = "BRM_phytoRUEN_CD_P_260604.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#.................

start_time <- Sys.time()

model <- brm(formula = log_RUE_B_N ~ pH + MaxDepth + daylight_hours +
               genus_richness*Temperature + genus_richness*Salinity + genus_richness*TotalNitrogen_unfiltered +
               FEve*Temperature + FEve*Salinity + FEve*TotalNitrogen_unfiltered + 
               FDis*Temperature + FDis*Salinity + FDis*TotalNitrogen_unfiltered +
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

save(model, file = "BRM_phytoRUEN_CD_LI_260604.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#................

start_time <- Sys.time()

model <- brm(formula = log_RUE_B_N ~ 
               #env variables
               Temperature + I(Temperature^2) + Salinity + I(Salinity^2) + 
               pH + I(pH^2) + MaxDepth + I(MaxDepth^2) + 
               TotalNitrogen_unfiltered + I(TotalNitrogen_unfiltered^2) + 
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

save(model, file = "BRM_phytoRUEN_CD_PI_260604.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)


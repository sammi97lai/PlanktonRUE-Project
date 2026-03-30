################################################################################
######################    Plankton RUE and div BRMS:  ##########################
######################    N:P ratio                   ##########################
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
load('plankton_RUE_data_scaled_260330.RDATA')

setwd("~/GitHub/PlanktonRUE-Project/BRM_model_test_260330/RDATA")

#make BRMS model- phytoplankton diversity #-----------------------------------
#genus richness.......................................................
start_time <- Sys.time()

model <- brm( formula= genus_richness ~ poly(Temperature,2,raw = TRUE) + 
                poly(pH,2,raw = TRUE) + 
                poly(log_N_P_ratio,2,raw = TRUE) + poly(MaxDepth,2,raw = TRUE) + 
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

save(model, file = "BRM_phyto_GD_P_260330_NP.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#Phytoplankton FEve .......................................................
start_time <- Sys.time()

model <- brm( formula= FEve ~ poly(Temperature,2,raw = TRUE) + 
                poly(pH,2,raw = TRUE) + 
                poly(log_N_P_ratio,2,raw = TRUE) + poly(MaxDepth,2,raw = TRUE) +
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

save(model, file = "BRM_phyto_FEve_P_260330_NP.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#Phytoplankton FDis .................................................
start_time <- Sys.time()

model <- brm( formula= FDis ~ poly(Temperature,2,raw = TRUE) + 
                poly(pH,2,raw = TRUE) + 
                poly(log_N_P_ratio,2,raw = TRUE) + poly(MaxDepth,2,raw = TRUE) +
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

save(model, file = "BRM_phyto_FDis_P_260330_NP.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#make BRMS model - Cyanobacteria dominance- zero inflated #--------------------
start_time <- Sys.time()

model <- brm( formula= cyano_prop ~ poly(Temperature,2,raw = TRUE) + 
                poly(pH,2,raw = TRUE) + 
                poly(log_N_P_ratio,2,raw = TRUE) + poly(MaxDepth,2,raw = TRUE) +
                s(DOY) + poly(daylight_hours,2,raw = TRUE) +
                (1 | Year) + (1 | Site_Id),
              data = plankton_RUE_data_scaled,
              family = zero_one_inflated_beta(), #changed family
              iter = 4000, warmup = 2000, cores = 4, chains = 4,
              save_pars = save_pars(all = TRUE),
              prior = c(
                set_prior("normal(0, 1)", class = "b"),
                set_prior("exponential(1)", class = "sd"),
                set_prior("exponential(1)", class = "phi"),
                set_prior("normal(-2, 1)", class = "zoi")  # prior for zero inflation
              ),
              control = list(adapt_delta = 0.99, max_treedepth = 15),
              seed = 2345
)

save(model, file = "BRM_cyanodom_P_260330_NP.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#make BRMS model- phytoplankton RUE with CYANO DOM# ----------------------------
#Phytoplankton RUE with phosphorus .............................................
start_time <- Sys.time()

model <- brm(formula = log_RUE_B_P ~ 
               #env variables
               Temperature + I(Temperature^2) + 
               pH + I(pH^2) + MaxDepth + I(MaxDepth^2) + 
               log_N_P_ratio + I(log_N_P_ratio^2) + s(DOY) +
               daylight_hours + I(daylight_hours^2) + 
               # phytoplankton biodiversity
               genus_richness * (Temperature + log_N_P_ratio) +
               FEve * (Temperature + log_N_P_ratio) +
               FDis * (Temperature + log_N_P_ratio) +
               cyano_prop * (Temperature + log_N_P_ratio) +
               # zooplankton biodiversity
               species_richness * (Temperature + log_N_P_ratio) +
               zoo_FEve * (Temperature + log_N_P_ratio) +
               zoo_FDis * (Temperature + log_N_P_ratio) +
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

save(model, file = "BRM_phytoRUEP_CD_PI_260330_NP.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#Phytoplankton RUE with nitrogen ...............................................
start_time <- Sys.time()

model <- brm(formula = log_RUE_B_N ~ 
               #env variables
               Temperature + I(Temperature^2) + 
               pH + I(pH^2) + MaxDepth + I(MaxDepth^2) + 
               log_N_P_ratio + I(log_N_P_ratio^2) + s(DOY) +
               daylight_hours + I(daylight_hours^2) + 
               # phytoplankton biodiversity
               genus_richness * (Temperature + log_N_P_ratio) +
               FEve * (Temperature + log_N_P_ratio) +
               FDis * (Temperature + log_N_P_ratio) +
               cyano_prop * (Temperature + log_N_P_ratio) +
               # zooplankton biodiversity
               species_richness * (Temperature + log_N_P_ratio) +
               zoo_FEve * (Temperature + log_N_P_ratio) +
               zoo_FDis * (Temperature + log_N_P_ratio) +
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

save(model, file = "BRM_phytoRUEN_CD_PI_260330_NP.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#make BRMS model - Zooplankton diversity with CYANO DOM # ----------------------
#Zooplankton Species Richness #..................................................
start_time <- Sys.time()

model <- brm(formula = species_richness ~ poly(Temperature,2,raw = TRUE) + 
               poly(pH,2,raw = TRUE) +
               poly(log_N_P_ratio,2,raw = TRUE) + poly(Zooplankton_MaxDepth,2,raw = TRUE) +
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

save(model, file = "BRM_zoo_SR_CD_P_260330_NP.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#Zooplankton Functional Evenness ...............................................
start_time <- Sys.time()

model <- brm(formula = zoo_FEve ~ poly(Temperature,2,raw = TRUE) + 
               poly(pH,2,raw = TRUE) +
               poly(log_N_P_ratio,2,raw = TRUE) + poly(Zooplankton_MaxDepth,2,raw = TRUE) +
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

save(model, file = "BRM_zoo_FE_CD_P_260330_NP.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#Zooplankton Functional Dispersion .............................................
start_time <- Sys.time()

model <- brm(formula = zoo_FDis ~ poly(Temperature,2,raw = TRUE) + 
               poly(pH,2,raw = TRUE) +
               poly(log_N_P_ratio,2,raw = TRUE) + poly(Zooplankton_MaxDepth,2,raw = TRUE) +
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

save(model, file = "BRM_zoo_FD_CD_P_260330_NP.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#make BRMS model- zooplankton RUE with CYANO DOM#-------------------------------
start_time <- Sys.time()

model <- brm(formula = log_RUE_zooplankton ~ 
               #env variables
               Temperature + I(Temperature^2) + 
               pH + I(pH^2) + Zooplankton_MaxDepth + I(Zooplankton_MaxDepth^2) + 
               log_N_P_ratio + I(log_N_P_ratio^2) + s(DOY) +
               daylight_hours + I(daylight_hours^2) + 
               # phytoplankton biodiversity
               genus_richness * (Temperature + log_N_P_ratio) +
               FEve_selected * (Temperature + log_N_P_ratio) +
               FDis_selected * (Temperature + log_N_P_ratio) +
               cyano_prop * (Temperature + log_N_P_ratio) +
               # zooplankton biodiversity
               species_richness * (Temperature + log_N_P_ratio) +
               zoo_FEve * (Temperature + log_N_P_ratio) +
               zoo_FDis * (Temperature + log_N_P_ratio) +
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

save(model, file = "BRM_zoo_RUE_CD_PI_260330_NP.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)



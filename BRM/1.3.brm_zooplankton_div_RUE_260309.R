################################################################################
######################    RUE Zooplankton BRMS   ###############################
################################################################################

library(dplyr)
library(readxl)
library(tidyr)
library(lubridate)
library(brms)

rm(list = ls())
#------------------------------------------------------------------------------
setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/RDATA")
load('plankton_RUE_data_scaled_260309.RDATA')

setwd("~/GitHub/PlanktonRUE-Project/BRM")

#make BRMS model - Zooplankton diversity # -------------------------------------
#Zooplankton Species Richness #..................................................
start_time <- Sys.time()

model <- brm(formula = species_richness ~ Temperature + Salinity + pH +
               genus_richness + FEve_selected + FDis_selected +
               DOY + 
               (1 | Year) + (1 | Site_Id) + (1 | Zooplankton_MaxDepth), 
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

save(model, file = "BRM_zoo_SR_L_260109.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#...................
start_time <- Sys.time()

model <- brm(formula = species_richness ~ poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) +
               genus_richness + FEve_selected + FDis_selected +
               poly(DOY,2,raw = TRUE) + 
               (1 | Year) + (1 | Site_Id) + (1 | Zooplankton_MaxDepth), 
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

save(model, file = "BRM_zoo_SR_P_260109.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#Zooplankton Species Richness- poisson #........................................
start_time <- Sys.time()

model <- brm(formula = species_richness ~ Temperature + Salinity + pH +
               genus_richness + FEve_selected + FDis_selected +
               DOY + 
               (1 | Year) + (1 | Site_Id) + (1 | Zooplankton_MaxDepth), 
             data = plankton_RUE_data_scaled,
             family = poisson(), #changed family
             iter = 4000, warmup = 2000, cores = 4, chains = 4,
             save_pars = save_pars(all = TRUE),
             prior = c(
               set_prior("normal(0, 1)", class = "b"),
               set_prior("exponential(1)", class = "sd")
             ),
             control = list(adapt_delta = 0.99, max_treedepth = 15),
             seed = 2345
)

save(model, file = "BRM_zoo_SR_L_p_260127.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#...................
start_time <- Sys.time()

model <- brm(formula = species_richness ~ poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) +
               genus_richness + FEve_selected + FDis_selected +
               poly(DOY,2,raw = TRUE) + 
               (1 | Year) + (1 | Site_Id) + (1 | Zooplankton_MaxDepth), 
             data = plankton_RUE_data_scaled,
             family = poisson(), #changed family
             iter = 4000, warmup = 2000, cores = 4, chains = 4,
             save_pars = save_pars(all = TRUE),
             prior = c(
               set_prior("normal(0, 1)", class = "b"),
               set_prior("exponential(1)", class = "sd")
             ),
             control = list(adapt_delta = 0.99, max_treedepth = 15),
             seed = 2345
)

save(model, file = "BRM_zoo_SR_P_p_260127.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)


#Functional Evenness ...........................................................
start_time <- Sys.time()

model <- brm(formula = zoo_FEve ~ Temperature + Salinity + pH +
               genus_richness + FEve_selected + FDis_selected +
               DOY + 
               (1 | Year) + (1 | Site_Id) + (1 | Zooplankton_MaxDepth), 
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

save(model, file = "BRM_zoo_FE_L_260109.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)


#...................
start_time <- Sys.time()

model <- brm(formula = zoo_FEve ~ poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) +
               genus_richness + FEve_selected + FDis_selected +
               poly(DOY,2,raw = TRUE) + 
               (1 | Year) + (1 | Site_Id) + (1 | Zooplankton_MaxDepth), 
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

save(model, file = "BRM_zoo_FE_P_260109.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#Functional Evenness -beta .....................................................
start_time <- Sys.time()

model <- brm(formula = zoo_FEve ~ Temperature + Salinity + pH +
               genus_richness + FEve_selected + FDis_selected +
               DOY + 
               (1 | Year) + (1 | Site_Id) + (1 | Zooplankton_MaxDepth), 
             data = plankton_RUE_data_scaled,
             family = Beta(), #changed family
             iter = 4000, warmup = 2000, cores = 4, chains = 4,
             save_pars = save_pars(all = TRUE),
             prior = c(
               set_prior("normal(0, 1)", class = "b"),
               set_prior("exponential(1)", class = "sd"),
               set_prior("exponential(1)", class = "phi") # prior for beta
             ),
             control = list(adapt_delta = 0.99, max_treedepth = 15),
             seed = 2345
)

save(model, file = "BRM_zoo_FE_L_b_260109.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)


#...................
start_time <- Sys.time()

model <- brm(formula = zoo_FEve ~ poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) +
               genus_richness + FEve_selected + FDis_selected +
               poly(DOY,2,raw = TRUE) + 
               (1 | Year) + (1 | Site_Id) + (1 | Zooplankton_MaxDepth), 
             data = plankton_RUE_data_scaled,
             family = Beta(), #changed family
             iter = 4000, warmup = 2000, cores = 4, chains = 4,
             save_pars = save_pars(all = TRUE),
             prior = c(
               set_prior("normal(0, 1)", class = "b"),
               set_prior("exponential(1)", class = "sd"),
               set_prior("exponential(1)", class = "phi") # prior for beta
             ),
             control = list(adapt_delta = 0.99, max_treedepth = 15),
             seed = 2345
)

save(model, file = "BRM_zoo_FE_P_b_260109.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#Functional Dispersion .........................................................
start_time <- Sys.time()

model <- brm(formula = zoo_FDis ~ Temperature + Salinity  + pH +
               genus_richness + FEve_selected + FDis_selected +
               DOY + 
               (1 | Year) + (1 | Site_Id) + (1 | Zooplankton_MaxDepth), 
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

save(model, file = "BRM_zoo_FD_L_260109.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#...................
start_time <- Sys.time()

model <- brm(formula = zoo_FDis ~ poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) +
               genus_richness + FEve_selected + FDis_selected +
               poly(DOY,2,raw = TRUE) + 
               (1 | Year) + (1 | Site_Id) + (1 | Zooplankton_MaxDepth), 
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

save(model, file = "BRM_zoo_FD_P_260109.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#Functional Dispersion -beta ...................................................
start_time <- Sys.time()

model <- brm(formula = zoo_FDis ~ Temperature + Salinity  + pH +
               genus_richness + FEve_selected + FDis_selected +
               DOY + 
               (1 | Year) + (1 | Site_Id) + (1 | Zooplankton_MaxDepth), 
             data = plankton_RUE_data_scaled,
             family = Beta(), #changed family
             iter = 4000, warmup = 2000, cores = 4, chains = 4,
             save_pars = save_pars(all = TRUE),
             prior = c(
               set_prior("normal(0, 1)", class = "b"),
               set_prior("exponential(1)", class = "sd"),
               set_prior("exponential(1)", class = "phi") # prior for beta
             ),
             control = list(adapt_delta = 0.99, max_treedepth = 15),
             seed = 2345
)

save(model, file = "BRM_zoo_FD_L_b_260127.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#...................
start_time <- Sys.time()

model <- brm(formula = zoo_FDis ~ poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) +
               genus_richness + FEve_selected + FDis_selected +
               poly(DOY,2,raw = TRUE) + 
               (1 | Year) + (1 | Site_Id) + (1 | Zooplankton_MaxDepth), 
             data = plankton_RUE_data_scaled,
             family = Beta(), #changed family
             iter = 4000, warmup = 2000, cores = 4, chains = 4,
             save_pars = save_pars(all = TRUE),
             prior = c(
               set_prior("normal(0, 1)", class = "b"),
               set_prior("exponential(1)", class = "sd"),
               set_prior("exponential(1)", class = "phi") # prior for beta
             ),
             control = list(adapt_delta = 0.99, max_treedepth = 15),
             seed = 2345
)

save(model, file = "BRM_zoo_FD_P_b_260127.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)


#make BRMS model - Zooplankton diversity with CYANO DOM # ----------------------
#Zooplankton Species Richness#..................................................
start_time <- Sys.time()

model <- brm(formula = species_richness ~ Temperature + Salinity  + pH +
               genus_richness + FEve_selected + FDis_selected + cyano_prop +
               DOY + 
               (1 | Year) + (1 | Site_Id) + (1 | Zooplankton_MaxDepth), 
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

save(model, file = "BRM_zoo_SR_CD_L_260109.RData")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#...................
start_time <- Sys.time()

model <- brm(formula = species_richness ~ poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) +
               genus_richness + FEve_selected + FDis_selected + cyano_prop +
               poly(DOY,2,raw = TRUE) + 
               (1 | Year) + (1 | Site_Id) + (1 | Zooplankton_MaxDepth), 
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

save(model, file = "BRM_zoo_SR_CD_P_260109.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#Zooplankton Species Richness- poisson#.........................................
start_time <- Sys.time()

model <- brm(formula = species_richness ~ Temperature + Salinity  + pH +
               genus_richness + FEve_selected + FDis_selected + cyano_prop +
               DOY + 
               (1 | Year) + (1 | Site_Id) + (1 | Zooplankton_MaxDepth), 
             data = plankton_RUE_data_scaled,
             family = poisson(), #changed family
             iter = 4000, warmup = 2000, cores = 4, chains = 4,
             save_pars = save_pars(all = TRUE),
             prior = c(
               set_prior("normal(0, 1)", class = "b"),
               set_prior("exponential(1)", class = "sd")
             ),
             control = list(adapt_delta = 0.99, max_treedepth = 15),
             seed = 2345
)

save(model, file = "BRM_zoo_SR_CD_L_p_260127.RData")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#...................
start_time <- Sys.time()

model <- brm(formula = species_richness ~ poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) +
               genus_richness + FEve_selected + FDis_selected + cyano_prop +
               poly(DOY,2,raw = TRUE) + 
               (1 | Year) + (1 | Site_Id) + (1 | Zooplankton_MaxDepth), 
             data = plankton_RUE_data_scaled,
             family = poisson(), #changed family
             iter = 4000, warmup = 2000, cores = 4, chains = 4,
             save_pars = save_pars(all = TRUE),
             prior = c(
               set_prior("normal(0, 1)", class = "b"),
               set_prior("exponential(1)", class = "sd")
             ),
             control = list(adapt_delta = 0.99, max_treedepth = 15),
             seed = 2345
)

save(model, file = "BRM_zoo_SR_CD_P_p_260127.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#Functional Evenness ...........................................................
start_time <- Sys.time()

model <- brm(formula = zoo_FEve ~ Temperature + Salinity  + pH +
               genus_richness + FEve_selected + FDis_selected + cyano_prop +
               DOY + 
               (1 | Year) + (1 | Site_Id) + (1 | Zooplankton_MaxDepth), 
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

save(model, file = "BRM_zoo_FE_CD_L_260109.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#...................
start_time <- Sys.time()

model <- brm(formula = zoo_FEve ~ poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) +
               genus_richness + FEve_selected + FDis_selected + cyano_prop +
               poly(DOY,2,raw = TRUE) + 
               (1 | Year) + (1 | Site_Id) + (1 | Zooplankton_MaxDepth), 
             data = plankton_RUE_data_scaled,
             iter = 6000, warmup = 3000, cores = 4, chains = 4, #increased inter
             save_pars = save_pars(all = TRUE),
             prior = c(
               set_prior("normal(0, 1)", class = "b"),
               set_prior("exponential(1)", class = "sd"),
               set_prior("exponential(1)", class = "sigma")
             ),
             control = list(adapt_delta = 0.995, max_treedepth = 15), #increased adapt_delta
             seed = 2345
)

save(model, file = "BRM_zoo_FE_CD_P_260109.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)
# Tail Effective Samples Size (ESS) was too low so model parameters were changed

#Functional Evenness -beta .....................................................
start_time <- Sys.time()

model <- brm(formula = zoo_FEve ~ Temperature + Salinity  + pH +
               genus_richness + FEve_selected + FDis_selected + cyano_prop +
               DOY + 
               (1 | Year) + (1 | Site_Id) + (1 | Zooplankton_MaxDepth), 
             data = plankton_RUE_data_scaled,
             family = Beta(), #changed family
             iter = 4000, warmup = 2000, cores = 4, chains = 4,
             save_pars = save_pars(all = TRUE),
             prior = c(
               set_prior("normal(0, 1)", class = "b"),
               set_prior("exponential(1)", class = "sd"),
               set_prior("exponential(1)", class = "phi") # prior for beta
             ),
             control = list(adapt_delta = 0.99, max_treedepth = 15),
             seed = 2345
)

save(model, file = "BRM_zoo_FE_CD_L_b_260127.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#...................
start_time <- Sys.time()

model <- brm(formula = zoo_FEve ~ poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) +
               genus_richness + FEve_selected + FDis_selected + cyano_prop +
               poly(DOY,2,raw = TRUE) + 
               (1 | Year) + (1 | Site_Id) + (1 | Zooplankton_MaxDepth), 
             data = plankton_RUE_data_scaled,
             family = Beta(), #changed family
             iter = 6000, warmup = 3000, cores = 4, chains = 4, #increased iter
             save_pars = save_pars(all = TRUE),
             prior = c(
               set_prior("normal(0, 1)", class = "b"),
               set_prior("exponential(1)", class = "sd"),
               set_prior("exponential(1)", class = "phi") # prior for beta
             ),
             control = list(adapt_delta = 0.995, max_treedepth = 15), #increased adapt_delta
             seed = 2345
)

save(model, file = "BRM_zoo_FE_CD_P_b_260127.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)
# Tail Effective Samples Size (ESS) was too low so model parameters were changed

#Functional Dispersion .........................................................
start_time <- Sys.time()

model <- brm(formula = zoo_FDis ~ Temperature + Salinity  + pH +
               genus_richness + FEve_selected + FDis_selected + cyano_prop +
               DOY + 
               (1 | Year) + (1 | Site_Id) + (1 | Zooplankton_MaxDepth), 
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

save(model, file = "BRM_zoo_FD_CD_L_260109.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#...................
start_time <- Sys.time()

model <- brm(formula = zoo_FDis ~ poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) +
               genus_richness + FEve_selected + FDis_selected + cyano_prop +
               poly(DOY,2,raw = TRUE) + 
               (1 | Year) + (1 | Site_Id) + (1 | Zooplankton_MaxDepth), 
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

save(model, file = "BRM_zoo_FD_CD_P_260109.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#Functional Dispersion -beta....................................................
start_time <- Sys.time()

model <- brm(formula = zoo_FDis ~ Temperature + Salinity  + pH +
               genus_richness + FEve_selected + FDis_selected + cyano_prop +
               DOY + 
               (1 | Year) + (1 | Site_Id) + (1 | Zooplankton_MaxDepth), 
             data = plankton_RUE_data_scaled,
             family = Beta(), #changed family
             iter = 4000, warmup = 2000, cores = 4, chains = 4,
             save_pars = save_pars(all = TRUE),
             prior = c(
               set_prior("normal(0, 1)", class = "b"),
               set_prior("exponential(1)", class = "sd"),
               set_prior("exponential(1)", class = "phi") # prior for beta
             ),
             control = list(adapt_delta = 0.99, max_treedepth = 15),
             seed = 2345
)

save(model, file = "BRM_zoo_FD_CD_L_b_260127.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#...................
start_time <- Sys.time()

model <- brm(formula = zoo_FDis ~ poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) +
               genus_richness + FEve_selected + FDis_selected + cyano_prop +
               poly(DOY,2,raw = TRUE) + 
               (1 | Year) + (1 | Site_Id) + (1 | Zooplankton_MaxDepth), 
             data = plankton_RUE_data_scaled,
             family = Beta(), #changed family
             iter = 4000, warmup = 2000, cores = 4, chains = 4,
             save_pars = save_pars(all = TRUE),
             prior = c(
               set_prior("normal(0, 1)", class = "b"),
               set_prior("exponential(1)", class = "sd"),
               set_prior("exponential(1)", class = "phi") # prior for beta
             ),
             control = list(adapt_delta = 0.99, max_treedepth = 15),
             seed = 2345
)

save(model, file = "BRM_zoo_FD_CD_P_b_260127.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#make BRMS model- zooplankton RUE#----------------------------------------------
start_time <- Sys.time()

model <- brm(formula = log_RUE_zooplankton ~ Temperature + Salinity  + pH +
               genus_richness + FEve_selected + FDis_selected + 
               species_richness + zoo_FEve + zoo_FDis +
               DOY + 
               (1 | Year) + (1 | Site_Id) + (1 | Zooplankton_MaxDepth), 
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

save(model, file = "BRM_zoo_RUE_L_260109.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#.........................

start_time <- Sys.time()

model <- brm(formula = log_RUE_zooplankton ~ poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) +
               genus_richness + FEve_selected + FDis_selected + 
               species_richness + zoo_FEve + zoo_FDis +
               poly(DOY,2,raw = TRUE) + 
               (1 | Year) + (1 | Site_Id) + (1 | Zooplankton_MaxDepth), 
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

save(model, file = "BRM_zoo_RUE_P_260109.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)


#make BRMS model- zooplankton RUE with CYANO DOM#-------------------------------
start_time <- Sys.time()

model <- brm(formula = log_RUE_zooplankton ~ Temperature + Salinity + pH +
               genus_richness + FEve_selected + FDis_selected + cyano_prop +
               species_richness + zoo_FEve + zoo_FDis +
               DOY + 
               (1 | Year) + (1 | Site_Id) + (1 | Zooplankton_MaxDepth), 
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

save(model, file = "BRM_zoo_RUE_CD_L_260109.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#....................
start_time <- Sys.time()

model <- brm(formula = log_RUE_zooplankton ~ poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) +
               genus_richness + FEve_selected + FDis_selected + cyano_prop +
               species_richness + zoo_FEve + zoo_FDis +
               poly(DOY,2,raw = TRUE) + 
               (1 | Year) + (1 | Site_Id) + (1 | Zooplankton_MaxDepth), 
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

save(model, file = "BRM_zoo_RUE_CD_P_260109.RData")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)


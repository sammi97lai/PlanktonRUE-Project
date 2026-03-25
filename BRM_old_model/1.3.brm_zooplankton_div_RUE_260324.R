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
load('plankton_RUE_data_scaled_old_260324.RDATA')

setwd("~/GitHub/PlanktonRUE-Project/BRM_old_model/RDATA")

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

save(model, file = "BRM_zoo_SR_L_260324.RDATA")

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

save(model, file = "BRM_zoo_SR_P_260324.RDATA")

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

save(model, file = "BRM_zoo_FE_L_260324.RDATA")

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

save(model, file = "BRM_zoo_FE_P_260324.RDATA")

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

save(model, file = "BRM_zoo_FD_L_260324.RDATA")

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

save(model, file = "BRM_zoo_FD_P_260324.RDATA")

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

save(model, file = "BRM_zoo_SR_CD_L_260324.RData")

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

save(model, file = "BRM_zoo_SR_CD_P_260324.RDATA")

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

save(model, file = "BRM_zoo_FE_CD_L_260324.RDATA")

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

save(model, file = "BRM_zoo_FE_CD_P_260324.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

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

save(model, file = "BRM_zoo_FD_CD_L_260324.RDATA")

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

save(model, file = "BRM_zoo_FD_CD_P_260324.RDATA")

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

save(model, file = "BRM_zoo_RUE_L_260324.RDATA")

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

save(model, file = "BRM_zoo_RUE_P_260324.RDATA")

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

save(model, file = "BRM_zoo_RUE_CD_L_260324.RDATA")

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

save(model, file = "BRM_zoo_RUE_CD_P_260324.RData")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)


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
setwd("~/GitHub/Phytoplankton-project/phyto_RUE/dataformatting")
load('plankton_RUE_data_scaled_260109.RDATA')

setwd("~/GitHub/Phytoplankton-project/phyto_RUE/BRM_260109")

#make BRMS model- phytoplankton RUE#--------------------------------------------
#Phytoplankton RUE with phosphorus .............................................
start_time <- Sys.time()

model <- brm(formula = log_RUE_A_P ~ genus_richness + FEve + FDis + 
               Temperature + Salinity + pH + log_N_P_ratio +
               DOY + 
               (1 | Year) + (1 | Site_Id) + (1 | MaxDepth), 
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

save(model, file = "BRM_phytoRUEP_L_260109.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#................

start_time <- Sys.time()

model <- brm(formula = log_RUE_A_P ~ genus_richness + FEve + FDis + 
               poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) + 
               poly(log_N_P_ratio,2,raw = TRUE) +
               poly(DOY,2,raw = TRUE) + 
               (1 | Year) + (1 | Site_Id) + (1 | MaxDepth), 
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

save(model, file = "BRM_phytoRUEP_P_260109.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#.................

start_time <- Sys.time()

model <- brm(formula = log_RUE_A_P ~ genus_richness + FEve + FDis + 
               Temperature + Salinity + pH + log_N_P_ratio*DOY + 
               (1 | Year) + (1 | Site_Id) + (1 | MaxDepth), 
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

save(model, file = "BRM_phytoRUEP_LI_260109.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#................

start_time <- Sys.time()

model <- brm(formula = log_RUE_A_P ~ genus_richness + FEve + FDis + 
               poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) + 
               log_N_P_ratio*DOY + I(log_N_P_ratio^2)*I(DOY^2) +
               (1 | Year) + (1 | Site_Id) + (1 | MaxDepth), 
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

save(model, file = "BRM_phytoRUEP_PI_260109.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#Phytoplankton RUE with nitrogen ...............................................

start_time <- Sys.time()

model <- brm(formula = log_RUE_A_N ~ genus_richness + FEve + FDis + 
               Temperature + Salinity + pH + log_N_P_ratio +
               DOY + 
               (1 | Year) + (1 | Site_Id) + (1 | MaxDepth), 
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

save(model, file = "BRM_phytoRUEN_L_260109.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#................

start_time <- Sys.time()

model <- brm(formula = log_RUE_A_N ~ genus_richness + FEve + FDis + 
               poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) + 
               poly(log_N_P_ratio,2,raw = TRUE) +
               poly(DOY,2,raw = TRUE) + 
               (1 | Year) + (1 | Site_Id) + (1 | MaxDepth), 
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

save(model, file = "BRM_phytoRUEN_P_260109.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#.................

start_time <- Sys.time()

model <- brm(formula = log_RUE_A_N ~ genus_richness + FEve + FDis + 
               Temperature + Salinity + pH + log_N_P_ratio*DOY + 
               (1 | Year) + (1 | Site_Id) + (1 | MaxDepth), 
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

save(model, file = "BRM_phytoRUEN_LI_260109.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#................

start_time <- Sys.time()

model <- brm(formula = log_RUE_A_N ~ genus_richness + FEve + FDis + 
               poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) + 
               log_N_P_ratio*DOY + I(log_N_P_ratio^2)*I(DOY^2) +
               (1 | Year) + (1 | Site_Id) + (1 | MaxDepth), 
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

save(model, file = "BRM_phytoRUEN_PI_260109.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)


#make BRMS model- phytoplankton RUE with CYANO DOM# ----------------------------
#Phytoplankton RUE with phosphorus .............................................
start_time <- Sys.time()

model <- brm(formula = log_RUE_A_P ~ genus_richness + FEve + FDis + cyano_prop +
               Temperature + Salinity + pH + log_N_P_ratio +
               DOY + 
               (1 | Year) + (1 | Site_Id) + (1 | MaxDepth), 
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

save(model, file = "BRM_phytoRUEP_CD_L_260109.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#................

start_time <- Sys.time()

model <- brm(formula = log_RUE_A_P ~ genus_richness + FEve + FDis + cyano_prop +
               poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) + 
               poly(log_N_P_ratio,2,raw = TRUE) +
               poly(DOY,2,raw = TRUE) + 
               (1 | Year) + (1 | Site_Id) + (1 | MaxDepth), 
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

save(model, file = "BRM_phytoRUEP_CD_P_260109.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#.................

start_time <- Sys.time()

model <- brm(formula = log_RUE_A_P ~ genus_richness + FEve + FDis + cyano_prop +
               Temperature + Salinity + pH + log_N_P_ratio*DOY +
               (1 | Year) + (1 | Site_Id) + (1 | MaxDepth), 
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

save(model, file = "BRM_phytoRUEP_CD_LI_260109.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#................

start_time <- Sys.time()

model <- brm(formula = log_RUE_A_P ~ genus_richness + FEve + FDis + cyano_prop +
               poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) + 
               log_N_P_ratio*DOY + I(log_N_P_ratio^2)*I(DOY^2) +
               (1 | Year) + (1 | Site_Id) + (1 | MaxDepth), 
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

save(model, file = "BRM_phytoRUEP_CD_PI_260109.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#Phytoplankton RUE with nitrogen ...............................................

start_time <- Sys.time()

model <- brm(formula = log_RUE_A_N ~ genus_richness + FEve + FDis + cyano_prop +
               Temperature + Salinity + pH + log_N_P_ratio +
               DOY + 
               (1 | Year) + (1 | Site_Id) + (1 | MaxDepth), 
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

save(model, file = "BRM_phytoRUEN_CD_L_260109.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#................

start_time <- Sys.time()

model <- brm(formula = log_RUE_A_N ~ genus_richness + FEve + FDis + cyano_prop +
               poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) +
               poly(log_N_P_ratio,2,raw = TRUE) +
               poly(DOY,2,raw = TRUE) + 
               (1 | Year) + (1 | Site_Id) + (1 | MaxDepth), 
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

save(model, file = "BRM_phytoRUEN_CD_P_260109.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#.................

start_time <- Sys.time()

model <- brm(formula = log_RUE_A_N ~ genus_richness + FEve + FDis + cyano_prop +
               Temperature + Salinity + pH + log_N_P_ratio*DOY + 
               (1 | Year) + (1 | Site_Id) + (1 | MaxDepth), 
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

save(model, file = "BRM_phytoRUEN_CD_LI_260109.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#................

start_time <- Sys.time()

model <- brm(formula = log_RUE_A_N ~ genus_richness + FEve + FDis + cyano_prop +
               poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) + 
               log_N_P_ratio*DOY + I(log_N_P_ratio^2)*I(DOY^2) +
               (1 | Year) + (1 | Site_Id) + (1 | MaxDepth), 
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

save(model, file = "BRM_phytoRUEN_CD_PI_260109.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)


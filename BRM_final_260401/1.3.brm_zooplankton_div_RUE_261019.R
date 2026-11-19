################################################################################
######################    RUE Zooplankton and div BRMS   #######################
################################################################################

library(dplyr)
library(tidyr)
library(brms)

rm(list = ls())
#------------------------------------------------------------------------------
setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/RDATA")
load('plankton_RUE_data_scaled_Zdiv_260916.RDATA')
load('plankton_RUE_data_scaled_RUE_260916.RDATA')

setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA")

project_dir <- "~/GitHub/PlanktonRUE-Project"
model_dir <- file.path(
  project_dir,
  "BRM_final_260401/RDATA"
)

#load this when using remote desktop
# Local directory for Stan compilation
stan_dir <- "C:/R_temp"
if (!dir.exists(stan_dir)) {
  dir.create(stan_dir, recursive = TRUE)
}
# Use local disk for Stan compilation
setwd(stan_dir)

#make BRMS model - Zooplankton diversity # -------------------------------------
#Zooplankton species richness #..................................................
start_time <- Sys.time()

model <- brm(formula = species_richness ~ Temperature + Salinity + pH +
               TotalNitrogen_unfiltered + Zooplankton_MaxDepth +
               genus_richness + phyto_evenness + phyto_FEve_ZD + phyto_FDis_ZD + 
               daylight_hours +
               (1 | Year) + (1 | Site_Id), 
             data = plankton_RUE_data_scaled_Zdiv,
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

save(model, file = file.path(model_dir, "BRM_zoo_SR_L_261013_All.RDATA"))

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#...................
start_time <- Sys.time()

model <- brm(formula = species_richness ~ poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) +
               poly(TotalNitrogen_unfiltered,2,raw = TRUE) + poly(Zooplankton_MaxDepth,2,raw = TRUE) +
               genus_richness + phyto_evenness + phyto_FEve_ZD + phyto_FDis_ZD + 
               poly(daylight_hours,2,raw = TRUE) +
               (1 | Year) + (1 | Site_Id), 
             data = plankton_RUE_data_scaled_Zdiv,
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

save(model, file = file.path(model_dir, "BRM_zoo_SR_P_261013_All.RDATA"))

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#Zooplankton shannon evenness #..................................................
start_time <- Sys.time()

model <- brm(formula = zoo_evenness ~ Temperature + Salinity + pH +
               TotalNitrogen_unfiltered + Zooplankton_MaxDepth +
               genus_richness + phyto_evenness + phyto_FEve_ZD + phyto_FDis_ZD + 
               daylight_hours +
               (1 | Year) + (1 | Site_Id), 
             data = plankton_RUE_data_scaled_Zdiv,
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

save(model, file = file.path(model_dir, "BRM_zoo_SE_L_261013_All.RDATA"))

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#...................
start_time <- Sys.time()

model <- brm(formula = zoo_evenness ~ poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) +
               poly(TotalNitrogen_unfiltered,2,raw = TRUE) + poly(Zooplankton_MaxDepth,2,raw = TRUE) +
               genus_richness + phyto_evenness + phyto_FEve_ZD + phyto_FDis_ZD + 
               poly(daylight_hours,2,raw = TRUE) +
               (1 | Year) + (1 | Site_Id), 
             data = plankton_RUE_data_scaled_Zdiv,
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

save(model, file = file.path(model_dir, "BRM_zoo_SE_P_261013_All.RDATA"))

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#Functional Evenness ...........................................................
start_time <- Sys.time()

model <- brm(formula = zoo_FEve ~ Temperature + Salinity + pH +
               TotalNitrogen_unfiltered + Zooplankton_MaxDepth +
               genus_richness + phyto_evenness + phyto_FEve_ZD + phyto_FDis_ZD + 
               daylight_hours +
               (1 | Year) + (1 | Site_Id), 
             data = plankton_RUE_data_scaled_Zdiv,
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

save(model, file = file.path(model_dir, "BRM_zoo_FE_L_261013_All.RDATA"))

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)


#...................
start_time <- Sys.time()

model <- brm(formula = zoo_FEve ~ poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) +
               poly(TotalNitrogen_unfiltered,2,raw = TRUE) + poly(Zooplankton_MaxDepth,2,raw = TRUE) +
               genus_richness + phyto_evenness + phyto_FEve_ZD + phyto_FDis_ZD + 
               poly(daylight_hours,2,raw = TRUE) +
               (1 | Year) + (1 | Site_Id), 
             data = plankton_RUE_data_scaled_Zdiv,
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

save(model, file = file.path(model_dir, "BRM_zoo_FE_P_261013_All.RDATA"))

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#Functional Dispersion .........................................................
start_time <- Sys.time()

model <- brm(formula = zoo_FDis ~ Temperature + Salinity + pH +
               TotalNitrogen_unfiltered + Zooplankton_MaxDepth +
               genus_richness + phyto_evenness + phyto_FEve_ZD + phyto_FDis_ZD + 
               daylight_hours +
               (1 | Year) + (1 | Site_Id), 
             data = plankton_RUE_data_scaled_Zdiv,
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

save(model, file = file.path(model_dir, "BRM_zoo_FD_L_261013_All.RDATA"))

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#...................
start_time <- Sys.time()

model <- brm(formula = zoo_FDis ~ poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) +
               poly(TotalNitrogen_unfiltered,2,raw = TRUE) + poly(Zooplankton_MaxDepth,2,raw = TRUE) +
               genus_richness + phyto_evenness + phyto_FEve_ZD + phyto_FDis_ZD + 
               poly(daylight_hours,2,raw = TRUE) +
               (1 | Year) + (1 | Site_Id), 
             data = plankton_RUE_data_scaled_Zdiv,
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

save(model, file = file.path(model_dir, "BRM_zoo_FD_P_261013_All.RDATA"))

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)


#make BRMS model - Zooplankton diversity with CYANO DOM # ----------------------
#Zooplankton species richness #..................................................
start_time <- Sys.time()

model <- brm(formula = species_richness ~ Temperature + Salinity  + pH +
               TotalNitrogen_unfiltered + Zooplankton_MaxDepth +
               genus_richness + phyto_evenness + phyto_FEve_ZD + phyto_FDis_ZD + cyano_prop +
               daylight_hours +
               (1 | Year) + (1 | Site_Id), 
             data = plankton_RUE_data_scaled_Zdiv,
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

save(model, file = file.path(model_dir, "BRM_zoo_SR_CD_L_261013_All.RDATA"))

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#...................
start_time <- Sys.time()

model <- brm(formula = species_richness ~ poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) +
               poly(TotalNitrogen_unfiltered,2,raw = TRUE) + poly(Zooplankton_MaxDepth,2,raw = TRUE) +
               genus_richness + phyto_evenness + phyto_FEve_ZD + phyto_FDis_ZD + cyano_prop +
               poly(daylight_hours,2,raw = TRUE) +
               (1 | Year) + (1 | Site_Id), 
             data = plankton_RUE_data_scaled_Zdiv,
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

save(model, file = file.path(model_dir, "BRM_zoo_SR_CD_P_261013_All.RDATA"))

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#Zooplankton Shannon evenness#..................................................
start_time <- Sys.time()

model <- brm(formula = zoo_evenness ~ Temperature + Salinity  + pH +
               TotalNitrogen_unfiltered + Zooplankton_MaxDepth +
               genus_richness + phyto_evenness + phyto_FEve_ZD + phyto_FDis_ZD + cyano_prop +
               daylight_hours +
               (1 | Year) + (1 | Site_Id), 
             data = plankton_RUE_data_scaled_Zdiv,
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

save(model, file = file.path(model_dir, "BRM_zoo_SE_CD_L_261013_All.RDATA"))

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#...................
start_time <- Sys.time()

model <- brm(formula = zoo_evenness ~ poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) +
               poly(TotalNitrogen_unfiltered,2,raw = TRUE) + poly(Zooplankton_MaxDepth,2,raw = TRUE) +
               genus_richness + phyto_evenness + phyto_FEve_ZD + phyto_FDis_ZD + cyano_prop +
               poly(daylight_hours,2,raw = TRUE) +
               (1 | Year) + (1 | Site_Id), 
             data = plankton_RUE_data_scaled_Zdiv,
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

save(model, file = file.path(model_dir, "BRM_zoo_SE_CD_P_261013_All.RDATA"))

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#Functional Evenness ...........................................................
start_time <- Sys.time()

model <- brm(formula = zoo_FEve ~ Temperature + Salinity  + pH +
               TotalNitrogen_unfiltered + Zooplankton_MaxDepth +
               genus_richness + phyto_evenness + phyto_FEve_ZD + phyto_FDis_ZD + cyano_prop +
               daylight_hours +
               (1 | Year) + (1 | Site_Id), 
             data = plankton_RUE_data_scaled_Zdiv,
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

save(model, file = file.path(model_dir, "BRM_zoo_FE_CD_L_261013_All.RDATA"))

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#...................
start_time <- Sys.time()

model <- brm(formula = zoo_FEve ~ poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) +
               poly(TotalNitrogen_unfiltered,2,raw = TRUE) + poly(Zooplankton_MaxDepth,2,raw = TRUE) +
               genus_richness + phyto_evenness + phyto_FEve_ZD + phyto_FDis_ZD + cyano_prop +
               poly(daylight_hours,2,raw = TRUE) +
               (1 | Year) + (1 | Site_Id), 
             data = plankton_RUE_data_scaled_Zdiv,
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

save(model, file = file.path(model_dir, "BRM_zoo_FE_CD_P_261013_All.RDATA"))

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#Functional Dispersion .........................................................
start_time <- Sys.time()

model <- brm(formula = zoo_FDis ~ Temperature + Salinity  + pH +
               TotalNitrogen_unfiltered + Zooplankton_MaxDepth +
               genus_richness + phyto_evenness + phyto_FEve_ZD + phyto_FDis_ZD + cyano_prop +
               daylight_hours +
               (1 | Year) + (1 | Site_Id), 
             data = plankton_RUE_data_scaled_Zdiv,
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

save(model, file = file.path(model_dir, "BRM_zoo_FD_CD_L_261013_All.RDATA"))

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#...................
start_time <- Sys.time()

model <- brm(formula = zoo_FDis ~ poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) +
               poly(TotalNitrogen_unfiltered,2,raw = TRUE) + poly(Zooplankton_MaxDepth,2,raw = TRUE) +
               genus_richness + phyto_evenness + phyto_FEve_ZD + phyto_FDis_ZD + cyano_prop +
               poly(daylight_hours,2,raw = TRUE) +
               (1 | Year) + (1 | Site_Id), 
             data = plankton_RUE_data_scaled_Zdiv,
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

save(model, file = file.path(model_dir, "BRM_zoo_FD_CD_P_261013_All.RDATA"))

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#make BRMS model- zooplankton RUE#----------------------------------------------
start_time <- Sys.time()

model <- brm(formula = log_RUE_zooplankton ~ Temperature + Salinity  + pH +
               TotalNitrogen_unfiltered + Zooplankton_MaxDepth +
               genus_richness + phyto_evenness + phyto_FEve_ZRUE + phyto_FDis_ZRUE +
               species_richness + zoo_evenness + zoo_FEve + zoo_FDis + 
               daylight_hours +
               (1 | Year) + (1 | Site_Id), 
             data = plankton_RUE_data_scaled_RUE,
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

save(model, file = file.path(model_dir, "BRM_zoo_RUE_L_261013_All.RDATA"))

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#.........................

start_time <- Sys.time()

model <- brm(formula = log_RUE_zooplankton ~ poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) +
               poly(TotalNitrogen_unfiltered,2,raw = TRUE) + poly(Zooplankton_MaxDepth,2,raw = TRUE) +
               genus_richness + phyto_evenness + phyto_FEve_ZRUE + phyto_FDis_ZRUE + 
               species_richness + zoo_evenness + zoo_FEve + zoo_FDis + 
               poly(daylight_hours,2,raw = TRUE) +
               (1 | Year) + (1 | Site_Id), 
             data = plankton_RUE_data_scaled_RUE,
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

save(model, file = file.path(model_dir, "BRM_zoo_RUE_P_261013_All.RDATA"))

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#...........................

start_time <- Sys.time()

model <- brm(formula = log_RUE_zooplankton ~ pH + Zooplankton_MaxDepth +
               daylight_hours + 
               genus_richness*Temperature + genus_richness*Salinity + genus_richness*TotalNitrogen_unfiltered +
               phyto_evenness*Temperature + phyto_evenness*Salinity + phyto_evenness*TotalNitrogen_unfiltered +
               phyto_FEve_ZRUE*Temperature + phyto_FEve_ZRUE*Salinity + phyto_FEve_ZRUE*TotalNitrogen_unfiltered + 
               phyto_FDis_ZRUE*Temperature + phyto_FDis_ZRUE*Salinity + phyto_FDis_ZRUE*TotalNitrogen_unfiltered + 
               species_richness*Temperature + species_richness*Salinity + species_richness*TotalNitrogen_unfiltered +
               zoo_evenness*Temperature + zoo_evenness*Salinity + zoo_evenness*TotalNitrogen_unfiltered +
               zoo_FEve*Temperature + zoo_FEve*Salinity + zoo_FEve*TotalNitrogen_unfiltered +
               zoo_FDis*Temperature + zoo_FDis*Salinity + zoo_FDis*TotalNitrogen_unfiltered +
               (1 | Year) + (1 | Site_Id),  
             data = plankton_RUE_data_scaled_RUE,
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

save(model, file = file.path(model_dir, "BRM_zoo_RUE_LI_261013_All.RDATA"))

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#.........................

start_time <- Sys.time()

model <- brm(formula = log_RUE_zooplankton ~ 
               #env variables
               pH + I(pH^2) + Zooplankton_MaxDepth + I(Zooplankton_MaxDepth^2) + 
               daylight_hours + I(daylight_hours^2) + Temperature + I(Temperature^2) +
               Salinity + I(Salinity^2) + TotalNitrogen_unfiltered + I(TotalNitrogen_unfiltered^2) +
               # phytoplankton biodiversity
               genus_richness * (Temperature + Salinity + TotalNitrogen_unfiltered) +
               phyto_evenness * (Temperature + Salinity + TotalNitrogen_unfiltered) +
               phyto_FEve_ZRUE * (Temperature + Salinity + TotalNitrogen_unfiltered) +
               phyto_FDis_ZRUE * (Temperature + Salinity + TotalNitrogen_unfiltered) +
               # zooplankton biodiversity
               species_richness * (Temperature + Salinity + TotalNitrogen_unfiltered) +
               zoo_evenness * (Temperature + Salinity + TotalNitrogen_unfiltered) +
               zoo_FEve * (Temperature + Salinity + TotalNitrogen_unfiltered) +
               zoo_FDis * (Temperature + Salinity + TotalNitrogen_unfiltered) +
               # random effects
               (1 | Year) + (1 | Site_Id),
             data = plankton_RUE_data_scaled_RUE,
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

save(model, file = file.path(model_dir, "BRM_zoo_RUE_PI_261013_All.RDATA"))

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)


#make BRMS model- zooplankton RUE with CYANO DOM#-------------------------------
start_time <- Sys.time()

model <- brm(formula = log_RUE_zooplankton ~ Temperature + Salinity  + pH +
               TotalNitrogen_unfiltered + Zooplankton_MaxDepth +
               genus_richness + phyto_evenness + phyto_FEve_ZRUE + phyto_FDis_ZRUE + cyano_prop +
               species_richness + zoo_evenness + zoo_FEve + zoo_FDis + 
               daylight_hours +
               (1 | Year) + (1 | Site_Id), 
             data = plankton_RUE_data_scaled_RUE,
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

save(model, file = file.path(model_dir, "BRM_zoo_RUE_CD_L_261013_All.RDATA"))

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#....................
start_time <- Sys.time()

model <- brm(formula = log_RUE_zooplankton ~ poly(Temperature,2,raw = TRUE) + 
               poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) +
               poly(TotalNitrogen_unfiltered,2,raw = TRUE) + poly(Zooplankton_MaxDepth,2,raw = TRUE) +
               genus_richness + phyto_evenness + phyto_FEve_ZRUE + phyto_FDis_ZRUE + cyano_prop +
               species_richness + zoo_evenness + zoo_FEve + zoo_FDis +
               poly(daylight_hours,2,raw = TRUE) +
               (1 | Year) + (1 | Site_Id), 
             data = plankton_RUE_data_scaled_RUE,
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

save(model, file = file.path(model_dir, "BRM_zoo_RUE_CD_P_261013_All.RDATA"))

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#.......................

start_time <- Sys.time()

model <- brm(formula = log_RUE_zooplankton ~ pH + Zooplankton_MaxDepth +
               daylight_hours +
               genus_richness*Temperature + genus_richness*Salinity + genus_richness*TotalNitrogen_unfiltered +
               phyto_evenness*Temperature + phyto_evenness*Salinity + phyto_evenness*TotalNitrogen_unfiltered +
               phyto_FEve_ZRUE*Temperature + phyto_FEve_ZRUE*Salinity + phyto_FEve_ZRUE*TotalNitrogen_unfiltered + 
               phyto_FDis_ZRUE*Temperature + phyto_FDis_ZRUE*Salinity + phyto_FDis_ZRUE*TotalNitrogen_unfiltered + 
               cyano_prop*Temperature + cyano_prop*Salinity + cyano_prop*TotalNitrogen_unfiltered +
               species_richness*Temperature + species_richness*Salinity + species_richness*TotalNitrogen_unfiltered +
               zoo_evenness*Temperature + zoo_evenness*Salinity + zoo_evenness*TotalNitrogen_unfiltered +
               zoo_FEve*Temperature + zoo_FEve*Salinity + zoo_FEve*TotalNitrogen_unfiltered +
               zoo_FDis*Temperature + zoo_FDis*Salinity + zoo_FDis*TotalNitrogen_unfiltered +
               (1 | Year) + (1 | Site_Id),  
             data = plankton_RUE_data_scaled_RUE,
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

save(model, file = file.path(model_dir, "BRM_zoo_RUE_CD_LI_261013_All.RDATA"))

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#....................
start_time <- Sys.time()

model <- brm(formula = log_RUE_zooplankton ~ 
               #env variables
               Temperature + I(Temperature^2) + Salinity + I(Salinity^2) + 
               pH + I(pH^2) + Zooplankton_MaxDepth + I(Zooplankton_MaxDepth^2) + 
               TotalNitrogen_unfiltered + I(TotalNitrogen_unfiltered^2) + 
               daylight_hours + I(daylight_hours^2) + 
               # phytoplankton biodiversity
               genus_richness * (Temperature + Salinity + TotalNitrogen_unfiltered) +
               phyto_evenness * (Temperature + Salinity + TotalNitrogen_unfiltered) +
               phyto_FEve_ZRUE * (Temperature + Salinity + TotalNitrogen_unfiltered) +
               phyto_FDis_ZRUE * (Temperature + Salinity + TotalNitrogen_unfiltered) +
               cyano_prop * (Temperature + Salinity + TotalNitrogen_unfiltered) +
               # zooplankton biodiversity
               species_richness * (Temperature + Salinity + TotalNitrogen_unfiltered) +
               zoo_evenness * (Temperature + Salinity + TotalNitrogen_unfiltered) +
               zoo_FEve * (Temperature + Salinity + TotalNitrogen_unfiltered) +
               zoo_FDis * (Temperature + Salinity + TotalNitrogen_unfiltered) +
               # random effects
               (1 | Year) + (1 | Site_Id),
             data = plankton_RUE_data_scaled_RUE,
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

save(model, file = file.path(model_dir, "BRM_zoo_RUE_CD_PI_261013_All.RDATA"))

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

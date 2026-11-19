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
setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/spatial_sensitivity_analysis/Larger_radius/RDATA")
load('plankton_RUE_data_scaled_RUE_260526.RDATA')

setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/spatial_sensitivity_analysis/Larger_radius/RDATA")

project_dir <- "~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/spatial_sensitivity_analysis"
model_dir <- file.path(
  project_dir,
  "Larger_radius/RDATA"
)

#load this when using remote desktop
# Local directory for Stan compilation
stan_dir <- "C:/R_temp"
if (!dir.exists(stan_dir)) {
  dir.create(stan_dir, recursive = TRUE)
}
# Use local disk for Stan compilation
setwd(stan_dir)

#make BRMS model- phytoplankton RUE with CYANO DOM# ----------------------------
#Phytoplankton RUE with phosphorus .............................................

start_time <- Sys.time()

model <- brm(formula = log_RUE_B_P ~ 
               #env variables
               Temperature + I(Temperature^2) + Salinity + I(Salinity^2) + 
               pH + I(pH^2) + MaxDepth + I(MaxDepth^2) + 
               TotalNitrogen_unfiltered + I(TotalNitrogen_unfiltered^2) + 
               daylight_hours + I(daylight_hours^2) + 
               # phytoplankton biodiversity
               genus_richness * (Temperature + Salinity + TotalNitrogen_unfiltered) +
               phyto_evenness * (Temperature + Salinity + TotalNitrogen_unfiltered) +
               phyto_FEve_PRUE * (Temperature + Salinity + TotalNitrogen_unfiltered) +
               phyto_FDis_PRUE * (Temperature + Salinity + TotalNitrogen_unfiltered) +
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

save(model, file = file.path(model_dir, "BRM_phytoRUEP_CD_PI_260526.RDATA"))

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#Phytoplankton RUE with nitrogen ...............................................

start_time <- Sys.time()

model <- brm(formula = log_RUE_B_N ~ 
               #env variables
               Temperature + I(Temperature^2) + Salinity + I(Salinity^2) + 
               pH + I(pH^2) + MaxDepth + I(MaxDepth^2) + 
               TotalNitrogen_unfiltered + I(TotalNitrogen_unfiltered^2) + 
               daylight_hours + I(daylight_hours^2) + 
               # phytoplankton biodiversity
               genus_richness * (Temperature + Salinity + TotalNitrogen_unfiltered) +
               phyto_evenness * (Temperature + Salinity + TotalNitrogen_unfiltered) +
               phyto_FEve_PRUE * (Temperature + Salinity + TotalNitrogen_unfiltered) +
               phyto_FDis_PRUE * (Temperature + Salinity + TotalNitrogen_unfiltered) +
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

save(model, file = file.path(model_dir, "BRM_phytoRUEN_CD_PI_260526.RDATA"))

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)


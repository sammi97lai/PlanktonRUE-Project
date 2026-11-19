################################################################################
#######################    Phytoplankton Div BRMS   ############################
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
load('plankton_RUE_data_scaled_260526.RDATA')

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

#make BRMS model- phytoplankton diversity #-----------------------------------
#genus richness...............................................................

start_time <- Sys.time()

model <- brm( formula= genus_richness ~ poly(Temperature,2,raw = TRUE) + 
                poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) + 
                poly(TotalNitrogen_unfiltered,2,raw = TRUE) + poly(MaxDepth,2,raw = TRUE) + 
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

save(model, file = file.path(model_dir, "BRM_phyto_GD_P_260526.RDATA"))

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#shannon evenness...............................................................

start_time <- Sys.time()

model <- brm( formula= phyto_evenness ~ poly(Temperature,2,raw = TRUE) + 
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

save(model, file = file.path(model_dir, "BRM_phyto_SE_P_260526.RDATA"))

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#Functional Evenness ...........................................................

start_time <- Sys.time()

model <- brm( formula= phyto_FEve_PD ~ poly(Temperature,2,raw = TRUE) + 
                poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) + 
                poly(TotalNitrogen_unfiltered,2,raw = TRUE) + poly(MaxDepth,2,raw = TRUE) +
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

save(model, file = file.path(model_dir, "BRM_phyto_FEve_P_260526.RDATA"))

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)


#Functional Dispersion .........................................................

start_time <- Sys.time()

model <- brm( formula= phyto_FDis_PD ~ poly(Temperature,2,raw = TRUE) + 
                poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) + 
                poly(TotalNitrogen_unfiltered,2,raw = TRUE) + poly(MaxDepth,2,raw = TRUE) +
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

save(model, file = file.path(model_dir, "BRM_phyto_FDis_P_260526.RDATA"))

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)


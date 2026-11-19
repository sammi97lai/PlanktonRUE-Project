################################################################################
###################    Cyanobacteria Dominance BRMS   ##########################
################################################################################

library(dplyr)
library(tidyr)
library(brms)

rm(list = ls())
#------------------------------------------------------------------------------
setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/RDATA")
load('plankton_RUE_data_scaled_260916.RDATA')

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

#make BRMS model - Cyanobacteria dominance- zero inflated #--------------------
start_time <- Sys.time()

model <- brm( formula= cyano_prop ~ Temperature + 
                Salinity + pH + TotalNitrogen_unfiltered + MaxDepth + 
                daylight_hours +
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

save(model, file = file.path(model_dir, "BRM_cyanodom_L_261009_z.RDATA"))

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#................

start_time <- Sys.time()

model <- brm( formula= cyano_prop ~ poly(Temperature,2,raw = TRUE) + 
                poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) + 
                poly(TotalNitrogen_unfiltered,2,raw = TRUE) + poly(MaxDepth,2,raw = TRUE) +
                poly(daylight_hours,2,raw = TRUE) +
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

save(model, file = file.path(model_dir, "BRM_cyanodom_P_261009_z.RDATA"))

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)




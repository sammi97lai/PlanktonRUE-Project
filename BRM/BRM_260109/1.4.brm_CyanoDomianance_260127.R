################################################################################
###################    Cyanobacteria Dominance BRMS   ##########################
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

#make BRMS model - Cyanobacteria dominance- zero inflated #--------------------
start_time <- Sys.time()

model <- brm( formula= cyano_prop ~ Temperature + 
                Salinity + pH + log_N_P_ratio + DOY +
                (1 | Year) + (1 | Site_Id) + (1 | MaxDepth), 
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

save(model, file = "BRM_cyanodom_L_260109_z.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#................

start_time <- Sys.time()

model <- brm( formula= cyano_prop ~ poly(Temperature,2,raw = TRUE) + 
                poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) + 
                poly(log_N_P_ratio,2,raw = TRUE) +
                poly(DOY,2,raw = TRUE) + 
                (1 | Year) + (1 | Site_Id) + (1 | MaxDepth), 
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

save(model, file = "BRM_cyanodom_P_260109_z.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#.............
start_time <- Sys.time()

model <- brm( formula= cyano_prop ~ Temperature + 
                Salinity + pH + log_N_P_ratio*DOY + 
                (1 | Year) + (1 | Site_Id) + (1 | MaxDepth), 
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

save(model, file = "BRM_cyanodom_LI_260109_z.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#.............
start_time <- Sys.time()

model <- brm( formula= cyano_prop ~ poly(Temperature,2,raw = TRUE) + 
                poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) + 
                poly(log_N_P_ratio, 2, raw = TRUE) * poly(DOY, 2, raw = TRUE) +
                (1 | Year) + (1 | Site_Id) + (1 | MaxDepth), 
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

save(model, file = "BRM_cyanodom_PI_260109_z.RData")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)




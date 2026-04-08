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
setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/RDATA")
load('plankton_RUE_data_scaled_260408.RDATA')

setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA")

#make BRMS model- phytoplankton diversity #-----------------------------------
#genus richness...............................................................
start_time <- Sys.time()

model <- brm( formula= genus_richness ~ Temperature + 
                Salinity + pH + TotalNitrogen_unfiltered + MaxDepth + 
                DOY + daylight_hours +
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

save(model, file = "BRM_phyto_GD_L_260408.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#................

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

save(model, file = "BRM_phyto_GD_P_260408.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#Functional Evenness ...........................................................
start_time <- Sys.time()

model <- brm( formula= FEve ~ Temperature + 
                Salinity + pH + TotalNitrogen_unfiltered + MaxDepth + 
                DOY + daylight_hours +
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

save(model, file = "BRM_phyto_FEve_L_260408.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)


#.................

start_time <- Sys.time()

model <- brm( formula= FEve ~ poly(Temperature,2,raw = TRUE) + 
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

save(model, file = "BRM_phyto_FEve_P_260408.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)


#Functional Dispersion .........................................................
start_time <- Sys.time()

model <- brm( formula= FDis ~ Temperature + 
                Salinity + pH + TotalNitrogen_unfiltered + MaxDepth + 
                DOY + daylight_hours +
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

save(model, file = "BRM_phyto_FDis_L_260408.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#.................

start_time <- Sys.time()

model <- brm( formula= FDis ~ poly(Temperature,2,raw = TRUE) + 
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

save(model, file = "BRM_phyto_FDis_P_260408.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)


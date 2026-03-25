################################################################################
#######################    RUE Phytoplankton BRMS   ############################
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
load('plankton_RUE_data6_260309.RDATA')
### calculate N:P ratio........................................................
plankton_RUE_data7 <- plankton_RUE_data6 %>% mutate(N_P_ratio = TotalNitrogen_unfiltered / TotalPhosphorus_unfiltered)
colnames(plankton_RUE_data7)

#### check distribution of variables ...........................................
#response
# hist(plankton_RUE_data7$genus_richness)
# hist(plankton_RUE_data7$FEve)
# hist(plankton_RUE_data7$FDis)
# hist(plankton_RUE_data7$cyano_prop)
plankton_RUE_data7$log_cyano_prop <- log(plankton_RUE_data7$cyano_prop+0.0001)
# hist(plankton_RUE_data7$log_cyano_prop)
# hist(plankton_RUE_data7$species_richness)
# hist(plankton_RUE_data7$zoo_FEve)
# hist(plankton_RUE_data7$zoo_FDis)
# hist(plankton_RUE_data7$log_RUE_B_P)
# hist(plankton_RUE_data7$log_RUE_B_N)
# hist(plankton_RUE_data7$log_RUE_zooplankton)


#predictors
# hist(plankton_RUE_data7$Temperature)
# hist(plankton_RUE_data7$Salinity)
# hist(plankton_RUE_data7$pH)

# hist(plankton_RUE_data7$CWM_log_Biovolume)
# hist(sqrt(plankton_RUE_data7$TotalPhosphorus_unfiltered))
# hist(log(plankton_RUE_data7$TotalNitrogen_unfiltered))
# hist(plankton_RUE_data7$N_P_ratio)
plankton_RUE_data7$log_N_P_ratio <- log(plankton_RUE_data7$N_P_ratio)
# hist(plankton_RUE_data7$log_N_P_ratio)

plankton_RUE_data7$MaxDepth <- as.factor(plankton_RUE_data7$MaxDepth)

# hist(plankton_RUE_data7$MinDepth)
# hist(plankton_RUE_data7$MaxDepth)
# hist(plankton_RUE_data7$DOY)
# hist(plankton_RUE_data7$Zooplankton_MaxDepth)

plankton_RUE_data7$Zooplankton_MaxDepth <- as.factor(plankton_RUE_data7$Zooplankton_MaxDepth)

#test skewness
# library(moments)
# skewness(log(plankton_RUE_data7$RUE_C_N))
# skewness(sqrt(plankton_RUE_data7$MaxDepth))

save(plankton_RUE_data7, file = 'plankton_RUE_data7_old_260324.RDATA')
# write.csv(plankton_RUE_data7, "plankton_RUE_data_260324.csv", row.names = FALSE)

#### make everything scaled and add poly terms .................................
load('plankton_RUE_data7_old_260324.RDATA')

plankton_RUE_data_scaled <- plankton_RUE_data7 %>%
  mutate(across(c(Temperature, Salinity, pH, log_N_P_ratio, DOY),
                ~ as.numeric(scale(.))))

save(plankton_RUE_data_scaled, file = 'plankton_RUE_data_scaled_old_260324.RDATA')
# write.csv(plankton_RUE_data_scaled, "plankton_RUE_data_scaled_260324.csv", row.names = FALSE)

setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/RDATA")
load('plankton_RUE_data_scaled_old_260324.RDATA')

setwd("~/GitHub/PlanktonRUE-Project/BRM_old_model/RDATA")

#make BRMS model- phytoplankton diversity #-----------------------------------
#genus richness...............................................................
start_time <- Sys.time()

model <- brm( formula= genus_richness ~ Temperature + 
                Salinity + pH + log_N_P_ratio + DOY +
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

save(model, file = "BRM_phyto_GD_L_260324.RData")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#................

start_time <- Sys.time()

model <- brm( formula= genus_richness ~ poly(Temperature,2,raw = TRUE) + 
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

save(model, file = "BRM_phyto_GD_P_260324.RData")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#.............
start_time <- Sys.time()

model <- brm( formula= genus_richness ~ Temperature + 
                Salinity + pH + log_N_P_ratio*DOY + 
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

save(model, file = "BRM_phyto_GD_LI_260324.RData")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#.............
start_time <- Sys.time()

model <- brm( formula= genus_richness ~ poly(Temperature,2,raw = TRUE) + 
                poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) + 
                poly(log_N_P_ratio, 2, raw = TRUE) * poly(DOY, 2, raw = TRUE) +
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

save(model, file = "BRM_phyto_GD_PI_260324.RData")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#Functional Evenness ...........................................................
start_time <- Sys.time()

model <- brm( formula= FEve ~ Temperature + 
                Salinity + pH + log_N_P_ratio + DOY + 
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

save(model, file = "BRM_phyto_FEve_L_260324.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)


#.................

start_time <- Sys.time()

model <- brm( formula= FEve ~ poly(Temperature,2,raw = TRUE) + 
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

save(model, file = "BRM_phyto_FEve_P_260324.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)
#1 divergent

#.............
start_time <- Sys.time()

model <- brm( formula= FEve ~ Temperature + 
                Salinity + pH + log_N_P_ratio*DOY + 
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

save(model, file = "BRM_phyto_FEve_LI_260324.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#3 divergence

#.............
start_time <- Sys.time()

model <- brm( formula= FEve ~ poly(Temperature,2,raw = TRUE) + 
                poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) +
                poly(log_N_P_ratio, 2, raw = TRUE) * poly(DOY, 2, raw = TRUE) +
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

save(model, file = "BRM_phyto_FEve_PI_260324.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)
#1 divergence

#Functional Dispersion .........................................................
start_time <- Sys.time()

model <- brm( formula= FDis ~ Temperature + 
                Salinity + pH + log_N_P_ratio + DOY + 
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

save(model, file = "BRM_phyto_FDis_L_260324.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#.................

start_time <- Sys.time()

model <- brm( formula= FDis ~ poly(Temperature,2,raw = TRUE) + 
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

save(model, file = "BRM_phyto_FDis_P_260324.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#.............
start_time <- Sys.time()

model <- brm( formula= FDis ~ Temperature + 
                Salinity + pH + log_N_P_ratio*DOY + 
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

save(model, file = "BRM_phyto_FDis_LI_260324.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

#.............
start_time <- Sys.time()

model <- brm( formula= FDis ~ poly(Temperature,2,raw = TRUE) + 
                poly(Salinity,2,raw = TRUE) + poly(pH,2,raw = TRUE) + 
                poly(log_N_P_ratio, 2, raw = TRUE) * poly(DOY, 2, raw = TRUE) +
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

save(model, file = "BRM_phyto_FDis_PI_260324.RDATA")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)


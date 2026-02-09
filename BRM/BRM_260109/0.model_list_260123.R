#############List of models ############################

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

save(model, file = "BRM_phyto_GD_PI_260109.RData")
#......................................
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

save(model, file = "BRM_phytoRUEP_CD_P_260109.RData")
#.....................................................
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
#.......................................................
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
#........................................................
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
# ===========================================================
# Salinity and NP ratio Model list 
# ===========================================================
model_list <- list(
  #PHYTOPLANKTON DIVERSITY
  # ---------------------------------------------------------
  # PHYTOPLANKTON GENUS RICHNESS MODELS
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260330/RDATA/BRM_phyto_GD_P_260330_S_NP.RDATA",
       name = "Phyto_Genus_Richness_P_S_NP",
       resp = "genus_richness",
       predictors = c("Temperature", "Salinity", "pH", "log_N_P_ratio", "MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic Variables for Phytoplankton Genus Richness (Salinity and NP ratio)"),
  
  # ---------------------------------------------------------
  # PHYTOPLANKTON FUNCTIONAL EVENNESS MODELS
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260330/RDATA/BRM_phyto_FEve_P_260330_S_NP.RDATA",
       name = "Phyto_FEve_P_S_NP",
       resp = "FEve",
       predictors = c("Temperature", "Salinity", "pH", "log_N_P_ratio", "MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic Variables for Phytoplankton Functional Evenness (Salinity and NP ratio)"),
  
  # ---------------------------------------------------------
  # PHYTOPLANKTON FUNCTIONAL DISPERSION MODELS
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260330/RDATA/BRM_phyto_FDis_P_260330_S_NP.RDATA",
       name = "Phyto_FDis_P_S_NP",
       resp = "FDis",
       predictors = c("Temperature", "Salinity", "pH", "log_N_P_ratio", "MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic Variables for Phytoplankton Functional Dispersion (Salinity and NP ratio)"),
  
  # ---------------------------------------------------------
  # CYANOBACTERIA DOMINANCE MODELS
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260330/RDATA/BRM_cyanodom_P_260330_S_NP.RDATA",
       name = "cyanodom_P_S_NP",
       resp = "cyano_prop",
       predictors = c("Temperature", "Salinity", "pH", "log_N_P_ratio", "MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic Variables for Cyanobacteria Dominance (Salinity and NP ratio)"),  
  
  # ---------------------------------------------------------
  # PHYTOPLANKTON RUE WITH PHOSPHORUS MODELS WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260330/RDATA/BRM_phytoRUEP_CD_PI_260330_S_NP.RDATA",
       name = "phytoRUEP_CD_PI_S_NP",
       resp = "log_RUE_B_P",
       predictors = c("Temperature", "Salinity", "pH", "log_N_P_ratio", "MaxDepth", "genus_richness",
                      "FEve", "FDis", "cyano_prop", "species_richness", "zoo_FEve", "zoo_FDis",
                      "genus_richness:Temperature", "genus_richness:Salinity", "genus_richness:log_N_P_ratio",
                      "FEve:Temperature", "FEve:Salinity", "FEve:log_N_P_ratio", 
                      "FDis:Temperature", "FDis:Salinity", "FDis:log_N_P_ratio",
                      "cyano_prop:Temperature", "cyano_prop:Salinity", "cyano_prop:log_N_P_ratio",
                      "species_richness:Temperature", "species_richness:Salinity", "species_richness:log_N_P_ratio",
                      "zoo_FEve:Temperature", "zoo_FEve:Salinity", "zoo_FEve:log_N_P_ratio",
                      "zoo_FDis:Temperature", "zoo_FDis:Salinity", "zoo_FDis:log_N_P_ratio",
                      "DOY", "daylight_hours"),
       ncol_ce = 6,
       ce_title = "Marginal Effects of Abiotic AND Biotic Variables for Phytoplankton RUE per P (Salinity and NP ratio)"),
  
  # ---------------------------------------------------------
  # PHYTOPLANKTON RUE WITH NITROGEN MODELS WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260330/RDATA/BRM_phytoRUEN_CD_PI_260330_S_NP.RDATA",
       name = "phytoRUEN_CD_PI_S_NP",
       resp = "log_RUE_B_N",
       predictors = c("Temperature", "Salinity", "pH", "log_N_P_ratio", "MaxDepth", "genus_richness",
                      "FEve", "FDis", "cyano_prop", "species_richness", "zoo_FEve", "zoo_FDis",
                      "genus_richness:Temperature", "genus_richness:Salinity", "genus_richness:log_N_P_ratio",
                      "FEve:Temperature", "FEve:Salinity", "FEve:log_N_P_ratio", 
                      "FDis:Temperature", "FDis:Salinity", "FDis:log_N_P_ratio",
                      "cyano_prop:Temperature", "cyano_prop:Salinity", "cyano_prop:log_N_P_ratio",
                      "species_richness:Temperature", "species_richness:Salinity", "species_richness:log_N_P_ratio",
                      "zoo_FEve:Temperature", "zoo_FEve:Salinity", "zoo_FEve:log_N_P_ratio",
                      "zoo_FDis:Temperature", "zoo_FDis:Salinity", "zoo_FDis:log_N_P_ratio",
                      "DOY", "daylight_hours"),
       ncol_ce = 6,
       ce_title = "Marginal Effects of Abiotic AND Biotic Variables for Phytoplankton RUE per N (Salinity and NP ratio)"),
  
  # ---------------------------------------------------------
  # ZOOPLANKTON SPECIES RICHNESS MODELS WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260330/RDATA/BRM_zoo_SR_CD_P_260330_S_NP.RDATA",
       name = "Zoo_SR_CD_P_S_NP",
       resp = "species_richness",
       predictors = c("Temperature", "Salinity", "pH", "log_N_P_ratio", "Zooplankton_MaxDepth", 
                      "genus_richness", "FEve_selected", "FDis_selected", "cyano_prop",
                      "DOY", "daylight_hours"),
       ncol_ce = 4,
       ce_title = "Marginal Effects of Abiotic and Biotic Variables for Zooplankton Species Richness (Salinity and NP ratio)"),
  
  # ---------------------------------------------------------
  # ZOOPLANKTON FUNCTIONAL EVENNESS MODELS WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260330/RDATA/BRM_zoo_FE_CD_P_260330_S_NP.RDATA",
       name = "Zoo_FEve_CD_P_S_NP",
       resp = "zoo_FEve",
       predictors = c("Temperature", "Salinity", "pH", "log_N_P_ratio", "Zooplankton_MaxDepth", 
                      "genus_richness", "FEve_selected", "FDis_selected", "cyano_prop",
                      "DOY", "daylight_hours"),
       ncol_ce = 4,
       ce_title = "Marginal Effects of Abiotic and Biotic Variables for Zooplankton Functional Evenness (Salinity and NP ratio)"),
  
  # ---------------------------------------------------------
  # ZOOPLANKTON FUNCTIONAL DISPERSION MODELS WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260330/RDATA/BRM_zoo_FD_CD_P_260330_S_NP.RDATA",
       name = "Zoo_FDis_CD_P_S_NP",
       resp = "zoo_FDis",
       predictors = c("Temperature", "Salinity", "pH", "log_N_P_ratio", "Zooplankton_MaxDepth", 
                      "genus_richness", "FEve_selected", "FDis_selected", "cyano_prop",
                      "DOY", "daylight_hours"),
       ncol_ce = 4,
       ce_title = "Marginal Effects of Abiotic and Biotic Variables for Zooplankton Functional Dispersion (Salinity and NP ratio)"),
  
  # ---------------------------------------------------------
  # ZOOPLANKTON RUE WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260330/RDATA/BRM_zoo_RUE_CD_PI_260330_S_NP.RDATA",
       name = "zoo_RUE_CD_PI_S_NP",
       resp = "log_RUE_zooplankton",
       predictors = c("Temperature", "Salinity", "pH", "log_N_P_ratio", "Zooplankton_MaxDepth", "genus_richness",
                      "FEve_selected", "FDis_selected",  "cyano_prop",
                      "species_richness", "zoo_FEve", "zoo_FDis",
                      "genus_richness:Temperature", "genus_richness:Salinity", "genus_richness:log_N_P_ratio",
                      "FEve_selected:Temperature", "FEve_selected:Salinity", "FEve_selected:log_N_P_ratio", 
                      "FDis_selected:Temperature", "FDis_selected:Salinity", "FDis_selected:log_N_P_ratio",
                      "cyano_prop:Temperature", "cyano_prop:Salinity", "cyano_prop:log_N_P_ratio",
                      "species_richness:Temperature", "species_richness:Salinity", "species_richness:log_N_P_ratio",
                      "zoo_FEve:Temperature", "zoo_FEve:Salinity", "zoo_FEve:log_N_P_ratio", 
                      "zoo_FDis:Temperature", "zoo_FDis:Salinity", "zoo_FDis:log_N_P_ratio",
                      "DOY", "daylight_hours"),
       ncol_ce = 6,
       ce_title = "Marginal Effects of Abiotic and Biotic Variables for Zooplankton RUE (Salinity and NP ratio)")
  
)
  


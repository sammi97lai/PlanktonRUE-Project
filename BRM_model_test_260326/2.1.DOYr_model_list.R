# ===========================================================
# DOY-random effect Model list 
# ===========================================================
model_list <- list(
  #PHYTOPLANKTON DIVERSITY
  # ---------------------------------------------------------
  # PHYTOPLANKTON GENUS RICHNESS MODELS
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260326/RDATA/BRM_phyto_GD_P_260326_DOYr.RData",
       name = "Phyto_Genus_Richness_P_DOYr",
       resp = "genus_richness",
       predictors = c("Temperature", "Salinity", "pH", "MaxDepth", "daylight_hours"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic Variables for Phytoplankton Genus Richness (DOY-random)"),
  
  # ---------------------------------------------------------
  # PHYTOPLANKTON FUNCTIONAL EVENNESS MODELS
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260326/RDATA/BRM_phyto_FEve_P_260326_DOYr.RDATA",
       name = "Phyto_FEve_P_DOYr",
       resp = "FEve",
       predictors = c("Temperature", "Salinity", "pH", "MaxDepth", "daylight_hours"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic Variables for Phytoplankton Functional Evenness (DOY-random)"),
  
  # ---------------------------------------------------------
  # PHYTOPLANKTON FUNCTIONAL DISPERSION MODELS
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260326/RDATA/BRM_phyto_FDis_P_260326_DOYr.RDATA",
       name = "Phyto_FDis_P_DOYr",
       resp = "FDis",
       predictors = c("Temperature", "Salinity", "pH", "MaxDepth", "daylight_hours"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic Variables for Phytoplankton Functional Dispersion (DOY-random)"),
  
  # ---------------------------------------------------------
  # CYANOBACTERIA DOMINANCE MODELS
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260326/RDATA/BRM_cyanodom_P_260326_DOYr.RDATA",
       name = "cyanodom_P_DOYr",
       resp = "cyano_prop",
       predictors = c("Temperature", "Salinity", "pH", "MaxDepth", "daylight_hours"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic Variables for Cyanobacteria Dominance (DOY-random)"),  
  
  # ---------------------------------------------------------
  # PHYTOPLANKTON RUE WITH PHOSPHORUS MODELS WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260326/RDATA/BRM_phytoRUEP_CD_PI_260326_DOYr.RDATA",
       name = "phytoRUEP_CD_PI_DOYr",
       resp = "log_RUE_B_P",
       predictors = c("Temperature", "Salinity", "pH", "MaxDepth", "genus_richness",
                      "FEve", "FDis", "cyano_prop", "species_richness", "zoo_FEve", "zoo_FDis",
                      "genus_richness:Temperature", "genus_richness:Salinity",
                      "FEve:Temperature", "FEve:Salinity", 
                      "FDis:Temperature", "FDis:Salinity",
                      "cyano_prop:Temperature", "cyano_prop:Salinity",
                      "species_richness:Temperature", "species_richness:Salinity",
                      "zoo_FEve:Temperature", "zoo_FEve:Salinity", 
                      "zoo_FDis:Temperature", "zoo_FDis:Salinity",
                      "daylight_hours"),
       ncol_ce = 6,
       ce_title = "Marginal Effects of Abiotic AND Biotic Variables for Phytoplankton RUE per P (DOY-random)"),
  
  # ---------------------------------------------------------
  # PHYTOPLANKTON RUE WITH NITROGEN MODELS WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260326/RDATA/BRM_phytoRUEN_CD_PI_260326_DOYr.RDATA",
       name = "phytoRUEN_CD_PI_DOYr",
       resp = "log_RUE_B_N",
       predictors = c("Temperature", "Salinity", "pH", "MaxDepth", "genus_richness",
                      "FEve", "FDis", "cyano_prop", "species_richness", "zoo_FEve", "zoo_FDis",
                      "genus_richness:Temperature", "genus_richness:Salinity",
                      "FEve:Temperature", "FEve:Salinity", 
                      "FDis:Temperature", "FDis:Salinity",
                      "cyano_prop:Temperature", "cyano_prop:Salinity",
                      "species_richness:Temperature", "species_richness:Salinity",
                      "zoo_FEve:Temperature", "zoo_FEve:Salinity", 
                      "zoo_FDis:Temperature", "zoo_FDis:Salinity",
                      "daylight_hours"),
       ncol_ce = 6,
       ce_title = "Marginal Effects of Abiotic AND Biotic Variables for Phytoplankton RUE per N (DOY-random)"),
  
  # ---------------------------------------------------------
  # ZOOPLANKTON SPECIES RICHNESS MODELS WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260326/RDATA/BRM_zoo_SR_CD_P_260326_DOYr.RDATA",
       name = "Zoo_SR_CD_P_DOYr",
       resp = "species_richness",
       predictors = c("Temperature", "Salinity", "pH", "Zooplankton_MaxDepth", 
                      "genus_richness", "FEve_selected", "FDis_selected", "cyano_prop",
                      "daylight_hours"),
       ncol_ce = 4,
       ce_title = "Marginal Effects of Abiotic and Biotic Variables for Zooplankton Species Richness (DOY-random)"),
  
  # ---------------------------------------------------------
  # ZOOPLANKTON FUNCTIONAL EVENNESS MODELS WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260326/RDATA/BRM_zoo_FE_CD_P_260326_DOYr.RDATA",
       name = "Zoo_FEve_CD_P_DOYr",
       resp = "zoo_FEve",
       predictors = c("Temperature", "Salinity", "pH", "Zooplankton_MaxDepth", 
                      "genus_richness", "FEve_selected", "FDis_selected", "cyano_prop",
                      "daylight_hours"),
       ncol_ce = 4,
       ce_title = "Marginal Effects of Abiotic and Biotic Variables for Zooplankton Functional Evenness (DOY-random)"),
  
  # ---------------------------------------------------------
  # ZOOPLANKTON FUNCTIONAL DISPERSION MODELS WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260326/RDATA/BRM_zoo_FD_CD_P_260326_DOYr.RDATA",
       name = "Zoo_FDis_CD_P_DOYr",
       resp = "zoo_FDis",
       predictors = c("Temperature", "Salinity", "pH", "Zooplankton_MaxDepth", 
                      "genus_richness", "FEve_selected", "FDis_selected", "cyano_prop",
                      "daylight_hours"),
       ncol_ce = 4,
       ce_title = "Marginal Effects of Abiotic and Biotic Variables for Zooplankton Functional Dispersion (DOY-random)"),
  
  # ---------------------------------------------------------
  # ZOOPLANKTON RUE WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260326/RDATA/BRM_zoo_RUE_CD_PI_260326_DOYr.RData",
       name = "zoo_RUE_CD_PI_DOYr",
       resp = "log_RUE_zooplankton",
       predictors = c("Temperature", "Salinity", "pH", "Zooplankton_MaxDepth", "genus_richness",
                      "FEve_selected", "FDis_selected",  "cyano_prop",
                      "species_richness", "zoo_FEve", "zoo_FDis",
                      "genus_richness:Temperature", "genus_richness:Salinity",
                      "FEve_selected:Temperature", "FEve_selected:Salinity", 
                      "FDis_selected:Temperature", "FDis_selected:Salinity",
                      "cyano_prop:Temperature", "cyano_prop:Salinity",
                      "species_richness:Temperature", "species_richness:Salinity",
                      "zoo_FEve:Temperature", "zoo_FEve:Salinity", 
                      "zoo_FDis:Temperature", "zoo_FDis:Salinity",
                      "daylight_hours"),
       ncol_ce = 6,
       ce_title = "Marginal Effects of Abiotic and Biotic Variables for Zooplankton RUE (DOY-random)")
  
  
)
  


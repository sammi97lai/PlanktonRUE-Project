# ===========================================================
# sDOY Model list 
# ===========================================================
model_list <- list(
  #PHYTOPLANKTON DIVERSITY
  # ---------------------------------------------------------
  # PHYTOPLANKTON GENUS RICHNESS MODELS
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260326/RDATA/BRM_phyto_GD_P_260326_sDOY.RData",
       name = "Phyto_Genus_Richness_P_sDOY",
       resp = "genus_richness",
       predictors = c("Temperature", "Salinity", "pH", "MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic Variables for Phytoplankton Genus Richness (spline DOY)"),
  
  # ---------------------------------------------------------
  # PHYTOPLANKTON FUNCTIONAL EVENNESS MODELS
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260326/RDATA/BRM_phyto_FEve_P_260326_sDOY.RDATA",
       name = "Phyto_FEve_P_sDOY",
       resp = "FEve",
       predictors = c("Temperature", "Salinity", "pH", "MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic Variables for Phytoplankton Functional Evenness (spline DOY)"),
  
  # ---------------------------------------------------------
  # PHYTOPLANKTON FUNCTIONAL DISPERSION MODELS
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260326/RDATA/BRM_phyto_FDis_P_260326_sDOY.RDATA",
       name = "Phyto_FDis_P_sDOY",
       resp = "FDis",
       predictors = c("Temperature", "Salinity", "pH", "MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic Variables for Phytoplankton Functional Dispersion (spline DOY)"),
  
  # ---------------------------------------------------------
  # CYANOBACTERIA DOMINANCE MODELS
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260326/RDATA/BRM_cyanodom_P_260326_sDOY_z.RDATA",
       name = "cyanodom_P_sDOY",
       resp = "cyano_prop",
       predictors = c("Temperature", "Salinity", "pH", "MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic Variables for Cyanobacteria Dominance (spline DOY)"),  
  
  # ---------------------------------------------------------
  # PHYTOPLANKTON RUE WITH PHOSPHORUS MODELS WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260326/RDATA/BRM_phytoRUEP_CD_PI_260326_sDOY.RDATA",
       name = "phytoRUEP_CD_PI_sDOY",
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
                      "DOY", "daylight_hours"),
       ncol_ce = 6,
       ce_title = "Marginal Effects of Abiotic AND Biotic Variables for Phytoplankton RUE per P (spline DOY)"),
  
  # ---------------------------------------------------------
  # PHYTOPLANKTON RUE WITH NITROGEN MODELS WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260326/RDATA/BRM_phytoRUEN_CD_PI_260326_sDOY.RDATA",
       name = "phytoRUEN_CD_PI_sDOY",
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
                      "DOY", "daylight_hours"),
       ncol_ce = 6,
       ce_title = "Marginal Effects of Abiotic AND Biotic Variables for Phytoplankton RUE per N (spline DOY)"),
  
  # ---------------------------------------------------------
  # ZOOPLANKTON SPECIES RICHNESS MODELS WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260326/RDATA/BRM_zoo_SR_CD_P_260326_sDOY.RDATA",
       name = "Zoo_SR_CD_P_sDOY",
       resp = "species_richness",
       predictors = c("Temperature", "Salinity", "pH", "Zooplankton_MaxDepth", 
                      "genus_richness", "FEve_selected", "FDis_selected", "cyano_prop",
                      "DOY", "daylight_hours"),
       ncol_ce = 4,
       ce_title = "Marginal Effects of Abiotic and Biotic Variables for Zooplankton Species Richness (spline DOY)"),
  
  # ---------------------------------------------------------
  # ZOOPLANKTON FUNCTIONAL EVENNESS MODELS WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260326/RDATA/BRM_zoo_FE_CD_P_260326_sDOY.RDATA",
       name = "Zoo_FEve_CD_P_sDOY",
       resp = "zoo_FEve",
       predictors = c("Temperature", "Salinity", "pH", "Zooplankton_MaxDepth", 
                      "genus_richness", "FEve_selected", "FDis_selected", "cyano_prop",
                      "DOY", "daylight_hours"),
       ncol_ce = 4,
       ce_title = "Marginal Effects of Abiotic and Biotic Variables for Zooplankton Functional Evenness (spline DOY)"),
  
  # ---------------------------------------------------------
  # ZOOPLANKTON FUNCTIONAL DISPERSION MODELS WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260326/RDATA/BRM_zoo_FD_CD_P_260326_sDOY.RDATA",
       name = "Zoo_FDis_CD_P_sDOY",
       resp = "zoo_FDis",
       predictors = c("Temperature", "Salinity", "pH", "Zooplankton_MaxDepth", 
                      "genus_richness", "FEve_selected", "FDis_selected", "cyano_prop",
                      "DOY", "daylight_hours"),
       ncol_ce = 4,
       ce_title = "Marginal Effects of Abiotic and Biotic Variables for Zooplankton Functional Dispersion (spline DOY)"),
  
  # ---------------------------------------------------------
  # ZOOPLANKTON RUE WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260326/RDATA/BRM_zoo_RUE_CD_PI_260326_sDOY.RData",
       name = "zoo_RUE_CD_PI_sDOY",
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
                      "DOY", "daylight_hours"),
       ncol_ce = 6,
       ce_title = "Marginal Effects of Abiotic and Biotic Variables for Zooplankton RUE (spline DOY)")
  
)
  


# ===========================================================
# N:P ratio Model list 
# ===========================================================
model_list <- list(
  #PHYTOPLANKTON DIVERSITY
  # ---------------------------------------------------------
  # PHYTOPLANKTON GENUS RICHNESS MODELS
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260326/RDATA/BRM_phyto_GD_P_260326_NP.RData",
       name = "Phyto_Genus_Richness_P_NP",
       resp = "genus_richness",
       predictors = c("Temperature", "Salinity", "pH", "log_N_P_ratio", "MaxDepth", "DOY"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic Variables for Phytoplankton Genus Richness (with N:P ratio)"),

  # ---------------------------------------------------------
  # PHYTOPLANKTON FUNCTIONAL EVENNESS MODELS
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260326/RDATA/BRM_phyto_FEve_P_260326_NP.RDATA",
       name = "Phyto_FEve_P_NP",
       resp = "FEve",
       predictors = c("Temperature", "Salinity", "pH", "log_N_P_ratio", "MaxDepth", "DOY"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic Variables for Phytoplankton Functional Evenness (with N:P ratio)"),

  # ---------------------------------------------------------
  # PHYTOPLANKTON FUNCTIONAL DISPERSION MODELS
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260326/RDATA/BRM_phyto_FDis_P_260326_NP.RDATA",
       name = "Phyto_FDis_P_NP",
       resp = "FDis",
       predictors = c("Temperature", "Salinity", "pH", "log_N_P_ratio", "MaxDepth", "DOY"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic Variables for Phytoplankton Functional Dispersion (with N:P ratio)"),
  
  # ---------------------------------------------------------
  # CYANOBACTERIA DOMINANCE MODELS
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260326/RDATA/BRM_cyanodom_P_260326_NP_z.RDATA",
       name = "cyanodom_P_NP",
       resp = "cyano_prop",
       predictors = c("Temperature", "Salinity", "pH", "log_N_P_ratio", "MaxDepth", "DOY"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic Variables for Cyanobacteria Dominance (with N:P ratio)"),  
  
  # ---------------------------------------------------------
  # PHYTOPLANKTON RUE WITH PHOSPHORUS MODELS WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260326/RDATA/BRM_phytoRUEP_CD_PI_260326_NP.RDATA",
       name = "phytoRUEP_CD_PI_NP",
       resp = "log_RUE_B_P",
       predictors = c("Temperature", "Salinity", "pH", "log_N_P_ratio", "MaxDepth", "genus_richness",
                      "FEve", "FDis", "cyano_prop", "species_richness", "zoo_FEve", "zoo_FDis",
                      "genus_richness:Temperature", "genus_richness:Salinity",
                      "FEve:Temperature", "FEve:Salinity",
                      "FDis:Temperature", "FDis:Salinity",
                      "cyano_prop:Temperature", "cyano_prop:Salinity",
                      "species_richness:Temperature", "species_richness:Salinity",
                      "zoo_FEve:Temperature", "zoo_FEve:Salinity",
                      "zoo_FDis:Temperature", "zoo_FDis:Salinity",
                      "DOY"),
       ncol_ce = 6,
       ce_title = "Marginal Effects of Abiotic AND Biotic Variables for Phytoplankton RUE per P (with N:P ratio)"),

  # ---------------------------------------------------------
  # PHYTOPLANKTON RUE WITH NITROGEN MODELS WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260326/RDATA/BRM_phytoRUEN_CD_PI_260326_NP.RDATA",
       name = "phytoRUEN_CD_PI_NP",
       resp = "log_RUE_B_N",
       predictors = c("Temperature", "Salinity", "pH", "log_N_P_ratio", "MaxDepth", "genus_richness",
                      "FEve", "FDis", "cyano_prop", "species_richness", "zoo_FEve", "zoo_FDis",
                      "genus_richness:Temperature", "genus_richness:Salinity",
                      "FEve:Temperature", "FEve:Salinity",
                      "FDis:Temperature", "FDis:Salinity",
                      "cyano_prop:Temperature", "cyano_prop:Salinity",
                      "species_richness:Temperature", "species_richness:Salinity",
                      "zoo_FEve:Temperature", "zoo_FEve:Salinity",
                      "zoo_FDis:Temperature", "zoo_FDis:Salinity",
                      "DOY"),
       ncol_ce = 6,
       ce_title = "Marginal Effects of Abiotic AND Biotic Variables for Phytoplankton RUE per N (with N:P ratio)"),

  # ---------------------------------------------------------
  # ZOOPLANKTON SPECIES RICHNESS MODELS WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260326/RDATA/BRM_zoo_SR_CD_P_260326_NP.RDATA",
       name = "Zoo_SR_CD_P_NP",
       resp = "species_richness",
       predictors = c("Temperature", "Salinity", "pH", "log_N_P_ratio", "Zooplankton_MaxDepth",
                      "genus_richness", "FEve_selected", "FDis_selected", "cyano_prop",
                      "DOY"),
       ncol_ce = 4,
       ce_title = "Marginal Effects of Abiotic and Biotic Variables for Zooplankton Species Richness (with N:P ratio)"),

  # ---------------------------------------------------------
  # ZOOPLANKTON FUNCTIONAL EVENNESS MODELS WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260326/RDATA/BRM_zoo_FE_CD_P_260326_NP.RDATA",
       name = "Zoo_FEve_CD_P_NP",
       resp = "zoo_FEve",
       predictors = c("Temperature", "Salinity", "pH", "log_N_P_ratio", "Zooplankton_MaxDepth",
                      "genus_richness", "FEve_selected", "FDis_selected", "cyano_prop",
                      "DOY"),
       ncol_ce = 4,
       ce_title = "Marginal Effects of Abiotic and Biotic Variables for Zooplankton Functional Evenness (with N:P ratio)"),

  # ---------------------------------------------------------
  # ZOOPLANKTON FUNCTIONAL DISPERSION MODELS WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260326/RDATA/BRM_zoo_FD_CD_P_260326_NP.RDATA",
       name = "Zoo_FDis_CD_P_NP",
       resp = "zoo_FDis",
       predictors = c("Temperature", "Salinity", "pH", "log_N_P_ratio", "Zooplankton_MaxDepth",
                      "genus_richness", "FEve_selected", "FDis_selected", "cyano_prop",
                      "DOY"),
       ncol_ce = 4,
       ce_title = "Marginal Effects of Abiotic and Biotic Variables for Zooplankton Functional Dispersion (with N:P ratio)"),

  # ---------------------------------------------------------
  # ZOOPLANKTON RUE WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260326/RDATA/BRM_zoo_RUE_CD_PI_260326_NP.RData",
       name = "zoo_RUE_CD_PI_NP",
       resp = "log_RUE_zooplankton",
       predictors = c("Temperature", "Salinity", "pH", "log_N_P_ratio", "Zooplankton_MaxDepth", "genus_richness",
                      "FEve_selected", "FDis_selected",  "cyano_prop",
                      "species_richness", "zoo_FEve", "zoo_FDis",
                      "genus_richness:Temperature", "genus_richness:Salinity",
                      "FEve_selected:Temperature", "FEve_selected:Salinity",
                      "FDis_selected:Temperature", "FDis_selected:Salinity",
                      "cyano_prop:Temperature", "cyano_prop:Salinity",
                      "species_richness:Temperature", "species_richness:Salinity",
                      "zoo_FEve:Temperature", "zoo_FEve:Salinity",
                      "zoo_FDis:Temperature", "zoo_FDis:Salinity",
                      "DOY"),
       ncol_ce = 6,
       ce_title = "Marginal Effects of Abiotic and Biotic Variables for Zooplankton RUE (with N:P ratio)")
  
)
  


# ===========================================================
# Zooplankton Model list 
# ===========================================================
model_list <- list(
  #ZOOPLANKTON DIVERSITY#
  # ---------------------------------------------------------
  # ZOOPLANKTON SPECIES RICHNESS MODELS
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM/RDATA/BRM_zoo_SR_L_260313.RDATA",
       name = "Zoo_SR_L",
       resp = "species_richness",
       predictors = c("Temperature", "Salinity", "pH", 
                      "genus_richness", "FEve_selected", "FDis_selected", 
                      "Zooplankton_MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic and Biotic Variables for Zooplankton Species Richness (Linear)"),

  list(file = "~/GitHub/PlanktonRUE-Project/BRM/RDATA/BRM_zoo_SR_P_260313.RDATA",
       name = "Zoo_SR_P",
       resp = "species_richness",
       predictors = c("Temperature", "Salinity", "pH", 
                      "genus_richness", "FEve_selected", "FDis_selected", 
                      "Zooplankton_MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic and Biotic Variables for Zooplankton Species Richness (Quadratic)"),

  # ---------------------------------------------------------
  # ZOOPLANKTON SPECIES RICHNESS MODELS - POISSON
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM/RDATA/BRM_zoo_SR_L_p_260313.RDATA",
       name = "Zoo_SR_L_p",
       resp = "species_richness",
       predictors = c("Temperature", "Salinity", "pH", 
                      "genus_richness", "FEve_selected", "FDis_selected", 
                      "Zooplankton_MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic and Biotic Variables for Zooplankton Species Richness (Linear) - poisson"),

  list(file = "~/GitHub/PlanktonRUE-Project/BRM/RDATA/BRM_zoo_SR_P_p_260313.RDATA",
       name = "Zoo_SR_P_p",
       resp = "species_richness",
       predictors = c("Temperature", "Salinity", "pH", 
                      "genus_richness", "FEve_selected", "FDis_selected", 
                      "Zooplankton_MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic and Biotic Variables for Zooplankton Species Richness (Quadratic)- poisson"),

  # ---------------------------------------------------------
  # ZOOPLANKTON FUNCTIONAL EVENNESS MODELS
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM/RDATA/BRM_zoo_FE_L_260313.RDATA",
       name = "Zoo_FEve_L",
       resp = "zoo_FEve",
       predictors = c("Temperature", "Salinity", "pH", 
                      "genus_richness", "FEve_selected", "FDis_selected", 
                      "Zooplankton_MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic and Biotic Variables for Zooplankton Functional Evenness (Linear)"),

  list(file = "~/GitHub/PlanktonRUE-Project/BRM/RDATA/BRM_zoo_FE_P_260313.RDATA",
       name = "Zoo_FEve_P",
       resp = "zoo_FEve",
       predictors = c("Temperature", "Salinity", "pH", 
                      "genus_richness", "FEve_selected", "FDis_selected", 
                      "Zooplankton_MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic and Biotic Variables for Zooplankton Functional Evenness (Quadratic)"),

  # ---------------------------------------------------------
  # ZOOPLANKTON FUNCTIONAL EVENNESS MODELS - BETA
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM/RDATA/BRM_zoo_FE_L_b_260313.RDATA",
       name = "Zoo_FEve_L_b",
       resp = "zoo_FEve",
       predictors = c("Temperature", "Salinity", "pH", 
                      "genus_richness", "FEve_selected", "FDis_selected", 
                      "Zooplankton_MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic and Biotic Variables for Zooplankton Functional Evenness (Linear) - beta"),

  list(file = "~/GitHub/PlanktonRUE-Project/BRM/RDATA/BRM_zoo_FE_P_b_260313.RDATA",
       name = "Zoo_FEve_P_b",
       resp = "zoo_FEve",
       predictors = c("Temperature", "Salinity", "pH", 
                      "genus_richness", "FEve_selected", "FDis_selected", 
                      "Zooplankton_MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic and Biotic Variables for Zooplankton Functional Evenness (Quadratic) - beta"),

  # ---------------------------------------------------------
  # ZOOPLANKTON FUNCTIONAL DISPERSION MODELS
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM/RDATA/BRM_zoo_FD_L_260313.RDATA",
       name = "Zoo_FDis_L",
       resp = "zoo_FDis",
       predictors = c("Temperature", "Salinity", "pH", 
                      "genus_richness", "FEve_selected", "FDis_selected", 
                      "Zooplankton_MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic and Biotic Variables for Zooplankton Functional Dispersion (Linear)"),

  list(file = "~/GitHub/PlanktonRUE-Project/BRM/RDATA/BRM_zoo_FD_P_260313.RDATA",
       name = "Zoo_FDis_P",
       resp = "zoo_FDis",
       predictors = c("Temperature", "Salinity", "pH", 
                      "genus_richness", "FEve_selected", "FDis_selected", 
                      "Zooplankton_MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic and Biotic Variables for Zooplankton Functional Dispersion (Quadratic)"),

  # ---------------------------------------------------------
  # ZOOPLANKTON FUNCTIONAL DISPERSION MODELS - BETA
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM/RDATA/BRM_zoo_FD_L_b_260313.RDATA",
       name = "Zoo_FDis_L_b",
       resp = "zoo_FDis",
       predictors = c("Temperature", "Salinity", "pH", 
                      "genus_richness", "FEve_selected", "FDis_selected", 
                      "Zooplankton_MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic and Biotic Variables for Zooplankton Functional Dispersion (Linear) - beta"),

  list(file = "~/GitHub/PlanktonRUE-Project/BRM/RDATA/BRM_zoo_FD_P_b_260313.RDATA",
       name = "Zoo_FDis_P_b",
       resp = "zoo_FDis",
       predictors = c("Temperature", "Salinity", "pH", 
                      "genus_richness", "FEve_selected", "FDis_selected", 
                      "Zooplankton_MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic and Biotic Variables for Zooplankton Functional Dispersion (Quadratic) -beta"),

  #ZOOPLANKTON DIVERSITY WITH CYANO DOM MODELS#
  # ---------------------------------------------------------
  # ZOOPLANKTON SPECIES RICHNESS MODELS WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM/RDATA/BRM_zoo_SR_CD_L_260313.RDATA",
       name = "Zoo_SR_CD_L",
       resp = "species_richness",
       predictors = c("Temperature", "Salinity", "pH", 
                      "genus_richness", "FEve_selected", "FDis_selected", "cyano_prop",
                      "Zooplankton_MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 4,
       ce_title = "Marginal Effects of Abiotic and Biotic (incl. CD) Variables for Zooplankton Species Richness (Linear)"),

  list(file = "~/GitHub/PlanktonRUE-Project/BRM/RDATA/BRM_zoo_SR_CD_P_260313.RDATA",
       name = "Zoo_SR_CD_P",
       resp = "species_richness",
       predictors = c("Temperature", "Salinity", "pH", 
                      "genus_richness", "FEve_selected", "FDis_selected", "cyano_prop",
                      "Zooplankton_MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 4,
       ce_title = "Marginal Effects of Abiotic and Biotic (incl. CD) Variables for Zooplankton Species Richness (Quadratic)"),

  # ---------------------------------------------------------
  # ZOOPLANKTON SPECIES RICHNESS MODELS WITH CYANO DOM - POISSON
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM/RDATA/BRM_zoo_SR_CD_L_p_260313.RData",
       name = "Zoo_SR_CD_L_p",
       resp = "species_richness",
       predictors = c("Temperature", "Salinity", "pH", 
                      "genus_richness", "FEve_selected", "FDis_selected", "cyano_prop",
                      "Zooplankton_MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 4,
       ce_title = "Marginal Effects of Abiotic and Biotic (incl. CD) Variables for Zooplankton Species Richness (Linear) - poisson"),

  list(file = "~/GitHub/PlanktonRUE-Project/BRM/RDATA/BRM_zoo_SR_CD_P_p_260313.RDATA",
       name = "Zoo_SR_CD_P_p",
       resp = "species_richness",
       predictors = c("Temperature", "Salinity", "pH", 
                      "genus_richness", "FEve_selected", "FDis_selected", "cyano_prop",
                      "Zooplankton_MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 4,
       ce_title = "Marginal Effects of Abiotic and Biotic (incl. CD) Variables for Zooplankton Species Richness (Quadratic) - poisson"),

  # ---------------------------------------------------------
  # ZOOPLANKTON FUNCTIONAL EVENNESS MODELS WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM/RDATA/BRM_zoo_FE_CD_L_260313.RDATA",
       name = "Zoo_FEve_CD_L",
       resp = "zoo_FEve",
       predictors = c("Temperature", "Salinity", "pH", 
                      "genus_richness", "FEve_selected", "FDis_selected", "cyano_prop",
                      "Zooplankton_MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 4,
       ce_title = "Marginal Effects of Abiotic and Biotic (incl. CD) Variables for Zooplankton Functional Evenness (Linear)"),

  list(file = "~/GitHub/PlanktonRUE-Project/BRM/RDATA/BRM_zoo_FE_CD_P_260313.RDATA",
       name = "Zoo_FEve_CD_P",
       resp = "zoo_FEve",
       predictors = c("Temperature", "Salinity", "pH", 
                      "genus_richness", "FEve_selected", "FDis_selected", "cyano_prop",
                      "Zooplankton_MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 4,
       ce_title = "Marginal Effects of Abiotic and Biotic (incl. CD) Variables for Zooplankton Functional Evenness (Quadratic)"),

  # ---------------------------------------------------------
  # ZOOPLANKTON FUNCTIONAL EVENNESS MODELS WITH CYANO DOM - BETA
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM/RDATA/BRM_zoo_FE_CD_L_b_260313.RDATA",
       name = "Zoo_FEve_CD_L_b",
       resp = "zoo_FEve",
       predictors = c("Temperature", "Salinity", "pH", 
                      "genus_richness", "FEve_selected", "FDis_selected", "cyano_prop",
                      "Zooplankton_MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 4,
       ce_title = "Marginal Effects of Abiotic and Biotic (incl. CD) Variables for Zooplankton Functional Evenness (Linear) - beta"),

  list(file = "~/GitHub/PlanktonRUE-Project/BRM/RDATA/BRM_zoo_FE_CD_P_b_260313.RDATA",
       name = "Zoo_FEve_CD_P_b",
       resp = "zoo_FEve",
       predictors = c("Temperature", "Salinity", "pH", 
                      "genus_richness", "FEve_selected", "FDis_selected", "cyano_prop",
                      "Zooplankton_MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 4,
       ce_title = "Marginal Effects of Abiotic and Biotic (incl. CD) Variables for Zooplankton Functional Evenness (Quadratic) - beta"),


  # ---------------------------------------------------------
  # ZOOPLANKTON FUNCTIONAL DISPERSION MODELS WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM/RDATA/BRM_zoo_FD_CD_L_260313.RDATA",
       name = "Zoo_FDis_CD_L",
       resp = "zoo_FDis",
       predictors = c("Temperature", "Salinity", "pH", 
                      "genus_richness", "FEve_selected", "FDis_selected", "cyano_prop",
                      "Zooplankton_MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 4,
       ce_title = "Marginal Effects of Abiotic and Biotic (incl. CD) Variables for Zooplankton Functional Dispersion (Linear)"),

  list(file = "~/GitHub/PlanktonRUE-Project/BRM/RDATA/BRM_zoo_FD_CD_P_260313.RDATA",
       name = "Zoo_FDis_CD_P",
       resp = "zoo_FDis",
       predictors = c("Temperature", "Salinity", "pH", 
                      "genus_richness", "FEve_selected", "FDis_selected", "cyano_prop",
                      "Zooplankton_MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 4,
       ce_title = "Marginal Effects of Abiotic and Biotic (incl. CD) Variables for Zooplankton Functional Dispersion (Quadratic)"),

  # ---------------------------------------------------------
  # ZOOPLANKTON FUNCTIONAL DISPERSION MODELS WITH CYANO DOM - BETA
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM/RDATA/BRM_zoo_FD_CD_L_b_260313.RDATA",
       name = "Zoo_FDis_CD_L_b",
       resp = "zoo_FDis",
       predictors = c("Temperature", "Salinity", "pH", 
                      "genus_richness", "FEve_selected", "FDis_selected", "cyano_prop",
                      "Zooplankton_MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 4,
       ce_title = "Marginal Effects of Abiotic and Biotic (incl. CD) Variables for Zooplankton Functional Dispersion (Linear) - beta"),

  list(file = "~/GitHub/PlanktonRUE-Project/BRM/RDATA/BRM_zoo_FD_CD_P_b_260313.RDATA",
       name = "Zoo_FDis_CD_P_b",
       resp = "zoo_FDis",
       predictors = c("Temperature", "Salinity", "pH", 
                      "genus_richness", "FEve_selected", "FDis_selected", "cyano_prop",
                      "Zooplankton_MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 4,
       ce_title = "Marginal Effects of Abiotic and Biotic (incl. CD) Variables for Zooplankton Functional Dispersion (Quadratic) - beta"),


  #ZOOPLANKTON RUE MODELS#
  # ---------------------------------------------------------
  # ZOOPLANKTON RUE
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM/RDATA/BRM_zoo_RUE_L_260313.RDATA",
       name = "zoo_RUE_L",
       resp = "log_RUE_zooplankton",
       predictors = c("Temperature", "Salinity", "pH", 
                      "genus_richness", "FEve_selected", "FDis_selected", 
                      "species_richness", "zoo_FEve", "zoo_FDis",
                      "Zooplankton_MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 4,
       ce_title = "Marginal Effects of Abiotic and Biotic Variables for Zooplankton RUE (Linear)"),

  list(file = "~/GitHub/PlanktonRUE-Project/BRM/RDATA/BRM_zoo_RUE_P_260313.RDATA",
       name = "zoo_RUE_P",
       resp = "log_RUE_zooplankton",
       predictors = c("Temperature", "Salinity", "pH", 
                      "genus_richness", "FEve_selected", "FDis_selected", 
                      "species_richness", "zoo_FEve", "zoo_FDis",
                      "Zooplankton_MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 4,
       ce_title = "Marginal Effects of Abiotic and Biotic Variables for Zooplankton RUE (Quadratic)"),

  list(file = "~/GitHub/PlanktonRUE-Project/BRM/RDATA/BRM_zoo_RUE_LI_260313.RDATA",
       name = "zoo_RUE_LI",
       resp = "log_RUE_zooplankton",
       predictors = c("Temperature", "Salinity", "pH", "genus_richness",
                      "FEve_selected", "FDis_selected", "species_richness", "zoo_FEve", "zoo_FDis",
                      "genus_richness:Temperature", "genus_richness:Salinity",
                      "FEve_selected:Temperature", "FEve_selected:Salinity", 
                      "FDis_selected:Temperature", "FDis_selected:Salinity",
                      "species_richness:Temperature", "species_richness:Salinity",
                      "zoo_FEve:Temperature", "zoo_FEve:Salinity", 
                      "zoo_FDis:Temperature", "zoo_FDis:Salinity",
                      "Zooplankton_MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 5,
       ce_title = "Marginal Effects of Abiotic and Biotic Variables with Interactions for Zooplankton RUE (Linear)"),
  
  list(file = "~/GitHub/PlanktonRUE-Project/BRM/RDATA/BRM_zoo_RUE_PI_260313.RDATA",
       name = "zoo_RUE_PI",
       resp = "log_RUE_zooplankton",
       predictors = c("Temperature", "Salinity", "pH", "genus_richness",
                      "FEve_selected", "FDis_selected", "species_richness", "zoo_FEve", "zoo_FDis",
                      "genus_richness:Temperature", "genus_richness:Salinity",
                      "FEve_selected:Temperature", "FEve_selected:Salinity", 
                      "FDis_selected:Temperature", "FDis_selected:Salinity",
                      "species_richness:Temperature", "species_richness:Salinity",
                      "zoo_FEve:Temperature", "zoo_FEve:Salinity", 
                      "zoo_FDis:Temperature", "zoo_FDis:Salinity",
                      "Zooplankton_MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 5,
       ce_title = "Marginal Effects of Abiotic and Biotic Variables with Interactions for Zooplankton RUE (Quadratic)"),
  
  # ---------------------------------------------------------
  # ZOOPLANKTON RUE WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM/RDATA/BRM_zoo_RUE_CD_L_260313.RDATA",
       name = "zoo_RUE_CD_L",
       resp = "log_RUE_zooplankton",
       predictors = c("Temperature", "Salinity", "pH", 
                      "genus_richness", "FEve_selected", "FDis_selected", "cyano_prop",
                      "species_richness", "zoo_FEve", "zoo_FDis",
                      "Zooplankton_MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 4,
       ce_title = "Marginal Effects of Abiotic and Biotic (incl. CD) Variables for Zooplankton RUE (Linear)"),

  list(file = "~/GitHub/PlanktonRUE-Project/BRM/RDATA/BRM_zoo_RUE_CD_P_260313.RDATA",
       name = "zoo_RUE_CD_P",
       resp = "log_RUE_zooplankton",
       predictors = c("Temperature", "Salinity", "pH", 
                      "genus_richness", "FEve_selected", "FDis_selected", "cyano_prop",
                      "species_richness", "zoo_FEve", "zoo_FDis",
                      "Zooplankton_MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 4,
       ce_title = "Marginal Effects of Abiotic and Biotic (incl. CD) Variables for Zooplankton RUE (Quadratic)"),
  
  list(file = "~/GitHub/PlanktonRUE-Project/BRM/RDATA/BRM_zoo_RUE_CD_LI_260313.RDATA",
       name = "zoo_RUE_CD_LI",
       resp = "log_RUE_zooplankton",
       predictors = c("Temperature", "Salinity", "pH", "genus_richness",
                      "FEve_selected", "FDis_selected", "cyano_prop",
                      "species_richness", "zoo_FEve", "zoo_FDis",
                      "genus_richness:Temperature", "genus_richness:Salinity",
                      "FEve_selected:Temperature", "FEve_selected:Salinity", 
                      "FDis_selected:Temperature", "FDis_selected:Salinity",
                      "cyano_prop:Temperature", "cyano_prop:Salinity",
                      "species_richness:Temperature", "species_richness:Salinity",
                      "zoo_FEve:Temperature", "zoo_FEve:Salinity", 
                      "zoo_FDis:Temperature", "zoo_FDis:Salinity",
                      "Zooplankton_MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 5,
       ce_title = "Marginal Effects of Abiotic and Biotic (incl. CD) Variables with Interactions for Zooplankton RUE (Linear)"),
  
  list(file = "~/GitHub/PlanktonRUE-Project/BRM/RDATA/BRM_zoo_RUE_CD_PI_260313.RDATA",
       name = "zoo_RUE_CD_PI",
       resp = "log_RUE_zooplankton",
       predictors = c("Temperature", "Salinity", "pH", "genus_richness",
                      "FEve_selected", "FDis_selected",  "cyano_prop",
                      "species_richness", "zoo_FEve", "zoo_FDis",
                      "genus_richness:Temperature", "genus_richness:Salinity",
                      "FEve_selected:Temperature", "FEve_selected:Salinity", 
                      "FDis_selected:Temperature", "FDis_selected:Salinity",
                      "cyano_prop:Temperature", "cyano_prop:Salinity",
                      "species_richness:Temperature", "species_richness:Salinity",
                      "zoo_FEve:Temperature", "zoo_FEve:Salinity", 
                      "zoo_FDis:Temperature", "zoo_FDis:Salinity",
                      "Zooplankton_MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 5,
       ce_title = "Marginal Effects of Abiotic and Biotic (incl. CD) Variables with Interactions for Zooplankton RUE (Quadratic)")
)




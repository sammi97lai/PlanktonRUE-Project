# ===========================================================
# Salinity and Nitrogen Model list 
# ===========================================================
model_list <- list(
  #PHYTOPLANKTON DIVERSITY
  # ---------------------------------------------------------
  # PHYTOPLANKTON GENUS RICHNESS MODELS
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260330/RDATA/BRM_phyto_GD_P_260330_S_Nit.RDATA",
       name = "Phyto_Genus_Richness_P_S_Nit",
       resp = "genus_richness",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered", "MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic Variables for Phytoplankton Genus Richness (Salinity and Nitrogen)"),
  
  # ---------------------------------------------------------
  # PHYTOPLANKTON FUNCTIONAL EVENNESS MODELS
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260330/RDATA/BRM_phyto_FEve_P_260330_S_Nit.RDATA",
       name = "Phyto_FEve_P_S_Nit",
       resp = "FEve",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered", "MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic Variables for Phytoplankton Functional Evenness (Salinity and Nitrogen)"),
  
  # ---------------------------------------------------------
  # PHYTOPLANKTON FUNCTIONAL DISPERSION MODELS
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260330/RDATA/BRM_phyto_FDis_P_260330_S_Nit.RDATA",
       name = "Phyto_FDis_P_S_Nit",
       resp = "FDis",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered", "MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic Variables for Phytoplankton Functional Dispersion (Salinity and Nitrogen)"),
  
  # ---------------------------------------------------------
  # CYANOBACTERIA DOMINANCE MODELS
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260330/RDATA/BRM_cyanodom_P_260330_S_Nit.RDATA",
       name = "cyanodom_P_S_Nit",
       resp = "cyano_prop",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered", "MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic Variables for Cyanobacteria Dominance (Salinity and Nitrogen)"),  
  
  # ---------------------------------------------------------
  # PHYTOPLANKTON RUE WITH PHOSPHORUS MODELS WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260330/RDATA/BRM_phytoRUEP_CD_PI_260330_S_Nit.RDATA",
       name = "phytoRUEP_CD_PI_S_Nit",
       resp = "log_RUE_B_P",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered", "MaxDepth", "genus_richness",
                      "FEve", "FDis", "cyano_prop", "species_richness", "zoo_FEve", "zoo_FDis",
                      "genus_richness:Temperature", "genus_richness:Salinity", "genus_richness:TotalNitrogen_unfiltered",
                      "FEve:Temperature", "FEve:Salinity", "FEve:TotalNitrogen_unfiltered", 
                      "FDis:Temperature", "FDis:Salinity", "FDis:TotalNitrogen_unfiltered",
                      "cyano_prop:Temperature", "cyano_prop:Salinity", "cyano_prop:TotalNitrogen_unfiltered",
                      "species_richness:Temperature", "species_richness:Salinity", "species_richness:TotalNitrogen_unfiltered",
                      "zoo_FEve:Temperature", "zoo_FEve:Salinity", "zoo_FEve:TotalNitrogen_unfiltered",
                      "zoo_FDis:Temperature", "zoo_FDis:Salinity", "zoo_FDis:TotalNitrogen_unfiltered",
                      "DOY", "daylight_hours"),
       ncol_ce = 6,
       ce_title = "Marginal Effects of Abiotic AND Biotic Variables for Phytoplankton RUE per P (Salinity and Nitrogen)"),
  
  # ---------------------------------------------------------
  # PHYTOPLANKTON RUE WITH NITROGEN MODELS WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260330/RDATA/BRM_phytoRUEN_CD_PI_260330_S_Nit.RDATA",
       name = "phytoRUEN_CD_PI_S_Nit",
       resp = "log_RUE_B_N",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered", "MaxDepth", "genus_richness",
                      "FEve", "FDis", "cyano_prop", "species_richness", "zoo_FEve", "zoo_FDis",
                      "genus_richness:Temperature", "genus_richness:Salinity", "genus_richness:TotalNitrogen_unfiltered",
                      "FEve:Temperature", "FEve:Salinity", "FEve:TotalNitrogen_unfiltered", 
                      "FDis:Temperature", "FDis:Salinity", "FDis:TotalNitrogen_unfiltered",
                      "cyano_prop:Temperature", "cyano_prop:Salinity", "cyano_prop:TotalNitrogen_unfiltered",
                      "species_richness:Temperature", "species_richness:Salinity", "species_richness:TotalNitrogen_unfiltered",
                      "zoo_FEve:Temperature", "zoo_FEve:Salinity", "zoo_FEve:TotalNitrogen_unfiltered",
                      "zoo_FDis:Temperature", "zoo_FDis:Salinity", "zoo_FDis:TotalNitrogen_unfiltered",
                      "DOY", "daylight_hours"),
       ncol_ce = 6,
       ce_title = "Marginal Effects of Abiotic AND Biotic Variables for Phytoplankton RUE per N (Salinity and Nitrogen)"),
  
  # ---------------------------------------------------------
  # ZOOPLANKTON SPECIES RICHNESS MODELS WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260330/RDATA/BRM_zoo_SR_CD_P_260330_S_Nit.RDATA",
       name = "Zoo_SR_CD_P_S_Nit",
       resp = "species_richness",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered", "Zooplankton_MaxDepth", 
                      "genus_richness", "FEve_selected", "FDis_selected", "cyano_prop",
                      "DOY", "daylight_hours"),
       ncol_ce = 4,
       ce_title = "Marginal Effects of Abiotic and Biotic Variables for Zooplankton Species Richness (Salinity and Nitrogen)"),
  
  # ---------------------------------------------------------
  # ZOOPLANKTON FUNCTIONAL EVENNESS MODELS WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260330/RDATA/BRM_zoo_FE_CD_P_260330_S_Nit.RDATA",
       name = "Zoo_FEve_CD_P_S_Nit",
       resp = "zoo_FEve",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered", "Zooplankton_MaxDepth", 
                      "genus_richness", "FEve_selected", "FDis_selected", "cyano_prop",
                      "DOY", "daylight_hours"),
       ncol_ce = 4,
       ce_title = "Marginal Effects of Abiotic and Biotic Variables for Zooplankton Functional Evenness (Salinity and Nitrogen)"),
  
  # ---------------------------------------------------------
  # ZOOPLANKTON FUNCTIONAL DISPERSION MODELS WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260330/RDATA/BRM_zoo_FD_CD_P_260330_S_Nit.RDATA",
       name = "Zoo_FDis_CD_P_S_Nit",
       resp = "zoo_FDis",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered", "Zooplankton_MaxDepth", 
                      "genus_richness", "FEve_selected", "FDis_selected", "cyano_prop",
                      "DOY", "daylight_hours"),
       ncol_ce = 4,
       ce_title = "Marginal Effects of Abiotic and Biotic Variables for Zooplankton Functional Dispersion (Salinity and Nitrogen)"),
  
  # ---------------------------------------------------------
  # ZOOPLANKTON RUE WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_model_test_260330/RDATA/BRM_zoo_RUE_CD_PI_260330_S_Nit.RDaTA",
       name = "zoo_RUE_CD_PI_S_Nit",
       resp = "log_RUE_zooplankton",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered", "Zooplankton_MaxDepth", "genus_richness",
                      "FEve_selected", "FDis_selected",  "cyano_prop",
                      "species_richness", "zoo_FEve", "zoo_FDis",
                      "genus_richness:Temperature", "genus_richness:Salinity", "genus_richness:TotalNitrogen_unfiltered",
                      "FEve_selected:Temperature", "FEve_selected:Salinity", "FEve_selected:TotalNitrogen_unfiltered", 
                      "FDis_selected:Temperature", "FDis_selected:Salinity", "FDis_selected:TotalNitrogen_unfiltered",
                      "cyano_prop:Temperature", "cyano_prop:Salinity", "cyano_prop:TotalNitrogen_unfiltered",
                      "species_richness:Temperature", "species_richness:Salinity", "species_richness:TotalNitrogen_unfiltered",
                      "zoo_FEve:Temperature", "zoo_FEve:Salinity", "zoo_FEve:TotalNitrogen_unfiltered", 
                      "zoo_FDis:Temperature", "zoo_FDis:Salinity", "zoo_FDis:TotalNitrogen_unfiltered",
                      "DOY", "daylight_hours"),
       ncol_ce = 6,
       ce_title = "Marginal Effects of Abiotic and Biotic Variables for Zooplankton RUE (Salinity and Nitrogen)")
  
)
  


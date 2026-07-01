# ===========================================================
# Phytoplankton Model list 
# ===========================================================
model_list <- list(
  #PHYTOPLANKTON DIVERSITY
  # ---------------------------------------------------------
  # PHYTOPLANKTON GENUS RICHNESS MODELS
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA/BRM_phyto_GD_L_260604.RDATA",
       name = "Phyto_Genus_Richness_L",
       resp = "genus_richness",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered", "MaxDepth", "daylight_hours"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic Variables for Phytoplankton Genus Richness (Linear)"),
  
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA/BRM_phyto_GD_P_260604.RDATA",
       name = "Phyto_Genus_Richness_P",
       resp = "genus_richness",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered", "MaxDepth", "daylight_hours"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic Variables for Phytoplankton Genus Richness (Quadratic)"),
  
  # ---------------------------------------------------------
  # PHYTOPLANKTON FUNCTIONAL EVENNESS MODELS
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA/BRM_phyto_FEve_L_260604.RDATA",
       name = "Phyto_FEve_L",
       resp = "FEve",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered", "MaxDepth", "daylight_hours"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic Variables for Phytoplankton Functional Evenness (Linear)"),

  list(file = "~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA/BRM_phyto_FEve_P_260604.RDATA",
       name = "Phyto_FEve_P",
       resp = "FEve",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered", "MaxDepth", "daylight_hours"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic Variables for Phytoplankton Functional Evenness (Quadratic)"),
  
  # ---------------------------------------------------------
  # PHYTOPLANKTON FUNCTIONAL DISPERSION MODELS
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA/BRM_phyto_FDis_L_260604.RDATA",
       name = "Phyto_FDis_L",
       resp = "FDis",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered", "MaxDepth", "daylight_hours"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic Variables for Phytoplankton Functional Dispersion (Linear)"),

  list(file = "~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA/BRM_phyto_FDis_P_260604.RDATA",
       name = "Phyto_FDis_P",
       resp = "FDis",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered", "MaxDepth", "daylight_hours"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic Variables for Phytoplankton Functional Dispersion (Quadratic)")
)


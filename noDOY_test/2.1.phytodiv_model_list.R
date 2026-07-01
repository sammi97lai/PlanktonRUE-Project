# ===========================================================
# Phytoplankton Model list 
# ===========================================================
model_list <- list(
  #PHYTOPLANKTON DIVERSITY
  # ---------------------------------------------------------
  # PHYTOPLANKTON GENUS RICHNESS MODELS
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/noDOY_test/RDATA/BRM_phyto_GD_P_260602.RDATA",
       name = "Phyto_Genus_Richness_P",
       resp = "genus_richness",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered", "MaxDepth", "daylight_hours"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic Variables for Phytoplankton Genus Richness (Quadratic)"),
  
  # ---------------------------------------------------------
  # PHYTOPLANKTON FUNCTIONAL EVENNESS MODELS
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/noDOY_test/RDATA/BRM_phyto_FEve_P_260602.RDATA",
       name = "Phyto_FEve_P",
       resp = "FEve",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered", "MaxDepth", "daylight_hours"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic Variables for Phytoplankton Functional Evenness (Quadratic)"),
  
  # ---------------------------------------------------------
  # PHYTOPLANKTON FUNCTIONAL DISPERSION MODELS
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/noDOY_test/RDATA/BRM_phyto_FDis_P_260602.RDATA",
       name = "Phyto_FDis_P",
       resp = "FDis",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered", "MaxDepth", "daylight_hours"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic Variables for Phytoplankton Functional Dispersion (Quadratic)")
)


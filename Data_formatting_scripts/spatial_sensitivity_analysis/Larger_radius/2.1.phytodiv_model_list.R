# ===========================================================
# Phytoplankton Model list 
# ===========================================================
model_list <- list(
  #PHYTOPLANKTON DIVERSITY
  # ---------------------------------------------------------
  # PHYTOPLANKTON GENUS RICHNESS MODELS
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/spatial_sensitivity_analysis/Larger_radius/RDATA/BRM_phyto_GD_P_260526.RDATA",
       name = "Phyto_Genus_Richness_P",
       resp = "genus_richness",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered", "MaxDepth", "daylight_hours")),
  
  # ---------------------------------------------------------
  # PHYTOPLANKTON SHANNON EVENNESS MODELS
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/spatial_sensitivity_analysis/Larger_radius/RDATA/BRM_phyto_SE_P_260526.RDATA",
       name = "Phyto_Evenness_P",
       resp = "phyto_evenness",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered", "MaxDepth", "daylight_hours")),
  
  # ---------------------------------------------------------
  # PHYTOPLANKTON FUNCTIONAL EVENNESS MODELS
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/spatial_sensitivity_analysis/Larger_radius/RDATA/BRM_phyto_FEve_P_260526.RDATA",
       name = "Phyto_FEve_P",
       resp = "phyto_FEve_PD",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered", "MaxDepth", "daylight_hours")),
  
  # ---------------------------------------------------------
  # PHYTOPLANKTON FUNCTIONAL DISPERSION MODELS
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/spatial_sensitivity_analysis/Larger_radius/RDATA/BRM_phyto_FDis_P_260526.RDATA",
       name = "Phyto_FDis_P",
       resp = "phyto_FDis_PD",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered", "MaxDepth", "daylight_hours"))
)


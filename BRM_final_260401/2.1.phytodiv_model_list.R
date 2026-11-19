# ===========================================================
# Phytoplankton Model list 
# ===========================================================
model_list <- list(
  #PHYTOPLANKTON DIVERSITY
  # ---------------------------------------------------------
  # PHYTOPLANKTON GENUS RICHNESS MODELS
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA/BRM_phyto_GD_L_260916.RDATA",
       name = "Phyto_Genus_Richness_L",
       resp = "genus_richness",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered", "MaxDepth", "daylight_hours")),
  
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA/BRM_phyto_GD_P_260916.RDATA",
       name = "Phyto_Genus_Richness_P",
       resp = "genus_richness",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered", "MaxDepth", "daylight_hours")),  
  
  # ---------------------------------------------------------
  # PHYTOPLANKTON SHANNON EVENNESS MODELS
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA/BRM_phyto_SE_L_260916.RDATA",
       name = "Phyto_Evenness_L",
       resp = "phyto_evenness",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered", "MaxDepth", "daylight_hours")),
  
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA/BRM_phyto_SE_P_260916.RDATA",
       name = "Phyto_Evenness_P",
       resp = "phyto_evenness",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered", "MaxDepth", "daylight_hours")),
  
  # ---------------------------------------------------------
  # PHYTOPLANKTON FUNCTIONAL EVENNESS MODELS
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA/BRM_phyto_FEve_L_260916.RDATA",
       name = "Phyto_FEve_L",
       resp = "phyto_FEve_PD",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered", "MaxDepth", "daylight_hours")),

  list(file = "~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA/BRM_phyto_FEve_P_260916.RDATA",
       name = "Phyto_FEve_P",
       resp = "phyto_FEve_PD",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered", "MaxDepth", "daylight_hours")),
  
  # ---------------------------------------------------------
  # PHYTOPLANKTON FUNCTIONAL DISPERSION MODELS
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA/BRM_phyto_FDis_L_260916.RDATA",
       name = "Phyto_FDis_L",
       resp = "phyto_FDis_PD",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered", "MaxDepth", "daylight_hours")),

  list(file = "~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA/BRM_phyto_FDis_P_260916.RDATA",
       name = "Phyto_FDis_P",
       resp = "phyto_FDis_PD",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered", "MaxDepth", "daylight_hours"))
)


# ===========================================================
# Cyano Dom Model list 
# ===========================================================
model_list <- list(
  
  # ---------------------------------------------------------
  # CYANOBACTERIA DOMINANCE MODELS - zero inflated
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/spatial_sensitivity_analysis/Larger_radius/RDATA/BRM_cyanodom_P_260526_Z.RDATA",
       name = "cyanodom_P_z",
       resp = "cyano_prop",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered", "MaxDepth", "daylight_hours"))
)


# ===========================================================
# Cyano Dom Model list 
# ===========================================================
model_list <- list(
  
  # ---------------------------------------------------------
  # CYANOBACTERIA DOMINANCE MODELS - zero inflated
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA/BRM_cyanodom_L_261009_z.RDATA",
       name = "cyanodom_L_z",
       resp = "cyano_prop",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered", "MaxDepth", "daylight_hours")),
  
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA/BRM_cyanodom_P_261009_z.RDATA",
       name = "cyanodom_P_z",
       resp = "cyano_prop",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered", "MaxDepth", "daylight_hours"))
)


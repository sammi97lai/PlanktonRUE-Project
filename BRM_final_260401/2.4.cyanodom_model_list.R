# ===========================================================
# Cyano Dom Model list 
# ===========================================================
model_list <- list(
  
  # ---------------------------------------------------------
  # CYANOBACTERIA DOMINANCE MODELS - zero inflated
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM/RDATA/BRM_cyanodom_L_260311_Z.RDATA",
       name = "cyanodom_L_z",
       resp = "cyano_prop",
       predictors = c("Temperature", "Salinity", "pH", "MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic Variables for Cyanobacteria Dominance (Linear)- zeroinflated"),
  
  list(file = "~/GitHub/PlanktonRUE-Project/BRM/RDATA/BRM_cyanodom_P_260311_Z.RDATA",
       name = "cyanodom_P_z",
       resp = "cyano_prop",
       predictors = c("Temperature", "Salinity", "pH", "MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic Variables for Cyanobacteria Dominance (Quadratic)- zeroinflated")
)


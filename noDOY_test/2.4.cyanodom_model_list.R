# ===========================================================
# Cyano Dom Model list 
# ===========================================================
model_list <- list(
  
  # ---------------------------------------------------------
  # CYANOBACTERIA DOMINANCE MODELS - zero inflated
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/noDOY_test/RDATA/BRM_cyanodom_P_260602_Z.RDATA",
       name = "cyanodom_P_z",
       resp = "cyano_prop",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered", "MaxDepth", "daylight_hours"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic Variables for Cyanobacteria Dominance (Quadratic)- zeroinflated")
)


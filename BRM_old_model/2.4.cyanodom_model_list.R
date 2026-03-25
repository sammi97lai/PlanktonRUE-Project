# ===========================================================
# Cyano Dom Model list 
# ===========================================================
model_list <- list(
  
  # ---------------------------------------------------------
  # CYANOBACTERIA DOMINANCE MODELS - zero inflated
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_old_model/RDATA/BRM_cyanodom_L_260324_Z.RDATA",
       name = "cyanodom_L_z",
       resp = "cyano_prop",
       predictors = c("Temperature", "Salinity", "pH", "log_N_P_ratio", "DOY"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic Variables for Cyanobacteria Dominance (Linear)- zeroinflated"),
  
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_old_model/RDATA/BRM_cyanodom_P_260324_Z.RDATA",
       name = "cyanodom_P_z",
       resp = "cyano_prop",
       predictors = c("Temperature", "Salinity", "pH", "log_N_P_ratio", "DOY"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic Variables for Cyanobacteria Dominance (Quadratic)- zeroinflated"),
  
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_old_model/RDATA/BRM_cyanodom_LI_260324_z.RDATA",
       name = "cyanodom_LI_z",
       resp = "cyano_prop",
       predictors = c("Temperature", "Salinity", "pH", "log_N_P_ratio", "DOY", "log_N_P_ratio:DOY"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic Variables with Interaction for Cyanobacteria Dominance (Linear)- zeroinflated"),
  
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_old_model/RDATA/BRM_cyanodom_PI_260324_z.RDATA",
       name = "cyanodom_PI_z",
       resp = "cyano_prop",
       predictors = c("Temperature", "Salinity", "pH", "log_N_P_ratio", "DOY", "log_N_P_ratio:DOY"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic Variables with Interaction for Cyanobacteria Dominance (Quadratic)- zeroinflated")
)


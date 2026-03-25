# ===========================================================
# Phytoplankton Model list 
# ===========================================================
model_list <- list(
  #PHYTOPLANKTON DIVERSITY
  # ---------------------------------------------------------
  # PHYTOPLANKTON GENUS RICHNESS MODELS
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_old_model/RDATA/BRM_phyto_GD_L_260324.RDATA",
       name = "Phyto_Genus_Richness_L",
       resp = "genus_richness",
       predictors = c("Temperature", "Salinity", "pH", "log_N_P_ratio", "DOY"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic Variables for Phytoplankton Genus Richness (Linear)"),
  
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_old_model/RDATA/BRM_phyto_GD_P_260324.RDATA",
       name = "Phyto_Genus_Richness_P",
       resp = "genus_richness",
       predictors = c("Temperature", "Salinity", "pH", "log_N_P_ratio", "DOY"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic Variables for Phytoplankton Genus Richness (Quadratic)"),
  
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_old_model/RDATA/BRM_phyto_GD_LI_260324.RDATA",
       name = "Phyto_Genus_Richness_LI",
       resp = "genus_richness",
       predictors = c("Temperature", "Salinity", "pH", "log_N_P_ratio", "DOY", "log_N_P_ratio:DOY"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic Variables with Interaction for Phytoplankton Genus Richness (Linear)"),

  list(file = "~/GitHub/PlanktonRUE-Project/BRM_old_model/RDATA/BRM_phyto_GD_PI_260324.RDATA",
       name = "Phyto_Genus_Richness_PI",
       resp = "genus_richness",
       predictors = c("Temperature", "Salinity", "pH", "log_N_P_ratio", "DOY", "log_N_P_ratio:DOY"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic Variables with Interaction for Phytoplankton Genus Richness (Quadratic)"),

  # ---------------------------------------------------------
  # PHYTOPLANKTON FUNCTIONAL EVENNESS MODELS
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_old_model/RDATA/BRM_phyto_FEve_L_260324.RDATA",
       name = "Phyto_FEve_L",
       resp = "FEve",
       predictors = c("Temperature", "Salinity", "pH", "log_N_P_ratio", "DOY"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic Variables for Phytoplankton Functional Evenness (Linear)"),

  list(file = "~/GitHub/PlanktonRUE-Project/BRM_old_model/RDATA/BRM_phyto_FEve_P_260324.RDATA",
       name = "Phyto_FEve_P",
       resp = "FEve",
       predictors = c("Temperature", "Salinity", "pH", "log_N_P_ratio", "DOY"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic Variables for Phytoplankton Functional Evenness (Quadratic)"),

  list(file = "~/GitHub/PlanktonRUE-Project/BRM_old_model/RDATA/BRM_phyto_FEve_LI_260324.RDATA",
       name = "Phyto_FEve_LI",
       resp = "FEve",
       predictors = c("Temperature", "Salinity", "pH", "log_N_P_ratio", "DOY", "log_N_P_ratio:DOY"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic Variables with Interaction for Phytoplankton Functional Evenness (Linear)"),

  list(file = "~/GitHub/PlanktonRUE-Project/BRM_old_model/RDATA/BRM_phyto_FEve_PI_260324.RDATA",
       name = "Phyto_FEve_PI",
       resp = "FEve",
       predictors = c("Temperature", "Salinity", "pH", "log_N_P_ratio", "DOY", "log_N_P_ratio:DOY"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic Variables with Interaction for Phytoplankton Functional Evenness (Quadratic)"),
  
  # ---------------------------------------------------------
  # PHYTOPLANKTON FUNCTIONAL DISPERSION MODELS
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_old_model/RDATA/BRM_phyto_FDis_L_260324.RDATA",
       name = "Phyto_FDis_L",
       resp = "FDis",
       predictors = c("Temperature", "Salinity", "pH", "log_N_P_ratio", "DOY"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic Variables for Phytoplankton Functional Dispersion (Linear)"),

  list(file = "~/GitHub/PlanktonRUE-Project/BRM_old_model/RDATA/BRM_phyto_FDis_P_260324.RDATA",
       name = "Phyto_FDis_P",
       resp = "FDis",
       predictors = c("Temperature", "Salinity", "pH", "log_N_P_ratio", "DOY"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic Variables for Phytoplankton Functional Dispersion (Quadratic)"),

  list(file = "~/GitHub/PlanktonRUE-Project/BRM_old_model/RDATA/BRM_phyto_FDis_LI_260324.RDATA",
       name = "Phyto_FDis_LI",
       resp = "FDis",
       predictors = c("Temperature", "Salinity", "pH", "log_N_P_ratio", "DOY", "log_N_P_ratio:DOY"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic Variables with Interaction for Phytoplankton Functional Dispersion (Linear)"),

  list(file = "~/GitHub/PlanktonRUE-Project/BRM_old_model/RDATA/BRM_phyto_FDis_PI_260324.RDATA",
       name = "Phyto_FDis_PI",
       resp = "FDis",
       predictors = c("Temperature", "Salinity", "pH", "log_N_P_ratio", "DOY", "log_N_P_ratio:DOY"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic Variables with Interaction for Phytoplankton Functional Dispersion (Quadratic)")
)




# ===========================================================
# Phytoplankton RUE Model list 
# ===========================================================
model_list <- list(
  #PHYTOPLANKTON RUE MODELS
  # ---------------------------------------------------------
  # PHYTOPLANKTON RUE WITH PHOSPHORUS MODELS
  # ---------------------------------------------------------
  list(file = "~/GitHub/Phytoplankton-project/phyto_RUE/BRM_260109/BRM_phytoRUEP_L_260109.RDATA",
       name = "phytoRUEP_L",
       resp = "log_RUE_A_P",
       predictors = c("Temperature", "Salinity", "pH", "genus_richness",
                      "FEve", "FDis", "log_N_P_ratio", "DOY"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic AND Biotic Variables for Phytoplankton RUE per P (Linear)"),

  list(file = "~/GitHub/Phytoplankton-project/phyto_RUE/BRM_260109/BRM_phytoRUEP_P_260109.RDATA",
       name = "phytoRUEP_P",
       resp = "log_RUE_A_P",
       predictors = c("Temperature", "Salinity", "pH", "genus_richness",
                      "FEve", "FDis", "log_N_P_ratio", "DOY"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic AND Biotic Variables for Phytoplankton RUE per P (Quadratic)"),

  list(file = "~/GitHub/Phytoplankton-project/phyto_RUE/BRM_260109/BRM_phytoRUEP_LI_260109.RDATA",
       name = "phytoRUEP_LI",
       resp = "log_RUE_A_P",
       predictors = c("Temperature", "Salinity", "pH", "genus_richness",
                      "FEve", "FDis", "log_N_P_ratio", "DOY",  "log_N_P_ratio:DOY"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic AND Biotic Variables with Interaction for Phytoplankton RUE per P (Linear)"),

  list(file = "~/GitHub/Phytoplankton-project/phyto_RUE/BRM_260109/BRM_phytoRUEP_PI_260109.RDATA",
       name = "phytoRUEP_PI",
       resp = "log_RUE_A_P",
       predictors = c("Temperature", "Salinity", "pH", "genus_richness",
                      "FEve", "FDis", "log_N_P_ratio", "DOY",  "log_N_P_ratio:DOY"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic AND Biotic Variables with Interaction for Phytoplankton RUE per P (Quadratic)"),

  # ---------------------------------------------------------
  # PHYTOPLANKTON RUE WITH NITROGEN MODELS
  # ---------------------------------------------------------
  list(file = "~/GitHub/Phytoplankton-project/phyto_RUE/BRM_260109/BRM_phytoRUEN_L_260109.RDATA",
       name = "phytoRUEN_L",
       resp = "log_RUE_A_N",
       predictors = c("Temperature", "Salinity", "pH", "genus_richness",
                      "FEve", "FDis", "log_N_P_ratio", "DOY"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic AND Biotic Variables for Phytoplankton RUE per N (Linear)"),

  list(file = "~/GitHub/Phytoplankton-project/phyto_RUE/BRM_260109/BRM_phytoRUEN_P_260109.RDATA",
       name = "phytoRUEN_P",
       resp = "log_RUE_A_N",
       predictors = c("Temperature", "Salinity", "pH", "genus_richness",
                      "FEve", "FDis", "log_N_P_ratio", "DOY"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic AND Biotic Variables for Phytoplankton RUE per N (Quadratic)"),

  list(file = "~/GitHub/Phytoplankton-project/phyto_RUE/BRM_260109/BRM_phytoRUEN_LI_260109.RDATA",
       name = "phytoRUEN_LI",
       resp = "log_RUE_A_N",
       predictors = c("Temperature", "Salinity", "pH", "genus_richness",
                      "FEve", "FDis", "log_N_P_ratio", "DOY",  "log_N_P_ratio:DOY"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic AND Biotic Variables with Interaction for Phytoplankton RUE per N (Linear)"),

  list(file = "~/GitHub/Phytoplankton-project/phyto_RUE/BRM_260109/BRM_phytoRUEN_PI_260109.RDATA",
       name = "phytoRUEN_PI",
       resp = "log_RUE_A_N",
       predictors = c("Temperature", "Salinity", "pH", "genus_richness",
                      "FEve", "FDis", "log_N_P_ratio", "DOY",  "log_N_P_ratio:DOY"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic AND Biotic Variables with Interaction for Phytoplankton RUE per N (Quadratic)"),

  #PHYTOPLANKTON RUE MODELS WITH CYANO DOM
  # ---------------------------------------------------------
  # PHYTOPLANKTON RUE WITH PHOSPHORUS MODELS WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/Phytoplankton-project/phyto_RUE/BRM_260109/BRM_phytoRUEP_CD_L_260109.RDATA",
       name = "phytoRUEP_CD_L",
       resp = "log_RUE_A_P",
       predictors = c("Temperature", "Salinity", "pH", "genus_richness",
                      "FEve", "FDis", "cyano_prop", "log_N_P_ratio", "DOY"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic AND Biotic (incl. CD) Variables for Phytoplankton RUE per P (Linear)"),

  list(file = "~/GitHub/Phytoplankton-project/phyto_RUE/BRM_260109/BRM_phytoRUEP_CD_P_260109.RDATA",
       name = "phytoRUEP_CD_P",
       resp = "log_RUE_A_P",
       predictors = c("Temperature", "Salinity", "pH", "genus_richness",
                      "FEve", "FDis", "cyano_prop", "log_N_P_ratio", "DOY"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic AND Biotic (incl. CD) Variables for Phytoplankton RUE per P (Quadratic)"),

  list(file = "~/GitHub/Phytoplankton-project/phyto_RUE/BRM_260109/BRM_phytoRUEP_CD_LI_260109.RDATA",
       name = "phytoRUEP_CD_LI",
       resp = "log_RUE_A_P",
       predictors = c("Temperature", "Salinity", "pH", "genus_richness",
                      "FEve", "FDis", "cyano_prop", "log_N_P_ratio", "DOY",
                      "log_N_P_ratio:DOY"),
       ncol_ce = 4,
       ce_title = "Marginal Effects of Abiotic AND Biotic (incl. CD) Variables with Interaction for Phytoplankton RUE per P (Linear)"),

  list(file = "~/GitHub/Phytoplankton-project/phyto_RUE/BRM_260109/BRM_phytoRUEP_CD_PI_260109.RDATA",
       name = "phytoRUEP_CD_PI",
       resp = "log_RUE_A_P",
       predictors = c("Temperature", "Salinity", "pH", "genus_richness",
                      "FEve", "FDis", "cyano_prop", "log_N_P_ratio", "DOY",
                      "log_N_P_ratio:DOY"),
       ncol_ce = 4,
       ce_title = "Marginal Effects of Abiotic AND Biotic (incl. CD) Variables with Interaction for Phytoplankton RUE per P (Quadratic)"),

  # ---------------------------------------------------------
  # PHYTOPLANKTON RUE WITH NITROGEN MODELS WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/Phytoplankton-project/phyto_RUE/BRM_260109/BRM_phytoRUEN_CD_L_260109.RDATA",
       name = "phytoRUEN_CD_L",
       resp = "log_RUE_A_N",
       predictors = c("Temperature", "Salinity", "pH", "genus_richness",
                      "FEve", "FDis", "cyano_prop", "log_N_P_ratio", "DOY"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic AND Biotic (incl. CD) Variables for Phytoplankton RUE per N (Linear)"),

  list(file = "~/GitHub/Phytoplankton-project/phyto_RUE/BRM_260109/BRM_phytoRUEN_CD_P_260109.RDATA",
       name = "phytoRUEN_CD_P",
       resp = "log_RUE_A_N",
       predictors = c("Temperature", "Salinity", "pH", "genus_richness",
                      "FEve", "FDis", "cyano_prop", "log_N_P_ratio", "DOY"),
       ncol_ce = 3,
       ce_title = "Marginal Effects of Abiotic AND Biotic (incl. CD) Variables for Phytoplankton RUE per N (Quadratic)"),

  list(file = "~/GitHub/Phytoplankton-project/phyto_RUE/BRM_260109/BRM_phytoRUEN_CD_LI_260109.RDATA",
       name = "phytoRUEN_CD_LI",
       resp = "log_RUE_A_N",
       predictors = c("Temperature", "Salinity", "pH", "genus_richness",
                      "FEve", "FDis", "cyano_prop", "log_N_P_ratio", "DOY",
                      "log_N_P_ratio:DOY"),
       ncol_ce = 4,
       ce_title = "Marginal Effects of Abiotic AND Biotic (incl. CD) Variables with Interaction for Phytoplankton RUE per N (Linear)"),

  list(file = "~/GitHub/Phytoplankton-project/phyto_RUE/BRM_260109/BRM_phytoRUEN_CD_PI_260109.RDATA",
       name = "phytoRUEN_CD_PI",
       resp = "log_RUE_A_N",
       predictors = c("Temperature", "Salinity", "pH", "genus_richness",
                      "FEve", "FDis", "cyano_prop", "log_N_P_ratio", "DOY",
                      "log_N_P_ratio:DOY"),
       ncol_ce = 4,
       ce_title = "Marginal Effects of Abiotic AND Biotic (incl. CD) Variables with Interaction for Phytoplankton RUE per N (Quadratic)")
  )

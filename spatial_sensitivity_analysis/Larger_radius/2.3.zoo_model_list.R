# ===========================================================
# Zooplankton Model list 
# ===========================================================
model_list <- list(
  #ZOOPLANKTON DIVERSITY WITH CYANO DOM MODELS#
  # ---------------------------------------------------------
  # ZOOPLANKTON SPECIES RICHNESS MODELS WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/spatial_sensitivity_analysis/Larger_radius/RDATA/BRM_zoo_SR_CD_P_260526.RDATA",
       name = "Zoo_SR_CD_P",
       resp = "species_richness",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered",
                      "genus_richness", "FEve", "FDis", "cyano_prop",
                      "Zooplankton_MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 4,
       ce_title = "Marginal Effects of Abiotic and Biotic (incl. CD) Variables for Zooplankton Species Richness (Quadratic)"),

  # ---------------------------------------------------------
  # ZOOPLANKTON FUNCTIONAL EVENNESS MODELS WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/spatial_sensitivity_analysis/Larger_radius/RDATA/BRM_zoo_FE_CD_P_260526.RDATA",
       name = "Zoo_FEve_CD_P",
       resp = "zoo_FEve",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered",
                      "genus_richness", "FEve", "FDis", "cyano_prop",
                      "Zooplankton_MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 4,
       ce_title = "Marginal Effects of Abiotic and Biotic (incl. CD) Variables for Zooplankton Functional Evenness (Quadratic)"),


  # ---------------------------------------------------------
  # ZOOPLANKTON FUNCTIONAL DISPERSION MODELS WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/spatial_sensitivity_analysis/Larger_radius/RDATA/BRM_zoo_FD_CD_P_260526.RDATA",
       name = "Zoo_FDis_CD_P",
       resp = "zoo_FDis",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered",
                      "genus_richness", "FEve", "FDis", "cyano_prop",
                      "Zooplankton_MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 4,
       ce_title = "Marginal Effects of Abiotic and Biotic (incl. CD) Variables for Zooplankton Functional Dispersion (Quadratic)"),


  #ZOOPLANKTON RUE MODELS#
  # ---------------------------------------------------------
  # ZOOPLANKTON RUE WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/spatial_sensitivity_analysis/Larger_radius/RDATA/BRM_zoo_RUE_CD_PI_260526.RDATA",
       name = "zoo_RUE_CD_PI",
       resp = "log_RUE_zooplankton",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered",
                      "genus_richness", "FEve", "FDis", "cyano_prop",
                      "species_richness", "zoo_FEve", "zoo_FDis",
                      "genus_richness:Temperature", "genus_richness:Salinity",
                      "genus_richness:TotalNitrogen_unfiltered",
                      "FEve:Temperature", "FEve:Salinity",
                      "FEve:TotalNitrogen_unfiltered",
                      "FDis:Temperature", "FDis:Salinity",
                      "FDis:TotalNitrogen_unfiltered",
                      "cyano_prop:Temperature", "cyano_prop:Salinity",
                      "cyano_prop:TotalNitrogen_unfiltered",
                      "species_richness:Temperature", "species_richness:Salinity",
                      "species_richness:TotalNitrogen_unfiltered",
                      "zoo_FEve:Temperature", "zoo_FEve:Salinity",
                      "zoo_FEve:TotalNitrogen_unfiltered",
                      "zoo_FDis:Temperature", "zoo_FDis:Salinity",
                      "zoo_FDis:TotalNitrogen_unfiltered",
                      "Zooplankton_MaxDepth", "DOY", "daylight_hours"),
       ncol_ce = 6,
       ce_title = "Marginal Effects of Abiotic and Biotic (incl. CD) Variables with Interactions for Zooplankton RUE (Quadratic)")
)




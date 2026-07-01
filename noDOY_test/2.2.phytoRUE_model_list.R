# ===========================================================
# Phytoplankton RUE Model list 
# ===========================================================
model_list <- list(
  #PHYTOPLANKTON RUE MODELS WITH CYANO DOM
  # ---------------------------------------------------------
  # PHYTOPLANKTON RUE WITH PHOSPHORUS MODELS WITH CYANO DOM
  # ---------------------------------------------------------

  list(file = "~/GitHub/PlanktonRUE-Project/noDOY_test/RDATA/BRM_phytoRUEP_CD_PI_260602.RDATA",
       name = "phytoRUEP_CD_PI",
       resp = "log_RUE_B_P",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered", "genus_richness",
                      "FEve", "FDis", "cyano_prop", "species_richness", "zoo_FEve", "zoo_FDis",
                      "genus_richness:Temperature", "genus_richness:Salinity",
                      "genus_richness:TotalNitrogen_unfiltered",
                      "FEve:Temperature", "FEve:Salinity", "FEve:TotalNitrogen_unfiltered", 
                      "FDis:Temperature", "FDis:Salinity", "FDis:TotalNitrogen_unfiltered",
                      "cyano_prop:Temperature", "cyano_prop:Salinity",
                      "cyano_prop:TotalNitrogen_unfiltered",
                      "species_richness:Temperature", "species_richness:Salinity",
                      "species_richness:TotalNitrogen_unfiltered",
                      "zoo_FEve:Temperature", "zoo_FEve:Salinity", "zoo_FEve:TotalNitrogen_unfiltered", 
                      "zoo_FDis:Temperature", "zoo_FDis:Salinity", "zoo_FDis:TotalNitrogen_unfiltered",
                      "MaxDepth", "daylight_hours"),
       ncol_ce = 6,
       ce_title = "Marginal Effects of Abiotic AND Biotic (incl. CD) Variables with Interaction for Phytoplankton RUE per P (Quadratic)"),

  # ---------------------------------------------------------
  # PHYTOPLANKTON RUE WITH NITROGEN MODELS WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/noDOY_test/RDATA/BRM_phytoRUEN_CD_PI_260602.RDATA",
       name = "phytoRUEN_CD_PI",
       resp = "log_RUE_B_N",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered", "genus_richness",
                      "FEve", "FDis", "cyano_prop", "species_richness", "zoo_FEve", "zoo_FDis",
                      "genus_richness:Temperature", "genus_richness:Salinity",
                      "genus_richness:TotalNitrogen_unfiltered",
                      "FEve:Temperature", "FEve:Salinity", "FEve:TotalNitrogen_unfiltered", 
                      "FDis:Temperature", "FDis:Salinity", "FDis:TotalNitrogen_unfiltered",
                      "cyano_prop:Temperature", "cyano_prop:Salinity",
                      "cyano_prop:TotalNitrogen_unfiltered",
                      "species_richness:Temperature", "species_richness:Salinity",
                      "species_richness:TotalNitrogen_unfiltered",
                      "zoo_FEve:Temperature", "zoo_FEve:Salinity", "zoo_FEve:TotalNitrogen_unfiltered", 
                      "zoo_FDis:Temperature", "zoo_FDis:Salinity", "zoo_FDis:TotalNitrogen_unfiltered",
                      "MaxDepth", "daylight_hours"),
       ncol_ce = 6,
       ce_title = "Marginal Effects of Abiotic AND Biotic (incl. CD) Variables with Interaction for Phytoplankton RUE per N (Quadratic)")
  )

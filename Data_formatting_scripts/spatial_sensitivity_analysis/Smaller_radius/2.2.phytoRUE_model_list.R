# ===========================================================
# Phytoplankton RUE Model list 
# ===========================================================
model_list <- list(
  #PHYTOPLANKTON RUE MODELS WITH CYANO DOM
  # ---------------------------------------------------------
  # PHYTOPLANKTON RUE WITH PHOSPHORUS MODELS WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/spatial_sensitivity_analysis/Smaller_radius/RDATA/BRM_phytoRUEP_CD_PI_260526.RDATA",
       name = "phytoRUEP_CD_PI",
       resp = "log_RUE_B_P",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered", 
                      "genus_richness", "phyto_evenness", "phyto_FEve_PRUE", "phyto_FDis_PRUE", "cyano_prop", 
                      "species_richness","zoo_evenness", "zoo_FEve", "zoo_FDis", 
                      "genus_richness:Temperature", "genus_richness:Salinity",
                      "genus_richness:TotalNitrogen_unfiltered",
                      "phyto_evenness:Temperature", "phyto_evenness:Salinity",
                      "phyto_evenness:TotalNitrogen_unfiltered",
                      "phyto_FEve_PRUE:Temperature", "phyto_FEve_PRUE:Salinity", "phyto_FEve_PRUE:TotalNitrogen_unfiltered", 
                      "phyto_FDis_PRUE:Temperature", "phyto_FDis_PRUE:Salinity", "phyto_FDis_PRUE:TotalNitrogen_unfiltered", 
                      "cyano_prop:Temperature", "cyano_prop:Salinity",
                      "cyano_prop:TotalNitrogen_unfiltered",
                      "species_richness:Temperature", "species_richness:Salinity",
                      "species_richness:TotalNitrogen_unfiltered",
                      "zoo_evenness:Temperature", "zoo_evenness:Salinity",
                      "zoo_evenness:TotalNitrogen_unfiltered",
                      "zoo_FEve:Temperature", "zoo_FEve:Salinity", 
                      "zoo_FEve:TotalNitrogen_unfiltered", 
                      "zoo_FDis:Temperature", "zoo_FDis:Salinity", 
                      "zoo_FDis:TotalNitrogen_unfiltered", 
                      "MaxDepth", "daylight_hours")),

  # ---------------------------------------------------------
  # PHYTOPLANKTON RUE WITH NITROGEN MODELS WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/spatial_sensitivity_analysis/Smaller_radius/RDATA/BRM_phytoRUEN_CD_PI_260526.RDATA",
       name = "phytoRUEN_CD_PI",
       resp = "log_RUE_B_N",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered", 
                      "genus_richness", "phyto_evenness", "phyto_FEve_PRUE", "phyto_FDis_PRUE", "cyano_prop", 
                      "species_richness","zoo_evenness", "zoo_FEve", "zoo_FDis", 
                      "genus_richness:Temperature", "genus_richness:Salinity",
                      "genus_richness:TotalNitrogen_unfiltered",
                      "phyto_evenness:Temperature", "phyto_evenness:Salinity",
                      "phyto_evenness:TotalNitrogen_unfiltered",
                      "phyto_FEve_PRUE:Temperature", "phyto_FEve_PRUE:Salinity", "phyto_FEve_PRUE:TotalNitrogen_unfiltered", 
                      "phyto_FDis_PRUE:Temperature", "phyto_FDis_PRUE:Salinity", "phyto_FDis_PRUE:TotalNitrogen_unfiltered", 
                      "cyano_prop:Temperature", "cyano_prop:Salinity",
                      "cyano_prop:TotalNitrogen_unfiltered",
                      "species_richness:Temperature", "species_richness:Salinity",
                      "species_richness:TotalNitrogen_unfiltered",
                      "zoo_evenness:Temperature", "zoo_evenness:Salinity",
                      "zoo_evenness:TotalNitrogen_unfiltered",
                      "zoo_FEve:Temperature", "zoo_FEve:Salinity", 
                      "zoo_FEve:TotalNitrogen_unfiltered", 
                      "zoo_FDis:Temperature", "zoo_FDis:Salinity", 
                      "zoo_FDis:TotalNitrogen_unfiltered", 
                      "MaxDepth", "daylight_hours"))
  )

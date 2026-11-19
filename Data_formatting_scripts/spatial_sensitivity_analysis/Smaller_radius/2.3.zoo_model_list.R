# ===========================================================
# Zooplankton Model list 
# ===========================================================
model_list <- list(
  #ZOOPLANKTON DIVERSITY WITH CYANO DOM MODELS#
  # ---------------------------------------------------------
  # ZOOPLANKTON SPECIES RICHNESS MODELS WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/spatial_sensitivity_analysis/Smaller_radius/RDATA/BRM_zoo_SR_CD_P_260526.RDATA",
       name = "Zoo_Species_Richness_CD_P",
       resp = "species_richness",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered",
                      "genus_richness", "phyto_evenness", "phyto_FEve_ZD", "phyto_FDis_ZD", "cyano_prop",
                      "Zooplankton_MaxDepth", "daylight_hours")),

  # ---------------------------------------------------------
  # ZOOPLANKTON SHANNON EVENNESS MODELS WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/spatial_sensitivity_analysis/Smaller_radius/RDATA/BRM_zoo_SE_CD_P_260526.RDATA",
       name = "Zoo_SE_CD_P",
       resp = "zoo_evenness",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered",
                      "genus_richness", "phyto_evenness", "phyto_FEve_ZD", "phyto_FDis_ZD", "cyano_prop",
                      "Zooplankton_MaxDepth", "daylight_hours")),
  
  # ---------------------------------------------------------
  # ZOOPLANKTON FUNCTIONAL EVENNESS MODELS WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/spatial_sensitivity_analysis/Smaller_radius/RDATA/BRM_zoo_FE_CD_P_260526.RDATA",
       name = "Zoo_FEve_CD_P",
       resp = "zoo_FEve",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered",
                      "genus_richness", "phyto_evenness", "phyto_FEve_ZD", "phyto_FDis_ZD", "cyano_prop",
                      "Zooplankton_MaxDepth", "daylight_hours")),


  # ---------------------------------------------------------
  # ZOOPLANKTON FUNCTIONAL DISPERSION MODELS WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/spatial_sensitivity_analysis/Smaller_radius/RDATA/BRM_zoo_FD_CD_P_260526.RDATA",
       name = "Zoo_FDis_CD_P",
       resp = "zoo_FDis",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered",
                      "genus_richness", "phyto_evenness", "phyto_FEve_ZD", "phyto_FDis_ZD", "cyano_prop",
                      "Zooplankton_MaxDepth", "daylight_hours")),


  #ZOOPLANKTON RUE MODELS#
  # ---------------------------------------------------------
  # ZOOPLANKTON RUE WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/spatial_sensitivity_analysis/Smaller_radius/RDATA/BRM_zoo_RUE_CD_PI_260526.RDATA",
       name = "zoo_RUE_CD_PI",
       resp = "log_RUE_zooplankton",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered",
                      "genus_richness", "phyto_evenness", "phyto_FEve_ZRUE", "phyto_FDis_ZRUE", "cyano_prop",
                      "species_richness","zoo_evenness", "zoo_FEve", "zoo_FDis",
                      "genus_richness:Temperature", "genus_richness:Salinity",
                      "genus_richness:TotalNitrogen_unfiltered",
                      "phyto_evenness:Temperature", "phyto_evenness:Salinity",
                      "phyto_evenness:TotalNitrogen_unfiltered",
                      "phyto_FEve_ZRUE:Temperature", "phyto_FEve_ZRUE:Salinity",
                      "phyto_FEve_ZRUE:TotalNitrogen_unfiltered",
                      "phyto_FDis_ZRUE:Temperature", "phyto_FDis_ZRUE:Salinity",
                      "phyto_FDis_ZRUE:TotalNitrogen_unfiltered",
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
                      "Zooplankton_MaxDepth", "daylight_hours"))
)




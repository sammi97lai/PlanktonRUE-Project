# ===========================================================
# Zooplankton Model list 
# ===========================================================
model_list <- list(
  #ZOOPLANKTON DIVERSITY#
  # ---------------------------------------------------------
  # ZOOPLANKTON SPECIES RICHNESS MODELS
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA/BRM_zoo_SR_L_261013_All.RDATA",
       name = "Zoo_Species_Richness_L",
       resp = "species_richness",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered",
                      "genus_richness", "phyto_evenness", "phyto_FEve_ZD", "phyto_FDis_ZD",
                      "Zooplankton_MaxDepth", "daylight_hours")),
  
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA/BRM_zoo_SR_P_261013_All.RDATA",
       name = "Zoo_Species_Richness_P",
       resp = "species_richness",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered",
                      "genus_richness", "phyto_evenness", "phyto_FEve_ZD", "phyto_FDis_ZD",
                      "Zooplankton_MaxDepth", "daylight_hours")),
  
  # ---------------------------------------------------------
  # ZOOPLANKTON SHANNON EVENNESS MODELS
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA/BRM_zoo_SE_L_261013_All.RDATA",
       name = "Zoo_SE_L",
       resp = "zoo_evenness",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered",
                      "genus_richness", "phyto_evenness", "phyto_FEve_ZD", "phyto_FDis_ZD",
                      "Zooplankton_MaxDepth", "daylight_hours")),
  
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA/BRM_zoo_SE_P_261013_All.RDATA",
       name = "Zoo_SE_P",
       resp = "zoo_evenness",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered",
                      "genus_richness", "phyto_evenness", "phyto_FEve_ZD", "phyto_FDis_ZD",
                      "Zooplankton_MaxDepth", "daylight_hours")),
  
  
  # ---------------------------------------------------------
  # ZOOPLANKTON FUNCTIONAL EVENNESS MODELS
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA/BRM_zoo_FE_L_261013_All.RDATA",
       name = "Zoo_FEve_L",
       resp = "zoo_FEve",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered",
                      "genus_richness", "phyto_evenness", "phyto_FEve_ZD", "phyto_FDis_ZD",
                      "Zooplankton_MaxDepth", "daylight_hours")),
  
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA/BRM_zoo_FE_P_261013_All.RDATA",
       name = "Zoo_FEve_P",
       resp = "zoo_FEve",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered",
                      "genus_richness", "phyto_evenness", "phyto_FEve_ZD", "phyto_FDis_ZD",
                      "Zooplankton_MaxDepth", "daylight_hours")),
  
  # ---------------------------------------------------------
  # ZOOPLANKTON FUNCTIONAL DISPERSION MODELS
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA/BRM_zoo_FD_L_261013_All.RDATA",
       name = "Zoo_FDis_L",
       resp = "zoo_FDis",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered",
                      "genus_richness", "phyto_evenness", "phyto_FEve_ZD", "phyto_FDis_ZD",
                      "Zooplankton_MaxDepth", "daylight_hours")),
  
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA/BRM_zoo_FD_P_261013_All.RDATA",
       name = "Zoo_FDis_P",
       resp = "zoo_FDis",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered",
                      "genus_richness", "phyto_evenness", "phyto_FEve_ZD", "phyto_FDis_ZD",
                      "Zooplankton_MaxDepth", "daylight_hours")),
  
  
  #ZOOPLANKTON DIVERSITY WITH CYANO DOM MODELS#
  # ---------------------------------------------------------
  # ZOOPLANKTON SPECIES RICHNESS MODELS WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA/BRM_zoo_SR_CD_L_261013_All.RDATA",
       name = "Zoo_Species_Richness_CD_L",
       resp = "species_richness",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered",
                      "genus_richness", "phyto_evenness", "phyto_FEve_ZD", "phyto_FDis_ZD", "cyano_prop",
                      "Zooplankton_MaxDepth", "daylight_hours")),
  
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA/BRM_zoo_SR_CD_P_261013_All.RDATA",
       name = "Zoo_Species_Richness_CD_P",
       resp = "species_richness",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered",
                      "genus_richness", "phyto_evenness", "phyto_FEve_ZD", "phyto_FDis_ZD", "cyano_prop",
                      "Zooplankton_MaxDepth", "daylight_hours")),
  
  # ---------------------------------------------------------
  # ZOOPLANKTON SHANNON EVENNESS MODELS WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA/BRM_zoo_SE_CD_L_261013_All.RDATA",
       name = "Zoo_SE_CD_L",
       resp = "zoo_evenness",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered",
                      "genus_richness", "phyto_evenness", "phyto_FEve_ZD", "phyto_FDis_ZD", "cyano_prop",
                      "Zooplankton_MaxDepth", "daylight_hours")),
  
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA/BRM_zoo_SE_CD_P_261013_All.RDATA",
       name = "Zoo_SE_CD_P",
       resp = "zoo_evenness",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered",
                      "genus_richness", "phyto_evenness", "phyto_FEve_ZD", "phyto_FDis_ZD", "cyano_prop",
                      "Zooplankton_MaxDepth", "daylight_hours")),
  
  # ---------------------------------------------------------
  # ZOOPLANKTON FUNCTIONAL EVENNESS MODELS WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA/BRM_zoo_FE_CD_L_261013_All.RDATA",
       name = "Zoo_FEve_CD_L",
       resp = "zoo_FEve",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered",
                      "genus_richness", "phyto_evenness", "phyto_FEve_ZD", "phyto_FDis_ZD", "cyano_prop",
                      "Zooplankton_MaxDepth", "daylight_hours")),
  
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA/BRM_zoo_FE_CD_P_261013_All.RDATA",
       name = "Zoo_FEve_CD_P",
       resp = "zoo_FEve",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered",
                      "genus_richness", "phyto_evenness", "phyto_FEve_ZD", "phyto_FDis_ZD", "cyano_prop",
                      "Zooplankton_MaxDepth", "daylight_hours")),
  
  
  # ---------------------------------------------------------
  # ZOOPLANKTON FUNCTIONAL DISPERSION MODELS WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA/BRM_zoo_FD_CD_L_261013_All.RDATA",
       name = "Zoo_FDis_CD_L",
       resp = "zoo_FDis",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered",
                      "genus_richness", "phyto_evenness", "phyto_FEve_ZD", "phyto_FDis_ZD", "cyano_prop",
                      "Zooplankton_MaxDepth", "daylight_hours")),
  
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA/BRM_zoo_FD_CD_P_261013_All.RDATA",
       name = "Zoo_FDis_CD_P",
       resp = "zoo_FDis",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered",
                      "genus_richness", "phyto_evenness", "phyto_FEve_ZD", "phyto_FDis_ZD", "cyano_prop",
                      "Zooplankton_MaxDepth", "daylight_hours")),
  
  #ZOOPLANKTON RUE MODELS#
  # ---------------------------------------------------------
  # ZOOPLANKTON RUE
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA/BRM_zoo_RUE_L_261013_All.RDATA",
       name = "zoo_RUE_L",
       resp = "log_RUE_zooplankton",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered",
                      "genus_richness", "phyto_evenness", "phyto_FEve_ZRUE", "phyto_FDis_ZRUE",
                      "species_richness", "zoo_evenness", "zoo_FEve", "zoo_FDis",
                      "Zooplankton_MaxDepth", "daylight_hours")),
  
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA/BRM_zoo_RUE_P_261013_All.RDATA",
       name = "zoo_RUE_P",
       resp = "log_RUE_zooplankton",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered",
                      "genus_richness", "phyto_evenness", "phyto_FEve_ZRUE", "phyto_FDis_ZRUE",
                      "species_richness", "zoo_evenness", "zoo_FEve", "zoo_FDis",
                      "Zooplankton_MaxDepth", "daylight_hours")),
  
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA/BRM_zoo_RUE_LI_261013_All.RDATA",
       name = "zoo_RUE_LI",
       resp = "log_RUE_zooplankton",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered",
                      "genus_richness", "phyto_evenness", "phyto_FEve_ZRUE", "phyto_FDis_ZRUE", 
                      "species_richness", "zoo_evenness", "zoo_FEve", "zoo_FDis",
                      "genus_richness:Temperature", "genus_richness:Salinity",
                      "genus_richness:TotalNitrogen_unfiltered",
                      "phyto_evenness:Temperature", "phyto_evenness:Salinity",
                      "phyto_evenness:TotalNitrogen_unfiltered",
                      "phyto_FEve_ZRUE:Temperature", "phyto_FEve_ZRUE:Salinity",
                      "phyto_FEve_ZRUE:TotalNitrogen_unfiltered",
                      "phyto_FDis_ZRUE:Temperature", "phyto_FDis_ZRUE:Salinity",
                      "phyto_FDis_ZRUE:TotalNitrogen_unfiltered",
                      "species_richness:Temperature", "species_richness:Salinity",
                      "species_richness:TotalNitrogen_unfiltered",
                      "zoo_evenness:Temperature", "zoo_evenness:Salinity",
                      "zoo_evenness:TotalNitrogen_unfiltered",
                      "zoo_FEve:Temperature", "zoo_FEve:Salinity",
                      "zoo_FEve:TotalNitrogen_unfiltered",
                      "zoo_FDis:Temperature", "zoo_FDis:Salinity",
                      "zoo_FDis:TotalNitrogen_unfiltered",
                      "Zooplankton_MaxDepth", "daylight_hours")),
  
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA/BRM_zoo_RUE_PI_261013_All.RDATA",
       name = "zoo_RUE_PI",
       resp = "log_RUE_zooplankton",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered",
                      "genus_richness", "phyto_evenness", "phyto_FEve_ZRUE", "phyto_FDis_ZRUE", 
                      "species_richness", "zoo_evenness", "zoo_FEve", "zoo_FDis",
                      "genus_richness:Temperature", "genus_richness:Salinity",
                      "genus_richness:TotalNitrogen_unfiltered",
                      "phyto_evenness:Temperature", "phyto_evenness:Salinity",
                      "phyto_evenness:TotalNitrogen_unfiltered",
                      "phyto_FEve_ZRUE:Temperature", "phyto_FEve_ZRUE:Salinity",
                      "phyto_FEve_ZRUE:TotalNitrogen_unfiltered",
                      "phyto_FDis_ZRUE:Temperature", "phyto_FDis_ZRUE:Salinity",
                      "phyto_FDis_ZRUE:TotalNitrogen_unfiltered",
                      "species_richness:Temperature", "species_richness:Salinity",
                      "species_richness:TotalNitrogen_unfiltered",
                      "zoo_evenness:Temperature", "zoo_evenness:Salinity",
                      "zoo_evenness:TotalNitrogen_unfiltered",
                      "zoo_FEve:Temperature", "zoo_FEve:Salinity",
                      "zoo_FEve:TotalNitrogen_unfiltered",
                      "zoo_FDis:Temperature", "zoo_FDis:Salinity",
                      "zoo_FDis:TotalNitrogen_unfiltered",
                      "Zooplankton_MaxDepth", "daylight_hours")),
  
  # ---------------------------------------------------------
  # ZOOPLANKTON RUE WITH CYANO DOM
  # ---------------------------------------------------------
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA/BRM_zoo_RUE_CD_L_261013_All.RDATA",
       name = "zoo_RUE_CD_L",
       resp = "log_RUE_zooplankton",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered",
                      "genus_richness", "phyto_evenness", "phyto_FEve_ZRUE", "phyto_FDis_ZRUE", "cyano_prop",
                      "species_richness", "zoo_evenness", "zoo_FEve", "zoo_FDis",
                      "Zooplankton_MaxDepth", "daylight_hours")),
  
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA/BRM_zoo_RUE_CD_P_261013_All.RDATA",
       name = "zoo_RUE_CD_P",
       resp = "log_RUE_zooplankton",
       predictors = c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered",
                      "genus_richness", "phyto_evenness", "phyto_FEve_ZRUE", "phyto_FDis_ZRUE", "cyano_prop",
                      "species_richness", "zoo_evenness", "zoo_FEve", "zoo_FDis",
                      "Zooplankton_MaxDepth", "daylight_hours")),
  
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA/BRM_zoo_RUE_CD_LI_261013_All.RDATA",
       name = "zoo_RUE_CD_LI",
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
                      "Zooplankton_MaxDepth", "daylight_hours")),
  
  list(file = "~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA/BRM_zoo_RUE_CD_PI_261013_All.RDATA",
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




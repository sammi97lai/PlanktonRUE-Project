####################### Calculate Edibility based off of phytoplankton ##################
library(dplyr)
library(readxl)

rm(list = ls())

setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/RDATA")
load('plankton_RUE_data6_260309.RDATA')

#load ediblity data
setwd("~/GitHub/PlanktonRUE-Project/Data/Phytoplankton")
Edibility_values <- read_excel("Sciecies list nutritional value or edibility_Sirpa.xlsx")
Edibility_values <- Edibility_values %>% select(species, Edibility)

setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/RDATA")
#calculate Edibility values ...................................................
# Extract phytoplankton columns
phyto <- plankton_RUE_data6[, 16:84]

# Match edibility values to columns
edibility_vec <- Edibility_values$Edibility[
  match(colnames(phyto), Edibility_values$species)
]

# Check for unmatched species
sum(is.na(edibility_vec))

# See which ones failed to match
setdiff(colnames(phyto), Edibility_values$species)

plankton_RUE_data7 <- plankton_RUE_data6
# ---- Weighted mean (abundance-weighted) ----
plankton_RUE_data7$mean_edibility_weighted <- apply(phyto, 1, function(x) {
  if (all(is.na(x)) || sum(x, na.rm = TRUE) == 0) return(NA)
  weighted.mean(edibility_vec, w = x, na.rm = TRUE)
})

# ---- Unweighted mean (presence/absence) ----
plankton_RUE_data7$mean_edibility_unweighted <- apply(phyto, 1, function(x) {
  present <- x > 0
  if (all(!present, na.rm = TRUE)) return(NA)
  mean(edibility_vec[present], na.rm = TRUE)
})

save(plankton_RUE_data7, file = 'plankton_RUE_data7_260408.RDATA')

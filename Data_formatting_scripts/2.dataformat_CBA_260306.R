############### Formatting carbon, biomass, and abundance dataset #############
library(dplyr)
library(readxl)
library(tidyr)
library(vegan)
library(openxlsx) # Save Excel file

setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/RDATA")
load('Phytoplankton_data4.RDATA')

########################### Carbon Dataset ####################################
carbon_data <- filtered_Phytoplankton_data %>%
  select(Site_Id, Sampling_Id, SettledSample_Id, Date, Time, Region, latitude, longitude, MinDepth, MaxDepth, Genus, 
         Temperature, Salinity, pH, Oxygen_soluble, TotalPhosphorus_unfiltered, TotalNitrogen_unfiltered, 
         `Chlorophyll-a`, Carbon_ugPerLiter)%>%
  pivot_wider(names_from = Genus, values_from = Carbon_ugPerLiter,values_fill = list(Carbon_ugPerLiter = 0))
carbon_data[] <- lapply(carbon_data, function(col) {
  col[is.nan(col)] <- NA
  col
})
#colSums(is.na(carbon_data)) #check NAs#
carbon_data[17:289] <- lapply(carbon_data[17:289], function(col) {
  col[is.na(col)] <- 0
  col
}) #change NA absence to zero#


#...................... remove low species presence .............................
#look at species presence#
#make presence absence data#
just_env<-carbon_data %>% ungroup() %>%
  select(c(Site_Id, Sampling_Id, SettledSample_Id, Date, Time, Region, latitude, longitude, MinDepth, MaxDepth, 
           Temperature, Salinity, pH, Oxygen_soluble, TotalPhosphorus_unfiltered, TotalNitrogen_unfiltered, 
           `Chlorophyll-a`)) #have dataset just of env factors
carbon_data_PA<-carbon_data %>% ungroup() %>%
  select(-c(Site_Id, Sampling_Id, SettledSample_Id, Date, Time, Region, latitude, longitude, MinDepth, MaxDepth, 
            Temperature, Salinity, pH, Oxygen_soluble, TotalPhosphorus_unfiltered, TotalNitrogen_unfiltered, 
            `Chlorophyll-a`)) #remove env factors
carbon_data_PA2 <- carbon_data_PA #create new dataset for PA
carbon_data_PA2[carbon_data_PA2>0] = 1

#calculate species presence#
species_percentage <- (colSums(carbon_data_PA2) / nrow(carbon_data_PA2)) * 100
print(species_percentage)
#remove species presence that is lower than 10%#
# filtered_phyto_data <- carbon_data_PA2[, species_percentage >= 10]
# print(filtered_phyto_data) #only 52 species remaining
#remove species presence that is lower than 5%#
filtered_phyto_data <- carbon_data_PA2[, species_percentage >= 5]
print(filtered_phyto_data) #only 69 species remaining

occurrence <- colSums(carbon_data_PA2)
summary(occurrence)
carbon_data_filtered <- carbon_data_PA2[, occurrence >= 3]

thresholds <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)

sapply(thresholds, function(t){
  sum(species_percentage >= t)
})

species_remaining <- sapply(thresholds, function(t){
  sum(species_percentage >= t)
})

#remove species with low presence#
# Remove species with presence < 5%
valid_species <- names(species_percentage[species_percentage >= 5])
filtered_carbon_data <- carbon_data_PA %>%
  select(any_of(valid_species))

#combine dataset#
filtered_carbon_env_data<-cbind(just_env,filtered_carbon_data)

save(filtered_carbon_env_data, file = 'filtered_carbon_env_data.RDATA')

load('filtered_carbon_env_data.RDATA')

########################### Biomass Dataset ###################################
biomass_data <- filtered_Phytoplankton_data %>%
  select(Site_Id, Sampling_Id, SettledSample_Id, Date, Time, Region, latitude, longitude, MinDepth, MaxDepth, Genus, 
         Temperature, Salinity, pH, Oxygen_soluble, TotalPhosphorus_unfiltered, TotalNitrogen_unfiltered, 
         `Chlorophyll-a`, Biomass_ugPerLiter)%>%
  pivot_wider(names_from = Genus, values_from = Biomass_ugPerLiter,values_fill = list(Biomass_ugPerLiter = 0))
biomass_data[] <- lapply(biomass_data, function(col) {
  col[is.nan(col)] <- NA
  col
})
#colSums(is.na(carbon_data)) #check NAs#
biomass_data[17:289] <- lapply(biomass_data[17:289], function(col) {
  col[is.na(col)] <- 0
  col
}) #change NA absence to zero#

#...................... remove low species presence .............................
#look at species presence#
#make presence absence data#
just_env<-biomass_data %>% ungroup() %>%
  select(c(Site_Id, Sampling_Id, SettledSample_Id, Date, Time, Region, latitude, longitude, MinDepth, MaxDepth, 
           Temperature, Salinity, pH, Oxygen_soluble, TotalPhosphorus_unfiltered, TotalNitrogen_unfiltered, 
           `Chlorophyll-a`)) #have dataset just of env factors

#....................... Total biomass per species .............................
# species_biomass: columns with biomass (same order as PA data)
species_cols <- 18:ncol(biomass_data)
species_biomass <- biomass_data[, species_cols]

# Sum biomass across all samples
species_total_biomass <- colSums(species_biomass, na.rm = TRUE)

biomass_retained <- sapply(thresholds, function(t){
  sum(species_total_biomass[species_percentage >= t])
})

#plot#
png("~/GitHub/PlanktonRUE-Project/Figures/occurrence_threshold.png", 
    width = 1400, height = 600, res = 150)
par(mfrow = c(1,2), mar = c(5,5,3,2))

plot(thresholds, species_remaining,
     type="b",
     xlab="Occurrence threshold (%)",
     ylab="Number of species retained")
abline(v = 5, col = "red", lty = 2, lwd = 2)
mtext("a", side = 3, line = 0.5, at = par("usr")[1], font = 2, cex = 1.5)

plot(thresholds, biomass_retained,
     type = "b",
     xlab = "Occurrence threshold (%)",
     ylab = "Total biomass of retained species")
abline(v = 5, col = "red", lty = 2, lwd = 2)
mtext("b", side = 3, line = 0.5, at = par("usr")[1], font = 2, cex = 1.5)

dev.off()

biomass_data_PA<-biomass_data %>% ungroup() %>%
  select(-c(Site_Id, Sampling_Id, SettledSample_Id, Date, Time, Region, latitude, longitude, MinDepth, MaxDepth, 
            Temperature, Salinity, pH, Oxygen_soluble, TotalPhosphorus_unfiltered, TotalNitrogen_unfiltered, 
            `Chlorophyll-a`)) #remove env factors
filtered_biomass_data <- biomass_data_PA %>%
  select(any_of(valid_species)) #filter out low presence#

filtered_biomass_env_data<-cbind(just_env,filtered_biomass_data) #combine dataset#

species_cols <- colnames(filtered_biomass_env_data)[18:min(289, ncol(filtered_biomass_env_data))]
biomass_matrix <- filtered_biomass_env_data[, species_cols]

# Check which rows sum to zero
zero_rows <- which(rowSums(biomass_matrix, na.rm = TRUE) == 0)
zero_rows #1555

#remove zero rows from dataset
# Keep only rows where at least one species has biomass
filtered_biomass_env_data_nonzero <- filtered_biomass_env_data[
  rowSums(filtered_biomass_env_data[, species_cols], na.rm = TRUE) > 0, 
]

# Optional: check that row 1555 is gone
zero_rows <- which(rowSums(filtered_biomass_env_data_nonzero[, species_cols], na.rm = TRUE) == 0)
zero_rows  # should return integer(0)

filtered_biomass_env_data<-filtered_biomass_env_data_nonzero

save(filtered_biomass_env_data, file = 'filtered_biomass_env_data.RDATA')
save(biomass_data, file = 'non_filtered_biomass_data.RDATA') #include rare species

########################### Abundance ###################################
abd_data <- filtered_Phytoplankton_data %>%
  select(Site_Id, Sampling_Id, SettledSample_Id, Date, Time, Region, latitude, longitude, MinDepth, MaxDepth, Genus, 
         Temperature, Salinity, pH, Oxygen_soluble, TotalPhosphorus_unfiltered, TotalNitrogen_unfiltered, 
         `Chlorophyll-a`, UnitsPerLiter)%>%
  pivot_wider(names_from = Genus, values_from = UnitsPerLiter,values_fill = list(UnitsPerLiter = 0))
abd_data[] <- lapply(abd_data, function(col) {
  col[is.nan(col)] <- NA
  col
})
#colSums(is.na(carbon_data)) #check NAs#
abd_data[17:289] <- lapply(abd_data[17:289], function(col) {
  col[is.na(col)] <- 0
  col
}) #change NA absence to zero#

#remove low presence species#
just_env<-abd_data %>% ungroup() %>%
  select(c(Site_Id, Sampling_Id, SettledSample_Id, Date, Time, Region, latitude, longitude, MinDepth, MaxDepth, 
           Temperature, Salinity, pH, Oxygen_soluble, TotalPhosphorus_unfiltered, TotalNitrogen_unfiltered, 
           `Chlorophyll-a`)) #have dataset just of env factors
abd_data_PA<-abd_data %>% ungroup() %>%
  select(-c(Site_Id, Sampling_Id, SettledSample_Id, Date, Time, Region, latitude, longitude, MinDepth, MaxDepth, 
            Temperature, Salinity, pH, Oxygen_soluble, TotalPhosphorus_unfiltered, TotalNitrogen_unfiltered, 
            `Chlorophyll-a`)) #remove env factors
filtered_abd_data <- abd_data_PA %>%
  select(any_of(valid_species)) #filter out low presence#
filtered_abd_env_data<-cbind(just_env,filtered_abd_data) #combine dataset#

save(filtered_abd_env_data, file = 'filtered_abd_env_data.RDATA')


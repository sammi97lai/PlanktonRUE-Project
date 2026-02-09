############### Formatting carbon, biomass, and abundance dataset #############
library(dplyr)
library(readxl)
library(tidyr)
library(vegan)
library(openxlsx) # Save Excel file

setwd("~/GitHub/Phytoplankton-project/phyto_RUE/dataformatting")
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
carbon_data[17:335] <- lapply(carbon_data[17:335], function(col) {
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
# print(filtered_phyto_data) #only 77 species remaining
#remove species presence that is lower than 5%#
filtered_phyto_data <- carbon_data_PA2[, species_percentage >= 5]
print(filtered_phyto_data) #only 103 species remaining

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
biomass_data[17:335] <- lapply(biomass_data[17:335], function(col) {
  col[is.na(col)] <- 0
  col
}) #change NA absence to zero#

#remove low presence species#
just_env<-biomass_data %>% ungroup() %>%
  select(c(Site_Id, Sampling_Id, SettledSample_Id, Date, Time, Region, latitude, longitude, MinDepth, MaxDepth, 
           Temperature, Salinity, pH, Oxygen_soluble, TotalPhosphorus_unfiltered, TotalNitrogen_unfiltered, 
           `Chlorophyll-a`)) #have dataset just of env factors
biomass_data_PA<-biomass_data %>% ungroup() %>%
  select(-c(Site_Id, Sampling_Id, SettledSample_Id, Date, Time, Region, latitude, longitude, MinDepth, MaxDepth, 
            Temperature, Salinity, pH, Oxygen_soluble, TotalPhosphorus_unfiltered, TotalNitrogen_unfiltered, 
            `Chlorophyll-a`)) #remove env factors
filtered_biomass_data <- biomass_data_PA %>%
  select(any_of(valid_species)) #filter out low presence#
filtered_biomass_env_data<-cbind(just_env,filtered_biomass_data) #combine dataset#

save(filtered_biomass_env_data, file = 'filtered_biomass_env_data.RDATA')

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
abd_data[17:335] <- lapply(abd_data[17:335], function(col) {
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


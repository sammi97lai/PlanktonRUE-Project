################################################################################
##################### Plankton data matching ###################################
################################################################################
library(dplyr)
library(readxl)
library(tidyr)
library(vegan)
library(openxlsx)
library(lubridate)
library(sf)
library(ggplot2)
library(maps)
library(viridis)
library(fuzzyjoin)
library(ggrepel)

rm(list = ls())

setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/RDATA")
load('RUE_biomass_com_260915.RDATA')
load('Zooplankton_REC_clean3_260307.RDATA')

setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/spatial_sensitivity_analysis/Smaller_radius/RDATA")

################################################################################
#clean Zooplankton dataset
Zooplankton_biomass_data <- Zooplankton_REC_clean3 %>%
  select(SampleEvent, Date, latitude, longitude, EnvironmentType, 
         Zooplankton_MinDepth, Zooplankton_MaxDepth, Species, WetWeight) %>%
  pivot_wider(
    names_from = Species,
    values_from = WetWeight,
    values_fill = 0
  ) # Fills missing combinations with 0

metadata_cols <- c("SampleEvent", "Date", "latitude", "longitude", "EnvironmentType", 
                   "Zooplankton_MinDepth", "Zooplankton_MaxDepth")

Zooplankton_biomass_data <- Zooplankton_biomass_data %>%
  rename_with(~ paste0(.x, "_z"), .cols = -all_of(metadata_cols))

#include environment type to RUE dataset
setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/RDATA")
load('Phytoplankton_data_format.RDATA')
Phytoplankton_data2<- Phytoplankton_data%>%
  select(Site_Id,`Coordinates_ EUREF_WGS84_lat`, `Coordinates_ EUREF_WGS84_long`, 
         EnvironmentalType, SurfacewaterType, WaterArea, WaterAreaNumber) %>%
  rename(latitude =`Coordinates_ EUREF_WGS84_lat`,
         longitude =`Coordinates_ EUREF_WGS84_long`,
         SiteID = Site_Id, EnvironmentType = EnvironmentalType, 
         Waterbody = WaterArea, WaterbodyCode = WaterAreaNumber)%>%
  distinct() %>%
  mutate(
    EnvironmentType = case_when(
      EnvironmentType %in% c("meri ulkos", "meri sisäs") ~ "Coastal",
      EnvironmentType == "avomeri" ~ "Open water",
      TRUE ~ EnvironmentType
    )
  )

RUE_with_env <- RUE_biomass_com %>%
  left_join(
    Phytoplankton_data2 %>%
      select(latitude, longitude, EnvironmentType, SurfacewaterType),
    by = c("latitude" = "latitude", "longitude" = "longitude")
  )

# 1. Prepare your environmental data
env_sf <- RUE_with_env %>%
  mutate(Date = as.Date(Date)) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs = 3857)

# 2. Prepare your zooplankton data
zoo_sf <- Zooplankton_biomass_data %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 3067) %>%
  st_transform(crs = 3857)

# 3. Spatial join:  SMALLER RADIUS  
# --- Conditional joins ---
## Coastal sites (5 km radius) #min distance between site with same habitat type
coastal_env <- env_sf %>% filter(EnvironmentType == "Coastal")
coastal_zoo <- zoo_sf %>% filter(EnvironmentType == "Coastal")
coastal_joined_env_to_zoo <- st_join(
  coastal_env,
  coastal_zoo,
  join = st_is_within_distance,
  dist = 5000,
  suffix = c("_env", "_zoo")
)
coastal_joined_zoo_to_env <- st_join(
  coastal_zoo,
  coastal_env,
  join = st_is_within_distance,
  dist = 5000,
  suffix = c("_zoo", "_env")
)

## Open water sites (30000 m radius) #min distance between each site
open_env <- env_sf %>% filter(EnvironmentType == "Open water")
open_zoo <- zoo_sf %>% filter(EnvironmentType == "Open water")
open_joined_env_to_zoo <- st_join(
  open_env,
  open_zoo,
  join = st_is_within_distance,
  dist = 30000,
  suffix = c("_env", "_zoo")
)
open_joined_zoo_to_env <- st_join(
  open_zoo,
  open_env,
  join = st_is_within_distance,
  dist = 30000,
  suffix = c("_zoo", "_env")
)

# Combine results
joined_all <- bind_rows(coastal_joined_env_to_zoo,
                        coastal_joined_zoo_to_env,
                        open_joined_env_to_zoo,
                        open_joined_zoo_to_env
                        ) %>%
  distinct()

# 4. Filter to get plus or minus 3 days date‐matches
merged <- joined_all %>%
  filter(abs(as.numeric(difftime(Date_env, Date_zoo, units = "days"))) <= 3)

# 5. Drop geometry to get plain data.frame and convert to degrees lat lon
merged_wgs <- st_transform(merged, crs = 4326)

final_df <- merged_wgs %>%
  mutate(
    longitude = st_coordinates(geometry)[, 1],
    latitude  = st_coordinates(geometry)[, 2]
  ) %>%
  st_drop_geometry() # 892 obs (930 obs manuscript final)
#less obs by splitting observations between coastal and open sea

final_df <- final_df %>%
  mutate(
    longitude = round(longitude, 5),
    latitude  = round(latitude, 5)
  )

#check number of sampling years ................................................
# Calculate sampling years count per site (lon, lat)
site_years <- final_df %>%
  mutate(year = format(Date_env, "%Y")) %>% #extract year from date
  group_by(longitude, latitude) %>%
  summarise(years_sampled = n_distinct(year), .groups = "drop") %>%
  mutate(years_cat = factor(years_sampled))  

site_years <- site_years %>%
  group_by(longitude, latitude) %>%
  mutate(site_id = cur_group_id()) %>%  # assigns a unique number per site
  ungroup()

#find number of sites
n_sites <- site_years %>%
  distinct(longitude, latitude)

n_sites

#calculate total Btot for zooplankton#........................................
filtered_plankton_data2<-final_df
filtered_plankton_data2$zoo_Btot <- rowSums(filtered_plankton_data2[, 112:162]) 

#RUE calculation
filtered_plankton_data2 <- filtered_plankton_data2 %>%
  mutate(RUE_zooplankton = zoo_Btot / Btot)
filtered_plankton_data2$log_RUE_zooplankton <- log(filtered_plankton_data2$RUE_zooplankton)
#check NAs
# sum(is.na(filtered_plankton_data2))
# filtered_plankton_data2[rowSums(is.na(filtered_plankton_data2)) > 0, ]
# colSums(is.na(filtered_plankton_data2))
# names(filtered_plankton_data2)[colSums(is.na(filtered_plankton_data2)) > 0]
# # only in oxygen soluble and surface water type

#clean up dataset
plankton_RUE_data <- filtered_plankton_data2 %>%
  select(-c(Date_zoo)) %>%
  rename(Date = Date_env)

setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/spatial_sensitivity_analysis/Smaller_radius/RDATA")
save(plankton_RUE_data, file = 'plankton_RUE_data_260526.RDATA')

#remove species with zero biomass #.............................................
load("plankton_RUE_data_260526.RDATA")
# Columns you want to always keep (metadata + Btot etc.)
keep_cols <- 1:110
extra_cols <- 162:166

# Define species columns (everything else except keep_cols and extra_cols)
species_cols <- setdiff(1:ncol(plankton_RUE_data), c(keep_cols, extra_cols))

# Keep only species with non-zero biomass
nonzero_species <- species_cols[colSums(plankton_RUE_data[, species_cols], na.rm = TRUE) > 0]

# Combine everything, keeping order: metadata, extra, non-zero species
plankton_RUE_data_filtered <- plankton_RUE_data[, c(keep_cols, extra_cols, nonzero_species)]

#identity species that were removed
removed_species <- species_cols[!(species_cols %in% nonzero_species)]
removed_species_names <- colnames(plankton_RUE_data)[removed_species]
removed_species_names

plankton_RUE_data2 <- plankton_RUE_data_filtered
setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/spatial_sensitivity_analysis/Smaller_radius/RDATA")
save(plankton_RUE_data2, file = 'plankton_RUE_data2_260526.RDATA')

### plot data ############......................................................
###plot points##

ggplot() +
  #geom_sf(data = biomass_sf, aes(color = "Phytoplankton biomass"), size = 2, alpha = 0.7) +  
  geom_sf(data = env_sf, aes(color = "Phytoplankton"), size = 3, shape= 15, alpha = 0.7) +
  geom_sf(data = zoo_sf, aes(color = "Zooplankton"), size = 2, alpha = 0.7) +
  geom_sf(data = joined_all, aes(color = "joined"), size = 1, alpha = 0.7) +
  scale_color_manual(values = c("Phytoplankton" = "green", "Zooplankton" = "blue",
                                "Phytoplankton biomass" = "green", "joined" = "black")) +
  theme_minimal() +
  labs(color = "Dataset", title = "Sampling points from datasets")

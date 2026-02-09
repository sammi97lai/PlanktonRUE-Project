#########################################################################################
## find the lowest distance between coastal sites that have similar characteristics #####
#########################################################################################
library(dplyr)
library(readxl)
library(tidyr)
library(vegan)
library(openxlsx)
library(lubridate)
library(sf)
library(ggplot2)

rm(list = ls())

setwd("~/GitHub/Phytoplankton-project/phyto_RUE/dataformatting")
#..............................................................................
#reformat zooplankton data
setwd("~/GitHub/Phytoplankton-project/Data/Zooplankton/Zooplankton_data_Syke")
Zooplankton_data <- read.csv("ZooplanktonMarine_V3EXT_11Sep2025.csv")
setwd("~/GitHub/Phytoplankton-project/phyto_RUE/dataformatting")

colnames(Zooplankton_data)
Zooplankton_data2 <- Zooplankton_data %>%
  filter(Year >= 2000 & Year <= 2024) %>%
  mutate(Date = make_date(Year, Month, Day)) %>%
  filter(!is.na(Latitude) & !is.na(Longitude)) %>%
  filter(EnvironmentType %in% c("Coast, outer archipelago", "Coast, inner archipelago")) %>%
  select(SiteID, Latitude, Longitude, EnvironmentType, Waterbody, WaterbodyCode) %>%
  distinct()
zoo_sf <- st_as_sf(Zooplankton_data2,
                   coords = c("Longitude", "Latitude"),
                   crs = 3067)   # Finnish CRS in meters
#..............................................................................
#reformat phytoplankton
load('Phytoplankton_data_format.RDATA')
Phytoplankton_data2<- Phytoplankton_data%>%
select(Site_Id,`Coordinates_ EUREF_WGS84_lat`, `Coordinates_ EUREF_WGS84_long`, 
       EnvironmentalType, SurfacewaterType, WaterArea, WaterAreaNumber) %>%
  rename(Latitude =`Coordinates_ EUREF_WGS84_lat`,
         Longitude =`Coordinates_ EUREF_WGS84_long`,
         SiteID = Site_Id, EnvironmentType = EnvironmentalType, 
         Waterbody = WaterArea, WaterbodyCode = WaterAreaNumber)%>%
  filter(EnvironmentType %in% c("meri ulkos", "meri sisäs")) %>%
  distinct() %>%
  mutate(
    EnvironmentType = recode(
      EnvironmentType,
      "meri ulkos" = "Coast, outer archipelago",
      "meri sisäs" = "Coast, inner archipelago"
    )
  )

phyto_sf <- st_as_sf(Phytoplankton_data2,
                     coords = c("Longitude", "Latitude"),
                     crs = 4326) %>% 
  st_transform(3067)   # ETRS-TM35FIN (EPSG:3067, meters)


#----------------------------------------------------------------------------------
#find smallest distance#

# Find nearest phytoplankton site for each zooplankton site
nearest_idx <- st_nearest_feature(zoo_sf, phyto_sf)

# Add SurfacewaterType (and optionally other attributes) from phyto
zoo_with_type <- zoo_sf %>%
  mutate(
    SurfacewaterType = phyto_sf$SurfacewaterType[nearest_idx],
    Phyto_SiteID     = phyto_sf$SiteID[nearest_idx],
    Distance_m       = st_distance(zoo_sf, phyto_sf[nearest_idx, ], by_element = TRUE) %>% as.numeric()
  )

# Pairwise distance matrix
dist_matrix <- st_distance(zoo_with_type, phyto_sf)

# Convert to data frame of unique pairs
pairs <- expand.grid(
  zoo = 1:nrow(zoo_with_type),
  phyto = 1:nrow(phyto_sf)
) %>%
  mutate(
    Env_match  = zoo_with_type$EnvironmentType[zoo] == phyto_sf$EnvironmentType[phyto],
    Surf_match = zoo_with_type$SurfacewaterType[zoo] == phyto_sf$SurfacewaterType[phyto],
    dist_m     = as.numeric(dist_matrix[cbind(zoo, phyto)])
  ) %>%
  filter(Env_match, Surf_match)   # keep only matching type pairs

min_radius <- min(pairs$dist_m, na.rm = TRUE)
min_radius
# 3.003839e-05

max_radius <- max(pairs$dist_m, na.rm = TRUE)
max_radius
# 243178.4

avg_radius <- mean(pairs$dist_m, na.rm = TRUE)
avg_radius
# 67149.8

#..............................................................................
#find smallest distance based on only phytoplankton habitat type
# Pairwise distance matrix
dist_matrix <- st_distance(phyto_sf)

# Convert to data frame of unique pairs
pairs <- expand.grid(i = 1:nrow(phyto_sf), j = 1:nrow(phyto_sf)) %>%
  filter(i < j) %>%
  mutate(dist = as.numeric(dist_matrix[cbind(i, j)]))

# Add attributes
pairs <- pairs %>%
  mutate(Env_i = phyto_sf$EnvironmentType[i],
         Env_j = phyto_sf$EnvironmentType[j],
         Water_i = phyto_sf$SurfacewaterType[i],
         Water_j = phyto_sf$SurfacewaterType[j]) %>%
  filter(Env_i == Env_j, Water_i == Water_j)

min_distance <- min(pairs$dist, na.rm = TRUE)
min_distance 
#3.734646

max_distance <- max(pairs$dist, na.rm = TRUE)
max_distance 
#273809.3

avg_distance <- mean(pairs$dist, na.rm = TRUE)
avg_distance
#64091.78

#.................................................................................

##########################################################################################
## find the lowest distance between open sea sites that have similar characteristics #####
##########################################################################################

Zooplankton_data2 <- Zooplankton_data %>%
  filter(Year >= 2000 & Year <= 2024) %>%
  mutate(Date = make_date(Year, Month, Day)) %>%
  filter(!is.na(Latitude) & !is.na(Longitude)) %>%
  filter(EnvironmentType %in% c("Open sea")) %>%
  select(SiteID, Latitude, Longitude, EnvironmentType, Waterbody, WaterbodyCode) %>%
  distinct()
zoo_sf <- st_as_sf(Zooplankton_data2,
                   coords = c("Longitude", "Latitude"),
                   crs = 3067)   # Finnish CRS in meters
#..............................................................................
#reformat phytoplankton
Phytoplankton_data2<- Phytoplankton_data%>%
  select(Site_Id,`Coordinates_ EUREF_WGS84_lat`, `Coordinates_ EUREF_WGS84_long`, 
         EnvironmentalType, SurfacewaterType, WaterArea, WaterAreaNumber) %>%
  rename(Latitude =`Coordinates_ EUREF_WGS84_lat`,
         Longitude =`Coordinates_ EUREF_WGS84_long`,
         SiteID = Site_Id, EnvironmentType = EnvironmentalType, 
         Waterbody = WaterArea, WaterbodyCode = WaterAreaNumber)%>%
  filter(EnvironmentType %in% c("avomeri")) %>%
  distinct() %>%
  mutate(
    EnvironmentType = recode(
      EnvironmentType,
      "avomeri" = "Open sea"
    )
  )

phyto_sf <- st_as_sf(Phytoplankton_data2,
                     coords = c("Longitude", "Latitude"),
                     crs = 4326) %>% 
  st_transform(3067)   # ETRS-TM35FIN (EPSG:3067, meters)


#----------------------------------------------------------------------------------
#find smallest distance#

# Find nearest phytoplankton site for each zooplankton site
nearest_idx <- st_nearest_feature(zoo_sf, phyto_sf)

# Add SurfacewaterType (and optionally other attributes) from phyto
zoo_with_type <- zoo_sf %>%
  mutate(
    SurfacewaterType = phyto_sf$SurfacewaterType[nearest_idx],
    Phyto_SiteID     = phyto_sf$SiteID[nearest_idx],
    Distance_m       = st_distance(zoo_sf, phyto_sf[nearest_idx, ], by_element = TRUE) %>% as.numeric()
  )

# Pairwise distance matrix
dist_matrix <- st_distance(zoo_with_type, phyto_sf)

# Convert to data frame of unique pairs
pairs <- expand.grid(
  zoo = 1:nrow(zoo_with_type),
  phyto = 1:nrow(phyto_sf)
) %>%
  mutate(
    Env_match  = zoo_with_type$EnvironmentType[zoo] == phyto_sf$EnvironmentType[phyto],
    dist_m     = as.numeric(dist_matrix[cbind(zoo, phyto)])
  ) %>%
  filter(Env_match)   # keep only matching type pairs

min_radius <- min(pairs$dist_m, na.rm = TRUE)
min_radius
# 0.0001026645

max_radius <- max(pairs$dist_m, na.rm = TRUE)
max_radius
# 1202978

avg_radius <- mean(pairs$dist_m, na.rm = TRUE)
avg_radius
# 357246.4

#..............................................................................
#find smallest distance based on only phytoplankton habitat type
# Pairwise distance matrix
dist_matrix <- st_distance(phyto_sf)

# Convert to data frame of unique pairs
pairs <- expand.grid(i = 1:nrow(phyto_sf), j = 1:nrow(phyto_sf)) %>%
  filter(i < j) %>%
  mutate(dist = as.numeric(dist_matrix[cbind(i, j)]))

# Add attributes
pairs <- pairs %>%
  mutate(
    Env_i = phyto_sf$EnvironmentType[i],
    Env_j = phyto_sf$EnvironmentType[j]
  ) %>%
  filter(Env_i == Env_j)

min_distance <- min(pairs$dist, na.rm = TRUE)
min_distance 
#592.3784

max_distance <- max(pairs$dist, na.rm = TRUE)
max_distance 
#972703.8

avg_distance <- mean(pairs$dist, na.rm = TRUE)
avg_distance
#359843.9
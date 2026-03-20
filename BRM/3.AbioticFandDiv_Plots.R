################################################################################
##################### Plots of Abiotic Factors #################################
################################################################################

library(dplyr)
library(tidyr)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(patchwork)
library(viridis)
library(ggpubr)
library(ggrepel)

setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/RDATA")
load('plankton_RUE_data_scaled_260309.RDATA')
load('plankton_RUE_data6_260309.RDATA')

setwd("~/GitHub/PlanktonRUE-Project/Figures")

#...............................................................................
# Convert your data to sf object
plankton_sf <- sf::st_as_sf(plankton_RUE_data6, 
                            coords = c("longitude", "latitude"),
                            crs = 4326)  # WGS84

# Get a simple map of the region (e.g., Baltic Sea area)
world <- ne_countries(scale = "medium", returnclass = "sf")

xlim <- c(17, 26)
ylim <- c(57, 66)


#.............................................................................
#Site map
coords <- sf::st_coordinates(plankton_sf)

plankton_sf$X <- coords[, 1]
plankton_sf$Y <- coords[, 2]

p1 <- ggplot() +
  geom_sf(data = world, fill = "gray95", color = "gray50") +
  geom_sf(data = plankton_sf,
          aes(color = Site_Id),
          size = 1.5, alpha = 0.7) +
  geom_text(
    data = plankton_sf,
    aes(x = X, y = Y + 0.1, label = Site_Id),
    size = 3
  ) +
  coord_sf(xlim = xlim, ylim = ylim) +
  theme_minimal() +
  labs(x = NULL, y = NULL) +
  ggtitle("Site ID")
p1
# p1 <- ggplot() +
#   geom_sf(data = world, fill = "gray95", color = "gray20") +
#   geom_sf(data = plankton_sf,
#           color = "black",
#           size = 1.5, alpha = 0.7) +
#   coord_sf(xlim = xlim, ylim = ylim) +
#   theme_minimal()

ggsave("Plankton_site_map.png", 
       plot = p1, width = 4.5, height = 7, dpi = 300)
#..................................
# #for poster
# world <- ne_countries(scale = "medium", returnclass = "sf")
# 
# # Reproject
# world_m <- st_transform(world, 3035)
# plankton_m <- st_transform(plankton_sf, 3035)
# 
# # Get bbox of sites
# bb <- st_bbox(plankton_m)
# 
# # Compute width & height
# width  <- bb["xmax"] - bb["xmin"]
# height <- bb["ymax"] - bb["ymin"]
# 
# # Take the larger one to make it square
# side <- max(width, height)
# 
# # Center of bbox
# cx <- (bb["xmax"] + bb["xmin"]) / 2
# cy <- (bb["ymax"] + bb["ymin"]) / 2
# 
# # Add some padding (e.g., 20%)
# pad <- 1.2
# half <- side * pad / 2
# 
# square_bb <- structure(
#   c(xmin = cx - half,
#     ymin = cy - half,
#     xmax = cx + half,
#     ymax = cy + half),
#   class = "bbox",
#   crs = st_crs(plankton_m)
# )
# 
# world_sq <- st_crop(world_m, square_bb)
# 
# p1 <- ggplot() +
#   geom_sf(data = world_sq, fill = "green4", color = "green4") +
#   geom_sf(data = plankton_m, color = "black", size = 3, alpha = 0.7) +
#   coord_sf(expand = FALSE) +
#   theme_minimal() +
#   theme(
#     panel.background = element_rect(fill = "lightblue", color = NA),
#     plot.background  = element_rect(fill = "lightblue", color = NA)
#   )
# 
# ggsave("Poster_plankton_site_map.png",
#        plot = p1, width = 5, height = 5, dpi = 300, bg = "lightblue")


#.................................................................................
plankton_mean <- plankton_RUE_data6 %>%
  group_by(longitude, latitude) %>%
  summarise(
    Temperature = mean(Temperature, na.rm = TRUE),
    Salinity = mean(Salinity, na.rm = TRUE),
    pH = mean(pH, na.rm = TRUE),
    TotalPhosphorus_unfiltered = mean(TotalPhosphorus_unfiltered, na.rm = TRUE),
    TotalNitrogen_unfiltered = mean(TotalNitrogen_unfiltered, na.rm = TRUE),
    .groups = "drop"
  )

plankton_sf <- st_as_sf(plankton_mean,
                        coords = c("longitude", "latitude"),
                        crs = 4326)

# Plot abiotic factors#
# Temperature map
p1 <- ggplot() +
  geom_sf(data = world, fill = "gray95", color = "gray20") +
  geom_sf(data = plankton_sf, aes(color = Temperature), size = 1.5, alpha = 0.7) +
  scale_color_viridis_c(option = "plasma") +
  coord_sf(xlim = xlim, ylim = ylim) +
  theme_minimal() +
  ggtitle("Temperature") +
  theme(legend.position = "right") +
  theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 5))

# Salinity map
p2 <- ggplot() +
  geom_sf(data = world, fill = "gray95", color = "gray20") +
  geom_sf(data = plankton_sf, aes(color = Salinity), size = 1.5, alpha = 0.7) +
  scale_color_viridis_c(option = "cividis") +
  coord_sf(xlim = xlim, ylim = ylim) +
  theme_minimal() +
  ggtitle("Salinity") +
  theme(legend.position = "right") +
  theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 5))

# pH map
p3 <- ggplot() +
  geom_sf(data = world, fill = "gray95", color = "gray20") +
  geom_sf(data = plankton_sf, aes(color = pH), size = 1.5, alpha = 0.7) +
  scale_color_viridis_c(option = "mako") +
  coord_sf(xlim = xlim, ylim = ylim) +
  theme_minimal() +
  ggtitle("pH") +
  theme(legend.position = "right") +
  theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 5))

# Phosphorus concentration map
p4 <- ggplot() +
  geom_sf(data = world, fill = "gray95", color = "gray20") +
  geom_sf(data = plankton_sf, aes(color = TotalPhosphorus_unfiltered), size = 1.5, alpha = 0.7) +
  scale_color_viridis_c(option = "rocket", name = "Total Phosphorus (µg/L)") +
  coord_sf(xlim = xlim, ylim = ylim) +
  theme_minimal() +
  ggtitle("Phosphorus concetrations") +
  theme(legend.position = "right") +
  theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 5))

# Nitrogen concentration map
p5 <- ggplot() +
  geom_sf(data = world, fill = "gray95", color = "gray20") +
  geom_sf(data = plankton_sf, aes(color = TotalNitrogen_unfiltered), size = 1.5, alpha = 0.7) +
  scale_color_viridis_c(option = "rocket", name = "Total Nitrogen (µg/L)") +
  coord_sf(xlim = xlim, ylim = ylim) +
  theme_minimal() +
  ggtitle("Nitrogen concentrations") +
  theme(legend.position = "right") +
  theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 5))

# N:P ratio map
# p6 <- ggplot() +
#   geom_sf(data = world, fill = "gray95", color = "gray20") +
#   geom_sf(data = plankton_sf, aes(color = N_P_ratio), size = 1.5, alpha = 0.7) +
#   scale_color_viridis_c(option = "rocket", name = "N:P ratio") +
#   coord_sf(xlim = xlim, ylim = ylim) +
#   theme_minimal() +
#   ggtitle("N:P ratio") +
#   theme(legend.position = "right") +
#   theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 5))


# Combine all four maps
combined_plot <- ggarrange(
  p1, p2, p3, 
  p4, p5,
  ncol = 3, nrow = 2,
  align = "hv"   # aligns horizontally and vertically
)

ggsave("Plankton_env_maps.png", 
       plot = combined_plot, width = 10, height = 6, dpi = 300)

#...............................................................................
#Plot Diversity Metrics and Cyanobacteria Dominance
plankton_sf <- sf::st_as_sf(plankton_RUE_data6, 
                            coords = c("longitude", "latitude"),
                            crs = 4326)  # WGS84

# Get a simple map of the region (e.g., Baltic Sea area)
world <- ne_countries(scale = "medium", returnclass = "sf")

xlim <- c(17, 26)
ylim <- c(57, 66)

coords <- sf::st_coordinates(plankton_sf)

plankton_sf$X <- coords[, 1]
plankton_sf$Y <- coords[, 2]
# Phytoplankton Genus Richness map
p1 <- ggplot() +
  geom_sf(data = world, fill = "gray95", color = "gray50") +
  geom_sf(data = plankton_sf, aes(color = genus_richness), size = 1.5, alpha = 0.7) +
  scale_color_gradient(low = "white", high = "darkgreen") +
  coord_sf(xlim = xlim, ylim = ylim) +
  theme_minimal() +
  ggtitle("Phyto Genus Richness") +
  theme(legend.position = "right") +
  theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 5))

# Phytoplankton FEve map
p2 <- ggplot() +
  geom_sf(data = world, fill = "gray95", color = "gray50") +
  geom_sf(data = plankton_sf, aes(color = FEve), size = 1.5, alpha = 0.7) +
  scale_color_gradient(low = "white", high = "darkgreen") +
  coord_sf(xlim = xlim, ylim = ylim) +
  theme_minimal() +
  ggtitle("Phyto Functional Evenness") +
  theme(legend.position = "right") +
  theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 5))

# Phytoplankton FDis map
p3 <- ggplot() +
  geom_sf(data = world, fill = "gray95", color = "gray50") +
  geom_sf(data = plankton_sf, aes(color = FDis), size = 1.5, alpha = 0.7) +
  scale_color_gradient(low = "white", high = "darkgreen") +
  coord_sf(xlim = xlim, ylim = ylim) +
  theme_minimal() +
  ggtitle("Phyto Functional Dispersion") +
  theme(legend.position = "right") +
  theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 5))

# Cyano Dom map
p4 <- ggplot() +
  geom_sf(data = world, fill = "gray95", color = "gray50") +
  geom_sf(data = plankton_sf, aes(color = cyano_prop), size = 1.5, alpha = 0.7) +
  scale_color_gradient(low = "white", high = "darkgreen") +
  coord_sf(xlim = xlim, ylim = ylim) +
  theme_minimal() +
  ggtitle("Proportion of Cyanobateria") +
  theme(legend.position = "right") +
  theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
        panel.border = element_rect(color = "green", fill = NA, linewidth = 1)
        ) 

# Zooplankton Species Richness map
p5 <- ggplot() +
  geom_sf(data = world, fill = "gray95", color = "gray50") +
  geom_sf(data = plankton_sf, aes(color = species_richness), size = 1.5, alpha = 0.7) +
  scale_color_gradient(low = "white", high = "purple") +
  coord_sf(xlim = xlim, ylim = ylim) +
  theme_minimal() +
  ggtitle("Zoo Species Richness") +
  theme(legend.position = "right") +
  theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 5))

# Zooplankton FEve map
p6 <- ggplot() +
  geom_sf(data = world, fill = "gray95", color = "gray50") +
  geom_sf(data = plankton_sf, aes(color = zoo_FEve), size = 1.5, alpha = 0.7) +
  scale_color_gradient(low = "white", high = "purple") +
  coord_sf(xlim = xlim, ylim = ylim) +
  theme_minimal() +
  ggtitle("Zoo Functional Evenness") +
  theme(legend.position = "right") +
  theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 5))

# Zooplankton FDis map
p7 <- ggplot() +
  geom_sf(data = world, fill = "gray95", color = "gray50") +
  geom_sf(data = plankton_sf, aes(color = zoo_FDis), size = 1.5, alpha = 0.7) +
  scale_color_gradient(low = "white", high = "purple") +
  coord_sf(xlim = xlim, ylim = ylim) +
  theme_minimal() +
  ggtitle("Zoo Functional Dispersion") +
  theme(legend.position = "right") +
  theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 5))

#combine maps
combined_plot <- (p1 | p2 | p3 | p4) / 
  (p5 | p6 | p7) + 
  plot_layout(guides = "keep") & # keeps separate legends 
  theme(legend.position = "right")

ggsave("Plankton_div_maps.png", 
       plot = combined_plot, width = 15, height = 8, dpi = 300)

####plot RUE ################################################################
#average RUE per site per year
avg_RUE <- plankton_sf %>%
  mutate(Year=as.numeric(format(Date, "%Y"))) %>%
  group_by(geometry,Year)%>%
  summarise(avg_RUE_B_P = mean(log_RUE_B_P, na.rm = TRUE),
            avg_RUE_B_N = mean(log_RUE_B_N, na.rm = TRUE),
            avg_zoo_RUE = mean(log_RUE_zooplankton, na.rm = TRUE),
            .groups = "drop")

# Phyto RUE P
RUE_P <- ggplot() +
  geom_sf(data = world, fill = "gray95", color = "gray50") +
  geom_sf(data = avg_RUE, aes(color = avg_RUE_B_P), size = 1.5, alpha = 0.7) +
  scale_color_gradient(low = "white", high = "darkgreen") +
  coord_sf(xlim = xlim, ylim = ylim) +
  theme_minimal() +
  ggtitle("Phyto RUE P across years") +
  theme(legend.position = "right",
        plot.margin = margin(t = 5, r = 5, b = 5, l = 5)) +
  facet_wrap(~Year, ncol = 10) 
RUE_P

# Phyto RUE N
RUE_N <- ggplot() +
  geom_sf(data = world, fill = "gray95", color = "gray50") +
  geom_sf(data = avg_RUE, aes(color = avg_RUE_B_N), size = 1.5, alpha = 0.7) +
  scale_color_gradient(low = "white", high = "darkgreen") +
  coord_sf(xlim = xlim, ylim = ylim) +
  theme_minimal() +
  ggtitle("Phyto RUE N across years") +
  theme(legend.position = "right",
        plot.margin = margin(t = 5, r = 5, b = 5, l = 5)) +
  facet_wrap(~Year, ncol = 10) 
RUE_N

# zooplankton RUE
zoo_RUE <- ggplot() +
  geom_sf(data = world, fill = "gray95", color = "gray50") +
  geom_sf(data = avg_RUE, aes(color = avg_zoo_RUE), size = 1.5, alpha = 0.7) +
  scale_color_gradient(low = "white", high = "darkgreen") +
  coord_sf(xlim = xlim, ylim = ylim) +
  theme_minimal() +
  ggtitle("Zooplankton RUE across years") +
  theme(legend.position = "right",
        plot.margin = margin(t = 5, r = 5, b = 5, l = 5)) +
  facet_wrap(~Year, ncol = 10) 
zoo_RUE

combined_plot <- ggarrange(
  RUE_P, RUE_N, zoo_RUE, 
  ncol = 1, nrow = 3,
  align = "hv"   # aligns horizontally and vertically
)

ggsave("Plankton_RUE_maps.png", 
       plot = combined_plot, width = 10, height = 18, dpi = 300)

# Plot sampling frequency for each site over DOY ###############################
site_order <- plankton_RUE_data6 %>%
  group_by(Site_Id) %>%
  summarise(latitude = mean(latitude, na.rm = TRUE)) %>%
  arrange(desc(latitude))   # north → south

site_order

site_freq <- plankton_RUE_data6 %>%
  mutate(Site_Id = factor(Site_Id, levels = site_order$Site_Id)) %>%
  count(Site_Id, DOY) %>%
  ggplot(aes(x = DOY, y = n)) +
  geom_col(width = 1) +
  facet_wrap(~ Site_Id) +
  labs(
    x = "Day of Year (DOY)",
    y = "Sampling frequency"
  ) +
  theme_bw()

ggsave("Site_freq.png", 
       plot = site_freq, width = 10, height = 6, dpi = 300)

#all sites
site_freq2 <- plankton_RUE_data6 %>%
  count(DOY) %>%
  ggplot(aes(x = DOY, y = n)) +
  geom_col(width = 1) +
  labs(
    x = "Day of Year (DOY)",
    y = "Sampling frequency (all sites)"
  ) +
  theme_bw()

ggsave("Site_freq2.png", 
       plot = site_freq2, width = 10, height = 6, dpi = 300)

# Count samples per site per year
plankton_RUE_data6 %>%
  group_by(Site_Id, Year) %>%
  summarise(n_samples = n(), .groups = "drop") %>%
  mutate(Site_Id = factor(Site_Id, levels = site_order$Site_Id)) %>%
  ggplot(aes(x = Year, y = n_samples, fill = Site_Id)) +
  geom_col() +
  facet_wrap(~ Site_Id) +
  labs(
    x = "Year",
    y = "Sampling frequency",
    title = "Sampling frequency per site across years"
  ) +
  theme_bw() +
  theme(legend.position = "none")

#Sampling over Year across all sites
plankton_RUE_data6 %>%
  group_by(Year) %>%
  summarise(n_samples = n(), .groups = "drop") %>%
  ggplot(aes(x = Year, y = n_samples)) +
  geom_col(fill = "steelblue") +
  labs(
    x = "Year",
    y = "Total sampling frequency (all sites)",
    title = "Overall sampling frequency across all sites per year"
  ) +
  theme_bw()

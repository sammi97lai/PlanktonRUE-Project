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
library(scales)

setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/RDATA")
load('plankton_RUE_data_scaled_260408.RDATA')
load('plankton_RUE_data7_260309.RDATA')

setwd("~/GitHub/PlanktonRUE-Project/Figures")

#...............................................................................
# Convert your data to sf object
plankton_sf <- sf::st_as_sf(plankton_RUE_data7, 
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
plankton_mean <- plankton_RUE_data7 %>%
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
p1 <- p1 + theme(plot.title = element_blank())
p2 <- p2 + theme(plot.title = element_blank())
p3 <- p3 + theme(plot.title = element_blank())
p5 <- p5 + theme(plot.title = element_blank())

combined_plot <- ggarrange(
  p1, p2, 
  p3, p5, 
  ncol = 2, nrow = 2,
  align = "hv",
  labels = c("a", "b", "c", "d"),
  label.x = 0.02,   # move slightly right
  label.y = 0.98,   # move slightly down
  hjust = 0, vjust = 1,
  font.label = list(size = 14, face = "bold")
)
combined_plot

ggsave("Plankton_env_maps.png", 
       plot = combined_plot, width = 9, height = 8, dpi = 300)

#...............................................................................
#Plot Diversity Metrics and Cyanobacteria Dominance
plankton_sf <- sf::st_as_sf(plankton_RUE_data7, 
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
  align = "hv"   
)

ggsave("Plankton_RUE_maps.png", 
       plot = combined_plot, width = 10, height = 18, dpi = 300)

#average RUE across DOY ..........................................................
avg_RUE_DOY <- plankton_sf %>%
  group_by(DOY)%>%
  summarise(avg_RUE_B_P = mean(log_RUE_B_P, na.rm = TRUE),
            avg_RUE_B_N = mean(log_RUE_B_N, na.rm = TRUE),
            avg_zoo_RUE = mean(log_RUE_zooplankton, na.rm = TRUE),
            .groups = "drop")

# Phyto RUE P
RUE_P_DOY <- ggplot(avg_RUE_DOY, aes(x = DOY, y = avg_RUE_B_P)) +
  geom_col(fill = "darkgreen") +
  labs(x = "Day of Year (DOY)", y = "Avg Phyto RUE P") +
  theme_minimal()
RUE_P_DOY
# Phyto RUE N
RUE_N_DOY <- ggplot(avg_RUE_DOY, aes(x = DOY, y = avg_RUE_B_N)) +
  geom_col(fill = "forestgreen") +
  labs(x = "Day of Year (DOY)", y = "Avg Phyto RUE N") +
  theme_minimal()
RUE_N_DOY

# zooplankton RUE
zoo_RUE_DOY <- ggplot(avg_RUE_DOY, aes(x = DOY, y = avg_zoo_RUE)) +
  geom_col(fill = "seagreen") +
  labs(x = "Day of Year (DOY)", y = "Avg Zooplankton RUE") +
  theme_minimal()
zoo_RUE_DOY

combined_bar <- ggarrange(
  RUE_P_bar, RUE_N_bar, zoo_RUE_bar,
  ncol = 1, nrow = 3,
  align = "v"
)

ggsave("Plankton_RUE_DOY_barplots.png", 
       plot = combined_bar, width = 12, height = 10, dpi = 300)

# Plot sampling frequency for each site over DOY ###############################
site_order <- plankton_RUE_data7 %>%
  group_by(Site_Id) %>%
  summarise(latitude = mean(latitude, na.rm = TRUE)) %>%
  arrange(desc(latitude))   # north → south

site_order

site_freq <- plankton_RUE_data7 %>%
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
# Precompute counts to avoid repeated calculations
doy_counts <- plankton_RUE_data7 %>% count(DOY)

# Max frequency for positioning labels
max_n <- max(doy_counts$n)

site_freq2 <- doy_counts %>%
  ggplot(aes(x = DOY, y = n)) +
  geom_col(width = 1, fill = "steelblue") +
  # vertical lines for July 1 and August 31
  geom_vline(xintercept = c(182, 243), linetype = "dashed", color = "red", linewidth = 0.8) +
  # labels
  annotate("text", x = 182, y = max_n * 1.05,
           label = "July 1", color = "red", hjust = 0.5) +
  annotate("text", x = 243, y = max_n * 1.05,
           label = "Aug 31", color = "red", hjust = 0.5) +
  labs(
    x = "Day of Year (DOY)",
    y = "Sampling frequency (all sites)"
  ) +
  theme_bw() +
  coord_cartesian(ylim = c(0, max_n * 1.2)) 
site_freq2

ggsave("Site_freq2.png", 
       plot = site_freq2, width = 10, height = 6, dpi = 300)

# Count samples per site per year
plankton_RUE_data7 %>%
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
plankton_RUE_data7 %>%
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


#DOY vs daylight hours and phytoplankton RUE plot ................................
scale_factor <- max(plankton_RUE_data7$RUE_B_P, na.rm = TRUE) /
  max(plankton_RUE_data7$daylight_hours, na.rm = TRUE)

p1 <- ggplot(plankton_RUE_data7, aes(x = DOY)) +
  geom_point(aes(y = daylight_hours, color = "Daylight hours")) +
  geom_point(aes(y = RUE_B_P / scale_factor, color = "RUE")) +
  scale_y_continuous(
    name = "Daylight hours",
    sec.axis = sec_axis(~ . * scale_factor, name = "RUE")
  ) +
  geom_smooth(aes(y = daylight_hours, color = "Daylight hours")) +
  geom_smooth(aes(y = RUE_B_P / scale_factor, color = "RUE")) +
  labs(x = "Day of Year", color = NULL) +
  theme_minimal()

scale_factor <- max(plankton_RUE_data7$genus_richness, na.rm = TRUE) /
  max(plankton_RUE_data7$daylight_hours, na.rm = TRUE)

p2 <- ggplot(plankton_RUE_data7, aes(x = DOY)) +
  geom_point(aes(y = daylight_hours, color = "Daylight hours")) +
  geom_point(aes(y = genus_richness / scale_factor, color = "genus richness")) +
  scale_y_continuous(
    name = "Daylight hours",
    sec.axis = sec_axis(~ . * scale_factor, name = "genus richness")
  ) +
  geom_smooth(aes(y = daylight_hours, color = "Daylight hours")) +
  geom_smooth(aes(y = genus_richness / scale_factor, color = "genus richness")) +
  labs(x = "Day of Year", color = NULL) +
  theme_minimal()

combined_plot <- ggarrange(
  p1, p2, 
  ncol = 2, nrow = 1,
  align = "hv"   
)

ggsave("DOY_DH_RUEdiv.png", 
       plot = combined_plot, width = 12, height = 4, dpi = 300)

#plot seasonal patterns of abiotic factors#...................................

m1 <- lm(TotalNitrogen_unfiltered ~ DOY, data = plankton_RUE_data7)
m2 <- lm(TotalNitrogen_unfiltered ~ poly(DOY, 2), data = plankton_RUE_data7)

AIC(m1,m2)

# Temperature
p_temp <- ggplot(plankton_RUE_data7, aes(x = DOY, y = Temperature)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = TRUE) +
  theme_bw() +
  labs(x = "Day of Year (DOY)", y = "Temperature (°C)")

# Salinity
p_sal <- ggplot(plankton_RUE_data7, aes(x = DOY, y = Salinity)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = TRUE) +
  theme_bw() +
  labs(x = "Day of Year (DOY)", y = "Salinity")

# pH
p_ph <- ggplot(plankton_RUE_data7, aes(x = DOY, y = pH)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_bw() +
  labs(x = "Day of Year (DOY)", y = "pH")

# Total Nitrogen
p_n <- ggplot(plankton_RUE_data7, aes(x = DOY, y = TotalNitrogen_unfiltered)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = TRUE) +
  theme_bw() +
  labs(x = "Day of Year (DOY)", y = "Nitrogen concentrations (µg/L)")

# Combine into one figure (2x2 layout)
combined_plot <- (p_temp | p_sal) / (p_ph | p_n)+
  plot_annotation(tag_levels = "a")

# Display
combined_plot

ggsave("environment_DOY_combined.png", combined_plot,
       width = 7, height = 6, dpi = 300)

#map of sites across seasons#.................................................
Data_prep <- plankton_sf
Data_prep$season <- cut(
  Data_prep$DOY,
  breaks = c(-Inf, 182, 243, Inf),
  labels = c("Early season", "Summer", "Late season")
)

early_data  <- subset(Data_prep, season == "Early season")
summer_data <- subset(Data_prep, season == "Summer")
late_data   <- subset(Data_prep, season == "Late season")

p_early <- ggplot() +
  geom_sf(data = world, fill = "gray95", color = "gray50") +
  geom_sf(data = early_data,
          color = "black",
          size = 1.5, alpha = 0.8) +
  coord_sf(xlim = xlim, ylim = ylim) +
  theme_minimal() +
  ggtitle("Early Season")

p_early

p_summer <- ggplot() +
  geom_sf(data = world, fill = "gray95", color = "gray50") +
  geom_sf(data = summer_data,
          color = "black",
          size = 1.5, alpha = 0.8) +
  coord_sf(xlim = xlim, ylim = ylim) +
  theme_minimal() +
  ggtitle("Summer")

p_summer

p_late <- ggplot() +
  geom_sf(data = world, fill = "gray95", color = "gray50") +
  geom_sf(data = late_data,
          color = "black",
          size = 1.5, alpha = 0.8) +
  coord_sf(xlim = xlim, ylim = ylim) +
  theme_minimal() +
  ggtitle("Late Season")

p_late

# Combine into one figure (1x3 layout)
combined_plot <- (
  (p_early | p_summer | p_late) /
    site_freq2
) + 
  plot_annotation(tag_levels = "a")

combined_plot

ggsave("season_sites.png", combined_plot,
       width = 8, height = 8, dpi = 300)

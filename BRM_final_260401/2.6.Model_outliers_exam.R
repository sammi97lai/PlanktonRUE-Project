################################################################################
####################### Examine model fitting and outliers #####################
################################################################################

library(brms)
library(dplyr)
library(tidyr)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(ggpubr)

rm(list = ls())

setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/RDATA")
load('plankton_RUE_data_scaled_260408.RDATA')

plankton_sf <- st_as_sf(
  plankton_RUE_data_scaled,
  coords = c("longitude", "latitude"),
  crs = 4326
)

world <- ne_countries(scale = "medium", returnclass = "sf")

xlim <- c(17, 26)
ylim <- c(57, 66)

#------------------------------------------------------------------------
# Phytoplankton FDis Model- polynomial

setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA")
load("BRM_phyto_FDis_P_260408.RDATA")

# posterior predictive draws
yrep <- posterior_predict(model)
# observed response
y <- model$data$FDis
# posterior predictive tail probability
pp_tail <- sapply(1:length(y), function(i) {
  mean(yrep[, i] <= y[i])
})
# flag extreme observations
outlier_idx <- which(pp_tail < 0.025 | pp_tail > 0.975)
length(outlier_idx)
outlier_idx

outliers <- model$data[outlier_idx, ]

outliers_phytoFDis <- plankton_RUE_data_scaled[outlier_idx, c(
  "FDis",
  "latitude",
  "longitude",
  "Site_Id",
  "Year",
  "DOY",
  "Temperature",
  "Salinity",
  "pH",
  "TotalNitrogen_unfiltered",
  "MaxDepth"
)]

head(outliers_phytoFDis)

outliers_phytoFDis$pp_tail <- pp_tail[outlier_idx]

outliers_phytoFDis <- outliers_phytoFDis[order(
  pmin(outliers_phytoFDis$pp_tail, 1 - outliers_phytoFDis$pp_tail)
), ]

head(outliers_phytoFDis, 20)

#compare distributions vs non outliers
plankton_RUE_data_scaled$outlier <- FALSE
plankton_RUE_data_scaled$outlier[outlier_idx] <- TRUE
ggplot(plankton_RUE_data_scaled, aes(x = TotalNitrogen_unfiltered, fill = outlier)) + #repeat for other abiotic var#
  geom_density(alpha = 0.4)

#model residual structure
yrep_mean <- apply(yrep, 2, mean)
residual <- y - yrep_mean
plankton_RUE_data_scaled$residual <- residual
ggplot(plankton_RUE_data_scaled, aes(Temperature, residual)) + #repeat for other abiotic var#
  geom_point(alpha = 0.4) +
  geom_smooth()

#plotting
outliers_phytoFDis_sf <- st_as_sf(
  outliers_phytoFDis,
  coords = c("longitude", "latitude"),
  crs = 4326
)

p1 <- ggplot() +
  geom_sf(data = world, fill = "gray95", color = "gray50") +
  geom_sf(
    data = plankton_sf,
    color = "gray60",
    alpha = 0.35,
    size = 1
  ) +
  # outliers highlighted
  geom_sf(
    data = outliers_phytoFDis_sf,
    color = "green",
    size = 3
  ) +
  coord_sf(xlim = xlim, ylim = ylim) +
  theme_minimal() +
  labs(
    title = "Posterior Predictive Outliers for Phyto FDis",
    x = NULL,
    y = NULL
  )
p1

p2 <- ggplot(plankton_RUE_data_scaled, aes(x = DOY, y = FDis)) +
  geom_point(color = "gray70", alpha = 0.4) +
  geom_point(
    data = outliers_phytoFDis,
    aes(x = DOY, y = FDis),
    color = "green",
    size = 2.5
  ) +
  theme_minimal() +
  labs(
    title = "FDis over Day of Year with Posterior Predictive Outliers",
    x = "Day of Year (DOY)",
    y = "Functional Diversity (FDis)"
  )
p2

combined_plot <- ggarrange(
  p1, p2,  ncol = 2, nrow = 1,
  align = "hv",
  labels = c("a", "b"))
combined_plot

#------------------------------------------------------------------------------
# Cyanobacteria dominance  Model- polynomial

setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA")
load("BRM_cyanodom_P_260408_z.RDATA")

# posterior predictive draws
yrep <- posterior_predict(model)
# observed response
y <- model$data$cyano_prop
# posterior predictive tail probability
pp_tail <- sapply(1:length(y), function(i) {
  mean(yrep[, i] <= y[i])
})
# flag extreme observations
outlier_idx <- which(pp_tail < 0.025 | pp_tail > 0.975)
length(outlier_idx)
outlier_idx

outliers <- model$data[outlier_idx, ]

outliers_cyano <- plankton_RUE_data_scaled[outlier_idx, c(
  "cyano_prop",
  "latitude",
  "longitude",
  "Site_Id",
  "Year",
  "DOY",
  "Temperature",
  "Salinity",
  "pH",
  "TotalNitrogen_unfiltered",
  "MaxDepth"
)]

head(outliers_cyano)

outliers_cyano$pp_tail <- pp_tail[outlier_idx]

outliers_cyano <- outliers_cyano[order(
  pmin(outliers_cyano$pp_tail, 1 - outliers_cyano$pp_tail)
), ]

head(outliers_cyano, 20)

#compare distributions vs non outliers
plankton_RUE_data_scaled$outlier <- FALSE
plankton_RUE_data_scaled$outlier[outlier_idx] <- TRUE
ggplot(plankton_RUE_data_scaled, aes(x = TotalNitrogen_unfiltered, fill = outlier)) + #repeat for other abiotic var#
  geom_density(alpha = 0.4)

#model residual structure
yrep_mean <- apply(yrep, 2, mean)
residual <- y - yrep_mean
plankton_RUE_data_scaled$residual <- residual
ggplot(plankton_RUE_data_scaled, aes(FDis, residual)) + #repeat for other abiotic var#
  geom_point(alpha = 0.4) +
  geom_smooth()

#plotting
outliers_cyano_sf <- st_as_sf(
  outliers_cyano,
  coords = c("longitude", "latitude"),
  crs = 4326
)

p1 <- ggplot() +
  geom_sf(data = world, fill = "gray95", color = "gray50") +
  geom_sf(
    data = plankton_sf,
    color = "gray60",
    alpha = 0.35,
    size = 1
  ) +
  # outliers highlighted
  geom_sf(
    data = outliers_cyano_sf,
    color = "green",
    size = 3
  ) +
  coord_sf(xlim = xlim, ylim = ylim) +
  theme_minimal() +
  labs(
    title = "Posterior Predictive Outliers for cyanobacteria",
    x = NULL,
    y = NULL
  )
p1

p2 <- ggplot(plankton_RUE_data_scaled, aes(x = DOY, y = cyano_prop)) +
  geom_point(color = "gray70", alpha = 0.4) +
  geom_point(
    data = outliers_cyano,
    aes(x = DOY, y = cyano_prop),
    color = "green",
    size = 2.5
  ) +
  theme_minimal() +
  labs(
    title = "Cyanobacteria over Day of Year with Posterior Predictive Outliers",
    x = "Day of Year (DOY)",
    y = "Cyanobacteria Dominance"
  )

p2

combined_plot <- ggarrange(
  p1, p2,  ncol = 2, nrow = 1,
  align = "hv",
  labels = c("a", "b"))
combined_plot

#-----------------------------------------------------------------------------
# Zooplankton FDis Model- polynomial

setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA")
load("BRM_zoo_FD_CD_P_260427.RDATA")

# posterior predictive draws
yrep <- posterior_predict(model)
# observed response
y <- model$data$zoo_FDis
# posterior predictive tail probability
pp_tail <- sapply(1:length(y), function(i) {
  mean(yrep[, i] <= y[i])
})
# flag extreme observations
outlier_idx <- which(pp_tail < 0.025 | pp_tail > 0.975)
length(outlier_idx)
outlier_idx

outliers <- model$data[outlier_idx, ]

outliers_zooFDis <- plankton_RUE_data_scaled[outlier_idx, c(
  "zoo_FDis",
  "latitude",
  "longitude",
  "Site_Id",
  "Year",
  "DOY",
  "Temperature",
  "Salinity",
  "pH",
  "TotalNitrogen_unfiltered",
  "genus_richness",
  "FEve", "FDis", "cyano_prop",
  "MaxDepth"
)]

head(outliers_zooFDis)

outliers_zooFDis$pp_tail <- pp_tail[outlier_idx]

outliers_zooFDis <- outliers_zooFDis[order(
  pmin(outliers_zooFDis$pp_tail, 1 - outliers_zooFDis$pp_tail)
), ]

head(outliers_zooFDis, 20)

#compare distributions vs non outliers
plankton_RUE_data_scaled$outlier <- FALSE
plankton_RUE_data_scaled$outlier[outlier_idx] <- TRUE
ggplot(plankton_RUE_data_scaled, aes(x = cyano_prop, fill = outlier)) + #repeat for other abiotic var#
  geom_density(alpha = 0.4)

#model residual structure
yrep_mean <- apply(yrep, 2, mean)
residual <- y - yrep_mean
plankton_RUE_data_scaled$residual <- residual
ggplot(plankton_RUE_data_scaled, aes(FDis, residual)) + #repeat for other abiotic var#
  geom_point(alpha = 0.4) +
  geom_smooth()

#plotting
outliers_zooFDis_sf <- st_as_sf(
  outliers_zooFDis,
  coords = c("longitude", "latitude"),
  crs = 4326
)

p1 <- ggplot() +
  geom_sf(data = world, fill = "gray95", color = "gray50") +
  geom_sf(
    data = plankton_sf,
    color = "gray60",
    alpha = 0.35,
    size = 1
  ) +
  # outliers highlighted
  geom_sf(
    data = outliers_zooFDis_sf,
    color = "green",
    size = 3
  ) +
  coord_sf(xlim = xlim, ylim = ylim) +
  theme_minimal() +
  labs(
    title = "Posterior Predictive Outliers for Zoo FDis",
    x = NULL,
    y = NULL
  )
p1

p2 <- ggplot(plankton_RUE_data_scaled, aes(x = DOY, y = zoo_FDis)) +
  geom_point(color = "gray70", alpha = 0.4) +
  geom_point(
    data = outliers_zooFDis,
    aes(x = DOY, y = zoo_FDis),
    color = "green",
    size = 2.5
  ) +
  theme_minimal() +
  labs(
    title = "FDis over Day of Year with Posterior Predictive Outliers",
    x = "Day of Year (DOY)",
    y = "Functional Diversity (FDis)"
  )

p2

combined_plot <- ggarrange(
  p1, p2,  ncol = 2, nrow = 1,
  align = "hv",
  labels = c("a", "b"))
combined_plot

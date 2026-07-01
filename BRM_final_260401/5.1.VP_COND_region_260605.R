##############################################################################
############## Conditional Variance partitioning by Region ###################
##############################################################################
# manipulating data
library(lubridate)
library(tidyverse)
library(dplyr)
library('Matrix')
library(reshape2)
library(data.table)
library(lme4)

# bayesian modeling
library(rstan) 
library("bayesplot")
library("rstanarm")
library(brms)

## plotting
library("ggplot2")
library(vioplot)
library(ggpubr)
library(maps)
library(patchwork)

rm(list = ls())
##...........................................................................................................

## Load variance functions modified from Torsti's
setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401")
source('Functions_for_brms_VP_260605.R')

# cleaner x axis labels
label_names <- c(
  Abiotic = "Abiotic",
  Phyto_genus_richness = "Phytoplankton\nGenus\nRichness",
  Phyto_FEve = "Phytoplankton\nFEve",
  Phyto_FDis = "Phytoplankton\nFDis",
  Cyano_dom = "Cyanobacteria\nDominance",
  Zoo_species_richness = "Zooplankton\nSpecies\nRichness",
  Zoo_FEve = "Zooplankton\nFEve",
  Zoo_FDis = "Zooplankton\nFDis",
  Site = "Site",
  Year = "Year"
)

#Phytoplankton RUE per P...........................................................
setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA")
load('BRM_phytoRUEP_CD_PI_260604.RDATA')

## this new function now allow you to extract the quantities in the proper format needed
modphytoRUEP_CD_PI = PrepbrmsVP(model)
n_iter <- nrow(modphytoRUEP_CD_PI$params)

params = modphytoRUEP_CD_PI$params
X = modphytoRUEP_CD_PI$X
y = modphytoRUEP_CD_PI$y

# select the names of the covariates
colnames(X)
model$formula

unique(model$data$Site_Id)

#define groups
groups <- list(
  Abiotic = c("Temperature", "ITemperatureE2",                     
              "Salinity", "ISalinityE2", 
              "pH", "IpHE2", 
              "TotalNitrogen_unfiltered", "ITotalNitrogen_unfilteredE2", 
              "MaxDepth",  "IMaxDepthE2",
              'daylight_hours','Idaylight_hoursE2'),  ## 
  Phyto_genus_richness = c("genus_richness", "Temperature:genus_richness", "Salinity:genus_richness", "TotalNitrogen_unfiltered:genus_richness"),
  Phyto_FEve = c("FEve", "Temperature:FEve", "Salinity:FEve", "TotalNitrogen_unfiltered:FEve"),
  Phyto_FDis = c("FDis", "Temperature:FDis", "Salinity:FDis", "TotalNitrogen_unfiltered:FDis"),
  Cyano_dom = c("cyano_prop", "Temperature:cyano_prop", "Salinity:cyano_prop", "TotalNitrogen_unfiltered:cyano_prop"),
  Zoo_species_richness = c("species_richness", "Temperature:species_richness", "Salinity:species_richness", "TotalNitrogen_unfiltered:species_richness"),
  Zoo_FEve = c("zoo_FEve", "Temperature:zoo_FEve", "Salinity:zoo_FEve", "TotalNitrogen_unfiltered:zoo_FEve"),
  Zoo_FDis = c("zoo_FDis", "Temperature:zoo_FDis", "Salinity:zoo_FDis", "TotalNitrogen_unfiltered:zoo_FDis"),
  Site = c('Site_Id'),
  Year = c('Year')#,
  # Residuals = 'residuals'
) #group variables into categories#

# define some needed parameters
linear_terms <- unname(unlist(groups)) #flattened vector of all variable names
n_lt <- length(linear_terms)
n_group <- length(groups)

# --- Define the conditional grouping (e.g., region) ---
#get site info (Helcome area and coordinates)
setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/RDATA")
load('Phytoplankton_data_format.RDATA')
site_info <- Phytoplankton_data %>%
  select(Site_Id,
         HELCOMarea,
         `Coordinates_ EUREF_WGS84_lat`,
         `Coordinates_ EUREF_WGS84_long`) %>%
  distinct()

#join site info with model data
Data_prep <- model$data

Data_prep$Site_Id <- as.character(Data_prep$Site_Id)
site_info$Site_Id <- as.character(site_info$Site_Id)

Data_prep <- Data_prep %>%
  left_join(site_info, by = "Site_Id")

# unique(Data_prep$HELCOMarea)
# unique(Data_prep$Site_Id)

Data_prep$region <- NA

lat <- Data_prep$`Coordinates_ EUREF_WGS84_lat`
lon <- Data_prep$`Coordinates_ EUREF_WGS84_long`

# 1. Bothnian Bay
Data_prep$region[lat >= 64] <- "Bothnian Bay"

# 2. Bothnian Sea
Data_prep$region[lat >= 61 & lat < 64] <- "Bothnian Sea"

# 3. Gulf of Finland
Data_prep$region[lon >= 24 & lat < 61.5] <- "Gulf of Finland"

# 4. Baltic Proper
Data_prep$region[lat < 61.5 & lon < 24] <- "Baltic Proper"

# Data_prep %>%
#   distinct(Site_Id, region) %>%
#   count(region)

g <- Data_prep$region
G <- unique(g)

n_cond <- length(G)

# plotting sites and regions
# get world map data
world <- map_data("world")

# prepare your site data
sites <- Data_prep %>%
  distinct(Site_Id, region, `Coordinates_ EUREF_WGS84_long`,
           `Coordinates_ EUREF_WGS84_lat`) %>%
  rename(lon = `Coordinates_ EUREF_WGS84_long`,
         lat = `Coordinates_ EUREF_WGS84_lat`)

#plot regions
regions <- ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),
               fill = "gray90", color = "gray60", size = 0.3) +
  geom_point(data = sites, aes(x = lon, y = lat, color = region),
             size = 2.5, alpha = 0.8) +
  coord_quickmap(xlim = c(17, 26), ylim = c(57, 66)) +
  scale_color_brewer(palette = "Set1") +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(legend.position = "bottom")
regions

setwd("~/GitHub/PlanktonRUE-Project/Figures")
ggsave("region_COND_VP.png",
       plot = regions, width = 6, height = 8, dpi = 300)

# Bothnian Bay
bb <- sites %>%
  filter(region == "Bothnian Bay")

bb_plot <- ggplot() +
  geom_polygon(
    data = world,
    aes(x = long, y = lat, group = group),
    fill = "gray90",
    color = "gray60",
    linewidth = 0.3
  ) +
  geom_point(
    data = bb,
    aes(x = lon, y = lat),
    color = "#377EB8FF",
    size = 2.5
  ) +
  geom_rect(
    aes(
      xmin = min(bb$lon) - 0.3,
      xmax = max(bb$lon) + 0.3,
      ymin = min(bb$lat) - 0.3,
      ymax = max(bb$lat) + 0.3
    ),
    fill = NA,
    color = "black",
    linewidth = 1
  ) +
  coord_quickmap(xlim = c(17, 26), ylim = c(57, 66)) +
  labs(title = "Bothnian Bay",
       x = "Longitude",
       y = "Latitude") +
  theme_minimal()

bb_plot

# Bothnian Sea
bs <- sites %>%
  filter(region == "Bothnian Sea")

bs_plot <- ggplot() +
  geom_polygon(data = world,
               aes(long, lat, group = group),
               fill = "gray90",
               color = "gray60") +
  geom_point(data = bs,
             aes(lon, lat),
             color = "#4DAF4AFF",
             size = 2.5) +
  geom_rect(
    aes(
      xmin = min(bs$lon)-0.3,
      xmax = max(bs$lon)+0.3,
      ymin = min(bs$lat)-0.3,
      ymax = max(bs$lat)+0.3
    ),
    fill = NA,
    color = "black"
  ) +
  coord_quickmap(xlim = c(17, 26), ylim = c(57, 66)) +
  labs(title = "Bothnian Sea",
       x = "Longitude",
       y = "Latitude") +
  theme_minimal()

bs_plot

# Gulf of Finland
gof <- sites %>%
  filter(region == "Gulf of Finland")

gof_plot <- ggplot() +
  geom_polygon(data = world,
               aes(long, lat, group = group),
               fill = "gray90",
               color = "gray60") +
  geom_point(data = gof,
             aes(lon, lat),
             color = "#984EA3FF",
             size = 2.5) +
  geom_rect(
    aes(
      xmin = min(gof$lon)-0.3,
      xmax = max(gof$lon)+0.3,
      ymin = min(gof$lat)-0.3,
      ymax = max(gof$lat)+0.3
    ),
    fill = NA,
    color = "black"
  ) +
  coord_quickmap(xlim = c(17, 26), ylim = c(57, 66)) +
  labs(title = "Gulf of Finland",
       x = "Longitude",
       y = "Latitude") +
  theme_minimal()

gof_plot

# Baltic Proper
bp <- sites %>%
  filter(region == "Baltic Proper")

bp_plot <- ggplot() +
  geom_polygon(data = world,
               aes(long, lat, group = group),
               fill = "gray90",
               color = "gray60") +
  geom_point(data = bp,
             aes(lon, lat),
             color = "#E41A1CFF",
             size = 2.5) +
  geom_rect(
    aes(
      xmin = min(bp$lon)-0.3,
      xmax = max(bp$lon)+0.3,
      ymin = min(bp$lat)-0.3,
      ymax = max(bp$lat)+0.3
    ),
    fill = NA,
    color = "black"
  ) +
  coord_quickmap(xlim = c(17, 26), ylim = c(57, 66)) +
  labs(title = "Baltic Proper",
       x = "Longitude",
       y = "Latitude") +
  theme_minimal()

bp_plot

setwd("~/GitHub/PlanktonRUE-Project/Figures")
ggsave("region_VP_bb.png",
       plot = bb_plot, width = 2, height = 3, dpi = 300)
setwd("~/GitHub/PlanktonRUE-Project/Figures")
ggsave("region_VP_bs.png",
       plot = bs_plot, width = 2, height = 3, dpi = 300)
setwd("~/GitHub/PlanktonRUE-Project/Figures")
ggsave("region_VP_gof.png",
       plot = gof_plot, width = 2, height = 3, dpi = 300)
setwd("~/GitHub/PlanktonRUE-Project/Figures")
ggsave("region_VP_bp.png",
       plot = bp_plot, width = 2, height = 3, dpi = 300)

# --- Conditional variance partitioning ---
setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/RDATA")
# Using Torsti's function (adapted for single-response model)
K_V <- apply(params, 'iterations', conditional_lin_term_cov, G, g, X, y, residuals = FALSE)

dim(K_V) <- c(n_cond, n_lt, n_lt, n_iter)
dimnames(K_V) <- list(
  condition = G,
  rows = linear_terms,
  cols = linear_terms,
  iterations = NULL
)

# Collapse by group
B <- sapply(groups, function(x) as.integer(linear_terms %in% x))
dimnames(B) <- list(linear_terms = linear_terms, groups = colnames(B))

K_V <- apply(K_V, c('condition', 'iterations'), function(x, B) { t(B) %*% x %*% B }, B)
dim(K_V) <- c(n_group, n_group, n_cond, n_iter)
dimnames(K_V) <- list(
  rows = names(groups),
  cols = names(groups),
  condition = G,
  iterations = NULL
)

# Compute conditional VP (sum(diag) normalization)
V_cond <- apply(K_V, c('condition', 'iterations'), function(x) diag(x)/sum(diag(x)))

#now weighted by R2
R2_model = bayes_R2(model)
R2_model[,1]  # we keep only the estimate although the uncertainy can be sometimes interesting it would be a bit complex and out of subject for you I guess

dim(V_cond) #  covar   10    cond  4     mcmc samples 8000

# create a weighting matrix of same dimension as you data (particularly good if you have also multiple species but here it is ok still)
varr <- aperm(array(R2_model[,1], dim = c(8000L, 4L, 10L)), perm = c(3L, 2L, 1L)) # check the dimensions carefully 10 covar (1L) 4 conditions (4L) for 8000 samples (2L)
dim(varr)  # [1]   10   4 8000

V_condweighted = V_cond*varr

# --- Prepare data for plotting ---
dim3W <- apply(V_condweighted, 3L, c)
dim12W <- expand.grid(dimnames(V_condweighted)[1:2])
dataW_wide <- data.frame(dim3W, dim12W)
dataW_long <- melt(setDT(dataW_wide), id.vars = c('rows','condition'), variable.name = 'iteration')

# Keep order of regions
dataW_long$condition <- factor(dataW_long$condition, levels = c("Bothnian Bay", "Bothnian Sea",
                                                                "Gulf of Finland", "Baltic Proper"))

#get vp summary 
vp_group_DOY_summary <- data.frame(
  expand.grid(
    Group = dimnames(V_condweighted)[[1]],
    Region = dimnames(V_condweighted)[[2]]
  ),
  
  MeanVP = apply(V_condweighted, c(1, 2), mean, na.rm = TRUE) |> c(),
  
  Lower95 = apply(V_condweighted, c(1, 2), quantile, probs = 0.025, na.rm = TRUE
  ) |> c(),
  
  Upper95 = apply(V_condweighted, c(1, 2), quantile, probs = 0.975, na.rm = TRUE
  ) |> c()
)

vp_group_DOY_summary

# Plotting
setwd("~/GitHub/PlanktonRUE-Project/Figures/VP_plots")

vp_RUEP_regionW <- ggplot(dataW_long, aes(x = rows, y = value, fill = rows)) +
  geom_violin(position="dodge", trim=F, scale="width", draw_quantiles=c(0.25,0.5,0.75), linewidth = 0.25) +
  scale_x_discrete(labels = label_names, guide = guide_axis(n.dodge = 2)) +
  scale_fill_manual(values=c('royalblue3','#278B9AFF','#44A57CFF','#819A7AFF','#58A449FF', '#403369FF','#5C5992FF','#AE93BEFF','salmon2','gold3')) +
  facet_grid(condition ~ .) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 8, lineheight = 0.85),
    axis.title.x = element_blank(),   # remove x-axis title
    axis.title.y = element_blank(),   # remove y-axis title
    legend.position = "none",          # remove colour legend
    plot.title = element_blank()
  )
vp_RUEP_regionW

ggsave(filename = "RUEP_vpcond_region_r2.png",
       plot = vp_RUEP_regionW,
       width = 5,
       height =7,
       dpi = 300)

#Phytoplankton RUE per N...........................................................
setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA")
load('BRM_phytoRUEN_CD_PI_260604.RDATA')

## this new function now allow you to extract the quantities in the proper format needed
modphytoRUEN_CD_PI = PrepbrmsVP(model)
n_iter <- nrow(modphytoRUEN_CD_PI$params)

params = modphytoRUEN_CD_PI$params
X = modphytoRUEN_CD_PI$X
y = modphytoRUEN_CD_PI$y

# select the names of the covariates
colnames(X)
model$formula

unique(model$data$Site_Id)

#define groups
groups <- list(
  Abiotic = c("Temperature", "ITemperatureE2",                     
              "Salinity", "ISalinityE2", 
              "pH", "IpHE2", 
              "TotalNitrogen_unfiltered", "ITotalNitrogen_unfilteredE2", 
              "MaxDepth",  "IMaxDepthE2",
              'daylight_hours','Idaylight_hoursE2'),  ## 
  Phyto_genus_richness = c("genus_richness", "Temperature:genus_richness", "Salinity:genus_richness", "TotalNitrogen_unfiltered:genus_richness"),
  Phyto_FEve = c("FEve", "Temperature:FEve", "Salinity:FEve", "TotalNitrogen_unfiltered:FEve"),
  Phyto_FDis = c("FDis", "Temperature:FDis", "Salinity:FDis", "TotalNitrogen_unfiltered:FDis"),
  Cyano_dom = c("cyano_prop", "Temperature:cyano_prop", "Salinity:cyano_prop", "TotalNitrogen_unfiltered:cyano_prop"),
  Zoo_species_richness = c("species_richness", "Temperature:species_richness", "Salinity:species_richness", "TotalNitrogen_unfiltered:species_richness"),
  Zoo_FEve = c("zoo_FEve", "Temperature:zoo_FEve", "Salinity:zoo_FEve", "TotalNitrogen_unfiltered:zoo_FEve"),
  Zoo_FDis = c("zoo_FDis", "Temperature:zoo_FDis", "Salinity:zoo_FDis", "TotalNitrogen_unfiltered:zoo_FDis"),
  Site = c('Site_Id'),
  Year = c('Year')#,
  # Residuals = 'residuals'
) #group variables into categories#

# define some needed parameters
linear_terms <- unname(unlist(groups)) #flattened vector of all variable names
n_lt <- length(linear_terms)
n_group <- length(groups)

# --- Define the conditional grouping (e.g., region) ---
#get site info (Helcome area and coordinates)
setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/RDATA")
load('Phytoplankton_data_format.RDATA')
site_info <- Phytoplankton_data %>%
  select(Site_Id,
         HELCOMarea,
         `Coordinates_ EUREF_WGS84_lat`,
         `Coordinates_ EUREF_WGS84_long`) %>%
  distinct()

#join site info with model data
Data_prep <- model$data

Data_prep$Site_Id <- as.character(Data_prep$Site_Id)
site_info$Site_Id <- as.character(site_info$Site_Id)

Data_prep <- Data_prep %>%
  left_join(site_info, by = "Site_Id")

# unique(Data_prep$HELCOMarea)
# unique(Data_prep$Site_Id)

Data_prep$region <- NA

lat <- Data_prep$`Coordinates_ EUREF_WGS84_lat`
lon <- Data_prep$`Coordinates_ EUREF_WGS84_long`

# 1. Bothnian Bay
Data_prep$region[lat >= 64] <- "Bothnian Bay"

# 2. Bothnian Sea
Data_prep$region[lat >= 61 & lat < 64] <- "Bothnian Sea"

# 3. Gulf of Finland
Data_prep$region[lon >= 24 & lat < 61.5] <- "Gulf of Finland"

# 4. Baltic Proper
Data_prep$region[lat < 61.5 & lon < 24] <- "Baltic Proper"

# Data_prep %>%
#   distinct(Site_Id, region) %>%
#   count(region)

g <- Data_prep$region
G <- unique(g)

n_cond <- length(G)

#plotting sites and regions
# get world map data
world <- map_data("world")

# prepare your site data
sites <- Data_prep %>%
  distinct(Site_Id, region, `Coordinates_ EUREF_WGS84_long`, 
           `Coordinates_ EUREF_WGS84_lat`) %>%
  rename(lon = `Coordinates_ EUREF_WGS84_long`,
         lat = `Coordinates_ EUREF_WGS84_lat`)

# --- Conditional variance partitioning ---
# Using Torsti's function (adapted for single-response model)
K_V <- apply(params, 'iterations', conditional_lin_term_cov, G, g, X, y, residuals = FALSE)

dim(K_V) <- c(n_cond, n_lt, n_lt, n_iter)
dimnames(K_V) <- list(
  condition = G,
  rows = linear_terms,
  cols = linear_terms,
  iterations = NULL
)

# Collapse by group
B <- sapply(groups, function(x) as.integer(linear_terms %in% x))
dimnames(B) <- list(linear_terms = linear_terms, groups = colnames(B))

K_V <- apply(K_V, c('condition', 'iterations'), function(x, B) { t(B) %*% x %*% B }, B)
dim(K_V) <- c(n_group, n_group, n_cond, n_iter)
dimnames(K_V) <- list(
  rows = names(groups),
  cols = names(groups),
  condition = G,
  iterations = NULL
)

# Compute conditional VP (sum(diag) normalization)
V_cond <- apply(K_V, c('condition', 'iterations'), function(x) diag(x)/sum(diag(x)))

#now weighted by R2
R2_model = bayes_R2(model)
R2_model[,1]  # we keep only the estimate although the uncertainy can be sometimes interesting it would be a bit complex and out of subject for you I guess

dim(V_cond) #  covar   10    cond  4     mcmc samples 8000

# create a weighting matrix of same dimension as you data (particularly good if you have also multiple species but here it is ok still)
varr <- aperm(array(R2_model[,1], dim = c(8000L, 4L, 10L)), perm = c(3L, 2L, 1L)) # check the dimensions carefully 10 covar (1L) 4 conditions (4L) for 8000 samples (2L)
dim(varr)  # [1]   10   4 8000

V_condweighted = V_cond*varr

# --- Prepare data for plotting ---
dim3W <- apply(V_condweighted, 3L, c)
dim12W <- expand.grid(dimnames(V_condweighted)[1:2])
dataW_wide <- data.frame(dim3W, dim12W)
dataW_long <- melt(setDT(dataW_wide), id.vars = c('rows','condition'), variable.name = 'iteration')

# Keep order of regions
dataW_long$condition <- factor(dataW_long$condition, levels = c("Bothnian Bay", "Bothnian Sea",
                                                              "Gulf of Finland", "Baltic Proper"))

#get vp summary 
vp_group_DOY_summary <- data.frame(
  expand.grid(
    Group = dimnames(V_condweighted)[[1]],
    Region = dimnames(V_condweighted)[[2]]
  ),
  
  MeanVP = apply(V_condweighted, c(1, 2), mean, na.rm = TRUE) |> c(),
  
  Lower95 = apply(V_condweighted, c(1, 2), quantile, probs = 0.025, na.rm = TRUE
  ) |> c(),
  
  Upper95 = apply(V_condweighted, c(1, 2), quantile, probs = 0.975, na.rm = TRUE
  ) |> c()
)

vp_group_DOY_summary
# Plotting
setwd("~/GitHub/PlanktonRUE-Project/Figures/VP_plots")

vp_RUEN_regionW <- ggplot(dataW_long, aes(x = rows, y = value, fill = rows)) +
  geom_violin(position="dodge", trim=F, scale="width", draw_quantiles=c(0.25,0.5,0.75), linewidth = 0.25) +
  scale_x_discrete(labels = label_names, guide = guide_axis(n.dodge = 2)) +
  scale_fill_manual(values=c('royalblue3','#278B9AFF','#44A57CFF','#819A7AFF','#58A449FF', '#403369FF','#5C5992FF','#AE93BEFF','salmon2','gold3')) +
  facet_grid(condition ~ .) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 8, lineheight = 0.85),
    axis.title.x = element_blank(),   # remove x-axis title
    axis.title.y = element_blank(),   # remove y-axis title
    legend.position = "none",          # remove colour legend
    plot.title = element_blank()
  )
vp_RUEN_regionW

ggsave(filename = "RUEN_vpcond_region_r2.png",
       plot = vp_RUEN_regionW,
       width = 5,
       height =7,
       dpi = 300)

#Zooplankton RUE ...........................................................
setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA")
load('BRM_zoo_RUE_CD_PI_260604.RDATA')

## this new function now allow you to extract the quantities in the proper format needed
modZoo_RUE_CD_PI = PrepbrmsVP(model)
n_iter <- nrow(modZoo_RUE_CD_PI$params)

params = modZoo_RUE_CD_PI$params
X = modZoo_RUE_CD_PI$X
y = modZoo_RUE_CD_PI$y

# select the names of the covariates
colnames(X)
model$formula

unique(model$data$Site_Id)

#define groups
groups <- list(
  Abiotic = c("Temperature", "ITemperatureE2",                     
              "Salinity", "ISalinityE2", 
              "pH", "IpHE2", 
              "TotalNitrogen_unfiltered", "ITotalNitrogen_unfilteredE2", 
              "Zooplankton_MaxDepth",  "IZooplankton_MaxDepthE2",
              'daylight_hours','Idaylight_hoursE2'),  ## 
  Phyto_genus_richness = c("genus_richness", "Temperature:genus_richness", "Salinity:genus_richness", "TotalNitrogen_unfiltered:genus_richness"),
  Phyto_FEve = c("FEve", "Temperature:FEve", "Salinity:FEve", "TotalNitrogen_unfiltered:FEve"),
  Phyto_FDis = c("FDis", "Temperature:FDis", "Salinity:FDis", "TotalNitrogen_unfiltered:FDis"),
  Cyano_dom = c("cyano_prop", "Temperature:cyano_prop", "Salinity:cyano_prop", "TotalNitrogen_unfiltered:cyano_prop"),
  Zoo_species_richness = c("species_richness", "Temperature:species_richness", "Salinity:species_richness", "TotalNitrogen_unfiltered:species_richness"),
  Zoo_FEve = c("zoo_FEve", "Temperature:zoo_FEve", "Salinity:zoo_FEve", "TotalNitrogen_unfiltered:zoo_FEve"),
  Zoo_FDis = c("zoo_FDis", "Temperature:zoo_FDis", "Salinity:zoo_FDis", "TotalNitrogen_unfiltered:zoo_FDis"),
  Site = c('Site_Id'),
  Year = c('Year')#,
  # Residuals = 'residuals'
) #group variables into categories#

# define some needed parameters
linear_terms <- unname(unlist(groups)) #flattened vector of all variable names
n_lt <- length(linear_terms)
n_group <- length(groups)

# --- Define the conditional grouping (e.g., region) ---
#get site info (Helcome area and coordinates)
setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/RDATA")
load('Phytoplankton_data_format.RDATA')
site_info <- Phytoplankton_data %>%
  select(Site_Id,
         HELCOMarea,
         `Coordinates_ EUREF_WGS84_lat`,
         `Coordinates_ EUREF_WGS84_long`) %>%
  distinct()

#join site info with model data
Data_prep <- model$data

Data_prep$Site_Id <- as.character(Data_prep$Site_Id)
site_info$Site_Id <- as.character(site_info$Site_Id)

Data_prep <- Data_prep %>%
  left_join(site_info, by = "Site_Id")

# unique(Data_prep$HELCOMarea)
# unique(Data_prep$Site_Id)

Data_prep$region <- NA

lat <- Data_prep$`Coordinates_ EUREF_WGS84_lat`
lon <- Data_prep$`Coordinates_ EUREF_WGS84_long`

# 1. Bothnian Bay
Data_prep$region[lat >= 64] <- "Bothnian Bay"

# 2. Bothnian Sea
Data_prep$region[lat >= 61 & lat < 64] <- "Bothnian Sea"

# 3. Gulf of Finland
Data_prep$region[lon >= 24 & lat < 61.5] <- "Gulf of Finland"

# 4. Baltic Proper
Data_prep$region[lat < 61.5 & lon < 24] <- "Baltic Proper"

# Data_prep %>%
#   distinct(Site_Id, region) %>%
#   count(region)

g <- Data_prep$region
G <- unique(g)

n_cond <- length(G)

#plotting sites and regions
# get world map data
world <- map_data("world")

# prepare your site data
sites <- Data_prep %>%
  distinct(Site_Id, region, `Coordinates_ EUREF_WGS84_long`, 
           `Coordinates_ EUREF_WGS84_lat`) %>%
  rename(lon = `Coordinates_ EUREF_WGS84_long`,
         lat = `Coordinates_ EUREF_WGS84_lat`)

# --- Conditional variance partitioning ---
# Using Torsti's function (adapted for single-response model)
K_V <- apply(params, 'iterations', conditional_lin_term_cov, G, g, X, y, residuals = FALSE)

dim(K_V) <- c(n_cond, n_lt, n_lt, n_iter)
dimnames(K_V) <- list(
  condition = G,
  rows = linear_terms,
  cols = linear_terms,
  iterations = NULL
)

# Collapse by group
B <- sapply(groups, function(x) as.integer(linear_terms %in% x))
dimnames(B) <- list(linear_terms = linear_terms, groups = colnames(B))

K_V <- apply(K_V, c('condition', 'iterations'), function(x, B) { t(B) %*% x %*% B }, B)
dim(K_V) <- c(n_group, n_group, n_cond, n_iter)
dimnames(K_V) <- list(
  rows = names(groups),
  cols = names(groups),
  condition = G,
  iterations = NULL
)

# Compute conditional VP (sum(diag) normalization)
V_cond <- apply(K_V, c('condition', 'iterations'), function(x) diag(x)/sum(diag(x)))

#now weighted by R2
R2_model = bayes_R2(model)
R2_model[,1]  # we keep only the estimate although the uncertainy can be sometimes interesting it would be a bit complex and out of subject for you I guess

dim(V_cond) #  covar   10    cond  4     mcmc samples 8000

# create a weighting matrix of same dimension as you data (particularly good if you have also multiple species but here it is ok still)
varr <- aperm(array(R2_model[,1], dim = c(8000L, 4L, 10L)), perm = c(3L, 2L, 1L)) # check the dimensions carefully 10 covar (1L) 4 conditions (4L) for 8000 samples (2L)
dim(varr)  # [1]   10   4 8000

V_condweighted = V_cond*varr

# --- Prepare data for plotting ---
dim3W <- apply(V_condweighted, 3L, c)
dim12W <- expand.grid(dimnames(V_condweighted)[1:2])
dataW_wide <- data.frame(dim3W, dim12W)
dataW_long <- melt(setDT(dataW_wide), id.vars = c('rows','condition'), variable.name = 'iteration')


# # Keep order of regions
dataW_long$condition <- factor(dataW_long$condition, levels = c("Bothnian Bay", "Bothnian Sea",
                                                                "Gulf of Finland","Baltic Proper"))

#get vp summary 
vp_group_DOY_summary <- data.frame(
  expand.grid(
    Group = dimnames(V_condweighted)[[1]],
    Region = dimnames(V_condweighted)[[2]]
  ),
  
  MeanVP = apply(V_condweighted, c(1, 2), mean, na.rm = TRUE) |> c(),
  
  Lower95 = apply(V_condweighted, c(1, 2), quantile, probs = 0.025, na.rm = TRUE
  ) |> c(),
  
  Upper95 = apply(V_condweighted, c(1, 2), quantile, probs = 0.975, na.rm = TRUE
  ) |> c()
)

vp_group_DOY_summary

#Plotting
setwd("~/GitHub/PlanktonRUE-Project/Figures/VP_plots")

vp_zooRUE_regionW <- ggplot(dataW_long, aes(x = rows, y = value, fill = rows)) +
  geom_violin(position="dodge", trim=F, scale="width", draw_quantiles=c(0.25,0.5,0.75), linewidth = 0.25) +
  scale_x_discrete(labels = label_names, guide = guide_axis(n.dodge = 2)) +
  scale_fill_manual(values=c('royalblue3','#278B9AFF','#44A57CFF','#819A7AFF','#58A449FF', '#403369FF','#5C5992FF','#AE93BEFF','#F7EABDFF','salmon2','gold3')) +
  facet_grid(condition ~ .) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 8, lineheight = 0.85),
    axis.title.x = element_blank(),   # remove x-axis title
    axis.title.y = element_blank(),   # remove y-axis title
    legend.position = "none",          # remove colour legend
    plot.title = element_blank()
  )
vp_zooRUE_regionW

ggsave(filename = "zooRUE_vpcond_region_r2.png",
       plot = vp_zooRUE_regionW,
       width = 5,
       height =7,
       dpi = 300)

#one big figure.............................................................
combined_plot <- (vp_RUEP_regionW | vp_RUEN_regionW | vp_zooRUE_regionW) +
  plot_annotation(tag_levels = "a")

combined_plot

ggsave("VP_region_combined_r2.png",
       plot = combined_plot,
       width = 16, height = 7, dpi = 300)

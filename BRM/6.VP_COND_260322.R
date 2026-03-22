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
#library(brms)

## plotting
library("ggplot2")
library(vioplot)
library(ggpubr)
library(maps)

##...........................................................................................................

## Load variance functions modified from Torsti's
setwd("~/GitHub/PlanktonRUE-Project/BRM")
source('Functions_for_brms_VP.R')

#Phytoplankton RUE per P...........................................................
setwd("~/GitHub/PlanktonRUE-Project/BRM/RDATA")
load('BRM_phytoRUEP_CD_PI_260311.RDATA')

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
              "MaxDepth",  "IMaxDepthE2"),  ## 
  Phyto_div = c("genus_richness", "FEve", "FDis"),
  Zoo_div = c("species_richness", "zoo_FEve", "zoo_FDis"),
  Cyano_dom = c("cyano_prop"),
  Abiotic_phyto_div = c("Temperature:genus_richness", "Salinity:genus_richness",
                        "Temperature:FEve", "Salinity:FEve",
                        "Temperature:FDis", "Salinity:FDis"),
  Abiotic_cyano = c("Temperature:cyano_prop", "Salinity:cyano_prop"),
  Abiotic_zoo_div = c("Temperature:species_richness", "Salinity:species_richness",
                      "Temperature:zoo_FEve", "Salinity:zoo_FEve",
                      "Temperature:zoo_FDis", "Salinity:zoo_FDis"),
  DOY = c('DOY', 'IDOYE2',
          'daylight_hours','Idaylight_hoursE2'),
  Site = c('Site_Id'),
  Year = c('Year')#,
  #Residuals = 'residuals'
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

# plot regions
# ggplot() +
#   geom_polygon(data = world, aes(x = long, y = lat, group = group),
#                fill = "gray90", color = "gray60", size = 0.3) +
#   geom_point(data = sites, aes(x = lon, y = lat, color = region),
#              size = 2.5, alpha = 0.8) +
#   coord_quickmap(xlim = c(19, 31), ylim = c(58, 66)) +
#   scale_color_brewer(palette = "Set1") +
#   labs(x = "Longitude", y = "Latitude",
#        title = "Sampling Sites in Baltic Sea Grouped by Region") +
#   theme_minimal() +
#   theme(legend.position = "bottom")

# --- Conditional variance partitioning ---
# Using Torsti's function (adapted for single-response model)
K_V <- apply(params, 'iterations', conditional_lin_term_cov, G, g, X, y)

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

# --- Prepare data for plotting ---
dim3 <- apply(V_cond, 3L, c)
dim12 <- expand.grid(dimnames(V_cond)[1:2])
data_wide <- data.frame(dim3, dim12)
data_long <- melt(setDT(data_wide), id.vars = c('rows','condition'), variable.name = 'iteration')

# Keep order of regions
data_long$condition <- factor(data_long$condition, levels = c("Bothnian Bay", "Bothnian Sea",
                                                              "Gulf of Finland", "Archipelago Sea",
                                                              "Baltic Proper"))

# Plotting
setwd("~/GitHub/PlanktonRUE-Project/Figures/VP_plots")

vp_RUEP_region <- ggplot(data_long, aes(x = rows, y = value, fill = rows)) +
  geom_violin(position="dodge", trim=F, scale="width", draw_quantiles=c(0.25,0.5,0.75), linewidth = 0.25) +
  scale_fill_manual(values=c('royalblue3','#278B9AFF','#1D2645FF','#44A57CFF','#819A7AFF','#58A449FF','#CEC917FF','violetred3','salmon2','gold3')) +
  facet_grid(condition ~ .) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  labs(title = "Phytoplankton RUE per P: Conditional Variance Partitioning by Region")
vp_RUEP_region

ggsave(filename = "RUEP_vpcond_region.png",
       plot = vp_RUEP_region,
       width = 7,
       height =10,
       dpi = 300)

#Phytoplankton RUE per N...........................................................
setwd("~/GitHub/PlanktonRUE-Project/BRM/RDATA")
load('BRM_phytoRUEN_CD_PI_260311.RDATA')

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
              "MaxDepth",  "IMaxDepthE2"),  ## 
  Phyto_div = c("genus_richness", "FEve", "FDis"),
  Zoo_div = c("species_richness", "zoo_FEve", "zoo_FDis"),
  Cyano_dom = c("cyano_prop"),
  Abiotic_phyto_div = c("Temperature:genus_richness", "Salinity:genus_richness",
                        "Temperature:FEve", "Salinity:FEve",
                        "Temperature:FDis", "Salinity:FDis"),
  Abiotic_cyano = c("Temperature:cyano_prop", "Salinity:cyano_prop"),
  Abiotic_zoo_div = c("Temperature:species_richness", "Salinity:species_richness",
                      "Temperature:zoo_FEve", "Salinity:zoo_FEve",
                      "Temperature:zoo_FDis", "Salinity:zoo_FDis"),
  DOY = c('DOY', 'IDOYE2',
          'daylight_hours','Idaylight_hoursE2'),
  Site = c('Site_Id'),
  Year = c('Year')#,
  #Residuals = 'residuals'
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
K_V <- apply(params, 'iterations', conditional_lin_term_cov, G, g, X, y)

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

# --- Prepare data for plotting ---
dim3 <- apply(V_cond, 3L, c)
dim12 <- expand.grid(dimnames(V_cond)[1:2])
data_wide <- data.frame(dim3, dim12)
data_long <- melt(setDT(data_wide), id.vars = c('rows','condition'), variable.name = 'iteration')

# Keep order of regions
data_long$condition <- factor(data_long$condition, levels = c("Bothnian Bay", "Bothnian Sea",
                                                              "Gulf of Finland", "Archipelago Sea",
                                                              "Baltic Proper"))

# Plotting
setwd("~/GitHub/PlanktonRUE-Project/Figures/VP_plots")

vp_RUEN_region <- ggplot(data_long, aes(x = rows, y = value, fill = rows)) +
  geom_violin(position="dodge", trim=F, scale="width", draw_quantiles=c(0.25,0.5,0.75), linewidth = 0.25) +
  scale_fill_manual(values=c('royalblue3','#278B9AFF','#1D2645FF','#44A57CFF','#819A7AFF','#58A449FF','#CEC917FF','violetred3','salmon2','gold3')) +
  facet_grid(condition ~ .) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  labs(title = "Phytoplankton RUE per N: Conditional Variance Partitioning by Region")
vp_RUEN_region

ggsave(filename = "RUEN_vpcond_region.png",
       plot = vp_RUEN_region,
       width = 7,
       height =10,
       dpi = 300)

#Zooplankton RUE ...........................................................
setwd("~/GitHub/PlanktonRUE-Project/BRM/RDATA")
load('BRM_zoo_RUE_CD_PI_260313.RDATA')

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
              "Zooplankton_MaxDepth",  "IZooplankton_MaxDepthE2"),  ## 
  Phyto_div = c("genus_richness", "FEve_selected", "FDis_selected"),
  Zoo_div = c("species_richness", "zoo_FEve", "zoo_FDis"),
  Cyano_dom = c("cyano_prop"),
  Abiotic_phyto_div = c("Temperature:genus_richness", "Salinity:genus_richness",
                        "Temperature:FEve_selected", "Salinity:FEve_selected",
                        "Temperature:FDis_selected", "Salinity:FDis_selected"),
  Abiotic_cyano = c("Temperature:cyano_prop", "Salinity:cyano_prop"),
  Abiotic_zoo_div = c("Temperature:species_richness", "Salinity:species_richness",
                      "Temperature:zoo_FEve", "Salinity:zoo_FEve",
                      "Temperature:zoo_FDis", "Salinity:zoo_FDis"),
  DOY = c('DOY', 'IDOYE2',
          'daylight_hours','Idaylight_hoursE2'),
  Site = c('Site_Id'),
  Year = c('Year')#,
  #Residuals = 'residuals'
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
K_V <- apply(params, 'iterations', conditional_lin_term_cov, G, g, X, y)

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

# --- Prepare data for plotting ---
dim3 <- apply(V_cond, 3L, c)
dim12 <- expand.grid(dimnames(V_cond)[1:2])
data_wide <- data.frame(dim3, dim12)
data_long <- melt(setDT(data_wide), id.vars = c('rows','condition'), variable.name = 'iteration')

# Keep order of regions
data_long$condition <- factor(data_long$condition, levels = c("Bothnian Bay", "Bothnian Sea",
                                                              "Gulf of Finland", "Archipelago Sea",
                                                              "Baltic Proper"))

# Plotting
setwd("~/GitHub/PlanktonRUE-Project/Figures/VP_plots")

vp_zooRUE_region <- ggplot(data_long, aes(x = rows, y = value, fill = rows)) +
  geom_violin(position="dodge", trim=F, scale="width", draw_quantiles=c(0.25,0.5,0.75), linewidth = 0.25) +
  scale_fill_manual(values=c('royalblue3','#278B9AFF','#1D2645FF','#44A57CFF','#819A7AFF','#58A449FF','#CEC917FF','violetred3','salmon2','gold3')) +
  facet_grid(condition ~ .) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  labs(title = "Zooplankton RUE: Conditional Variance Partitioning by Region")
vp_zooRUE_region

ggsave(filename = "zooRUE_vpcond_region.png",
       plot = vp_zooRUE_region,
       width = 7,
       height =10,
       dpi = 300)

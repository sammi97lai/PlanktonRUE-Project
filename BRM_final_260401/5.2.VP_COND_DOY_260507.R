##############################################################################
################ Conditional Variance partitioning by DOY ####################
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
library(patchwork)
library(vioplot)
library(ggpubr)
library(maps)

rm(list = ls())

##...........................................................................................................

## Load variance functions modified from Torsti's
setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401")
source('Functions_for_brms_VP_Ap26.R')

#Phytoplankton RUE per P...........................................................
setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA")
load('BRM_phytoRUEP_CD_PI_260408.RDATA')

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
  DOY = c('sDOY_1', 'sDOY'),
  Site = c('Site_Id'),
  Year = c('Year')#,
  # Residuals = 'residuals'
) #group variables into categories#

# define some needed parameters
linear_terms <- unname(unlist(groups)) #flattened vector of all variable names
n_lt <- length(linear_terms)
n_group <- length(groups)

# --- Define the conditional grouping (DOY) ---
setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/RDATA")
load('plankton_RUE_data7_260408.RDATA')

Data_prep <- model$data
Data_prep$season <- NA
mu <- mean(plankton_RUE_data7$DOY, na.rm = TRUE)
sd <- sd(plankton_RUE_data7$DOY, na.rm = TRUE)

doy_182_scaled <- (182 - mu) / sd
doy_243_scaled <- (243 - mu) / sd

Data_prep$season <- cut(
  Data_prep$DOY,
  breaks = c(-Inf, doy_182_scaled, doy_243_scaled, Inf),
  labels = c("Early season", "Summer", "Late season")
)

g <- Data_prep$season
G <- unique(g)

n_cond <- length(G)

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

dim(V_cond) #  covar   11    cond  3     mcmc samples 8000

# create a weighting matrix of same dimension as you data (particularly good if you have also multiple species but here it is ok still)
varr <- aperm(array(R2_model[,1], dim = c(8000L, 3L, 11L)), perm = c(3L, 2L, 1L)) # check the dimensions carefully 11 covar (1L) 3 conditions (2L) for 8000 samples (3L)
dim(varr)  # [1]   11   3 8000

V_condweighted = V_cond*varr

# --- Prepare data for plotting ---
dim3W <- apply(V_condweighted, 3L, c)
dim12W <- expand.grid(dimnames(V_condweighted)[1:2])
dataW_wide <- data.frame(dim3W, dim12W)
dataW_long <- melt(setDT(dataW_wide), id.vars = c('rows','condition'), variable.name = 'iteration')

# Keep order of regions
dataW_long$condition <- factor(dataW_long$condition, levels = c("Early season", "Summer", "Late season"))

#get vp summary 
vp_group_DOY_mean <- apply(V_condweighted, c(1, 2), mean)
RUEp_vp_group_DOY_df <- as.data.frame(as.table(vp_group_DOY_mean))
colnames(RUEp_vp_group_DOY_df) <- c("Group", "season", "MeanVP")
RUEp_vp_group_DOY_df

# Plotting
setwd("~/GitHub/PlanktonRUE-Project/Figures/VP_plots")

vp_RUEP_DOYW <- ggplot(dataW_long, aes(x = rows, y = value, fill = rows)) +
  geom_violin(position="dodge", trim=F, scale="width", draw_quantiles=c(0.25,0.5,0.75), linewidth = 0.25) +
  scale_fill_manual(values=c('royalblue3','#278B9AFF','#44A57CFF','#819A7AFF','#58A449FF', '#403369FF','#5C5992FF','#AE93BEFF','violetred3','salmon2','gold3')) +
  facet_grid(condition ~ .) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        axis.title.x = element_blank(),   # remove x-axis title
        axis.title.y = element_blank(),   # remove y-axis title
        legend.position = "none",          # remove colour legend
        plot.title = element_blank()
        )+
  labs(title = "Phytoplankton RUE per P: Conditional Variance Partitioning by Season")
vp_RUEP_DOYW

ggsave(filename = "RUEP_vpcond_season_r2.png",
       plot = vp_RUEP_DOYW,
       width = 5,
       height =6,
       dpi = 300)

#Phytoplankton RUE per N ...........................................................
setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA")
load('BRM_phytoRUEN_CD_PI_260408.RDATA')

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
  DOY = c('sDOY_1', 'sDOY'),
  Site = c('Site_Id'),
  Year = c('Year')#,
  # Residuals = 'residuals'
) #group variables into categories#

# define some needed parameters
linear_terms <- unname(unlist(groups)) #flattened vector of all variable names
n_lt <- length(linear_terms)
n_group <- length(groups)

# --- Define the conditional grouping (DOY) ---
setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/RDATA")
load('plankton_RUE_data7_260408.RDATA')

Data_prep <- model$data
Data_prep$season <- NA
mu <- mean(plankton_RUE_data7$DOY, na.rm = TRUE)
sd <- sd(plankton_RUE_data7$DOY, na.rm = TRUE)

doy_182_scaled <- (182 - mu) / sd
doy_243_scaled <- (243 - mu) / sd

Data_prep$season <- cut(
  Data_prep$DOY,
  breaks = c(-Inf, doy_182_scaled, doy_243_scaled, Inf),
  labels = c("Early season", "Summer", "Late season")
)

g <- Data_prep$season
G <- unique(g)

n_cond <- length(G)

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

dim(V_cond) #  covar   11    cond  3     mcmc samples 8000

# create a weighting matrix of same dimension as you data (particularly good if you have also multiple species but here it is ok still)
varr <- aperm(array(R2_model[,1], dim = c(8000L, 3L, 11L)), perm = c(3L, 2L, 1L)) # check the dimensions carefully 11 covar (1L) 3 conditions (2L) for 8000 samples (3L)
dim(varr)  # [1]   11   3 8000

V_condweighted = V_cond*varr

# --- Prepare data for plotting ---
dim3W <- apply(V_condweighted, 3L, c)
dim12W <- expand.grid(dimnames(V_condweighted)[1:2])
dataW_wide <- data.frame(dim3W, dim12W)
dataW_long <- melt(setDT(dataW_wide), id.vars = c('rows','condition'), variable.name = 'iteration')

# Keep order of regions
dataW_long$condition <- factor(dataW_long$condition, levels = c("Early season", "Summer", "Late season"))

#get vp summary 
vp_group_DOY_mean <- apply(V_condweighted, c(1, 2), mean)
RUEn_vp_group_DOY_df <- as.data.frame(as.table(vp_group_DOY_mean))
colnames(RUEn_vp_group_DOY_df) <- c("Group", "season", "MeanVP")
RUEn_vp_group_DOY_df

# Plotting
setwd("~/GitHub/PlanktonRUE-Project/Figures/VP_plots")

vp_RUEN_DOYW <- ggplot(dataW_long, aes(x = rows, y = value, fill = rows)) +
  geom_violin(position="dodge", trim=F, scale="width", draw_quantiles=c(0.25,0.5,0.75), linewidth = 0.25) +
  scale_fill_manual(values=c('royalblue3','#278B9AFF','#44A57CFF','#819A7AFF','#58A449FF', '#403369FF','#5C5992FF','#AE93BEFF','violetred3','salmon2','gold3')) +
  facet_grid(condition ~ .) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        axis.title.x = element_blank(),   # remove x-axis title
        axis.title.y = element_blank(),   # remove y-axis title
        legend.position = "none",          # remove colour legend
        plot.title = element_blank()
        )+
  labs(title = "Phytoplankton RUE per N: Conditional Variance Partitioning by Season")
vp_RUEN_DOYW

ggsave(filename = "RUEN_vpcond_season_r2.png",
       plot = vp_RUEN_DOYW,
       width = 5,
       height =6,
       dpi = 300)

# Zooplankton RUE ...............................................................
setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA")
load('BRM_zoo_RUE_CD_PI_260427.RDATA')

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
  DOY = c('sDOY_1', 'sDOY'),
  Site = c('Site_Id'),
  Year = c('Year')#,
  # Residuals = 'residuals'
) #group variables into categories#


# define some needed parameters
linear_terms <- unname(unlist(groups)) #flattened vector of all variable names
n_lt <- length(linear_terms)
n_group <- length(groups)

# --- Define the conditional grouping (DOY) ---
setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/RDATA")
load('plankton_RUE_data7_260408.RDATA')

Data_prep <- model$data
Data_prep$season <- NA
mu <- mean(plankton_RUE_data7$DOY, na.rm = TRUE)
sd <- sd(plankton_RUE_data7$DOY, na.rm = TRUE)

doy_182_scaled <- (182 - mu) / sd
doy_243_scaled <- (243 - mu) / sd

Data_prep$season <- cut(
  Data_prep$DOY,
  breaks = c(-Inf, doy_182_scaled, doy_243_scaled, Inf),
  labels = c("Early season", "Summer", "Late season")
)

g <- Data_prep$season
G <- unique(g)

n_cond <- length(G)

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

dim(V_cond) #  covar   11    cond  3     mcmc samples 8000

# create a weighting matrix of same dimension as you data (particularly good if you have also multiple species but here it is ok still)
varr <- aperm(array(R2_model[,1], dim = c(8000L, 3L, 11L)), perm = c(3L, 2L, 1L)) # check the dimensions carefully 11 covar (1L) 3 conditions (2L) for 8000 samples (3L)
dim(varr)  # [1]   11   3 8000

V_condweighted = V_cond*varr

# --- Prepare data for plotting ---
dim3W <- apply(V_condweighted, 3L, c)
dim12W <- expand.grid(dimnames(V_condweighted)[1:2])
dataW_wide <- data.frame(dim3W, dim12W)
dataW_long <- melt(setDT(dataW_wide), id.vars = c('rows','condition'), variable.name = 'iteration')

# Keep order of regions
dataW_long$condition <- factor(dataW_long$condition, levels = c("Early season", "Summer", "Late season"))

#get vp summary 
vp_group_DOY_mean <- apply(V_condweighted, c(1, 2), mean)
zoo_vp_group_DOY_df <- as.data.frame(as.table(vp_group_DOY_mean))
colnames(zoo_vp_group_DOY_df) <- c("Group", "season", "MeanVP")
zoo_vp_group_DOY_df

# Plotting
setwd("~/GitHub/PlanktonRUE-Project/Figures/VP_plots")

vp_zooRUE_DOYW <- ggplot(dataW_long, aes(x = rows, y = value, fill = rows)) +
  geom_violin(position="dodge", trim=F, scale="width", draw_quantiles=c(0.25,0.5,0.75), linewidth = 0.25) +
  scale_fill_manual(values=c('royalblue3','#278B9AFF','#44A57CFF','#819A7AFF','#58A449FF', '#403369FF','#5C5992FF','#AE93BEFF', 'violetred3','salmon2','gold3')) +
  facet_grid(condition ~ .) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        axis.title.x = element_blank(),   # remove x-axis title
        axis.title.y = element_blank(),   # remove y-axis title
        legend.position = "none",          # remove colour legend
        plot.title = element_blank()
        )+
  labs(title = "Zooplankton RUE: Conditional Variance Partitioning by Season")
vp_zooRUE_DOYW

ggsave(filename = "zooRUE_vpcond_season_r2.png",
       plot = vp_zooRUE_DOYW,
       width = 5,
       height =6,
       dpi = 300)

#one big figure.............................................................
combined_plot <- (vp_RUEP_DOYW | vp_RUEN_DOYW | vp_zooRUE_DOYW) +
  plot_annotation(tag_levels = "a")

combined_plot

ggsave("VP_DOY_combined_r2.png",
       plot = combined_plot,
       width = 13, height = 6, dpi = 300)


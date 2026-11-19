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
source('Functions_for_brms_VP_260605.R')

# cleaner x axis labels
label_names <- c(
  Abiotic = "Abiotic",
  Phyto_genus_richness = "Phytoplankton\nGenus\nRichness",
  Phyto_evenness = "Phytoplankton\nEvenness",
  Phyto_FEve = "Phytoplankton\nFEve",
  Phyto_FDis = "Phytoplankton\nFDis",
  Cyano_dom = "Cyanobacteria\nDominance",
  Zoo_species_richness = "Zooplankton\nSpecies\nRichness",
  Zoo_evenness = "Zooplankton\nEvenness",
  Zoo_FEve = "Zooplankton\nFEve",
  Zoo_FDis = "Zooplankton\nFDis",
  Random = "Random Effect"
)

#Phytoplankton RUE per P...........................................................
setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA")
load('BRM_phytoRUEP_CD_PI_261013_All.RDATA')

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
  Phyto_evenness = c("phyto_evenness", "Temperature:phyto_evenness", "Salinity:phyto_evenness", "TotalNitrogen_unfiltered:phyto_evenness"),
  Phyto_FEve = c("phyto_FEve_PRUE", "Temperature:phyto_FEve_PRUE", "Salinity:phyto_FEve_PRUE", "TotalNitrogen_unfiltered:phyto_FEve_PRUE"),
  Phyto_FDis = c("phyto_FDis_PRUE", "Temperature:phyto_FDis_PRUE", "Salinity:phyto_FDis_PRUE", "TotalNitrogen_unfiltered:phyto_FDis_PRUE"),
  Cyano_dom = c("cyano_prop", "Temperature:cyano_prop", "Salinity:cyano_prop", "TotalNitrogen_unfiltered:cyano_prop"),
  Zoo_species_richness = c("species_richness", "Temperature:species_richness", "Salinity:species_richness", "TotalNitrogen_unfiltered:species_richness"),
  Zoo_evenness = c("zoo_evenness", "Temperature:zoo_evenness", "Salinity:zoo_evenness", "TotalNitrogen_unfiltered:zoo_evenness"),
  Zoo_FEve = c("zoo_FEve", "Temperature:zoo_FEve", "Salinity:zoo_FEve", "TotalNitrogen_unfiltered:zoo_FEve"),
  Zoo_FDis = c("zoo_FDis", "Temperature:zoo_FDis", "Salinity:zoo_FDis", "TotalNitrogen_unfiltered:zoo_FDis"),
  Random = c('Site_Id', 'Year')#,
  # Residuals = 'residuals'
) #group variables into categories#

# Groups used ONLY for the detailed table
groups_table <- list(
  Temperature = c("Temperature", "ITemperatureE2"),
  Salinity = c("Salinity", "ISalinityE2"),
  pH = c("pH", "IpHE2"),
  Total_Nitrogen = c("TotalNitrogen_unfiltered", "ITotalNitrogen_unfilteredE2"),
  MaxDepth = c("MaxDepth", "IMaxDepthE2"),
  Daylight_hours = c("daylight_hours", "Idaylight_hoursE2"),
  Phyto_genus_richness = c("genus_richness", "Temperature:genus_richness", "Salinity:genus_richness", "TotalNitrogen_unfiltered:genus_richness"),
  Phyto_evenness = c("phyto_evenness", "Temperature:phyto_evenness", "Salinity:phyto_evenness", "TotalNitrogen_unfiltered:phyto_evenness" ),
  Phyto_FEve = c("phyto_FEve_PRUE", "Temperature:phyto_FEve_PRUE",  "Salinity:phyto_FEve_PRUE", "TotalNitrogen_unfiltered:phyto_FEve_PRUE"),
  Phyto_FDis = c("phyto_FDis_PRUE", "Temperature:phyto_FDis_PRUE", "Salinity:phyto_FDis_PRUE", "TotalNitrogen_unfiltered:phyto_FDis_PRUE"),
  Cyano_dom = c("cyano_prop", "Temperature:cyano_prop", "Salinity:cyano_prop", "TotalNitrogen_unfiltered:cyano_prop"),
  Zoo_species_richness = c("species_richness", "Temperature:species_richness", "Salinity:species_richness", "TotalNitrogen_unfiltered:species_richness"),
  oo_evenness = c("zoo_evenness", "Temperature:zoo_evenness", "Salinity:zoo_evenness", "TotalNitrogen_unfiltered:zoo_evenness"),
  Zoo_FEve = c("zoo_FEve", "Temperature:zoo_FEve", "Salinity:zoo_FEve", "TotalNitrogen_unfiltered:zoo_FEve"),
  Zoo_FDis = c("zoo_FDis", "Temperature:zoo_FDis", "Salinity:zoo_FDis", "TotalNitrogen_unfiltered:zoo_FDis"),
  Site = "Site_Id",
  Year = "Year"
)

# define some needed parameters
linear_terms <- unname(unlist(groups)) #flattened vector of all variable names
n_lt <- length(linear_terms)
n_group <- length(groups)

# --- Define the conditional grouping (DOY) ---
setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/RDATA")
load('plankton_RUE_data8_260916.RDATA')

Data_prep <- model$data
Data_prep$season <- NA

Data_prep$DOY <- plankton_RUE_data8$DOY

Data_prep$season <- cut(
  Data_prep$DOY,
  breaks = c(-Inf, 182, 243, Inf),
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

# Detailed VP table
# Map individual terms to the detailed table groups
B_table <- sapply(
  groups_table,
  function(x) as.integer(linear_terms %in% x)
)

dimnames(B_table) <- list(
  linear_terms = linear_terms,
  groups = names(groups_table)
)

n_table_groups <- length(groups_table)

# Grouped variance proportions:
# group x region x iteration
V_table <- array(
  NA_real_,
  dim = c(n_table_groups, n_cond, n_iter),
  dimnames = list(
    Group = names(groups_table),
    Region = G,
    Iteration = NULL
  )
)

for (r in seq_len(n_cond)) {
  for (i in seq_len(n_iter)) {
    
    K <- K_V[r, , , i]
    
    # Combine terms according to the TABLE grouping
    K_group <- t(B_table) %*% K %*% B_table
    
    V_table[, r, i] <- diag(K_group) / sum(diag(K_group))
  }
}

# Weight by model R², following your existing approach
R2_model <- bayes_R2(model)
R2_draws <- R2_model[, 1]

V_table_weighted <- sweep(
  V_table,
  MARGIN = 3,
  STATS = R2_draws,
  FUN = "*"
)

# Summarise posterior draws
vp_detailed_table <- expand.grid(
  Group = names(groups_table),
  Region = G,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

vp_detailed_table$MeanVP <- c(
  apply(V_table_weighted, c(1, 2), mean, na.rm = TRUE)
)

vp_detailed_table$PosteriorSD <- c(
  apply(V_table_weighted, c(1, 2), sd, na.rm = TRUE)
)

vp_detailed_table$Lower95 <- c(
  apply(
    V_table_weighted, c(1, 2),
    quantile, probs = 0.025, na.rm = TRUE
  )
)

vp_detailed_table$Upper95 <- c(
  apply(
    V_table_weighted, c(1, 2),
    quantile, probs = 0.975, na.rm = TRUE
  )
)

# Preserve the order specified in groups_table
vp_detailed_table <- vp_detailed_table %>%
  mutate(
    Group = factor(Group, levels = names(groups_table)),
    MeanVP = round(MeanVP, 4),
    PosteriorSD = round(PosteriorSD, 4),
    Lower95 = round(Lower95, 4),
    Upper95 = round(Upper95, 4)
  ) %>%
  arrange(Region, Group)

# Print table
options(na.print = "NA")
print(as.data.frame(vp_detailed_table), row.names = FALSE)

# Export table
setwd("~/GitHub/PlanktonRUE-Project/Results")
write.csv(
  vp_detailed_table,
  "DOY_VP_RUEp.csv",
  row.names = FALSE
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
vp_group_DOY_summary <- data.frame(
  expand.grid(
    Group = dimnames(V_condweighted)[[1]],
    season = dimnames(V_condweighted)[[2]]
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

vp_RUEP_DOYW <- ggplot(dataW_long, aes(x = rows, y = value, fill = rows)) +
  geom_violin(position="dodge", trim=F, scale="width", draw_quantiles=c(0.25,0.5,0.75), linewidth = 0.25) +
  scale_x_discrete(labels = label_names, guide = guide_axis(n.dodge = 2)) +
  scale_fill_manual(values=c("royalblue3",
                             "#278B9AFF", "#274637FF", "#44A57CFF",
                             "#819A7AFF", "#58A449FF", "#403369FF",
                             "#9F99B5FF", "#5C5992FF", "#AE93BEFF",
                             "gold3")) +
  facet_grid(condition ~ .) +
  coord_cartesian(ylim = c(0, 0.8)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 8, lineheight = 0.85),
    axis.title.x = element_blank(),   # remove x-axis title
    axis.title.y = element_blank(),   # remove y-axis title
    legend.position = "none",          # remove colour legend
    plot.title = element_blank()
  )
vp_RUEP_DOYW

ggsave(filename = "RUEP_vpcond_season_r2.png",
       plot = vp_RUEP_DOYW,
       width = 5,
       height =6,
       dpi = 300)

#Phytoplankton RUE per N ...........................................................
setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA")
load('BRM_phytoRUEN_CD_PI_261013_All.RDATA')

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
  Phyto_evenness = c("phyto_evenness", "Temperature:phyto_evenness", "Salinity:phyto_evenness", "TotalNitrogen_unfiltered:phyto_evenness"),
  Phyto_FEve = c("phyto_FEve_PRUE", "Temperature:phyto_FEve_PRUE", "Salinity:phyto_FEve_PRUE", "TotalNitrogen_unfiltered:phyto_FEve_PRUE"),
  Phyto_FDis = c("phyto_FDis_PRUE", "Temperature:phyto_FDis_PRUE", "Salinity:phyto_FDis_PRUE", "TotalNitrogen_unfiltered:phyto_FDis_PRUE"),
  Cyano_dom = c("cyano_prop", "Temperature:cyano_prop", "Salinity:cyano_prop", "TotalNitrogen_unfiltered:cyano_prop"),
  Zoo_species_richness = c("species_richness", "Temperature:species_richness", "Salinity:species_richness", "TotalNitrogen_unfiltered:species_richness"),
  Zoo_evenness = c("zoo_evenness", "Temperature:zoo_evenness", "Salinity:zoo_evenness", "TotalNitrogen_unfiltered:zoo_evenness"),
  Zoo_FEve = c("zoo_FEve", "Temperature:zoo_FEve", "Salinity:zoo_FEve", "TotalNitrogen_unfiltered:zoo_FEve"),
  Zoo_FDis = c("zoo_FDis", "Temperature:zoo_FDis", "Salinity:zoo_FDis", "TotalNitrogen_unfiltered:zoo_FDis"),
  Random = c('Site_Id', 'Year')#,
  # Residuals = 'residuals'
) #group variables into categories#
# Groups used ONLY for the detailed table

groups_table <- list(
  Temperature = c("Temperature", "ITemperatureE2"),
  Salinity = c("Salinity", "ISalinityE2"),
  pH = c("pH", "IpHE2"),
  Total_Nitrogen = c("TotalNitrogen_unfiltered", "ITotalNitrogen_unfilteredE2"),
  MaxDepth = c("MaxDepth", "IMaxDepthE2"),
  Daylight_hours = c("daylight_hours", "Idaylight_hoursE2"),
  Phyto_genus_richness = c("genus_richness", "Temperature:genus_richness", "Salinity:genus_richness", "TotalNitrogen_unfiltered:genus_richness"),
  Phyto_evenness = c("phyto_evenness", "Temperature:phyto_evenness", "Salinity:phyto_evenness", "TotalNitrogen_unfiltered:phyto_evenness" ),
  Phyto_FEve = c("phyto_FEve_PRUE", "Temperature:phyto_FEve_PRUE",  "Salinity:phyto_FEve_PRUE", "TotalNitrogen_unfiltered:phyto_FEve_PRUE"),
  Phyto_FDis = c("phyto_FDis_PRUE", "Temperature:phyto_FDis_PRUE", "Salinity:phyto_FDis_PRUE", "TotalNitrogen_unfiltered:phyto_FDis_PRUE"),
  Cyano_dom = c("cyano_prop", "Temperature:cyano_prop", "Salinity:cyano_prop", "TotalNitrogen_unfiltered:cyano_prop"),
  Zoo_species_richness = c("species_richness", "Temperature:species_richness", "Salinity:species_richness", "TotalNitrogen_unfiltered:species_richness"),
  Zoo_evenness = c("zoo_evenness", "Temperature:zoo_evenness", "Salinity:zoo_evenness", "TotalNitrogen_unfiltered:zoo_evenness"),
  Zoo_FEve = c("zoo_FEve", "Temperature:zoo_FEve", "Salinity:zoo_FEve", "TotalNitrogen_unfiltered:zoo_FEve"),
  Zoo_FDis = c("zoo_FDis", "Temperature:zoo_FDis", "Salinity:zoo_FDis", "TotalNitrogen_unfiltered:zoo_FDis"),
  Site = "Site_Id",
  Year = "Year"
)

# define some needed parameters
linear_terms <- unname(unlist(groups)) #flattened vector of all variable names
n_lt <- length(linear_terms)
n_group <- length(groups)

# --- Define the conditional grouping (DOY) ---
setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/RDATA")
load('plankton_RUE_data8_260916.RDATA')

Data_prep <- model$data
Data_prep$season <- NA

Data_prep$DOY <- plankton_RUE_data8$DOY

Data_prep$season <- cut(
  Data_prep$DOY,
  breaks = c(-Inf, 182, 243, Inf),
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

# Detailed VP table
# Map individual terms to the detailed table groups
B_table <- sapply(
  groups_table,
  function(x) as.integer(linear_terms %in% x)
)

dimnames(B_table) <- list(
  linear_terms = linear_terms,
  groups = names(groups_table)
)

n_table_groups <- length(groups_table)

# Grouped variance proportions:
# group x region x iteration
V_table <- array(
  NA_real_,
  dim = c(n_table_groups, n_cond, n_iter),
  dimnames = list(
    Group = names(groups_table),
    Region = G,
    Iteration = NULL
  )
)

for (r in seq_len(n_cond)) {
  for (i in seq_len(n_iter)) {
    
    K <- K_V[r, , , i]
    
    # Combine terms according to the TABLE grouping
    K_group <- t(B_table) %*% K %*% B_table
    
    V_table[, r, i] <- diag(K_group) / sum(diag(K_group))
  }
}

# Weight by model R², following your existing approach
R2_model <- bayes_R2(model)
R2_draws <- R2_model[, 1]

V_table_weighted <- sweep(
  V_table,
  MARGIN = 3,
  STATS = R2_draws,
  FUN = "*"
)

# Summarise posterior draws
vp_detailed_table <- expand.grid(
  Group = names(groups_table),
  Region = G,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

vp_detailed_table$MeanVP <- c(
  apply(V_table_weighted, c(1, 2), mean, na.rm = TRUE)
)

vp_detailed_table$PosteriorSD <- c(
  apply(V_table_weighted, c(1, 2), sd, na.rm = TRUE)
)

vp_detailed_table$Lower95 <- c(
  apply(
    V_table_weighted, c(1, 2),
    quantile, probs = 0.025, na.rm = TRUE
  )
)

vp_detailed_table$Upper95 <- c(
  apply(
    V_table_weighted, c(1, 2),
    quantile, probs = 0.975, na.rm = TRUE
  )
)

# Preserve the order specified in groups_table
vp_detailed_table <- vp_detailed_table %>%
  mutate(
    Group = factor(Group, levels = names(groups_table)),
    MeanVP = round(MeanVP, 4),
    PosteriorSD = round(PosteriorSD, 4),
    Lower95 = round(Lower95, 4),
    Upper95 = round(Upper95, 4)
  ) %>%
  arrange(Region, Group)

# Print table
options(na.print = "NA")
print(as.data.frame(vp_detailed_table), row.names = FALSE)

# Export table
setwd("~/GitHub/PlanktonRUE-Project/Results")
write.csv(
  vp_detailed_table,
  "DOY_VP_RUEn.csv",
  row.names = FALSE
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
vp_group_DOY_summary <- data.frame(
  expand.grid(
    Group = dimnames(V_condweighted)[[1]],
    season = dimnames(V_condweighted)[[2]]
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

vp_RUEN_DOYW <- ggplot(dataW_long, aes(x = rows, y = value, fill = rows)) +
  geom_violin(position="dodge", trim=F, scale="width", draw_quantiles=c(0.25,0.5,0.75), linewidth = 0.25) +
  scale_x_discrete(labels = label_names, guide = guide_axis(n.dodge = 2)) +
  scale_fill_manual(values=c("royalblue3", "#278B9AFF", "#274637FF", "#44A57CFF", "#819A7AFF", "#58A449FF", "#403369FF", "#9F99B5FF", "#5C5992FF", "#AE93BEFF", "gold3")) +
  facet_grid(condition ~ .) +
  coord_cartesian(ylim = c(0, 0.8)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 8, lineheight = 0.85),
    axis.title.x = element_blank(),   # remove x-axis title
    axis.title.y = element_blank(),   # remove y-axis title
    legend.position = "none",          # remove colour legend
    plot.title = element_blank()
  )
vp_RUEN_DOYW

ggsave(filename = "RUEN_vpcond_season_r2.png",
       plot = vp_RUEN_DOYW,
       width = 5,
       height =6,
       dpi = 300)

# Zooplankton RUE ...............................................................
setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA")
load('BRM_zoo_RUE_CD_PI_261013_All.RDATA')

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
  Phyto_evenness = c("phyto_evenness", "Temperature:phyto_evenness", "Salinity:phyto_evenness", "TotalNitrogen_unfiltered:phyto_evenness"),
  Phyto_FEve = c("phyto_FEve_ZRUE", "Temperature:phyto_FEve_ZRUE", "Salinity:phyto_FEve_ZRUE", "TotalNitrogen_unfiltered:phyto_FEve_ZRUE"),
  Phyto_FDis = c("phyto_FDis_ZRUE", "Temperature:phyto_FDis_ZRUE", "Salinity:phyto_FDis_ZRUE", "TotalNitrogen_unfiltered:phyto_FDis_ZRUE"),
  Cyano_dom = c("cyano_prop", "Temperature:cyano_prop", "Salinity:cyano_prop", "TotalNitrogen_unfiltered:cyano_prop"),
  Zoo_species_richness = c("species_richness", "Temperature:species_richness", "Salinity:species_richness", "TotalNitrogen_unfiltered:species_richness"),
  Zoo_evenness = c("zoo_evenness", "Temperature:zoo_evenness", "Salinity:zoo_evenness", "TotalNitrogen_unfiltered:zoo_evenness"),
  Zoo_FEve = c("zoo_FEve", "Temperature:zoo_FEve", "Salinity:zoo_FEve", "TotalNitrogen_unfiltered:zoo_FEve"),
  Zoo_FDis = c("zoo_FDis", "Temperature:zoo_FDis", "Salinity:zoo_FDis", "TotalNitrogen_unfiltered:zoo_FDis"),
  Random = c('Site_Id', 'Year')#,
  # Residuals = 'residuals'
) #group variables into categories#

groups_table <- list(
  Temperature = c("Temperature", "ITemperatureE2"),
  Salinity = c("Salinity", "ISalinityE2"),
  pH = c("pH", "IpHE2"),
  Total_Nitrogen = c("TotalNitrogen_unfiltered", "ITotalNitrogen_unfilteredE2"),
  MaxDepth = c("MaxDepth", "IMaxDepthE2"),
  Daylight_hours = c("daylight_hours", "Idaylight_hoursE2"),
  Phyto_genus_richness = c("genus_richness", "Temperature:genus_richness", "Salinity:genus_richness", "TotalNitrogen_unfiltered:genus_richness"),
  Phyto_evenness = c("phyto_evenness", "Temperature:phyto_evenness", "Salinity:phyto_evenness", "TotalNitrogen_unfiltered:phyto_evenness" ),
  Phyto_FEve = c("phyto_FEve_ZRUE", "Temperature:phyto_FEve_ZRUE", "Salinity:phyto_FEve_ZRUE", "TotalNitrogen_unfiltered:phyto_FEve_ZRUE"),
  Phyto_FDis = c("phyto_FDis_ZRUE", "Temperature:phyto_FDis_ZRUE", "Salinity:phyto_FDis_ZRUE", "TotalNitrogen_unfiltered:phyto_FDis_ZRUE"),
  Cyano_dom = c("cyano_prop", "Temperature:cyano_prop", "Salinity:cyano_prop", "TotalNitrogen_unfiltered:cyano_prop"),
  Zoo_species_richness = c("species_richness", "Temperature:species_richness", "Salinity:species_richness", "TotalNitrogen_unfiltered:species_richness"),
  Zoo_evenness = c("zoo_evenness", "Temperature:zoo_evenness", "Salinity:zoo_evenness", "TotalNitrogen_unfiltered:zoo_evenness"),
  Zoo_FEve = c("zoo_FEve", "Temperature:zoo_FEve", "Salinity:zoo_FEve", "TotalNitrogen_unfiltered:zoo_FEve"),
  Zoo_FDis = c("zoo_FDis", "Temperature:zoo_FDis", "Salinity:zoo_FDis", "TotalNitrogen_unfiltered:zoo_FDis"),
  Site = "Site_Id",
  Year = "Year"
)


# define some needed parameters
linear_terms <- unname(unlist(groups)) #flattened vector of all variable names
n_lt <- length(linear_terms)
n_group <- length(groups)

# --- Define the conditional grouping (DOY) ---
setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/RDATA")
load('plankton_RUE_data8_260916.RDATA')

Data_prep <- model$data
Data_prep$season <- NA

Data_prep$DOY <- plankton_RUE_data8$DOY

Data_prep$season <- cut(
  Data_prep$DOY,
  breaks = c(-Inf, 182, 243, Inf),
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

# Detailed VP table
# Map individual terms to the detailed table groups
B_table <- sapply(
  groups_table,
  function(x) as.integer(linear_terms %in% x)
)

dimnames(B_table) <- list(
  linear_terms = linear_terms,
  groups = names(groups_table)
)

n_table_groups <- length(groups_table)

# Grouped variance proportions:
# group x region x iteration
V_table <- array(
  NA_real_,
  dim = c(n_table_groups, n_cond, n_iter),
  dimnames = list(
    Group = names(groups_table),
    Region = G,
    Iteration = NULL
  )
)

for (r in seq_len(n_cond)) {
  for (i in seq_len(n_iter)) {
    
    K <- K_V[r, , , i]
    
    # Combine terms according to the TABLE grouping
    K_group <- t(B_table) %*% K %*% B_table
    
    V_table[, r, i] <- diag(K_group) / sum(diag(K_group))
  }
}

# Weight by model R², following your existing approach
R2_model <- bayes_R2(model)
R2_draws <- R2_model[, 1]

V_table_weighted <- sweep(
  V_table,
  MARGIN = 3,
  STATS = R2_draws,
  FUN = "*"
)

# Summarise posterior draws
vp_detailed_table <- expand.grid(
  Group = names(groups_table),
  Region = G,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

vp_detailed_table$MeanVP <- c(
  apply(V_table_weighted, c(1, 2), mean, na.rm = TRUE)
)

vp_detailed_table$PosteriorSD <- c(
  apply(V_table_weighted, c(1, 2), sd, na.rm = TRUE)
)

vp_detailed_table$Lower95 <- c(
  apply(
    V_table_weighted, c(1, 2),
    quantile, probs = 0.025, na.rm = TRUE
  )
)

vp_detailed_table$Upper95 <- c(
  apply(
    V_table_weighted, c(1, 2),
    quantile, probs = 0.975, na.rm = TRUE
  )
)

# Preserve the order specified in groups_table
vp_detailed_table <- vp_detailed_table %>%
  mutate(
    Group = factor(Group, levels = names(groups_table)),
    MeanVP = round(MeanVP, 4),
    PosteriorSD = round(PosteriorSD, 4),
    Lower95 = round(Lower95, 4),
    Upper95 = round(Upper95, 4)
  ) %>%
  arrange(Region, Group)

# Print table
options(na.print = "NA")
print(as.data.frame(vp_detailed_table), row.names = FALSE)

# Export table
setwd("~/GitHub/PlanktonRUE-Project/Results")
write.csv(
  vp_detailed_table,
  "DOY_VP_zooRUE.csv",
  row.names = FALSE
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
vp_group_DOY_summary <- data.frame(
  expand.grid(
    Group = dimnames(V_condweighted)[[1]],
    season = dimnames(V_condweighted)[[2]]
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

vp_zooRUE_DOYW <- ggplot(dataW_long, aes(x = rows, y = value, fill = rows)) +
  geom_violin(position="dodge", trim=F, scale="width", draw_quantiles=c(0.25,0.5,0.75), linewidth = 0.25) +
  scale_x_discrete(labels = label_names, guide = guide_axis(n.dodge = 2)) +
  scale_fill_manual(values=c("royalblue3", "#278B9AFF", "#274637FF", "#44A57CFF", "#819A7AFF", "#58A449FF", "#403369FF", "#9F99B5FF", "#5C5992FF", "#AE93BEFF", "gold3")) +
  facet_grid(condition ~ .) +
  coord_cartesian(ylim = c(0, 0.8)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 8, lineheight = 0.85),
    axis.title.x = element_blank(),   # remove x-axis title
    axis.title.y = element_blank(),   # remove y-axis title
    legend.position = "none",          # remove colour legend
    plot.title = element_blank()
  )
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
       width = 16, height = 6, dpi = 300)


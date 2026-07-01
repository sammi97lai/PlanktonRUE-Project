#
#               General script to plot results for brms results (variance partition)
#
#............................................................................................

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
library(ghibli)

rm(list = ls())

##...........................................................................................................

##.........................................................

## Get the functions
## ..........................................................................
## ............................................................................

## new and adapted functions
setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401")
source('Functions_for_brms_VP_260605.R')

##specify abiotic factors
abiotic_vars <- c(
  "Temperature",
  "Salinity",
  "pH",
  "Total_Nitrogen",
  "MaxDepth",
  "Daylight_hours"
)

#clear plot labels#
plot_labels <- c(
  Temperature = "Temperature",
  Salinity = "Salinity",
  pH = "pH",
  Total_Nitrogen = "Total Nitrogen\nConcentrations",
  MaxDepth = "Depth",
  Daylight_hours = "Daylight Hours",
  Site = "Site",
  Year = "Year"
)

plot_labels2 <- c(
  Temperature = "Temperature",
  Salinity = "Salinity",
  pH = "pH",
  Total_Nitrogen = "Total Nitrogen\nConcentrations",
  MaxDepth = "Depth",
  Daylight_hours = "Daylight\nHours",
  Phyto_genus_richness = "Phytoplankton\nGenus Richness",
  Phyto_FEve = "Phytoplankton FEve",
  Phyto_FDis = "Phytoplankton FDis",
  Cyano_dom = "Cyanobacteria\nDominance",
  Site = "Site",
  Year = "Year"
)

RUE_plot_labels <- c(
  Temperature = "Temperature",
  Salinity = "Salinity",
  pH = "pH",
  Total_Nitrogen = "Total Nitrogen\nConcentrations",
  MaxDepth = "Depth",
  Daylight_hours = "Daylight Hours",
  Phyto_genus_richness = "Phytoplankton\nGenus Richness",
  Phyto_FEve = "Phytoplankton FEve",
  Phyto_FDis = "Phytoplankton FDis",
  Cyano_dom = "Cyanobacteria\nDominance",
  Zoo_species_richness = "Zooplankton\nSpecies Richness",
  Zoo_FEve = "Zooplankton FEve",
  Zoo_FDis = "Zooplankton FDis",
  Site = "Site",
  Year = "Year"
)

##.....................................................................................................
##.....................................................................................................

##  Cyanobacteria dominance
setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA")
load('BRM_cyanodom_P_260604_Z.RDATA')

## this new function now allow you to extract the quantities in the proper format needed
modcyanodomP = PrepbrmsVP(model)
n_iter <- nrow(modcyanodomP$params)

params = modcyanodomP$params
X = modcyanodomP$X
y = modcyanodomP$y

# select the names of the covariates
colnames(X)
model$formula

## define some covariates to group # this you still need to choose and 
# the name need to match the names in X

groups <- list(
  Temperature = c("polyTemperature2rawEQTRUE1", "polyTemperature2rawEQTRUE2"), 
  Salinity = c("polySalinity2rawEQTRUE1", "polySalinity2rawEQTRUE2"), 
  pH = c("polypH2rawEQTRUE1", "polypH2rawEQTRUE2"),
  Total_Nitrogen = c("polyTotalNitrogen_unfiltered2rawEQTRUE1", "polyTotalNitrogen_unfiltered2rawEQTRUE2"),
  MaxDepth = c("polyMaxDepth2rawEQTRUE1",  "polyMaxDepth2rawEQTRUE2"),
  Daylight_hours = c("polydaylight_hours2rawEQTRUE1","polydaylight_hours2rawEQTRUE2"),  ## 
  Site = c('Site_Id'),
  Year = c('Year')#,
  #Residuals = 'residuals'
) #group variables into categories#

# define some needed parameters
linear_terms <- unname(unlist(groups)) #flattened vector of all variable names
n_lt <- length(linear_terms)
n_group <- length(groups)


### Variance partition $V_x$
#quantify how much variance each term explains relative to total variance

# A # No grouping
K <- apply(params, 'iterations', lin_term_cov, X, y, residuals=FALSE) 
dim(K) <- c(n_lt, n_lt, n_iter)
dimnames(K) <- list(rows = linear_terms, cols = linear_terms, iterations = NULL)

P_x <- apply(K, 'iterations', function(x) {diag(x) / sum(diag(x))})
dimnames(P_x) <- list(variance = linear_terms, iterations = NULL)

#summary
vp_summary <- data.frame(
  Variable = rownames(P_x),
  Mean = rowMeans(P_x)
)

vp_summary

setwd("~/GitHub/PlanktonRUE-Project/Figures/VP_plots")

# png("cyanobacteria_dominance_ng.png", width = 4000, height = 4000, res = 600)
# 
# par(mai = c(3.5, 0.8, 0.8, 0.4), # bottom, left, top, right inches
#     las = 2,
#     cex.axis = 1.1,
#     xpd = NA)
# 
# vp_plot(P_x, nm = linear_terms, main = 'Cyanobacteria Dominance (ng)', las = 2)
# 
# dev.off()

# B # Variance partitioning for groups of linear terms
### Variance partition $V_b$

B <- sapply(groups, function(x) { as.integer(linear_terms %in% x) })
dimnames(B) <- list(linear_terms = linear_terms, groups = colnames(B))

K_B <- apply(K, 'iterations', function(x, B) { t(B) %*% x %*% B }, B)
dim(K_B) <- c(n_group, n_group, n_iter)
dimnames(K_B) <- list(rows = names(groups), cols = names(groups), iterations = NULL)

P <- apply(K_B, 'iterations', function(x) {diag(x) / sum(diag(x))})
dimnames(P) <- list(variance = names(groups), iterations = NULL)

dim(P)
Pmean = apply(P, c(1), mean)
sum(Pmean) # = 1

# now weighted by R2
R2_model = bayes_R2(model)
R2_model[,1]  # we keep only the estimate although the uncertainy can be soemtimes interesting it would be a bit complex and out of subject for you I guess

# create a weighting matrix of same dimension as you data (particularly good if you have also multiple species but here it is ok still)
varr <- aperm(array(R2_model[,1], dim = c(8000L, 8L)), perm = c(2L, 1L)) # check the dimensions carefully 8 groups (1L) for 8000 samples (2L) change number of groups and samples for each model
dim(varr)  # [1]    8 8000

Pweighted = P*varr

Pm_weighted = apply(P*varr, c(1), mean)
dim(Pm_weighted)

P_all = as.matrix(cbind(P1 = Pmean, Pw = Pm_weighted)) # if combined with multiple species

vp_summary_grouped <- data.frame(
  Variable = rownames(Pweighted),
  Mean = rowMeans(Pweighted, na.rm = TRUE),
  Lower95 = apply(Pweighted, 1, quantile, probs = 0.025, na.rm = TRUE),
  Upper95 = apply(Pweighted, 1, quantile, probs = 0.975, na.rm = TRUE)
)

vp_summary_grouped

sum(vp_summary_grouped$Mean[vp_summary_grouped$Variable %in% abiotic_vars]) #get sum of abiotic factors#

#95%CI
abiotic_mat <- Pweighted[
  grepl(
    paste(abiotic_vars, collapse = "|"),
    rownames(Pweighted)
  ),
]
abiotic_sum <- colSums(abiotic_mat)
mean_total <- mean(abiotic_sum)

ci <- quantile(
  abiotic_sum,
  probs = c(0.025, 0.975)
)

data.frame(
  mean = mean_total,
  lower95 = ci[1],
  upper95 = ci[2])

png("cyanobacteria_dominance_r2.png", width = 4000, height = 2800, res = 600)

par(mai = c(1.5, 0.8, 0.5, 0.4), # bottom, left, top, right inches
    las = 2,
    cex.axis = 1.1,
    xpd = NA)

vp_plot(Pweighted, col = c('royalblue3','royalblue3','royalblue3','royalblue3','royalblue3','royalblue3','salmon2','gold3'),
        nm = plot_labels, main = 'Cyanobacteria Dominance', las=2)

dev.off()

##.....................................................................................................
##.....................................................................................................

##  Genus richness
setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA")
load('BRM_phyto_GD_P_260604.RDATA')

## this new function now allow you to extract the quantities in the proper format needed
modphyto_GD_P = PrepbrmsVP(model)
n_iter <- nrow(modphyto_GD_P$params)

params = modphyto_GD_P$params
X = modphyto_GD_P$X
y = modphyto_GD_P$y

# select the names of the covariates
colnames(X)
model$formula

## define some covariates to group # this you still need to choose and 
# the name need to match the names in X

groups <- list(
  Temperature = c("polyTemperature2rawEQTRUE1", "polyTemperature2rawEQTRUE2"), 
  Salinity = c("polySalinity2rawEQTRUE1", "polySalinity2rawEQTRUE2"), 
  pH = c("polypH2rawEQTRUE1", "polypH2rawEQTRUE2"),
  Total_Nitrogen = c("polyTotalNitrogen_unfiltered2rawEQTRUE1", "polyTotalNitrogen_unfiltered2rawEQTRUE2"),
  MaxDepth = c("polyMaxDepth2rawEQTRUE1",  "polyMaxDepth2rawEQTRUE2"),
  Daylight_hours = c("polydaylight_hours2rawEQTRUE1","polydaylight_hours2rawEQTRUE2"),  ## 
  Site = c('Site_Id'),
  Year = c('Year')#,
  #Residuals = 'residuals'
) #group variables into categories#

# define some needed parameters
linear_terms <- unname(unlist(groups)) #flattened vector of all variable names
n_lt <- length(linear_terms)
n_group <- length(groups)


### Variance partition $V_x$
#quantify how much variance each term explains relative to total variance

# A # No grouping
K <- apply(params, 'iterations', lin_term_cov, X, y, residuals=FALSE) 
dim(K) <- c(n_lt, n_lt, n_iter)
dimnames(K) <- list(rows = linear_terms, cols = linear_terms, iterations = NULL)

P_x <- apply(K, 'iterations', function(x) {diag(x) / sum(diag(x))})
dimnames(P_x) <- list(variance = linear_terms, iterations = NULL)

#summary
vp_summary <- data.frame(
  Variable = rownames(P_x),
  Mean = rowMeans(P_x)
)

vp_summary

setwd("~/GitHub/PlanktonRUE-Project/Figures/VP_plots")

# png("Phyto_GD_ng.png", width = 4000, height = 4000, res = 600)
# 
# par(mai = c(3.5, 0.8, 0.8, 0.4), # bottom, left, top, right inches
#     las = 2,
#     cex.axis = 1.1,
#     xpd = NA)
# 
# vp_plot(P_x, nm = linear_terms, main = 'Phytoplankton Genus Richness (ng)', las = 2)
# 
# dev.off()


# B # Variance partitioning for groups of linear terms
### Variance partition $V_b$

B <- sapply(groups, function(x) { as.integer(linear_terms %in% x) })
dimnames(B) <- list(linear_terms = linear_terms, groups = colnames(B))

K_B <- apply(K, 'iterations', function(x, B) { t(B) %*% x %*% B }, B)
dim(K_B) <- c(n_group, n_group, n_iter)
dimnames(K_B) <- list(rows = names(groups), cols = names(groups), iterations = NULL)

P <- apply(K_B, 'iterations', function(x) {diag(x) / sum(diag(x))})
dimnames(P) <- list(variance = names(groups), iterations = NULL)

dim(P)
Pmean = apply(P, c(1), mean)
sum(Pmean) # = 1

# now weighted by R2
R2_model = bayes_R2(model)
R2_model[,1]  # we keep only the estimate although the uncertainy can be soemtimes interesting it would be a bit complex and out of subject for you I guess

# create a weighting matrix of same dimension as you data (particularly good if you have also multiple species but here it is ok still)
varr <- aperm(array(R2_model[,1], dim = c(8000L, 8L)), perm = c(2L, 1L)) # check the dimensions carefully 8 groups (1L) for 8000 samples (2L) change number of groups and samples for each model
dim(varr)  # [1]    8 8000

Pweighted = P*varr

vp_summary_grouped <- data.frame(
  Variable = rownames(Pweighted),
  Mean = rowMeans(Pweighted, na.rm = TRUE),
  Lower95 = apply(Pweighted, 1, quantile, probs = 0.025, na.rm = TRUE),
  Upper95 = apply(Pweighted, 1, quantile, probs = 0.975, na.rm = TRUE)
)

vp_summary_grouped

sum(vp_summary_grouped$Mean[vp_summary_grouped$Variable %in% abiotic_vars]) #get sum of abiotic factors#

#95%CI
abiotic_mat <- Pweighted[
  grepl(
    paste(abiotic_vars, collapse = "|"),
    rownames(Pweighted)
  ),
]
abiotic_sum <- colSums(abiotic_mat)
mean_total <- mean(abiotic_sum)

ci <- quantile(
  abiotic_sum,
  probs = c(0.025, 0.975)
)

data.frame(
  mean = mean_total,
  lower95 = ci[1],
  upper95 = ci[2])

png("Phyto_GD_r2.png", width = 4000, height = 2800, res = 600)

par(mai = c(1.5, 0.8, 0.5, 0.4), # bottom, left, top, right inches
    las = 2,
    cex.axis = 1.1,
    xpd = NA)

vp_plot(Pweighted, col = c('royalblue3','royalblue3','royalblue3','royalblue3','royalblue3','royalblue3','salmon2','gold3'),
        nm = plot_labels, main = 'Phytoplankton Genus Richness', las=2)

dev.off()

##.....................................................................................................
##.....................................................................................................
##  phyto FEve
setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA")
load('BRM_phyto_FEve_P_260604.RDATA')

## this new function now allow you to extract the quantities in the proper format needed
modphyto_FEve_P = PrepbrmsVP(model)
n_iter <- nrow(modphyto_FEve_P$params)

params = modphyto_FEve_P$params
X = modphyto_FEve_P$X
y = modphyto_FEve_P$y

# select the names of the covariates
colnames(X)
model$formula

## define some covariates to group # this you still need to choose and 
# the name need to match the names in X

groups <- list(
  Temperature = c("polyTemperature2rawEQTRUE1", "polyTemperature2rawEQTRUE2"), 
  Salinity = c("polySalinity2rawEQTRUE1", "polySalinity2rawEQTRUE2"), 
  pH = c("polypH2rawEQTRUE1", "polypH2rawEQTRUE2"),
  Total_Nitrogen = c("polyTotalNitrogen_unfiltered2rawEQTRUE1", "polyTotalNitrogen_unfiltered2rawEQTRUE2"),
  MaxDepth = c("polyMaxDepth2rawEQTRUE1",  "polyMaxDepth2rawEQTRUE2"),
  Daylight_hours = c("polydaylight_hours2rawEQTRUE1","polydaylight_hours2rawEQTRUE2"),  ## 
  Site = c('Site_Id'),
  Year = c('Year')#,
  #Residuals = 'residuals'
) #group variables into categories#

# define some needed parameters
linear_terms <- unname(unlist(groups)) #flattened vector of all variable names
n_lt <- length(linear_terms)
n_group <- length(groups)

### Variance partition $V_x$
#quantify how much variance each term explains relative to total variance

# A # No grouping
K <- apply(params, 'iterations', lin_term_cov, X, y, residuals=FALSE) 
dim(K) <- c(n_lt, n_lt, n_iter)
dimnames(K) <- list(rows = linear_terms, cols = linear_terms, iterations = NULL)

P_x <- apply(K, 'iterations', function(x) {diag(x) / sum(diag(x))})
dimnames(P_x) <- list(variance = linear_terms, iterations = NULL)

#summary
vp_summary <- data.frame(
  Variable = rownames(P_x),
  Mean = rowMeans(P_x)
)

vp_summary

setwd("~/GitHub/PlanktonRUE-Project/Figures/VP_plots")

# png("Phyto_FEve_ng.png", width = 4000, height = 4000, res = 600)
# 
# par(mai = c(3.5, 0.8, 0.8, 0.4), # bottom, left, top, right inches
#     las = 2,
#     cex.axis = 1.1,
#     xpd = NA)
# 
# vp_plot(P_x, nm = linear_terms, main = 'Phytoplankton Functional Evenness (ng)', las = 2)
# 
# dev.off()


# B # Variance partitioning for groups of linear terms
### Variance partition $V_b$

B <- sapply(groups, function(x) { as.integer(linear_terms %in% x) })
dimnames(B) <- list(linear_terms = linear_terms, groups = colnames(B))

K_B <- apply(K, 'iterations', function(x, B) { t(B) %*% x %*% B }, B)
dim(K_B) <- c(n_group, n_group, n_iter)
dimnames(K_B) <- list(rows = names(groups), cols = names(groups), iterations = NULL)

P <- apply(K_B, 'iterations', function(x) {diag(x) / sum(diag(x))})
dimnames(P) <- list(variance = names(groups), iterations = NULL)

dim(P)
Pmean = apply(P, c(1), mean)
sum(Pmean) # = 1

# now weighted by R2
R2_model = bayes_R2(model)
R2_model[,1]  # we keep only the estimate although the uncertainy can be soemtimes interesting it would be a bit complex and out of subject for you I guess

# create a weighting matrix of same dimension as you data (particularly good if you have also multiple species but here it is ok still)
varr <- aperm(array(R2_model[,1], dim = c(8000L, 8L)), perm = c(2L, 1L)) # check the dimensions carefully 8 groups (1L) for 8000 samples (2L) change number of groups and samples for each model
dim(varr)  # [1]    8 8000

Pweighted = P*varr

vp_summary_grouped <- data.frame(
  Variable = rownames(Pweighted),
  Mean = rowMeans(Pweighted, na.rm = TRUE),
  Lower95 = apply(Pweighted, 1, quantile, probs = 0.025, na.rm = TRUE),
  Upper95 = apply(Pweighted, 1, quantile, probs = 0.975, na.rm = TRUE)
)

vp_summary_grouped

sum(vp_summary_grouped$Mean[vp_summary_grouped$Variable %in% abiotic_vars]) #get sum of abiotic factors#

#95%CI
abiotic_mat <- Pweighted[
  grepl(
    paste(abiotic_vars, collapse = "|"),
    rownames(Pweighted)
  ),
]
abiotic_sum <- colSums(abiotic_mat)
mean_total <- mean(abiotic_sum)

ci <- quantile(
  abiotic_sum,
  probs = c(0.025, 0.975)
)

data.frame(
  mean = mean_total,
  lower95 = ci[1],
  upper95 = ci[2])

png("Phyto_FEve_r2.png", width = 4000, height = 2800, res = 600)

par(mai = c(1.5, 0.8, 0.5, 0.4), # bottom, left, top, right inches
    las = 2,
    cex.axis = 1.1,
    xpd = NA)

vp_plot(Pweighted, col = c('royalblue3','royalblue3','royalblue3','royalblue3','royalblue3','royalblue3','salmon2','gold3'),
        nm = plot_labels, main = 'Phytoplankton Functional Evenness', las=2)

dev.off()

##.....................................................................................................
##.....................................................................................................
##  phyto FDis
setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA")
load('BRM_phyto_FDis_P_260604.RDATA')

## this new function now allow you to extract the quantities in the proper format needed
modphyto_FDis_P = PrepbrmsVP(model)
n_iter <- nrow(modphyto_FDis_P$params)

params = modphyto_FDis_P$params
X = modphyto_FDis_P$X
y = modphyto_FDis_P$y

# select the names of the covariates
colnames(X)
model$formula

## define some covariates to group # this you still need to choose and 
# the name need to match the names in X

groups <- list(
  Temperature = c("polyTemperature2rawEQTRUE1", "polyTemperature2rawEQTRUE2"), 
  Salinity = c("polySalinity2rawEQTRUE1", "polySalinity2rawEQTRUE2"), 
  pH = c("polypH2rawEQTRUE1", "polypH2rawEQTRUE2"),
  Total_Nitrogen = c("polyTotalNitrogen_unfiltered2rawEQTRUE1", "polyTotalNitrogen_unfiltered2rawEQTRUE2"),
  MaxDepth = c("polyMaxDepth2rawEQTRUE1",  "polyMaxDepth2rawEQTRUE2"),
  Daylight_hours = c("polydaylight_hours2rawEQTRUE1","polydaylight_hours2rawEQTRUE2"),  ## 
  Site = c('Site_Id'),
  Year = c('Year')#,
  #Residuals = 'residuals'
) #group variables into categories#

# define some needed parameters
linear_terms <- unname(unlist(groups)) #flattened vector of all variable names
n_lt <- length(linear_terms)
n_group <- length(groups)

### Variance partition $V_x$
#quantify how much variance each term explains relative to total variance

# A # No grouping
K <- apply(params, 'iterations', lin_term_cov, X, y, residuals=FALSE) 
dim(K) <- c(n_lt, n_lt, n_iter)
dimnames(K) <- list(rows = linear_terms, cols = linear_terms, iterations = NULL)

P_x <- apply(K, 'iterations', function(x) {diag(x) / sum(diag(x))})
dimnames(P_x) <- list(variance = linear_terms, iterations = NULL)

#summary
vp_summary <- data.frame(
  Variable = rownames(P_x),
  Mean = rowMeans(P_x)
)

vp_summary

setwd("~/GitHub/PlanktonRUE-Project/Figures/VP_plots")

# png("Phyto_FDis_ng.png", width = 4000, height = 4000, res = 600)
# 
# par(mai = c(3.5, 0.8, 0.8, 0.4), # bottom, left, top, right inches
#     las = 2,
#     cex.axis = 1.1,
#     xpd = NA)
# 
# vp_plot(P_x, nm = linear_terms, main = 'Phytoplankton Functional Dispersion (ng)', las = 2)
# 
# dev.off()


# B # Variance partitioning for groups of linear terms
### Variance partition $V_b$

B <- sapply(groups, function(x) { as.integer(linear_terms %in% x) })
dimnames(B) <- list(linear_terms = linear_terms, groups = colnames(B))

K_B <- apply(K, 'iterations', function(x, B) { t(B) %*% x %*% B }, B)
dim(K_B) <- c(n_group, n_group, n_iter)
dimnames(K_B) <- list(rows = names(groups), cols = names(groups), iterations = NULL)

P <- apply(K_B, 'iterations', function(x) {diag(x) / sum(diag(x))})
dimnames(P) <- list(variance = names(groups), iterations = NULL)

dim(P)
Pmean = apply(P, c(1), mean)
sum(Pmean) # = 1

# now weighted by R2
R2_model = bayes_R2(model)
R2_model[,1]  # we keep only the estimate although the uncertainy can be soemtimes interesting it would be a bit complex and out of subject for you I guess

# create a weighting matrix of same dimension as you data (particularly good if you have also multiple species but here it is ok still)
varr <- aperm(array(R2_model[,1], dim = c(8000L, 8L)), perm = c(2L, 1L)) # check the dimensions carefully 8 groups (1L) for 8000 samples (2L) change number of groups and samples for each model
dim(varr)  # [1]    8 8000

Pweighted = P*varr

vp_summary_grouped <- data.frame(
  Variable = rownames(Pweighted),
  Mean = rowMeans(Pweighted, na.rm = TRUE),
  Lower95 = apply(Pweighted, 1, quantile, probs = 0.025, na.rm = TRUE),
  Upper95 = apply(Pweighted, 1, quantile, probs = 0.975, na.rm = TRUE)
)

vp_summary_grouped

sum(vp_summary_grouped$Mean[vp_summary_grouped$Variable %in% abiotic_vars]) #get sum of abiotic factors#

#95%CI
abiotic_mat <- Pweighted[
  grepl(
    paste(abiotic_vars, collapse = "|"),
    rownames(Pweighted)
  ),
]
abiotic_sum <- colSums(abiotic_mat)
mean_total <- mean(abiotic_sum)

ci <- quantile(
  abiotic_sum,
  probs = c(0.025, 0.975)
)

data.frame(
  mean = mean_total,
  lower95 = ci[1],
  upper95 = ci[2])

png("Phyto_FDis_r2.png", width = 4000, height = 2800, res = 600)

par(mai = c(1.5, 0.8, 0.5, 0.4), # bottom, left, top, right inches
    las = 2,
    cex.axis = 1.1,
    xpd = NA)

vp_plot(Pweighted, col = c('royalblue3','royalblue3','royalblue3','royalblue3','royalblue3','royalblue3','salmon2','gold3'),
        nm = plot_labels, main = 'Phytoplankton Functional Dispersion', las=2)

dev.off()

##.....................................................................................................
##.....................................................................................................

##  Phyto RUE per P model
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

## define some covariates to group # this you still need to choose and 
# the name need to match the names in X

groups <- list(
  Temperature = c("Temperature", "ITemperatureE2"),  
  Salinity = c("Salinity", "ISalinityE2"), 
  pH = c("pH", "IpHE2"),
  Total_Nitrogen = c("TotalNitrogen_unfiltered", "ITotalNitrogen_unfilteredE2"), 
  MaxDepth = c("MaxDepth",  "IMaxDepthE2"),
  Daylight_hours = c('daylight_hours','Idaylight_hoursE2'),  ## 
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

### Variance partition $V_x$
#quantify how much variance each term explains relative to total variance

# A # No grouping
K <- apply(params, 'iterations', lin_term_cov, X, y, residuals=FALSE) 
dim(K) <- c(n_lt, n_lt, n_iter)
dimnames(K) <- list(rows = linear_terms, cols = linear_terms, iterations = NULL)

P_x <- apply(K, 'iterations', function(x) {diag(x) / sum(diag(x))})
dimnames(P_x) <- list(variance = linear_terms, iterations = NULL)

#summary
vp_summary <- data.frame(
  Variable = rownames(P_x),
  Mean = rowMeans(P_x)
)

vp_summary

setwd("~/GitHub/PlanktonRUE-Project/Figures/VP_plots")

# png("Phyto_RUE_P_ng.png", width = 5000, height = 4000, res = 600)
# 
# par(mai = c(2.5, 0.8, 0.8, 0.4), # bottom, left, top, right inches
#     las = 2,
#     cex.axis = 1.1,
#     xpd = NA)
# 
# vp_plot(P_x, nm = linear_terms, main = 'Phytoplankton RUE per P (ng)', las = 2)
# 
# dev.off()


# B # Variance partitioning for groups of linear terms
### Variance partition $V_b$

B <- sapply(groups, function(x) { as.integer(linear_terms %in% x) })
dimnames(B) <- list(linear_terms = linear_terms, groups = colnames(B))

K_B <- apply(K, 'iterations', function(x, B) { t(B) %*% x %*% B }, B)
dim(K_B) <- c(n_group, n_group, n_iter)
dimnames(K_B) <- list(rows = names(groups), cols = names(groups), iterations = NULL)

P <- apply(K_B, 'iterations', function(x) {diag(x) / sum(diag(x))})
dimnames(P) <- list(variance = names(groups), iterations = NULL)

dim(P)
Pmean = apply(P, c(1), mean)
sum(Pmean) # = 1

# now weighted by R2
R2_model = bayes_R2(model)
R2_model[,1]  # we keep only the estimate although the uncertainy can be soemtimes interesting it would be a bit complex and out of subject for you I guess

# create a weighting matrix of same dimension as you data (particularly good if you have also multiple species but here it is ok still)
varr <- aperm(array(R2_model[,1], dim = c(8000L, 15L)), perm = c(2L, 1L)) # check the dimensions carefully 15 groups (1L) for 8000 samples (2L) change number of groups and samples for each model
dim(varr)  # [1]    15 8000

Pweighted = P*varr

vp_summary_grouped <- data.frame(
  Variable = rownames(Pweighted),
  Mean = rowMeans(Pweighted, na.rm = TRUE),
  Lower95 = apply(Pweighted, 1, quantile, probs = 0.025, na.rm = TRUE),
  Upper95 = apply(Pweighted, 1, quantile, probs = 0.975, na.rm = TRUE)
)

vp_summary_grouped

sum(vp_summary_grouped$Mean[vp_summary_grouped$Variable %in% abiotic_vars]) #get sum of abiotic factors#

#95%CI
abiotic_mat <- Pweighted[
  grepl(
    paste(abiotic_vars, collapse = "|"),
    rownames(Pweighted)
  ),
]
abiotic_sum <- colSums(abiotic_mat)
mean_total <- mean(abiotic_sum)

ci <- quantile(
  abiotic_sum,
  probs = c(0.025, 0.975)
)

data.frame(
  mean = mean_total,
  lower95 = ci[1],
  upper95 = ci[2])

png("Phyto_RUE_P_r2.png", width = 5300, height = 3500, res = 600)

par(mai = c(2.0, 0.8, 0.5, 0.4), # bottom, left, top, right inches
    las = 2,
    cex.axis = 1.1,
    xpd = NA)

vp_plot(Pweighted, col = c('royalblue3','royalblue3','royalblue3','royalblue3','royalblue3','royalblue3','#278B9AFF',
                   '#44A57CFF','#819A7AFF','#58A449FF', '#403369FF','#5C5992FF','#AE93BEFF','salmon2','gold3'),
        nm = RUE_plot_labels, main = 'Phytoplankton RUE per P', las=2)

dev.off()

##.....................................................................................................
##.....................................................................................................


##  Phyto RUE per N model
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

## define some covariates to group # this you still need to choose and 
# the name need to match the names in X

groups <- list(
  Temperature = c("Temperature", "ITemperatureE2"),  
  Salinity = c("Salinity", "ISalinityE2"), 
  pH = c("pH", "IpHE2"),
  Total_Nitrogen = c("TotalNitrogen_unfiltered", "ITotalNitrogen_unfilteredE2"), 
  MaxDepth = c("MaxDepth",  "IMaxDepthE2"),
  Daylight_hours = c('daylight_hours','Idaylight_hoursE2'),  ## 
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

### Variance partition $V_x$
#quantify how much variance each term explains relative to total variance

# A # No grouping
K <- apply(params, 'iterations', lin_term_cov, X, y, residuals=FALSE) 
dim(K) <- c(n_lt, n_lt, n_iter)
dimnames(K) <- list(rows = linear_terms, cols = linear_terms, iterations = NULL)

P_x <- apply(K, 'iterations', function(x) {diag(x) / sum(diag(x))})
dimnames(P_x) <- list(variance = linear_terms, iterations = NULL)

#summary
vp_summary <- data.frame(
  Variable = rownames(P_x),
  Mean = rowMeans(P_x)
)

vp_summary

setwd("~/GitHub/PlanktonRUE-Project/Figures/VP_plots")

# png("Phyto_RUE_N_ng.png", width = 5000, height = 4000, res = 600)
# 
# par(mai = c(2.5, 0.8, 0.8, 0.4), # bottom, left, top, right inches
#     las = 2,
#     cex.axis = 1.1,
#     xpd = NA)
# 
# vp_plot(P_x, nm = linear_terms, main = 'Phytoplankton RUE per N (ng)', las = 2)
# 
# dev.off()

# B # Variance partitioning for groups of linear terms
### Variance partition $V_b$

B <- sapply(groups, function(x) { as.integer(linear_terms %in% x) })
dimnames(B) <- list(linear_terms = linear_terms, groups = colnames(B))

K_B <- apply(K, 'iterations', function(x, B) { t(B) %*% x %*% B }, B)
dim(K_B) <- c(n_group, n_group, n_iter)
dimnames(K_B) <- list(rows = names(groups), cols = names(groups), iterations = NULL)

P <- apply(K_B, 'iterations', function(x) {diag(x) / sum(diag(x))})
dimnames(P) <- list(variance = names(groups), iterations = NULL)

dim(P)
Pmean = apply(P, c(1), mean)
sum(Pmean) # = 1

# now weighted by R2
R2_model = bayes_R2(model)
R2_model[,1]  # we keep only the estimate although the uncertainy can be soemtimes interesting it would be a bit complex and out of subject for you I guess

# create a weighting matrix of same dimension as you data (particularly good if you have also multiple species but here it is ok still)
varr <- aperm(array(R2_model[,1], dim = c(8000L, 15L)), perm = c(2L, 1L)) # check the dimensions carefully 15 groups (1L) for 8000 samples (2L) change number of groups and samples for each model
dim(varr)  # [1]    15 8000

Pweighted = P*varr

vp_summary_grouped <- data.frame(
  Variable = rownames(Pweighted),
  Mean = rowMeans(Pweighted, na.rm = TRUE),
  Lower95 = apply(Pweighted, 1, quantile, probs = 0.025, na.rm = TRUE),
  Upper95 = apply(Pweighted, 1, quantile, probs = 0.975, na.rm = TRUE)
)

vp_summary_grouped

sum(vp_summary_grouped$Mean[vp_summary_grouped$Variable %in% abiotic_vars]) #get sum of abiotic factors#

#95%CI
abiotic_mat <- Pweighted[
  grepl(
    paste(abiotic_vars, collapse = "|"),
    rownames(Pweighted)
  ),
]
abiotic_sum <- colSums(abiotic_mat)
mean_total <- mean(abiotic_sum)

ci <- quantile(
  abiotic_sum,
  probs = c(0.025, 0.975)
)

data.frame(
  mean = mean_total,
  lower95 = ci[1],
  upper95 = ci[2])

png("Phyto_RUE_N_r2.png", width = 5300, height = 3500, res = 600)

par(mai = c(2.0, 0.8, 0.5, 0.4), # bottom, left, top, right inches
    las = 2,
    cex.axis = 1.1,
    xpd = NA)

vp_plot(Pweighted, col = c('royalblue3','royalblue3','royalblue3','royalblue3','royalblue3','royalblue3','#278B9AFF',
                   '#44A57CFF','#819A7AFF','#58A449FF', '#403369FF','#5C5992FF','#AE93BEFF','salmon2','gold3'),
        nm = RUE_plot_labels, main = 'Phytoplankton RUE per N', las=2)

dev.off()

##.....................................................................................................
##.....................................................................................................

##  Zoo Species Richness model
setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA")
load('BRM_zoo_SR_CD_P_260604.RDATA')

## this new function now allow you to extract the quantities in the proper format needed
modZoo_SR_CD = PrepbrmsVP(model)
n_iter <- nrow(modZoo_SR_CD$params)

params = modZoo_SR_CD$params
X = modZoo_SR_CD$X
y = modZoo_SR_CD$y

# select the names of the covariates
colnames(X)
model$formula

## define some covariates to group # this you still need to choose and 
# the name need to match the names in X

groups <- list(
  Temperature = c("polyTemperature2rawEQTRUE1", "polyTemperature2rawEQTRUE2"),
  Salinity = c("polySalinity2rawEQTRUE1", "polySalinity2rawEQTRUE2"),
  pH = c("polypH2rawEQTRUE1", "polypH2rawEQTRUE2"),
  Total_Nitrogen = c("polyTotalNitrogen_unfiltered2rawEQTRUE1", "polyTotalNitrogen_unfiltered2rawEQTRUE2"),
  MaxDepth = c("polyZooplankton_MaxDepth2rawEQTRUE1",  "polyZooplankton_MaxDepth2rawEQTRUE2"),
  Daylight_hours = c('polydaylight_hours2rawEQTRUE1','polydaylight_hours2rawEQTRUE2'),  ## 
  Phyto_genus_richness = c("genus_richness"),
  Phyto_FEve = c("FEve"),
  Phyto_FDis = c("FDis"),
  Cyano_dom = c("cyano_prop"),
  Site = c('Site_Id'),
  Year = c('Year')#,
  # Residuals = 'residuals'
) #group variables into categories#

# define some needed parameters
linear_terms <- unname(unlist(groups)) #flattened vector of all variable names
n_lt <- length(linear_terms)
n_group <- length(groups)

### Variance partition $V_x$
#quantify how much variance each term explains relative to total variance

# A # No grouping
K <- apply(params, 'iterations', lin_term_cov, X, y, residuals=FALSE) 
dim(K) <- c(n_lt, n_lt, n_iter)
dimnames(K) <- list(rows = linear_terms, cols = linear_terms, iterations = NULL)

P_x <- apply(K, 'iterations', function(x) {diag(x) / sum(diag(x))})
dimnames(P_x) <- list(variance = linear_terms, iterations = NULL)

#summary
vp_summary <- data.frame(
  Variable = rownames(P_x),
  Mean = rowMeans(P_x)
)

vp_summary

setwd("~/GitHub/PlanktonRUE-Project/Figures/VP_plots")

# png("Zoo_SR_ng.png", width = 4000, height = 4000, res = 600)
# 
# par(mai = c(3.5, 0.8, 0.8, 0.4), # bottom, left, top, right inches
#     las = 2,
#     cex.axis = 1.1,
#     xpd = NA)
# 
# vp_plot(P_x, nm = linear_terms, main = 'Zooplankton Species Richness (ng)', las = 2)
# 
# dev.off()


# B # Variance partitioning for groups of linear terms
### Variance partition $V_b$

B <- sapply(groups, function(x) { as.integer(linear_terms %in% x) })
dimnames(B) <- list(linear_terms = linear_terms, groups = colnames(B))

K_B <- apply(K, 'iterations', function(x, B) { t(B) %*% x %*% B }, B)
dim(K_B) <- c(n_group, n_group, n_iter)
dimnames(K_B) <- list(rows = names(groups), cols = names(groups), iterations = NULL)

P <- apply(K_B, 'iterations', function(x) {diag(x) / sum(diag(x))})
dimnames(P) <- list(variance = names(groups), iterations = NULL)

dim(P)
Pmean = apply(P, c(1), mean)
sum(Pmean) # = 1

# now weighted by R2
R2_model = bayes_R2(model)
R2_model[,1]  # we keep only the estimate although the uncertainy can be soemtimes interesting it would be a bit complex and out of subject for you I guess

# create a weighting matrix of same dimension as you data (particularly good if you have also multiple species but here it is ok still)
varr <- aperm(array(R2_model[,1], dim = c(8000L, 12L)), perm = c(2L, 1L)) # check the dimensions carefully 12 groups (1L) for 8000 samples (2L) change number of groups and samples for each model
dim(varr)  # [1]    12 8000

Pweighted = P*varr

vp_summary_grouped <- data.frame(
  Variable = rownames(Pweighted),
  Mean = rowMeans(Pweighted, na.rm = TRUE),
  Lower95 = apply(Pweighted, 1, quantile, probs = 0.025, na.rm = TRUE),
  Upper95 = apply(Pweighted, 1, quantile, probs = 0.975, na.rm = TRUE)
)

vp_summary_grouped

sum(vp_summary_grouped$Mean[vp_summary_grouped$Variable %in% abiotic_vars]) #get sum of abiotic factors#

#95%CI
abiotic_mat <- Pweighted[
  grepl(
    paste(abiotic_vars, collapse = "|"),
    rownames(Pweighted)
  ),
]
abiotic_sum <- colSums(abiotic_mat)
mean_total <- mean(abiotic_sum)

ci <- quantile(
  abiotic_sum,
  probs = c(0.025, 0.975)
)

data.frame(
  mean = mean_total,
  lower95 = ci[1],
  upper95 = ci[2])

png("Zoo_SR_r2.png", width = 4000, height = 3000, res = 600)

par(mai = c(2.0, 0.8, 0.5, 0.4), # bottom, left, top, right inches
    las = 2,
    cex.axis = 1.1,
    xpd = NA)

vp_plot(Pweighted, col = c('royalblue3','royalblue3','royalblue3','royalblue3','royalblue3','royalblue3','#278B9AFF',
                   '#44A57CFF','#819A7AFF','#58A449FF','salmon2','gold3'),
        nm = plot_labels2, main = 'Zooplankton Species Richness', las=2)

dev.off()

##.....................................................................................................
##.....................................................................................................


##  Zoo Functional Evenness model
setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA")
load('BRM_zoo_FE_CD_P_260604.RDATA')

## this new function now allow you to extract the quantities in the proper format needed
modZoo_FEve_CD = PrepbrmsVP(model)
n_iter <- nrow(modZoo_FEve_CD$params)

params = modZoo_FEve_CD$params
X = modZoo_FEve_CD$X
y = modZoo_FEve_CD$y

# select the names of the covariates
colnames(X)
model$formula

## define some covariates to group # this you still need to choose and 
# the name need to match the names in X

groups <- list(
  Temperature = c("polyTemperature2rawEQTRUE1", "polyTemperature2rawEQTRUE2"),
  Salinity = c("polySalinity2rawEQTRUE1", "polySalinity2rawEQTRUE2"),
  pH = c("polypH2rawEQTRUE1", "polypH2rawEQTRUE2"),
  Total_Nitrogen = c("polyTotalNitrogen_unfiltered2rawEQTRUE1", "polyTotalNitrogen_unfiltered2rawEQTRUE2"),
  MaxDepth = c("polyZooplankton_MaxDepth2rawEQTRUE1",  "polyZooplankton_MaxDepth2rawEQTRUE2"),
  Daylight_hours = c('polydaylight_hours2rawEQTRUE1','polydaylight_hours2rawEQTRUE2'),  ## 
  Phyto_genus_richness = c("genus_richness"),
  Phyto_FEve = c("FEve"),
  Phyto_FDis = c("FDis"),
  Cyano_dom = c("cyano_prop"),
  Site = c('Site_Id'),
  Year = c('Year')#,
  # Residuals = 'residuals'
) #group variables into categories#

# define some needed parameters
linear_terms <- unname(unlist(groups)) #flattened vector of all variable names
n_lt <- length(linear_terms)
n_group <- length(groups)

### Variance partition $V_x$
#quantify how much variance each term explains relative to total variance

# A # No grouping
K <- apply(params, 'iterations', lin_term_cov, X, y, residuals=FALSE) 
dim(K) <- c(n_lt, n_lt, n_iter)
dimnames(K) <- list(rows = linear_terms, cols = linear_terms, iterations = NULL)

P_x <- apply(K, 'iterations', function(x) {diag(x) / sum(diag(x))})
dimnames(P_x) <- list(variance = linear_terms, iterations = NULL)

#summary
vp_summary <- data.frame(
  Variable = rownames(P_x),
  Mean = rowMeans(P_x)
)

vp_summary

setwd("~/GitHub/PlanktonRUE-Project/Figures/VP_plots")

# png("Zoo_FEve_ng.png", width = 4000, height = 4000, res = 600)
# 
# par(mai = c(3.5, 0.8, 0.8, 0.4), # bottom, left, top, right inches
#     las = 2,
#     cex.axis = 1.1,
#     xpd = NA)
# 
# vp_plot(P_x, nm = linear_terms, main = 'Zooplankton Functional Evenness (ng)', las = 2)
# 
# dev.off()


# B # Variance partitioning for groups of linear terms
### Variance partition $V_b$

B <- sapply(groups, function(x) { as.integer(linear_terms %in% x) })
dimnames(B) <- list(linear_terms = linear_terms, groups = colnames(B))

K_B <- apply(K, 'iterations', function(x, B) { t(B) %*% x %*% B }, B)
dim(K_B) <- c(n_group, n_group, n_iter)
dimnames(K_B) <- list(rows = names(groups), cols = names(groups), iterations = NULL)

P <- apply(K_B, 'iterations', function(x) {diag(x) / sum(diag(x))})
dimnames(P) <- list(variance = names(groups), iterations = NULL)

dim(P)
Pmean = apply(P, c(1), mean)
sum(Pmean) # = 1

# now weighted by R2
R2_model = bayes_R2(model)
R2_model[,1]  # we keep only the estimate although the uncertainy can be soemtimes interesting it would be a bit complex and out of subject for you I guess

# create a weighting matrix of same dimension as you data (particularly good if you have also multiple species but here it is ok still)
varr <- aperm(array(R2_model[,1], dim = c(8000L, 12L)), perm = c(2L, 1L)) # check the dimensions carefully 12 groups (1L) for 8000 samples (2L) change number of groups and samples for each model
dim(varr)  # [1]    12 8000

Pweighted = P*varr

vp_summary_grouped <- data.frame(
  Variable = rownames(Pweighted),
  Mean = rowMeans(Pweighted, na.rm = TRUE),
  Lower95 = apply(Pweighted, 1, quantile, probs = 0.025, na.rm = TRUE),
  Upper95 = apply(Pweighted, 1, quantile, probs = 0.975, na.rm = TRUE)
)

vp_summary_grouped

sum(vp_summary_grouped$Mean[vp_summary_grouped$Variable %in% abiotic_vars]) #get sum of abiotic factors#

#95%CI
abiotic_mat <- Pweighted[
  grepl(
    paste(abiotic_vars, collapse = "|"),
    rownames(Pweighted)
  ),
]
abiotic_sum <- colSums(abiotic_mat)
mean_total <- mean(abiotic_sum)

ci <- quantile(
  abiotic_sum,
  probs = c(0.025, 0.975)
)

data.frame(
  mean = mean_total,
  lower95 = ci[1],
  upper95 = ci[2])

png("Zoo_FEve_r2.png", width = 4000, height = 3000, res = 600)

par(mai = c(2.0, 0.8, 0.5, 0.4), # bottom, left, top, right inches
    las = 2,
    cex.axis = 1.1,
    xpd = NA)

vp_plot(Pweighted, col = c('royalblue3','royalblue3','royalblue3','royalblue3','royalblue3','royalblue3','#278B9AFF',
                   '#44A57CFF','#819A7AFF','#58A449FF','salmon2','gold3'),
        nm = plot_labels2, main = 'Zooplankton Functional Evenness', las=2)

dev.off()

##.....................................................................................................
##.....................................................................................................


##  Zoo Functional Dispersion model
setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA")
load('BRM_zoo_FD_CD_P_260604.RDATA')

## this new function now allow you to extract the quantities in the proper format needed
modZoo_FDis_CD = PrepbrmsVP(model)
n_iter <- nrow(modZoo_FDis_CD$params)

params = modZoo_FDis_CD$params
X = modZoo_FDis_CD$X
y = modZoo_FDis_CD$y

# select the names of the covariates
colnames(X)
model$formula

## define some covariates to group # this you still need to choose and 
# the name need to match the names in X

groups <- list(
  Temperature = c("polyTemperature2rawEQTRUE1", "polyTemperature2rawEQTRUE2"),
  Salinity = c("polySalinity2rawEQTRUE1", "polySalinity2rawEQTRUE2"),
  pH = c("polypH2rawEQTRUE1", "polypH2rawEQTRUE2"),
  Total_Nitrogen = c("polyTotalNitrogen_unfiltered2rawEQTRUE1", "polyTotalNitrogen_unfiltered2rawEQTRUE2"),
  MaxDepth = c("polyZooplankton_MaxDepth2rawEQTRUE1",  "polyZooplankton_MaxDepth2rawEQTRUE2"),
  Daylight_hours = c('polydaylight_hours2rawEQTRUE1','polydaylight_hours2rawEQTRUE2'),  ## 
  Phyto_genus_richness = c("genus_richness"),
  Phyto_FEve = c("FEve"),
  Phyto_FDis = c("FDis"),
  Cyano_dom = c("cyano_prop"),
  Site = c('Site_Id'),
  Year = c('Year')#,
  # Residuals = 'residuals'
) #group variables into categories#

# define some needed parameters
linear_terms <- unname(unlist(groups)) #flattened vector of all variable names
n_lt <- length(linear_terms)
n_group <- length(groups)

### Variance partition $V_x$
#quantify how much variance each term explains relative to total variance

# A # No grouping
K <- apply(params, 'iterations', lin_term_cov, X, y, residuals=FALSE) 
dim(K) <- c(n_lt, n_lt, n_iter)
dimnames(K) <- list(rows = linear_terms, cols = linear_terms, iterations = NULL)

P_x <- apply(K, 'iterations', function(x) {diag(x) / sum(diag(x))})
dimnames(P_x) <- list(variance = linear_terms, iterations = NULL)

#summary
vp_summary <- data.frame(
  Variable = rownames(P_x),
  Mean = rowMeans(P_x)
)

vp_summary

setwd("~/GitHub/PlanktonRUE-Project/Figures/VP_plots")

# png("Zoo_FDis_ng.png", width = 4000, height = 4000, res = 600)
# 
# par(mai = c(3.5, 0.8, 0.8, 0.4), # bottom, left, top, right inches
#     las = 2,
#     cex.axis = 1.1,
#     xpd = NA)
# 
# vp_plot(P_x, nm = linear_terms, main = 'Zooplankton Dispersion (ng)', las = 2)
# 
# dev.off()


# B # Variance partitioning for groups of linear terms
### Variance partition $V_b$

B <- sapply(groups, function(x) { as.integer(linear_terms %in% x) })
dimnames(B) <- list(linear_terms = linear_terms, groups = colnames(B))

K_B <- apply(K, 'iterations', function(x, B) { t(B) %*% x %*% B }, B)
dim(K_B) <- c(n_group, n_group, n_iter)
dimnames(K_B) <- list(rows = names(groups), cols = names(groups), iterations = NULL)

P <- apply(K_B, 'iterations', function(x) {diag(x) / sum(diag(x))})
dimnames(P) <- list(variance = names(groups), iterations = NULL)

dim(P)
Pmean = apply(P, c(1), mean)
sum(Pmean) # = 1

# now weighted by R2
R2_model = bayes_R2(model)
R2_model[,1]  # we keep only the estimate although the uncertainy can be soemtimes interesting it would be a bit complex and out of subject for you I guess

# create a weighting matrix of same dimension as you data (particularly good if you have also multiple species but here it is ok still)
varr <- aperm(array(R2_model[,1], dim = c(8000L, 12L)), perm = c(2L, 1L)) # check the dimensions carefully 12 groups (1L) for 8000 samples (2L) change number of groups and samples for each model
dim(varr)  # [1]    12 8000

Pweighted = P*varr

vp_summary_grouped <- data.frame(
  Variable = rownames(Pweighted),
  Mean = rowMeans(Pweighted, na.rm = TRUE),
  Lower95 = apply(Pweighted, 1, quantile, probs = 0.025, na.rm = TRUE),
  Upper95 = apply(Pweighted, 1, quantile, probs = 0.975, na.rm = TRUE)
)

vp_summary_grouped

sum(vp_summary_grouped$Mean[vp_summary_grouped$Variable %in% abiotic_vars]) #get sum of abiotic factors#

#95%CI
abiotic_mat <- Pweighted[
  grepl(
    paste(abiotic_vars, collapse = "|"),
    rownames(Pweighted)
  ),
]
abiotic_sum <- colSums(abiotic_mat)
mean_total <- mean(abiotic_sum)

ci <- quantile(
  abiotic_sum,
  probs = c(0.025, 0.975)
)

data.frame(
  mean = mean_total,
  lower95 = ci[1],
  upper95 = ci[2])

png("Zoo_FDis_r2.png", width = 4000, height = 3000, res = 600)

par(mai = c(2.0, 0.8, 0.5, 0.4), # bottom, left, top, right inches
    las = 2,
    cex.axis = 1.1,
    xpd = NA)

vp_plot(Pweighted, col = c('royalblue3','royalblue3','royalblue3','royalblue3','royalblue3','royalblue3','#278B9AFF',
                   '#44A57CFF','#819A7AFF','#58A449FF','salmon2','gold3'),
        nm = plot_labels2, main = 'Zooplankton Functional Dispersion', las=2)

dev.off()

##.....................................................................................................
##.....................................................................................................

##  Zoo RUE model
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

## define some covariates to group # this you still need to choose and 
# the name need to match the names in X

groups <- list(
  Temperature = c("Temperature", "ITemperatureE2"),
  Salinity = c("Salinity", "ISalinityE2"),
  pH = c("pH", "IpHE2"),
  Total_Nitrogen = c("TotalNitrogen_unfiltered", "ITotalNitrogen_unfilteredE2"),
  MaxDepth = c("Zooplankton_MaxDepth",  "IZooplankton_MaxDepthE2"),
  Daylight_hours = c('daylight_hours','Idaylight_hoursE2'),  ## 
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

### Variance partition $V_x$
#quantify how much variance each term explains relative to total variance

# A # No grouping
K <- apply(params, 'iterations', lin_term_cov, X, y, residuals=FALSE) 
dim(K) <- c(n_lt, n_lt, n_iter)
dimnames(K) <- list(rows = linear_terms, cols = linear_terms, iterations = NULL)

P_x <- apply(K, 'iterations', function(x) {diag(x) / sum(diag(x))})
dimnames(P_x) <- list(variance = linear_terms, iterations = NULL)

#summary
vp_summary <- data.frame(
  Variable = rownames(P_x),
  Mean = rowMeans(P_x)
)

vp_summary

setwd("~/GitHub/PlanktonRUE-Project/Figures/VP_plots")

# png("Zoo_RUE_ng.png", width = 6000, height = 4500, res = 600)
# 
# par(mai = c(4, 0.8, 0.8, 0.4), # bottom, left, top, right inches
#     las = 2,
#     cex.axis = 1.1,
#     xpd = NA)
# 
# vp_plot(P_x, nm = linear_terms, main = 'Zooplankton RUE (ng)', las = 2)
# 
# dev.off()


# B # Variance partitioning for groups of linear terms
### Variance partition $V_b$

B <- sapply(groups, function(x) { as.integer(linear_terms %in% x) })
dimnames(B) <- list(linear_terms = linear_terms, groups = colnames(B))

K_B <- apply(K, 'iterations', function(x, B) { t(B) %*% x %*% B }, B)
dim(K_B) <- c(n_group, n_group, n_iter)
dimnames(K_B) <- list(rows = names(groups), cols = names(groups), iterations = NULL)

P <- apply(K_B, 'iterations', function(x) {diag(x) / sum(diag(x))})
dimnames(P) <- list(variance = names(groups), iterations = NULL)

dim(P)
Pmean = apply(P, c(1), mean)
sum(Pmean) # = 1

# now weighted by R2
R2_model = bayes_R2(model)
R2_model[,1]  # we keep only the estimate although the uncertainy can be soemtimes interesting it would be a bit complex and out of subject for you I guess

# create a weighting matrix of same dimension as you data (particularly good if you have also multiple species but here it is ok still)
varr <- aperm(array(R2_model[,1], dim = c(8000L, 15L)), perm = c(2L, 1L)) # check the dimensions carefully 15 groups (1L) for 8000 samples (2L) change number of groups and samples for each model
dim(varr)  # [1]    15 8000

Pweighted = P*varr

vp_summary_grouped <- data.frame(
  Variable = rownames(Pweighted),
  Mean = rowMeans(Pweighted, na.rm = TRUE),
  Lower95 = apply(Pweighted, 1, quantile, probs = 0.025, na.rm = TRUE),
  Upper95 = apply(Pweighted, 1, quantile, probs = 0.975, na.rm = TRUE)
)

vp_summary_grouped

sum(vp_summary_grouped$Mean[vp_summary_grouped$Variable %in% abiotic_vars]) #get sum of abiotic factors#

#95%CI
abiotic_mat <- Pweighted[
  grepl(
    paste(abiotic_vars, collapse = "|"),
    rownames(Pweighted)
  ),
]
abiotic_sum <- colSums(abiotic_mat)
mean_total <- mean(abiotic_sum)

ci <- quantile(
  abiotic_sum,
  probs = c(0.025, 0.975)
)

data.frame(
  mean = mean_total,
  lower95 = ci[1],
  upper95 = ci[2])

png("Zoo_RUE_r2.png", width = 5300, height = 3500, res = 600)

par(mai = c(2.0, 0.8, 0.5, 0.4), # bottom, left, top, right inches
    las = 2,
    cex.axis = 1.1,
    xpd = NA)

vp_plot(Pweighted, col = c('royalblue3','royalblue3','royalblue3','royalblue3','royalblue3','royalblue3','#278B9AFF',
                   '#44A57CFF','#819A7AFF','#58A449FF', '#403369FF','#5C5992FF','#AE93BEFF', 'salmon2','gold3'),
        nm = RUE_plot_labels, main = 'Zooplankton RUE', las=2)

dev.off()

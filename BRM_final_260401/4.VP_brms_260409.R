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
source('Functions_for_brms_VP_Ap26.R')

##.....................................................................................................
##.....................................................................................................

##  Cyanobacteria dominance
setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA")
load('BRM_cyanodom_P_260408_Z.RDATA')

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
  Abiotic = c("polyTemperature2rawEQTRUE1", "polyTemperature2rawEQTRUE2",                     
              "polySalinity2rawEQTRUE1", "polySalinity2rawEQTRUE2", 
              "polypH2rawEQTRUE1", "polypH2rawEQTRUE2", 
              "polyTotalNitrogen_unfiltered2rawEQTRUE1", "polyTotalNitrogen_unfiltered2rawEQTRUE2",
              "polyMaxDepth2rawEQTRUE1",  "polyMaxDepth2rawEQTRUE2",
              "polydaylight_hours2rawEQTRUE1","polydaylight_hours2rawEQTRUE2"),  ## 
  DOY = c('sDOY_1', 'sDOY'), 
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

png("cyanobacteria_dominance_ng.png", width = 4000, height = 4000, res = 600)

par(mai = c(3.5, 0.8, 0.8, 0.4), # bottom, left, top, right inches
    las = 2,
    cex.axis = 1.1,
    xpd = NA)

vp_plot(P_x, nm = linear_terms, main = 'Cyanobacteria Dominance (ng)', las = 2)

dev.off()

# B # Variance partitioning for groups of linear terms
### Variance partition $V_b$

B <- sapply(groups, function(x) { as.integer(linear_terms %in% x) })
dimnames(B) <- list(linear_terms = linear_terms, groups = colnames(B))

K_B <- apply(K, 'iterations', function(x, B) { t(B) %*% x %*% B }, B)
dim(K_B) <- c(n_group, n_group, n_iter)
dimnames(K_B) <- list(rows = names(groups), cols = names(groups), iterations = NULL)

P <- apply(K_B, 'iterations', function(x) {diag(x) / sum(diag(x))})
dimnames(P) <- list(variance = names(groups), iterations = NULL)

vp_summary_grouped <- data.frame(
  Variable = rownames(P),
  Mean = rowMeans(P)
)

vp_summary_grouped

png("cyanobacteria_dominance.png", width = 4000, height = 3000, res = 600)

par(mai = c(2.0, 0.8, 0.8, 0.4), # bottom, left, top, right inches
    las = 2,
    cex.axis = 1.1,
    xpd = NA)

vp_plot(P, col = c('royalblue3','violetred3','salmon2','gold3'),
        nm = names(groups), main = 'Cyanobacteria Dominance', las=2)

dev.off()

##.....................................................................................................
##.....................................................................................................

##  Genus richness
setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA")
load('BRM_phyto_GD_P_260408.RDATA')

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
  Abiotic = c("polyTemperature2rawEQTRUE1", "polyTemperature2rawEQTRUE2",                     
              "polySalinity2rawEQTRUE1", "polySalinity2rawEQTRUE2", 
              "polypH2rawEQTRUE1", "polypH2rawEQTRUE2", 
              "polyTotalNitrogen_unfiltered2rawEQTRUE1", "polyTotalNitrogen_unfiltered2rawEQTRUE2",
              "polyMaxDepth2rawEQTRUE1",  "polyMaxDepth2rawEQTRUE2",
              "polydaylight_hours2rawEQTRUE1","polydaylight_hours2rawEQTRUE2"),  ## 
  DOY = c('sDOY_1', 'sDOY'), 
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

png("Phyto_GD_ng.png", width = 4000, height = 4000, res = 600)

par(mai = c(3.5, 0.8, 0.8, 0.4), # bottom, left, top, right inches
    las = 2,
    cex.axis = 1.1,
    xpd = NA)

vp_plot(P_x, nm = linear_terms, main = 'Phytoplankton Genus Richness (ng)', las = 2)

dev.off()


# B # Variance partitioning for groups of linear terms
### Variance partition $V_b$

B <- sapply(groups, function(x) { as.integer(linear_terms %in% x) })
dimnames(B) <- list(linear_terms = linear_terms, groups = colnames(B))

K_B <- apply(K, 'iterations', function(x, B) { t(B) %*% x %*% B }, B)
dim(K_B) <- c(n_group, n_group, n_iter)
dimnames(K_B) <- list(rows = names(groups), cols = names(groups), iterations = NULL)

P <- apply(K_B, 'iterations', function(x) {diag(x) / sum(diag(x))})
dimnames(P) <- list(variance = names(groups), iterations = NULL)

vp_summary_grouped <- data.frame(
  Variable = rownames(P),
  Mean = rowMeans(P)
)

vp_summary_grouped

png("Phyto_GD.png", width = 4000, height = 3000, res = 600)

par(mai = c(2.0, 0.8, 0.8, 0.4), # bottom, left, top, right inches
    las = 2,
    cex.axis = 1.1,
    xpd = NA)

vp_plot(P, col = c('royalblue3','violetred3','salmon2','gold3'),
        nm = names(groups), main = 'Phytoplankton Genus Richness', las=2)

dev.off()

##.....................................................................................................
##.....................................................................................................
##  phyto FEve
setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA")
load('BRM_phyto_FEve_P_260408.RDATA')

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
  Abiotic = c("polyTemperature2rawEQTRUE1", "polyTemperature2rawEQTRUE2",                     
              "polySalinity2rawEQTRUE1", "polySalinity2rawEQTRUE2", 
              "polypH2rawEQTRUE1", "polypH2rawEQTRUE2", 
              "polyTotalNitrogen_unfiltered2rawEQTRUE1", "polyTotalNitrogen_unfiltered2rawEQTRUE2",
              "polyMaxDepth2rawEQTRUE1",  "polyMaxDepth2rawEQTRUE2",
              "polydaylight_hours2rawEQTRUE1","polydaylight_hours2rawEQTRUE2"),  ## 
  DOY = c('sDOY_1', 'sDOY'), 
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

png("Phyto_FEve_ng.png", width = 4000, height = 4000, res = 600)

par(mai = c(3.5, 0.8, 0.8, 0.4), # bottom, left, top, right inches
    las = 2,
    cex.axis = 1.1,
    xpd = NA)

vp_plot(P_x, nm = linear_terms, main = 'Phytoplankton Functional Evenness (ng)', las = 2)

dev.off()


# B # Variance partitioning for groups of linear terms
### Variance partition $V_b$

B <- sapply(groups, function(x) { as.integer(linear_terms %in% x) })
dimnames(B) <- list(linear_terms = linear_terms, groups = colnames(B))

K_B <- apply(K, 'iterations', function(x, B) { t(B) %*% x %*% B }, B)
dim(K_B) <- c(n_group, n_group, n_iter)
dimnames(K_B) <- list(rows = names(groups), cols = names(groups), iterations = NULL)

P <- apply(K_B, 'iterations', function(x) {diag(x) / sum(diag(x))})
dimnames(P) <- list(variance = names(groups), iterations = NULL)

vp_summary_grouped <- data.frame(
  Variable = rownames(P),
  Mean = rowMeans(P)
)

vp_summary_grouped

png("Phyto_FEve.png", width = 4000, height = 3000, res = 600)

par(mai = c(2.0, 0.8, 0.8, 0.4), # bottom, left, top, right inches
    las = 2,
    cex.axis = 1.1,
    xpd = NA)

vp_plot(P, col = c('royalblue3','violetred3','salmon2','gold3'),
        nm = names(groups), main = 'Phytoplankton Functional Evenness', las=2)

dev.off()

##.....................................................................................................
##.....................................................................................................
##  phyto FDis
setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA")
load('BRM_phyto_FDis_P_260408.RDATA')

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
  Abiotic = c("polyTemperature2rawEQTRUE1", "polyTemperature2rawEQTRUE2",                     
              "polySalinity2rawEQTRUE1", "polySalinity2rawEQTRUE2", 
              "polypH2rawEQTRUE1", "polypH2rawEQTRUE2", 
              "polyTotalNitrogen_unfiltered2rawEQTRUE1", "polyTotalNitrogen_unfiltered2rawEQTRUE2",
              "polyMaxDepth2rawEQTRUE1",  "polyMaxDepth2rawEQTRUE2",
              "polydaylight_hours2rawEQTRUE1","polydaylight_hours2rawEQTRUE2"),  ## 
  DOY = c('sDOY_1', 'sDOY'), 
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

png("Phyto_FDis_ng.png", width = 4000, height = 4000, res = 600)

par(mai = c(3.5, 0.8, 0.8, 0.4), # bottom, left, top, right inches
    las = 2,
    cex.axis = 1.1,
    xpd = NA)

vp_plot(P_x, nm = linear_terms, main = 'Phytoplankton Functional Dispersion (ng)', las = 2)

dev.off()


# B # Variance partitioning for groups of linear terms
### Variance partition $V_b$

B <- sapply(groups, function(x) { as.integer(linear_terms %in% x) })
dimnames(B) <- list(linear_terms = linear_terms, groups = colnames(B))

K_B <- apply(K, 'iterations', function(x, B) { t(B) %*% x %*% B }, B)
dim(K_B) <- c(n_group, n_group, n_iter)
dimnames(K_B) <- list(rows = names(groups), cols = names(groups), iterations = NULL)

P <- apply(K_B, 'iterations', function(x) {diag(x) / sum(diag(x))})
dimnames(P) <- list(variance = names(groups), iterations = NULL)

vp_summary_grouped <- data.frame(
  Variable = rownames(P),
  Mean = rowMeans(P)
)

vp_summary_grouped

png("Phyto_FDis.png", width = 4000, height = 3000, res = 600)

par(mai = c(2.0, 0.8, 0.8, 0.4), # bottom, left, top, right inches
    las = 2,
    cex.axis = 1.1,
    xpd = NA)

vp_plot(P, col = c('royalblue3','violetred3','salmon2','gold3'),
        nm = names(groups), main = 'Phytoplankton Functional Dispersion', las=2)

dev.off()

##.....................................................................................................
##.....................................................................................................

##  Phyto RUE per P model
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

## define some covariates to group # this you still need to choose and 
# the name need to match the names in X

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

png("Phyto_RUE_P_ng.png", width = 5000, height = 4000, res = 600)

par(mai = c(2.5, 0.8, 0.8, 0.4), # bottom, left, top, right inches
    las = 2,
    cex.axis = 1.1,
    xpd = NA)

vp_plot(P_x, nm = linear_terms, main = 'Phytoplankton RUE per P (ng)', las = 2)

dev.off()


# B # Variance partitioning for groups of linear terms
### Variance partition $V_b$

B <- sapply(groups, function(x) { as.integer(linear_terms %in% x) })
dimnames(B) <- list(linear_terms = linear_terms, groups = colnames(B))

K_B <- apply(K, 'iterations', function(x, B) { t(B) %*% x %*% B }, B)
dim(K_B) <- c(n_group, n_group, n_iter)
dimnames(K_B) <- list(rows = names(groups), cols = names(groups), iterations = NULL)

P <- apply(K_B, 'iterations', function(x) {diag(x) / sum(diag(x))})
dimnames(P) <- list(variance = names(groups), iterations = NULL)

vp_summary_grouped <- data.frame(
  Variable = rownames(P),
  Mean = rowMeans(P)
)

vp_summary_grouped

png("Phyto_RUE_P.png", width = 4000, height = 3000, res = 600)

par(mai = c(2.0, 0.8, 0.8, 0.4), # bottom, left, top, right inches
    las = 2,
    cex.axis = 1.1,
    xpd = NA)

vp_plot(P, col = c('royalblue3','#278B9AFF','#44A57CFF','#819A7AFF','#58A449FF', '#403369FF','#5C5992FF','#AE93BEFF','violetred3','salmon2','gold3'),
        nm = names(groups), main = 'Phytoplankton RUE per P', las=2)

dev.off()

##.....................................................................................................
##.....................................................................................................


##  Phyto RUE per N model
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

## define some covariates to group # this you still need to choose and 
# the name need to match the names in X

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

png("Phyto_RUE_N_ng.png", width = 5000, height = 4000, res = 600)

par(mai = c(2.5, 0.8, 0.8, 0.4), # bottom, left, top, right inches
    las = 2,
    cex.axis = 1.1,
    xpd = NA)

vp_plot(P_x, nm = linear_terms, main = 'Phytoplankton RUE per N (ng)', las = 2)

dev.off()

# B # Variance partitioning for groups of linear terms
### Variance partition $V_b$

B <- sapply(groups, function(x) { as.integer(linear_terms %in% x) })
dimnames(B) <- list(linear_terms = linear_terms, groups = colnames(B))

K_B <- apply(K, 'iterations', function(x, B) { t(B) %*% x %*% B }, B)
dim(K_B) <- c(n_group, n_group, n_iter)
dimnames(K_B) <- list(rows = names(groups), cols = names(groups), iterations = NULL)

P <- apply(K_B, 'iterations', function(x) {diag(x) / sum(diag(x))})
dimnames(P) <- list(variance = names(groups), iterations = NULL)

vp_summary_grouped <- data.frame(
  Variable = rownames(P),
  Mean = rowMeans(P)
)

vp_summary_grouped

png("Phyto_RUE_N.png", width = 4000, height = 3000, res = 600)

par(mai = c(2.0, 0.8, 0.8, 0.4), # bottom, left, top, right inches
    las = 2,
    cex.axis = 1.1,
    xpd = NA)

vp_plot(P, col = c('royalblue3','#278B9AFF','#44A57CFF','#819A7AFF','#58A449FF', '#403369FF','#5C5992FF','#AE93BEFF','violetred3','salmon2','gold3'),
        nm = names(groups), main = 'Phytoplankton RUE per N', las=2)

dev.off()

##.....................................................................................................
##.....................................................................................................

##  Zoo Species Richness model
setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA")
load('BRM_zoo_SR_CD_P_260408.RDATA')

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
  Abiotic = c("polyTemperature2rawEQTRUE1", "polyTemperature2rawEQTRUE2",                     
              "polySalinity2rawEQTRUE1", "polySalinity2rawEQTRUE2", 
              "polypH2rawEQTRUE1", "polypH2rawEQTRUE2", 
              "polyTotalNitrogen_unfiltered2rawEQTRUE1", "polyTotalNitrogen_unfiltered2rawEQTRUE2", 
              "polyZooplankton_MaxDepth2rawEQTRUE1",  "polyZooplankton_MaxDepth2rawEQTRUE2",
              'polydaylight_hours2rawEQTRUE1','polydaylight_hours2rawEQTRUE2'),  ## 
  Phyto_genus_richness = c("genus_richness"),
  Phyto_FEve = c("FEve_selected"),
  Phyto_FDis = c("FDis_selected"),
  Cyano_dom = c("cyano_prop"),
  DOY = c('sDOY_1', 'sDOY'),
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

png("Zoo_SR_ng.png", width = 4000, height = 4000, res = 600)

par(mai = c(3.5, 0.8, 0.8, 0.4), # bottom, left, top, right inches
    las = 2,
    cex.axis = 1.1,
    xpd = NA)

vp_plot(P_x, nm = linear_terms, main = 'Zooplankton Species Richness (ng)', las = 2)

dev.off()


# B # Variance partitioning for groups of linear terms
### Variance partition $V_b$

B <- sapply(groups, function(x) { as.integer(linear_terms %in% x) })
dimnames(B) <- list(linear_terms = linear_terms, groups = colnames(B))

K_B <- apply(K, 'iterations', function(x, B) { t(B) %*% x %*% B }, B)
dim(K_B) <- c(n_group, n_group, n_iter)
dimnames(K_B) <- list(rows = names(groups), cols = names(groups), iterations = NULL)

P <- apply(K_B, 'iterations', function(x) {diag(x) / sum(diag(x))})
dimnames(P) <- list(variance = names(groups), iterations = NULL)

vp_summary_grouped <- data.frame(
  Variable = rownames(P),
  Mean = rowMeans(P)
)

vp_summary_grouped

png("Zoo_SR.png", width = 4000, height = 3000, res = 600)

par(mai = c(2.0, 0.8, 0.8, 0.4), # bottom, left, top, right inches
    las = 2,
    cex.axis = 1.1,
    xpd = NA)

vp_plot(P, col = c('royalblue3','#278B9AFF','#44A57CFF','#819A7AFF','#58A449FF','violetred3','salmon2','gold3'),
        nm = names(groups), main = 'Zooplankton Species Richness', las=2)

dev.off()

##.....................................................................................................
##.....................................................................................................


##  Zoo Functional Evenness model
setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA")
load('BRM_zoo_FE_CD_P_260408.RDATA')

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
  Abiotic = c("polyTemperature2rawEQTRUE1", "polyTemperature2rawEQTRUE2",                     
              "polySalinity2rawEQTRUE1", "polySalinity2rawEQTRUE2", 
              "polypH2rawEQTRUE1", "polypH2rawEQTRUE2", 
              "polyTotalNitrogen_unfiltered2rawEQTRUE1", "polyTotalNitrogen_unfiltered2rawEQTRUE2", 
              "polyZooplankton_MaxDepth2rawEQTRUE1",  "polyZooplankton_MaxDepth2rawEQTRUE2",
              'polydaylight_hours2rawEQTRUE1','polydaylight_hours2rawEQTRUE2'),  ## 
  Phyto_genus_richness = c("genus_richness"),
  Phyto_FEve = c("FEve_selected"),
  Phyto_FDis = c("FDis_selected"),
  Cyano_dom = c("cyano_prop"),
  DOY = c('sDOY_1', 'sDOY'),
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

png("Zoo_FEve_ng.png", width = 4000, height = 4000, res = 600)

par(mai = c(3.5, 0.8, 0.8, 0.4), # bottom, left, top, right inches
    las = 2,
    cex.axis = 1.1,
    xpd = NA)

vp_plot(P_x, nm = linear_terms, main = 'Zooplankton Functional Evenness (ng)', las = 2)

dev.off()


# B # Variance partitioning for groups of linear terms
### Variance partition $V_b$

B <- sapply(groups, function(x) { as.integer(linear_terms %in% x) })
dimnames(B) <- list(linear_terms = linear_terms, groups = colnames(B))

K_B <- apply(K, 'iterations', function(x, B) { t(B) %*% x %*% B }, B)
dim(K_B) <- c(n_group, n_group, n_iter)
dimnames(K_B) <- list(rows = names(groups), cols = names(groups), iterations = NULL)

P <- apply(K_B, 'iterations', function(x) {diag(x) / sum(diag(x))})
dimnames(P) <- list(variance = names(groups), iterations = NULL)

vp_summary_grouped <- data.frame(
  Variable = rownames(P),
  Mean = rowMeans(P)
)

vp_summary_grouped

png("Zoo_FEve.png", width = 4000, height = 3000, res = 600)

par(mai = c(2.0, 0.8, 0.8, 0.4), # bottom, left, top, right inches
    las = 2,
    cex.axis = 1.1,
    xpd = NA)

vp_plot(P, col = c('royalblue3','#278B9AFF','#44A57CFF','#819A7AFF','#58A449FF','violetred3','salmon2','gold3'),
        nm = names(groups), main = 'Zooplankton Functional Evenness', las=2)

dev.off()

##.....................................................................................................
##.....................................................................................................


##  Zoo Functional Dispersion model
setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA")
load('BRM_zoo_FD_CD_P_260408.RDATA')

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
  Abiotic = c("polyTemperature2rawEQTRUE1", "polyTemperature2rawEQTRUE2",                     
              "polySalinity2rawEQTRUE1", "polySalinity2rawEQTRUE2", 
              "polypH2rawEQTRUE1", "polypH2rawEQTRUE2", 
              "polyTotalNitrogen_unfiltered2rawEQTRUE1", "polyTotalNitrogen_unfiltered2rawEQTRUE2", 
              "polyZooplankton_MaxDepth2rawEQTRUE1",  "polyZooplankton_MaxDepth2rawEQTRUE2",
              'polydaylight_hours2rawEQTRUE1','polydaylight_hours2rawEQTRUE2'),  ## 
  Phyto_genus_richness = c("genus_richness"),
  Phyto_FEve = c("FEve_selected"),
  Phyto_FDis = c("FDis_selected"),
  Cyano_dom = c("cyano_prop"),
  DOY = c('sDOY_1', 'sDOY'),
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

png("Zoo_FDis_ng.png", width = 4000, height = 4000, res = 600)

par(mai = c(3.5, 0.8, 0.8, 0.4), # bottom, left, top, right inches
    las = 2,
    cex.axis = 1.1,
    xpd = NA)

vp_plot(P_x, nm = linear_terms, main = 'Zooplankton Dispersion (ng)', las = 2)

dev.off()


# B # Variance partitioning for groups of linear terms
### Variance partition $V_b$

B <- sapply(groups, function(x) { as.integer(linear_terms %in% x) })
dimnames(B) <- list(linear_terms = linear_terms, groups = colnames(B))

K_B <- apply(K, 'iterations', function(x, B) { t(B) %*% x %*% B }, B)
dim(K_B) <- c(n_group, n_group, n_iter)
dimnames(K_B) <- list(rows = names(groups), cols = names(groups), iterations = NULL)

P <- apply(K_B, 'iterations', function(x) {diag(x) / sum(diag(x))})
dimnames(P) <- list(variance = names(groups), iterations = NULL)

vp_summary_grouped <- data.frame(
  Variable = rownames(P),
  Mean = rowMeans(P)
)

vp_summary_grouped

png("Zoo_FDis.png", width = 4000, height = 3000, res = 600)

par(mai = c(2.0, 0.8, 0.8, 0.4), # bottom, left, top, right inches
    las = 2,
    cex.axis = 1.1,
    xpd = NA)

vp_plot(P, col = c('royalblue3','#278B9AFF','#44A57CFF','#819A7AFF','#58A449FF','violetred3','salmon2','gold3'),
        nm = names(groups), main = 'Zooplankton Functional Dispersion', las=2)

dev.off()

##.....................................................................................................
##.....................................................................................................

##  Zoo RUE model
setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA")
load('BRM_zoo_RUE_CD_PI_260408.RDATA')

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
  Abiotic = c("Temperature", "ITemperatureE2",                     
              "Salinity", "ISalinityE2", 
              "pH", "IpHE2", 
              "TotalNitrogen_unfiltered", "ITotalNitrogen_unfilteredE2", 
              "Zooplankton_MaxDepth",  "IZooplankton_MaxDepthE2",
              'daylight_hours','Idaylight_hoursE2'),  ## 
  Phyto_genus_richness = c("genus_richness", "Temperature:genus_richness", "Salinity:genus_richness", "TotalNitrogen_unfiltered:genus_richness"),
  Phyto_FEve = c("FEve_selected", "Temperature:FEve_selected", "Salinity:FEve_selected", "TotalNitrogen_unfiltered:FEve_selected"),
  Phyto_FDis = c("FDis_selected", "Temperature:FDis_selected", "Salinity:FDis_selected", "TotalNitrogen_unfiltered:FDis_selected"),
  Cyano_dom = c("cyano_prop", "Temperature:cyano_prop", "Salinity:cyano_prop", "TotalNitrogen_unfiltered:cyano_prop"),
  Zoo_species_richness = c("species_richness", "Temperature:species_richness", "Salinity:species_richness", "TotalNitrogen_unfiltered:species_richness"),
  Zoo_FEve = c("zoo_FEve", "Temperature:zoo_FEve", "Salinity:zoo_FEve", "TotalNitrogen_unfiltered:zoo_FEve"),
  Zoo_FDis = c("zoo_FDis", "Temperature:zoo_FDis", "Salinity:zoo_FDis", "TotalNitrogen_unfiltered:zoo_FDis"),
  Edibility= c("mean_edibility_unweighted"),
  DOY = c('sDOY_1', 'sDOY'),
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

png("Zoo_RUE_ng.png", width = 6000, height = 4500, res = 600)

par(mai = c(2.5, 0.8, 0.8, 0.4), # bottom, left, top, right inches
    las = 2,
    cex.axis = 1.1,
    xpd = NA)

vp_plot(P_x, nm = linear_terms, main = 'Zooplankton RUE (ng)', las = 2)

dev.off()


# B # Variance partitioning for groups of linear terms
### Variance partition $V_b$

B <- sapply(groups, function(x) { as.integer(linear_terms %in% x) })
dimnames(B) <- list(linear_terms = linear_terms, groups = colnames(B))

K_B <- apply(K, 'iterations', function(x, B) { t(B) %*% x %*% B }, B)
dim(K_B) <- c(n_group, n_group, n_iter)
dimnames(K_B) <- list(rows = names(groups), cols = names(groups), iterations = NULL)

P <- apply(K_B, 'iterations', function(x) {diag(x) / sum(diag(x))})
dimnames(P) <- list(variance = names(groups), iterations = NULL)

vp_summary_grouped <- data.frame(
  Variable = rownames(P),
  Mean = rowMeans(P)
)

vp_summary_grouped

png("Zoo_RUE.png", width = 4000, height = 3000, res = 600)

par(mai = c(2.0, 0.8, 0.8, 0.4), # bottom, left, top, right inches
    las = 2,
    cex.axis = 1.1,
    xpd = NA)

vp_plot(P, col = c('royalblue3','#278B9AFF','#44A57CFF','#819A7AFF','#58A449FF', '#403369FF','#5C5992FF','#AE93BEFF','#F7EABDFF', 'violetred3','salmon2','gold3'),
        nm = names(groups), main = 'Zooplankton RUE', las=2)

dev.off()

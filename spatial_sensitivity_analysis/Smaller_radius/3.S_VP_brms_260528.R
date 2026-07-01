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
  DOY = "DOY",
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
  DOY = "DOY",
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
  DOY = "DOY",
  Site = "Site",
  Year = "Year"
)

##.....................................................................................................
##.....................................................................................................

##  Cyanobacteria dominance
setwd("~/GitHub/PlanktonRUE-Project/spatial_sensitivity_analysis/Smaller_radius/RDATA")
load('BRM_cyanodom_P_260526_Z.RDATA')

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
varr <- aperm(array(R2_model[,1], dim = c(8000L, 9L)), perm = c(2L, 1L)) # check the dimensions carefully 9 groups (1L) for 8000 samples (2L) change number of groups and samples for each model
dim(varr)  # [1]    9 8000

Pweighted = P*varr

Pm_weighted = apply(P*varr, c(1), mean)
dim(Pm_weighted)

P_all = as.matrix(cbind(P1 = Pmean, Pw = Pm_weighted)) # if combined with multiple species

vp_summary_grouped <- data.frame(
  Variable = rownames(Pweighted),
  Mean = rowMeans(Pweighted)
)

vp_summary_grouped

sum(vp_summary_grouped$Mean[vp_summary_grouped$Variable %in% abiotic_vars]) #get sum of abiotic factors#


##.....................................................................................................
##.....................................................................................................

##  Genus richness
setwd("~/GitHub/PlanktonRUE-Project/spatial_sensitivity_analysis/Smaller_radius/RDATA")
load('BRM_phyto_GD_P_260526.RDATA')

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
varr <- aperm(array(R2_model[,1], dim = c(8000L, 9L)), perm = c(2L, 1L)) # check the dimensions carefully 9 groups (1L) for 8000 samples (2L) change number of groups and samples for each model
dim(varr)  # [1]    9 8000

Pweighted = P*varr

vp_summary_grouped <- data.frame(
  Variable = rownames(Pweighted),
  Mean = rowMeans(Pweighted)
)

vp_summary_grouped

sum(vp_summary_grouped$Mean[vp_summary_grouped$Variable %in% abiotic_vars]) #get sum of abiotic factors#

##.....................................................................................................
##.....................................................................................................
##  phyto FEve
setwd("~/GitHub/PlanktonRUE-Project/spatial_sensitivity_analysis/Smaller_radius/RDATA")
load('BRM_phyto_FEve_P_260526.RDATA')

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
varr <- aperm(array(R2_model[,1], dim = c(8000L, 9L)), perm = c(2L, 1L)) # check the dimensions carefully 9 groups (1L) for 8000 samples (2L) change number of groups and samples for each model
dim(varr)  # [1]    9 8000

Pweighted = P*varr

vp_summary_grouped <- data.frame(
  Variable = rownames(Pweighted),
  Mean = rowMeans(Pweighted)
)

vp_summary_grouped

sum(vp_summary_grouped$Mean[vp_summary_grouped$Variable %in% abiotic_vars]) #get sum of abiotic factors#

##.....................................................................................................
##.....................................................................................................
##  phyto FDis
setwd("~/GitHub/PlanktonRUE-Project/spatial_sensitivity_analysis/Smaller_radius/RDATA")
load('BRM_phyto_FDis_P_260526.RDATA')

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
varr <- aperm(array(R2_model[,1], dim = c(8000L, 9L)), perm = c(2L, 1L)) # check the dimensions carefully 9 groups (1L) for 8000 samples (2L) change number of groups and samples for each model
dim(varr)  # [1]    9 8000

Pweighted = P*varr

vp_summary_grouped <- data.frame(
  Variable = rownames(Pweighted),
  Mean = rowMeans(Pweighted)
)

vp_summary_grouped

sum(vp_summary_grouped$Mean[vp_summary_grouped$Variable %in% abiotic_vars]) #get sum of abiotic factors#


##.....................................................................................................
##.....................................................................................................

##  Phyto RUE per P model
setwd("~/GitHub/PlanktonRUE-Project/spatial_sensitivity_analysis/Smaller_radius/RDATA")
load('BRM_phytoRUEP_CD_PI_260526.RDATA')

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
varr <- aperm(array(R2_model[,1], dim = c(8000L, 16L)), perm = c(2L, 1L)) # check the dimensions carefully 16 groups (1L) for 8000 samples (2L) change number of groups and samples for each model
dim(varr)  # [1]    9 8000

Pweighted = P*varr

vp_summary_grouped <- data.frame(
  Variable = rownames(Pweighted),
  Mean = rowMeans(Pweighted)
)

vp_summary_grouped

sum(vp_summary_grouped$Mean[vp_summary_grouped$Variable %in% abiotic_vars]) #get sum of abiotic factors#

##.....................................................................................................
##.....................................................................................................


##  Phyto RUE per N model
setwd("~/GitHub/PlanktonRUE-Project/spatial_sensitivity_analysis/Smaller_radius/RDATA")
load('BRM_phytoRUEN_CD_PI_260526.RDATA')

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
varr <- aperm(array(R2_model[,1], dim = c(8000L, 16L)), perm = c(2L, 1L)) # check the dimensions carefully 16 groups (1L) for 8000 samples (2L) change number of groups and samples for each model
dim(varr)  # [1]    9 8000

Pweighted = P*varr

vp_summary_grouped <- data.frame(
  Variable = rownames(Pweighted),
  Mean = rowMeans(Pweighted)
)

vp_summary_grouped

sum(vp_summary_grouped$Mean[vp_summary_grouped$Variable %in% abiotic_vars]) #get sum of abiotic factors#

##.....................................................................................................
##.....................................................................................................

##  Zoo Species Richness model
setwd("~/GitHub/PlanktonRUE-Project/spatial_sensitivity_analysis/Smaller_radius/RDATA")
load('BRM_zoo_SR_CD_P_260526.RDATA')

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
varr <- aperm(array(R2_model[,1], dim = c(8000L, 13L)), perm = c(2L, 1L)) # check the dimensions carefully 13 groups (1L) for 8000 samples (2L) change number of groups and samples for each model
dim(varr)  # [1]    9 8000

Pweighted = P*varr

vp_summary_grouped <- data.frame(
  Variable = rownames(Pweighted),
  Mean = rowMeans(Pweighted)
)

vp_summary_grouped

sum(vp_summary_grouped$Mean[vp_summary_grouped$Variable %in% abiotic_vars]) #get sum of abiotic factors#

##.....................................................................................................
##.....................................................................................................


##  Zoo Functional Evenness model
setwd("~/GitHub/PlanktonRUE-Project/spatial_sensitivity_analysis/Smaller_radius/RDATA")
load('BRM_zoo_FE_CD_P_260526.RDATA')

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
varr <- aperm(array(R2_model[,1], dim = c(8000L, 13L)), perm = c(2L, 1L)) # check the dimensions carefully 13 groups (1L) for 8000 samples (2L) change number of groups and samples for each model
dim(varr)  # [1]    13 8000

Pweighted = P*varr

vp_summary_grouped <- data.frame(
  Variable = rownames(Pweighted),
  Mean = rowMeans(Pweighted)
)

vp_summary_grouped

sum(vp_summary_grouped$Mean[vp_summary_grouped$Variable %in% abiotic_vars]) #get sum of abiotic factors#

##.....................................................................................................
##.....................................................................................................


##  Zoo Functional Dispersion model
setwd("~/GitHub/PlanktonRUE-Project/spatial_sensitivity_analysis/Smaller_radius/RDATA")
load('BRM_zoo_FD_CD_P_260526.RDATA')

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
varr <- aperm(array(R2_model[,1], dim = c(8000L, 13L)), perm = c(2L, 1L)) # check the dimensions carefully 13 groups (1L) for 8000 samples (2L) change number of groups and samples for each model
dim(varr)  # [1]    9 8000

Pweighted = P*varr

vp_summary_grouped <- data.frame(
  Variable = rownames(Pweighted),
  Mean = rowMeans(Pweighted)
)

vp_summary_grouped

sum(vp_summary_grouped$Mean[vp_summary_grouped$Variable %in% abiotic_vars]) #get sum of abiotic factors#

##.....................................................................................................
##.....................................................................................................

##  Zoo RUE model
setwd("~/GitHub/PlanktonRUE-Project/spatial_sensitivity_analysis/Smaller_radius/RDATA")
load('BRM_zoo_RUE_CD_PI_260526.RDATA')

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
varr <- aperm(array(R2_model[,1], dim = c(8000L, 16L)), perm = c(2L, 1L)) # check the dimensions carefully 16 groups (1L) for 8000 samples (2L) change number of groups and samples for each model
dim(varr)  # [1]    9 8000

Pweighted = P*varr

vp_summary_grouped <- data.frame(
  Variable = rownames(Pweighted),
  Mean = rowMeans(Pweighted)
)

vp_summary_grouped

sum(vp_summary_grouped$Mean[vp_summary_grouped$Variable %in% abiotic_vars]) #get sum of abiotic factors#

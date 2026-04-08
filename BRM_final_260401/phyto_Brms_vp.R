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

rm(list = ls())

##...........................................................................................................

##.........................................................

## Get the functions
## ..........................................................................
## ............................................................................

## new and adapted functions
setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401")
source('Functions_for_brms_VP_Ap26.R') ## you need the new one here

##  Zoo RUE model
setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA")
load('BRM_zoo_RUE_CD_PI_260330_S_Nit.RDaTA')

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
  Phyto_div = c("genus_richness", "FEve_selected", "FDis_selected"),
  Zoo_div = c("species_richness", "zoo_FEve", "zoo_FDis"),
  Cyano_dom = c("cyano_prop"),
  Abiotic_phyto_div = c("Temperature:genus_richness", "Salinity:genus_richness", "TotalNitrogen_unfiltered:genus_richness",
                        "Temperature:FEve_selected", "Salinity:FEve_selected", "TotalNitrogen_unfiltered:FEve_selected", 
                        "Temperature:FDis_selected", "Salinity:FDis_selected", "TotalNitrogen_unfiltered:FDis_selected"),
  Abiotic_cyano = c("Temperature:cyano_prop", "Salinity:cyano_prop", "TotalNitrogen_unfiltered:cyano_prop"),
  Abiotic_zoo_div = c("Temperature:species_richness", "Salinity:species_richness", "TotalNitrogen_unfiltered:species_richness",
                      "Temperature:zoo_FEve", "Salinity:zoo_FEve", "TotalNitrogen_unfiltered:zoo_FEve",
                      "Temperature:zoo_FDis", "Salinity:zoo_FDis", "TotalNitrogen_unfiltered:zoo_FDis"),
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
K <- apply(params, 'iterations', lin_term_cov, X, y)
dim(K) <- c(n_lt, n_lt, n_iter)
dimnames(K) <- list(rows = linear_terms, cols = linear_terms, iterations = NULL)

P_x <- apply(K, 'iterations', function(x) {diag(x) / sum(diag(x))})
dimnames(P_x) <- list(variance = linear_terms, iterations = NULL)

setwd("~/GitHub/PlanktonRUE-Project/Figures/VP_plots")

png("Zoo_RUE_ng.png", width = 5000, height = 4000, res = 600)

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

png("Zoo_RUE.png", width = 4000, height = 3000, res = 600)

par(mai = c(2.0, 0.8, 0.8, 0.4), # bottom, left, top, right inches
    las = 2,
    cex.axis = 1.1,
    xpd = NA)

vp_plot(P, c('royalblue3','#278B9AFF','#1D2645FF','#44A57CFF','#819A7AFF','#58A449FF','#CEC917FF','violetred3','salmon2','gold3'),
        nm = names(groups), main = 'Zooplankton RUE', las=2)

dev.off()

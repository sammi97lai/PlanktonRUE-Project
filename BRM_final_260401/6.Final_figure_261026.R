#
#               Script to plot Final Figure
#
#............................................................................................

# manipulating data
library(tidyverse)

# bayesian modeling
library("bayesplot")
library(brms)
library(posterior)

## plotting
library("ggplot2")
library(ggplotify)
library(patchwork)
library(see)
library(vioplot)
library(ggnewscale)
library(cowplot)

rm(list = ls())

##...........................................................................................................

############################# Variance Partitioning plot #############################
## Get the functions
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
RUE_plot_labels <- c(
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

##.....................................................................................................
##.....................................................................................................

##  Phyto RUE per P model
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
varr <- aperm(array(R2_model[,1], dim = c(8000L, 11L)), perm = c(2L, 1L)) # check the dimensions carefully 11 groups (1L) for 8000 samples (2L) change number of groups and samples for each model
dim(varr)  # [1]    11 8000

Pweighted = P*varr

RUE_p_vp <- as.ggplot(~{
  par(
    mai = c(2.0, 0.8, 0.5, 0.4),
    las = 1,
    cex.axis = 1.1,
    xpd = NA
  )
  
  vp_plot(
    Pweighted,
    col = c(
      "royalblue3",
      "#278B9AFF", "#274637FF", "#44A57CFF",
      "#819A7AFF", "#58A449FF", "#403369FF",
      "#9F99B5FF", "#5C5992FF", "#AE93BEFF",
      "gold3"
    ),
    nm = RUE_plot_labels,
    ylim = c(0, 0.8),
    main = "Phytoplankton RUE per P",
    ylab = "Proportion of explained variance",
    las = 1
  )
})

setwd("~/GitHub/PlanktonRUE-Project/Figures/Final_Figure")
png("Phyto_RUE_P.png", width = 5300, height = 3000, res = 600)

par(mai = c(0.5, 0.7, 0.7, 0.5), # bottom, left, top, right inches
    las = 1,
    cex.axis = 1.0,
    cex.main = 1.5,
    mgp = c(4.5, 1.2, 0),
    xpd = NA,
    xaxt = "n")

vp_plot(Pweighted, 
        col = c(
          "royalblue3",
          "#278B9AFF", "#274637FF", "#44A57CFF",
          "#819A7AFF", "#58A449FF", "#403369FF",
          "#9F99B5FF", "#5C5992FF", "#AE93BEFF",
          "gold3"
        ),
        nm = "", main = 'Phytoplankton RUE per P')

dev.off()

##.....................................................................................................
##.....................................................................................................

##  Phyto RUE per N model
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
varr <- aperm(array(R2_model[,1], dim = c(8000L, 11L)), perm = c(2L, 1L)) # check the dimensions carefully 11 groups (1L) for 8000 samples (2L) change number of groups and samples for each model
dim(varr)  # [1]    11 8000

Pweighted = P*varr

RUE_n_vp <- as.ggplot(~{
  par(
    mai = c(2.0, 0.8, 0.5, 0.4),
    las = 2,
    cex.axis = 1.1,
    xpd = NA
  )
  
  vp_plot(
    Pweighted,
    col = c(
      "royalblue3",
      "#278B9AFF", "#274637FF", "#44A57CFF",
      "#819A7AFF", "#58A449FF", "#403369FF",
      "#9F99B5FF", "#5C5992FF", "#AE93BEFF",
      "gold3"
    ),
    nm = RUE_plot_labels,
    ylim = c(0, 0.8),
    main = "Phytoplankton RUE per N",
    ylab = "",
    las = 2
  )
})

setwd("~/GitHub/PlanktonRUE-Project/Figures/Final_Figure")
png("Phyto_RUE_N.png", width = 5300, height = 3000, res = 600)

par(mai = c(0.5, 0.7, 0.7, 0.5), # bottom, left, top, right inches
    las = 1,
    cex.axis = 1.0,
    cex.main = 1.5,
    mgp = c(4.5, 1.2, 0),
    xpd = NA,
    xaxt = "n")

vp_plot(Pweighted, 
        col = c(
          "royalblue3",
          "#278B9AFF", "#274637FF", "#44A57CFF",
          "#819A7AFF", "#58A449FF", "#403369FF",
          "#9F99B5FF", "#5C5992FF", "#AE93BEFF",
          "gold3"
        ),
        nm = "", main = 'Phytoplankton RUE per N')

dev.off()

##.....................................................................................................
##.....................................................................................................

##  Zoo RUE model
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
varr <- aperm(array(R2_model[,1], dim = c(8000L, 11L)), perm = c(2L, 1L)) # check the dimensions carefully 11 groups (1L) for 8000 samples (2L) change number of groups and samples for each model
dim(varr)  # [1]    11 8000

Pweighted = P*varr

RUE_zoo_vp <- as.ggplot(~{
  par(
    mai = c(2.0, 0.8, 0.5, 0.4),
    las = 2,
    cex.axis = 1.1,
    xpd = NA
  )
  
  vp_plot(
    Pweighted,
    col = c(
      "royalblue3",
      "#278B9AFF", "#274637FF", "#44A57CFF",
      "#819A7AFF", "#58A449FF", "#403369FF",
      "#9F99B5FF", "#5C5992FF", "#AE93BEFF",
      "gold3"
    ),
    nm = RUE_plot_labels,
    ylim = c(0, 0.8),
    main = "Zooplankton RUE",
    ylab = "",
    las = 2
  )
})

setwd("~/GitHub/PlanktonRUE-Project/Figures/Final_Figure")
png("Zoo_RUE.png", width = 5300, height = 3000, res = 600)

par(mai = c(0.5, 0.7, 0.7, 0.5), # bottom, left, top, right inches
    las = 1,
    cex.axis = 1.0,
    cex.main = 1.5,
    mgp = c(4.5, 1.2, 0),
    xpd = NA,
    xaxt = "n")

vp_plot(Pweighted, 
        col = c(
          "royalblue3",
          "#278B9AFF", "#274637FF", "#44A57CFF",
          "#819A7AFF", "#58A449FF", "#403369FF",
          "#9F99B5FF", "#5C5992FF", "#AE93BEFF",
          "gold3"
        ),
        nm = "", main = 'Zooplankton RUE')

dev.off()

#############################################################################
############ Posterior distribution Violin plots ############################

# Biodiversity metrics and their display labels
biotic_effects <- c(
  "genus_richness",
  "phyto_evenness",
  "phyto_FEve_PRUE",
  "phyto_FDis_PRUE",
  "cyano_prop",
  "species_richness",
  "zoo_evenness",
  "zoo_FEve",
  "zoo_FDis"
)

biotic_labels <- c(
  "genus_richness"   = "Phytoplankton\nGenus\nRichness",
  "phyto_evenness"   = "Phytoplankton\nEvenness",
  "phyto_FEve_PRUE"  = "Phytoplankton\nFEve",
  "phyto_FDis_PRUE"  = "Phytoplankton\nFDis",
  "cyano_prop"       = "Cyanobacteria\nDominance",
  "species_richness" = "Zooplankton\nSpecies\nRichness",
  "zoo_evenness"     = "Zooplankton\nEvenness",
  "zoo_FEve"         = "Zooplankton\nFEve",
  "zoo_FDis"         = "Zooplankton\nFDis"
)

biotic_colours <- c(
  "genus_richness"   = "#278B9AFF",
  "phyto_evenness"   = "#274637FF",
  "phyto_FEve_PRUE"  = "#44A57CFF",
  "phyto_FDis_PRUE"  = "#819A7AFF",
  "cyano_prop"       = "#58A449FF",
  "species_richness" = "#403369FF",
  "zoo_evenness"     = "#9F99B5FF",
  "zoo_FEve"         = "#5C5992FF",
  "zoo_FDis"         = "#AE93BEFF"
)

interaction_env <- c(
  "Temperature",
  "Salinity",
  "TotalNitrogen_unfiltered"
)

env_levels <- c("Main effect", interaction_env)

env_colours <- c(
  "Main effect" = "grey25",
  "Temperature" = "red",
  "Salinity" = "orange",
  "TotalNitrogen_unfiltered" = "yellow"
)

##  Phyto RUE per P model ..........................................
setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA")
load('BRM_phytoRUEP_CD_PI_261013_All.RDATA')

# Significance stars for biodiversity main effects
RUEp_star_labels <- data.frame(
  biotic = biotic_effects,
  stars = c("***", "***", "***", "", "***", "*", "", "", "*"),
  y_pos = rep(0.5, length(biotic_effects))
)

RUEp_star_labels$x <- match(
  RUEp_star_labels$biotic,
  biotic_effects
)

# Extract posterior draws
draws <- posterior::as_draws_df(model)

# Find coefficient names safely, regardless of interaction order
find_coef <- function(a, b = NULL) {
  
  candidates <- grep("^b_", names(draws), value = TRUE)
  
  matches <- vapply(candidates, function(nm) {
    
    parts <- strsplit(
      sub("^b_", "", nm),
      ":",
      fixed = TRUE
    )[[1]]
    
    if (is.null(b)) {
      length(parts) == 1 && parts[1] == a
    } else {
      length(parts) == 2 &&
        setequal(parts, c(a, b))
    }
    
  }, logical(1))
  
  found <- candidates[matches]
  
  if (length(found) != 1) {
    stop(
      "Could not uniquely identify coefficient: ",
      a,
      if (!is.null(b)) paste0(":", b)
    )
  }
  
  found
}

# Horizontal offsets around each biodiversity main effect
# Negative offsets = low conditions; positive = high conditions
offsets <- data.frame(
  env = c(
    "Temperature", "Temperature",
    "Salinity", "Salinity",
    "TotalNitrogen_unfiltered",
    "TotalNitrogen_unfiltered"
  ),
  condition = c(
    "Low", "High",
    "Low", "High",
    "Low", "High"
  ),
  offset = c(-0.30, 0.30, -0.20, 0.20, -0.10, 0.10)
)

# Calculate posterior distributions
plot_list <- list()

for (i in seq_along(biotic_effects)) {
  
  bio <- biotic_effects[i]
  
  main_col <- find_coef(bio)
  
  # Main effect
  plot_list[[length(plot_list) + 1]] <- data.frame(
    x = rep(i, nrow(draws)),
    estimate = draws[[main_col]],
    env = "Main effect",
    condition = "Main effect",
    biotic = bio
  )
  
  # Conditional effects at low and high environmental values
  for (env in interaction_env) {
    
    int_col <- find_coef(bio, env)
    
    for (level in c("Low", "High")) {
      
      env_value <- if (level == "Low") -1 else 1
      
      # Select exactly one horizontal offset
      offset_index <- which(
        offsets$env == env &
          offsets$condition == level
      )
      
      if (length(offset_index) != 1) {
        stop("Expected one offset for ", env, " / ", level)
      }
      
      pos <- offsets$offset[offset_index]
      
      # Conditional posterior draws
      conditional_draws <-
        draws[[main_col]] +
        env_value * draws[[int_col]]
      
      plot_list[[length(plot_list) + 1]] <- data.frame(
        x = rep(i + pos, length(conditional_draws)),
        estimate = conditional_draws,
        env = env,
        condition = level,
        biotic = bio
      )
    }
  }
}

# Combine all posterior draws
interaction_plot_data <- dplyr::bind_rows(plot_list) %>%
  dplyr::mutate(
    env = factor(
      env,
      levels = env_levels
    ),
    condition = factor(
      condition,
      levels = c("Main effect", "Low", "High")
    ),
    
    # Biodiversity-specific colour for main effects;
    # environmental colour for conditional effects
    plot_colour = dplyr::if_else(
      as.character(env) == "Main effect",
      biotic,
      as.character(env)
    ),
    
    group_id = interaction(biotic, env, condition)
  )

# Calculate posterior medians and 95% credible intervals
interaction_summary <- interaction_plot_data %>%
  dplyr::filter(condition %in% c("Low", "High")) %>%
  dplyr::group_by(biotic, env, condition, x, plot_colour) %>%
  dplyr::summarise(
    posterior_median = median(estimate, na.rm = TRUE),
    posterior_lower = quantile(estimate, 0.025, na.rm = TRUE),
    posterior_upper = quantile(estimate, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::rename(
    estimate = posterior_median,
    lower = posterior_lower,
    upper = posterior_upper
  )

# Combine biodiversity and environmental colour palettes
all_colours <- c(
  biotic_colours,
  env_colours[interaction_env]
)

# Legend order
legend_breaks <- c(
  biotic_effects,
  interaction_env
)

legend_labels <- c(
  unname(biotic_labels[biotic_effects]),
  "Temperature",
  "Salinity",
  "Total Nitrogen"
)

# Alternating background shading for biodiversity metrics
background_data <- data.frame(
  xmin = seq(0.5, length(biotic_effects) - 0.5, by = 1),
  xmax = seq(1.5, length(biotic_effects) + 0.5, by = 1),
  ymin = -Inf,
  ymax = Inf,
  shade = rep(c("white", "grey90"), length.out = length(biotic_effects))
)

# Plot all main and conditional effects together
RUE_p_biotic <- ggplot(
  interaction_plot_data,
  aes(
    x = x,
    y = estimate,
    fill = plot_colour,
    colour = plot_colour,
    group = group_id,
    linetype = condition
  )
) +
  # Alternating background shading
  geom_rect(
    data = background_data,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = shade
    ),
    inherit.aes = FALSE,
    alpha = 1
  ) +
  scale_fill_identity() +
  ggnewscale::new_scale_fill() +
  # Zero reference line
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey60",
    linewidth = 0.4
  ) +
  # Main effects: half-violin distributions only
  see::geom_violinhalf(
    data = subset(
      interaction_plot_data,
      condition == "Main effect"
    ),
    aes(
      fill = plot_colour,
      colour = plot_colour,
      alpha = condition
    ),
    position = "identity",
    flip = FALSE,
    scale = "width",
    trim = FALSE,
    width = 0.15,
    linewidth = 0.3
  ) +
  # Transparency for main effects
  scale_alpha_manual(
    values = c(
      "Main effect" = 1,
      "Low" = 0.30,
      "High" = 0.30
    ),
    guide = "none"
  ) +
  # Environmental interactions:
  # Posterior median and 95% credible intervals
  geom_pointrange(
    data = interaction_summary,
    aes(
      x = x,
      y = estimate,
      ymin = lower,
      ymax = upper,
      colour = plot_colour,
      linetype = condition
    ),
    inherit.aes = FALSE,
    linewidth = 0.7,
    fatten = 2
  ) +
  # Significance stars for main effects
  geom_text(
    data = subset(
      RUEp_star_labels,
      stars != ""
    ),
    aes(
      x = x,
      y = y_pos,
      label = stars
    ),
    inherit.aes = FALSE,
    colour = "red",
    size = 5,
    vjust = 0
  ) +
  # Fill colours
  scale_fill_manual(
    values = all_colours,
    breaks = interaction_env,
    labels = c(
      "Temperature",
      "Salinity",
      "Total Nitrogen"
    ),
    name = "Environmental interactions"
  ) +
  # Line and point colours
  scale_colour_manual(
    values = all_colours,
    breaks = interaction_env,
    labels = c(
      "Temperature",
      "Salinity",
      "Total Nitrogen"
    ),
    name = "Environmental interactions"
  ) +
  # Y-axis
  scale_y_continuous(
    limits = c(-0.7, 0.6)
  ) +
  # Environmental condition line types
  scale_linetype_manual(
    values = c(
      "Main effect" = "solid",
      "Low" = "dashed",
      "High" = "solid"
    ),
    name = "Environmental level"
  ) +
  # X-axis
  scale_x_continuous(
    breaks = seq_along(biotic_effects),
    labels = unname(biotic_labels[biotic_effects]),
    expand = expansion(add = c(0.4, 0.4))
  ) +
  # Axis labels
  labs(
    x = NULL,
    y = NULL
  ) +
  # Theme
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    plot.title = element_text(size = 12, hjust = 0.5),
    legend.position = "none",
    plot.margin = margin(3, 5, 3, 5)
  )


# Display
RUE_p_biotic

setwd("~/GitHub/PlanktonRUE-Project/Figures/Final_Figure")
svg("PD_Phyto_RUE_P.svg", width = 8.83, height = 5, pointsize = 10)

RUE_p_biotic

dev.off()

# Extract legend from original plot
RUE_legend <- get_legend(
  RUE_p_biotic +
    theme(
      legend.position = "top",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 10)
    )
)

# Display separate legend
grid::grid.newpage()
grid::grid.draw(RUE_legend)

setwd("~/GitHub/PlanktonRUE-Project/Figures/Final_Figure")
svg("PD_legend.svg", width = 10, height = 1.5, pointsize = 10)

grid::grid.newpage()
grid::grid.draw(RUE_legend)

dev.off()

##  Phyto RUE per N model ..........................................
setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA")
load('BRM_phytoRUEN_CD_PI_261013_All.RDATA')

#set CI stars
RUEn_star_labels <- data.frame(
  biotic = biotic_effects,
  stars = c("***", "***", "***", "", "***", "", "", "***", ""),
  y_pos = rep(0.5, length(biotic_effects))
)

RUEn_star_labels$x <- match(
  RUEn_star_labels$biotic,
  biotic_effects
)

# Extract posterior draws
draws <- posterior::as_draws_df(model)

# Find coefficient names safely, regardless of interaction order
find_coef <- function(a, b = NULL) {
  
  candidates <- grep("^b_", names(draws), value = TRUE)
  
  matches <- vapply(candidates, function(nm) {
    
    parts <- strsplit(
      sub("^b_", "", nm),
      ":",
      fixed = TRUE
    )[[1]]
    
    if (is.null(b)) {
      length(parts) == 1 && parts[1] == a
    } else {
      length(parts) == 2 &&
        setequal(parts, c(a, b))
    }
    
  }, logical(1))
  
  found <- candidates[matches]
  
  if (length(found) != 1) {
    stop(
      "Could not uniquely identify coefficient: ",
      a,
      if (!is.null(b)) paste0(":", b)
    )
  }
  
  found
}

# Horizontal offsets around each biodiversity main effect
# Negative offsets = low conditions; positive = high conditions
offsets <- data.frame(
  env = c(
    "Temperature", "Temperature",
    "Salinity", "Salinity",
    "TotalNitrogen_unfiltered",
    "TotalNitrogen_unfiltered"
  ),
  condition = c(
    "Low", "High",
    "Low", "High",
    "Low", "High"
  ),
  offset = c(-0.30, 0.30, -0.20, 0.20, -0.10, 0.10)
)

# Calculate posterior distributions
plot_list <- list()

for (i in seq_along(biotic_effects)) {
  
  bio <- biotic_effects[i]
  
  main_col <- find_coef(bio)
  
  # Main effect
  plot_list[[length(plot_list) + 1]] <- data.frame(
    x = rep(i, nrow(draws)),
    estimate = draws[[main_col]],
    env = "Main effect",
    condition = "Main effect",
    biotic = bio
  )
  
  # Conditional effects at low and high environmental values
  for (env in interaction_env) {
    
    int_col <- find_coef(bio, env)
    
    for (level in c("Low", "High")) {
      
      env_value <- if (level == "Low") -1 else 1
      
      # Select exactly one horizontal offset
      offset_index <- which(
        offsets$env == env &
          offsets$condition == level
      )
      
      if (length(offset_index) != 1) {
        stop("Expected one offset for ", env, " / ", level)
      }
      
      pos <- offsets$offset[offset_index]
      
      # Conditional posterior draws
      conditional_draws <-
        draws[[main_col]] +
        env_value * draws[[int_col]]
      
      plot_list[[length(plot_list) + 1]] <- data.frame(
        x = rep(i + pos, length(conditional_draws)),
        estimate = conditional_draws,
        env = env,
        condition = level,
        biotic = bio
      )
    }
  }
}

# Combine all posterior draws
interaction_plot_data <- dplyr::bind_rows(plot_list) %>%
  dplyr::mutate(
    env = factor(
      env,
      levels = env_levels
    ),
    condition = factor(
      condition,
      levels = c("Main effect", "Low", "High")
    ),
    
    # Biodiversity-specific colour for main effects;
    # environmental colour for conditional effects
    plot_colour = dplyr::if_else(
      as.character(env) == "Main effect",
      biotic,
      as.character(env)
    ),
    
    group_id = interaction(biotic, env, condition)
  )

# Calculate posterior medians and 95% credible intervals
interaction_summary <- interaction_plot_data %>%
  dplyr::filter(condition %in% c("Low", "High")) %>%
  dplyr::group_by(biotic, env, condition, x, plot_colour) %>%
  dplyr::summarise(
    posterior_median = median(estimate, na.rm = TRUE),
    posterior_lower = quantile(estimate, 0.025, na.rm = TRUE),
    posterior_upper = quantile(estimate, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::rename(
    estimate = posterior_median,
    lower = posterior_lower,
    upper = posterior_upper
  )

# Combine biodiversity and environmental colour palettes
all_colours <- c(
  biotic_colours,
  env_colours[interaction_env]
)

# Legend order
legend_breaks <- c(
  biotic_effects,
  interaction_env
)

legend_labels <- c(
  unname(biotic_labels[biotic_effects]),
  "Temperature",
  "Salinity",
  "Total Nitrogen"
)

# Alternating background shading for biodiversity metrics
background_data <- data.frame(
  xmin = seq(0.5, length(biotic_effects) - 0.5, by = 1),
  xmax = seq(1.5, length(biotic_effects) + 0.5, by = 1),
  ymin = -Inf,
  ymax = Inf,
  shade = rep(c("white", "grey90"), length.out = length(biotic_effects))
)

# Plot all main and conditional effects together
RUE_n_biotic <- ggplot(
  interaction_plot_data,
  aes(
    x = x,
    y = estimate,
    fill = plot_colour,
    colour = plot_colour,
    group = group_id,
    linetype = condition
  )
) +
  # Alternating background shading
  geom_rect(
    data = background_data,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = shade
    ),
    inherit.aes = FALSE,
    alpha = 1
  ) +
  scale_fill_identity() +
  ggnewscale::new_scale_fill() +
  # Zero reference line
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey60",
    linewidth = 0.4
  ) +
  # Main effects: half-violin distributions only
  see::geom_violinhalf(
    data = subset(
      interaction_plot_data,
      condition == "Main effect"
    ),
    aes(
      fill = plot_colour,
      colour = plot_colour,
      alpha = condition
    ),
    position = "identity",
    flip = FALSE,
    scale = "width",
    trim = FALSE,
    width = 0.15,
    linewidth = 0.3
  ) +
  # Transparency for main effects
  scale_alpha_manual(
    values = c(
      "Main effect" = 1,
      "Low" = 0.30,
      "High" = 0.30
    ),
    guide = "none"
  ) +
  # Environmental interactions:
  # Posterior median and 95% credible intervals
  geom_pointrange(
    data = interaction_summary,
    aes(
      x = x,
      y = estimate,
      ymin = lower,
      ymax = upper,
      colour = plot_colour,
      linetype = condition
    ),
    inherit.aes = FALSE,
    linewidth = 0.7,
    fatten = 2
  ) +
  # Significance stars for main effects
  geom_text(
    data = subset(
      RUEn_star_labels,
      stars != ""
    ),
    aes(
      x = x,
      y = y_pos,
      label = stars
    ),
    inherit.aes = FALSE,
    colour = "red",
    size = 5,
    vjust = 0
  ) +
  # Fill colours
  scale_fill_manual(
    values = all_colours,
    breaks = interaction_env,
    labels = c(
      "Temperature",
      "Salinity",
      "Total Nitrogen"
    ),
    name = "Environmental interactions"
  ) +
  # Line and point colours
  scale_colour_manual(
    values = all_colours,
    breaks = interaction_env,
    labels = c(
      "Temperature",
      "Salinity",
      "Total Nitrogen"
    ),
    name = "Environmental interactions"
  ) +
  # Y-axis
  scale_y_continuous(
    limits = c(-0.7, 0.6)
  ) +
  # Environmental condition line types
  scale_linetype_manual(
    values = c(
      "Main effect" = "solid",
      "Low" = "dashed",
      "High" = "solid"
    ),
    name = "Environmental level"
  ) +
  # X-axis
  scale_x_continuous(
    breaks = seq_along(biotic_effects),
    labels = unname(biotic_labels[biotic_effects]),
    expand = expansion(add = c(0.4, 0.4))
  ) +
  # Axis labels
  labs(
    x = NULL,
    y = NULL
  ) +
  # Theme
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    plot.title = element_text(size = 12, hjust = 0.5),
    legend.position = "none",
    plot.margin = margin(3, 5, 3, 5)
  )

# Display
RUE_n_biotic

setwd("~/GitHub/PlanktonRUE-Project/Figures/Final_Figure")
svg("PD_Phyto_RUE_N.svg", width = 8.83, height = 5, pointsize = 10)

RUE_n_biotic

dev.off()

##  Zoo RUE model .......................................................
setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA")
load('BRM_zoo_RUE_CD_PI_261013_All.RDATA')

# Biotic main effects only
biotic_effects <- c(
  "genus_richness",
  "phyto_evenness",
  "phyto_FEve_ZRUE",
  "phyto_FDis_ZRUE",
  "cyano_prop",
  "species_richness",
  "zoo_evenness",
  "zoo_FEve",
  "zoo_FDis"
)

biotic_labels <- c(
  "genus_richness"   = "Phytoplankton\nGenus\nRichness",
  "phyto_evenness"   = "Phytoplankton\nEvenness",
  "phyto_FEve_ZRUE"  = "Phytoplankton\nFEve",
  "phyto_FDis_ZRUE"  = "Phytoplankton\nFDis",
  "cyano_prop"       = "Cyanobacteria\nDominance",
  "species_richness" = "Zooplankton\nSpecies\nRichness",
  "zoo_evenness"     = "Zooplankton\nEvenness",
  "zoo_FEve"         = "Zooplankton\nFEve",
  "zoo_FDis"         = "Zooplankton\nFDis"
)

# Colours in the same order
biotic_colours <- c(
  "genus_richness"   = "#278B9AFF",
  "phyto_evenness"   = "#274637FF",
  "phyto_FEve_ZRUE"  = "#44A57CFF",
  "phyto_FDis_ZRUE"  = "#819A7AFF",
  "cyano_prop"       = "#58A449FF",
  "species_richness" = "#403369FF",
  "zoo_evenness"     = "#9F99B5FF",
  "zoo_FEve"         = "#5C5992FF",
  "zoo_FDis"         = "#AE93BEFF"
)

#set CI stars
zRUE_star_labels <- data.frame(
  biotic = biotic_effects,
  stars = c("**", "***", "", "", "", "", "***", "", "***"),
  y_pos = rep(1.2, length(biotic_effects))
)

zRUE_star_labels$x <- match(
  zRUE_star_labels$biotic,
  biotic_effects
)

interaction_env <- c(
  "Temperature",
  "Salinity",
  "TotalNitrogen_unfiltered"
)

env_levels <- c("Main effect", interaction_env)

env_colours <- c(
  "Main effect" = "grey25",
  "Temperature" = "red",
  "Salinity" = "orange",
  "TotalNitrogen_unfiltered" = "yellow"
)

# Extract posterior draws
draws <- posterior::as_draws_df(model)

# Find coefficient names safely, regardless of interaction order
find_coef <- function(a, b = NULL) {
  
  candidates <- grep("^b_", names(draws), value = TRUE)
  
  matches <- vapply(candidates, function(nm) {
    
    parts <- strsplit(
      sub("^b_", "", nm),
      ":",
      fixed = TRUE
    )[[1]]
    
    if (is.null(b)) {
      length(parts) == 1 && parts[1] == a
    } else {
      length(parts) == 2 &&
        setequal(parts, c(a, b))
    }
    
  }, logical(1))
  
  found <- candidates[matches]
  
  if (length(found) != 1) {
    stop(
      "Could not uniquely identify coefficient: ",
      a,
      if (!is.null(b)) paste0(":", b)
    )
  }
  
  found
}

# Horizontal offsets around each biodiversity main effect
# Negative offsets = low conditions; positive = high conditions
offsets <- data.frame(
  env = c(
    "Temperature", "Temperature",
    "Salinity", "Salinity",
    "TotalNitrogen_unfiltered",
    "TotalNitrogen_unfiltered"
  ),
  condition = c(
    "Low", "High",
    "Low", "High",
    "Low", "High"
  ),
  offset = c(-0.30, 0.30, -0.20, 0.20, -0.10, 0.10)
)

# Calculate posterior distributions
plot_list <- list()

for (i in seq_along(biotic_effects)) {
  
  bio <- biotic_effects[i]
  
  main_col <- find_coef(bio)
  
  # Main effect
  plot_list[[length(plot_list) + 1]] <- data.frame(
    x = rep(i, nrow(draws)),
    estimate = draws[[main_col]],
    env = "Main effect",
    condition = "Main effect",
    biotic = bio
  )
  
  # Conditional effects at low and high environmental values
  for (env in interaction_env) {
    
    int_col <- find_coef(bio, env)
    
    for (level in c("Low", "High")) {
      
      env_value <- if (level == "Low") -1 else 1
      
      # Select exactly one horizontal offset
      offset_index <- which(
        offsets$env == env &
          offsets$condition == level
      )
      
      if (length(offset_index) != 1) {
        stop("Expected one offset for ", env, " / ", level)
      }
      
      pos <- offsets$offset[offset_index]
      
      # Conditional posterior draws
      conditional_draws <-
        draws[[main_col]] +
        env_value * draws[[int_col]]
      
      plot_list[[length(plot_list) + 1]] <- data.frame(
        x = rep(i + pos, length(conditional_draws)),
        estimate = conditional_draws,
        env = env,
        condition = level,
        biotic = bio
      )
    }
  }
}

# Combine all posterior draws
interaction_plot_data <- dplyr::bind_rows(plot_list) %>%
  dplyr::mutate(
    env = factor(
      env,
      levels = env_levels
    ),
    condition = factor(
      condition,
      levels = c("Main effect", "Low", "High")
    ),
    
    # Biodiversity-specific colour for main effects;
    # environmental colour for conditional effects
    plot_colour = dplyr::if_else(
      as.character(env) == "Main effect",
      biotic,
      as.character(env)
    ),
    
    group_id = interaction(biotic, env, condition)
  )

# Calculate posterior medians and 95% credible intervals
interaction_summary <- interaction_plot_data %>%
  dplyr::filter(condition %in% c("Low", "High")) %>%
  dplyr::group_by(biotic, env, condition, x, plot_colour) %>%
  dplyr::summarise(
    posterior_median = median(estimate, na.rm = TRUE),
    posterior_lower = quantile(estimate, 0.025, na.rm = TRUE),
    posterior_upper = quantile(estimate, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::rename(
    estimate = posterior_median,
    lower = posterior_lower,
    upper = posterior_upper
  )

# Combine biodiversity and environmental colour palettes
all_colours <- c(
  biotic_colours,
  env_colours[interaction_env]
)

# Legend order
legend_breaks <- c(
  biotic_effects,
  interaction_env
)

legend_labels <- c(
  unname(biotic_labels[biotic_effects]),
  "Temperature",
  "Salinity",
  "Total Nitrogen"
)

# Alternating background shading for biodiversity metrics
background_data <- data.frame(
  xmin = seq(0.5, length(biotic_effects) - 0.5, by = 1),
  xmax = seq(1.5, length(biotic_effects) + 0.5, by = 1),
  ymin = -Inf,
  ymax = Inf,
  shade = rep(c("white", "grey90"), length.out = length(biotic_effects))
)

# Plot all main and conditional effects together
ZRUE_biotic <- ggplot(
  interaction_plot_data,
  aes(x = x, y = estimate,
      fill = plot_colour, colour = plot_colour, group = group_id, linetype = condition
  )
) +
  geom_rect(
    data = background_data,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = shade
    ),
    inherit.aes = FALSE,
    alpha = 1
  ) +
  scale_fill_identity() +
  ggnewscale::new_scale_fill() +
  geom_hline(
    yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4
  ) +
  see::geom_violinhalf(
    data = subset(
      interaction_plot_data,
      condition == "Main effect"
    ),
    aes(
      fill = plot_colour,
      colour = plot_colour,
      alpha = condition
    ),
    position = "identity",
    flip = FALSE,
    scale = "width",
    trim = FALSE,
    width = 0.15,
    linewidth = 0.3
  ) +
  scale_alpha_manual(
    values = c(
      "Main effect" = 1,"Low" = 0.30,"High" = 0.30
    ),
    guide = "none"
  ) +
  geom_pointrange(
    data = interaction_summary,
    aes(
      x = x,
      y = estimate,
      ymin = lower,
      ymax = upper,
      colour = plot_colour,
      linetype = condition
    ),
    inherit.aes = FALSE,
    linewidth = 0.7,
    fatten = 2
  ) +
  geom_text(
    data = subset(zRUE_star_labels, stars != ""),
    aes(
      x = x,
      y = y_pos,
      label = stars
    ),
    inherit.aes = FALSE,colour = "red",size = 5,vjust = 0
  ) +
  scale_fill_manual(
    values = all_colours,
    breaks = interaction_env,
    labels = c("Temperature", "Salinity", "Total Nitrogen"),
    name = "Environmental interactions"
  ) +
  scale_colour_manual(
    values = all_colours,
    breaks = interaction_env,
    labels = c("Temperature", "Salinity", "Total Nitrogen"),
    name = "Environmental interactions"
  ) +
  scale_y_continuous(
    limits = c(-1.3, 1.3)
  ) +
  scale_linetype_manual(
    values = c(
      "Main effect" = "solid", "Low" = "dashed", "High" = "solid"
    ),
    name = "Environmental level"
  ) +
  scale_x_continuous(
    breaks = seq_along(biotic_effects),
    labels = unname(biotic_labels[biotic_effects]),
    expand = expansion(add = c(0.4, 0.4))
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    plot.title = element_text(size = 12, hjust = 0.5),
    legend.position = "none",
    plot.margin = margin(3, 5, 3, 5)
  )

# Display
ZRUE_biotic

setwd("~/GitHub/PlanktonRUE-Project/Figures/Final_Figure")
svg("PD_Zoo_RUE.svg", width = 8.83, height = 5, pointsize = 10)

ZRUE_biotic

dev.off()

################################################################################
### Plotting final figure ##################################################

# First row: three variance-partitioning plots
vp_row <- wrap_plots(
  RUE_p_vp, RUE_n_vp, RUE_zoo_vp,
  nrow = 1,
  ncol = 3
)

# Second row: three posterior plots
posterior_row <- wrap_plots(
  RUE_p_biotic + theme(legend.position = "bottom"),
  RUE_n_biotic + theme(legend.position = "bottom"),
  ZRUE_biotic + theme(legend.position = "bottom"),
  nrow = 1,
  ncol = 3,
  guides = "collect"
) +
  plot_layout(guides = "collect") +
  plot_annotation(
    theme = theme(legend.position = "bottom")
  )

# Combine the two rows
combined_plot <- vp_row / posterior_row +
  plot_layout(heights = c(1, 0.8))

setwd("~/GitHub/PlanktonRUE-Project/Figures")
ggsave(
  filename = "Final_Fig2.png",
  plot = combined_plot,
  width = 34,
  height = 15,
  units = "in",
  dpi = 300,
  bg = "white"
)


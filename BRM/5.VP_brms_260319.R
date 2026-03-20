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


##...........................................................................................................

##.........................................................

## Get the functions
## ..........................................................................
## ............................................................................

## new and adapted functions
setwd("~/GitHub/PlanktonRUE-Project/BRM")
source('Functions_for_brms_VP.R')

## VP function
run_vp <- function(model, groups, out_prefix, out_dir,
                   main_ng = "VP (no grouping)",
                   main_g  = "VP (grouped)",
                   cols_g = NULL,
                   base_png_ng = 4000,
                   base_png_g  = 4000,
                   res = 600) {
  
  # Prepare model
  mod <- PrepbrmsVP(model)
  colnames(mod$X)
  params <- mod$params
  X_full <- mod$X
  
  linear_terms <- unname(unlist(groups))
  linear_terms <- linear_terms[linear_terms %in% colnames(X_full)]
  
  if (length(linear_terms) == 0) {
    stop("No matching linear terms found.")
  }
  
  X <- X_full[, linear_terms, drop = FALSE]
  n_lt <- length(linear_terms)
  y <- mod$y
  n_iter <- nrow(params)
  
  # Flatten linear terms
  linear_terms <- unname(unlist(groups))
  n_lt <- length(linear_terms)
  n_group <- length(groups)
  
  ## --- No grouping ---
  # Adjust PNG size based on number of terms
  width_ng <- max(base_png_ng, n_lt * 300)
  height_ng <- max(base_png_ng, n_lt * 300)
  
  # Compute covariance matrix (lin_term_cov uses linear_terms from this environment)
  K <- apply(params, "iterations", function(par) lin_term_cov(par, X, y))
  dim(K) <- c(n_lt, n_lt, n_iter)
  dimnames(K) <- list(rows = linear_terms, cols = linear_terms, iterations = NULL)
  
  P_x <- apply(K, "iterations", function(x) diag(x) / sum(diag(x)))
  dimnames(P_x) <- list(variance = linear_terms, iterations = NULL)
  
  # Save no-grouping plot
  png(file.path(out_dir, paste0(out_prefix, "_ng.png")),
      width = width_ng, height = height_ng, res = res)
  op <- par(mai = c(3.5, 0.8, 0.8, 0.4), las = 2, cex.axis = 1.1, xpd = NA)
  vp_plot(P_x, nm = linear_terms, main = main_ng, las = 2)
  par(op)
  dev.off()
  
  ## --- Grouped ---
  B <- sapply(groups, function(x) as.integer(linear_terms %in% x))
  dimnames(B) <- list(linear_terms = linear_terms, groups = colnames(B))
  
  K_B <- apply(K, "iterations", function(x, B) t(B) %*% x %*% B, B)
  dim(K_B) <- c(n_group, n_group, n_iter)
  dimnames(K_B) <- list(rows = names(groups), cols = names(groups), iterations = NULL)
  
  P <- apply(K_B, "iterations", function(x) diag(x) / sum(diag(x)))
  dimnames(P) <- list(variance = names(groups), iterations = NULL)
  
  # Adjust PNG size for grouped plot
  width_g <- max(base_png_g, n_group * 300)
  height_g <- base_png_g
  
  png(file.path(out_dir, paste0(out_prefix, ".png")),
      width = width_g, height = height_g, res = res)
  op <- par(mai = c(1.5, 0.8, 0.8, 0.4), las = 2, cex.axis = 1.1, xpd = NA)
  if (is.null(cols_g)) {
    vp_plot(P, nm = names(groups), main = main_g, las = 2)
  } else {
    vp_plot(P, col = cols_g, nm = names(groups), main = main_g, las = 2)
  }
  par(op)
  dev.off()
  
  invisible(list(P_x = P_x, P = P))
}



#.............................................................................
#Model list:

base_dir <- "~/GitHub/PlanktonRUE-Project/BRM/RDATA"
out_dir  <- "~/GitHub/PlanktonRUE-Project/Figures/VP_plots"

configs <- list(
  
  cyano_dom = list(
    file = "BRM_cyanodom_P_260311_Z.RDATA",
    groups = list(
      Abiotic = c("polyTemperature2rawEQTRUE1","polyTemperature2rawEQTRUE2",
                  "polySalinity2rawEQTRUE1","polySalinity2rawEQTRUE2",
                  "polypH2rawEQTRUE1","polypH2rawEQTRUE2",
                  "polyMaxDepth2rawEQTRUE1", "polyMaxDepth2rawEQTRUE2"),
      DOY = c("polyDOY2rawEQTRUE1","polyDOY2rawEQTRUE2",
              "polydaylight_hours2rawEQTRUE1", "polydaylight_hours2rawEQTRUE2"),
      Site = "Site_Id",
      Year = "Year",
      Residuals = "residuals"
    ),
    out_prefix = "cyanobacteria_dominance",
    main_ng = "Cyanobacteria Dominance (ng)",
    main_g  = "Cyanobacteria Dominance",
    cols_g = c('royalblue3','violetred3','salmon2','gold3','black')
  ),
  
  phyto_GD = list(
    file = "BRM_phyto_GD_P_260310.RDATA",
    groups = list(
      Abiotic = c("polyTemperature2rawEQTRUE1","polyTemperature2rawEQTRUE2",
                  "polySalinity2rawEQTRUE1","polySalinity2rawEQTRUE2",
                  "polypH2rawEQTRUE1","polypH2rawEQTRUE2",
                  "polyMaxDepth2rawEQTRUE1", "polyMaxDepth2rawEQTRUE2"),
      DOY = c("polyDOY2rawEQTRUE1","polyDOY2rawEQTRUE2",
              "polydaylight_hours2rawEQTRUE1", "polydaylight_hours2rawEQTRUE2"),
      Site = "Site_Id",
      Year = "Year",
      Residuals = "residuals"
    ),
    out_prefix = "Phyto_GD",
    main_ng = "Phytoplankton Genus Richness (ng)",
    main_g  = "Phytoplankton Genus Richness",
    cols_g = c('royalblue3','violetred3','salmon2','gold3','black')
  ),
  
  phyto_FEve = list(
    file = "BRM_phyto_FE_P_260310.RDATA",
    groups = list(
      Abiotic = c("polyTemperature2rawEQTRUE1","polyTemperature2rawEQTRUE2",
                  "polySalinity2rawEQTRUE1","polySalinity2rawEQTRUE2",
                  "polypH2rawEQTRUE1","polypH2rawEQTRUE2",
                  "polyMaxDepth2rawEQTRUE1", "polyMaxDepth2rawEQTRUE2"),
      DOY = c("polyDOY2rawEQTRUE1","polyDOY2rawEQTRUE2",
              "polydaylight_hours2rawEQTRUE1", "polydaylight_hours2rawEQTRUE2"),
      Site = "Site_Id",
      Year = "Year",
      Residuals = "residuals"
    ),
    out_prefix = "Phyto_FEve",
    main_ng = "Phytoplankton Functional Evenness (ng)",
    main_g  = "Phytoplankton Functional Evenness",
    cols_g = c('royalblue3','violetred3','salmon2','gold3','black')
  ),
  
  phyto_FDis = list(
    file = "BRM_phyto_FD_P_260310.RDATA",
    groups = list(
      Abiotic = c("polyTemperature2rawEQTRUE1","polyTemperature2rawEQTRUE2",
                  "polySalinity2rawEQTRUE1","polySalinity2rawEQTRUE2",
                  "polypH2rawEQTRUE1","polypH2rawEQTRUE2",
                  "polyMaxDepth2rawEQTRUE1", "polyMaxDepth2rawEQTRUE2"),
      DOY = c("polyDOY2rawEQTRUE1","polyDOY2rawEQTRUE2",
              "polydaylight_hours2rawEQTRUE1", "polydaylight_hours2rawEQTRUE2"),
      Site = "Site_Id",
      Year = "Year",
      Residuals = "residuals"
    ),
    out_prefix = "Phyto_FDis",
    main_ng = "Phytoplankton Dispersion (ng)",
    main_g  = "Phytoplankton Dispersion",
    cols_g = c('royalblue3','violetred3','salmon2','gold3','black')
  ),
  
  phyto_RUE_P = list(
    file = "BRM_phytoRUEP_CD_PI_260311.RDATA",
    groups = list(
      Abiotic = c("Temperature","ITemperatureE2",
                  "Salinity","ISalinityE2",
                  "pH","IpHE2",
                  "MaxDepth","IMaxDepthE2"),
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
      DOY = c("DOY","IDOYE2",
              "daylight_hours", "Idaylight_hoursE2"),
      Site = "Site_Id",
      Year = "Year",
      Residuals = "residuals"
    ),
    out_prefix = "Phyto_RUEP_CD",
    main_ng = "Phytoplankton RUE per P (ng)",
    main_g  = "Phytoplankton RUE per P",
    cols_g = c('royalblue3','#278B9AFF','#1D2645FF','#44A57CFF','#819A7AFF','#58A449FF','#CEC917FF','violetred3','salmon2','gold3','black')
  ),
  
  phyto_RUE_N = list(
    file = "BRM_phytoRUEN_CD_PI_260311.RDATA",
    groups = list(
      Abiotic = c("Temperature","ITemperatureE2",
                  "Salinity","ISalinityE2",
                  "pH","IpHE2",
                  "MaxDepth","IMaxDepthE2"),
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
      DOY = c("DOY","IDOYE2",
              "daylight_hours", "Idaylight_hoursE2"),
      Site = "Site_Id",
      Year = "Year",
      Residuals = "residuals"
    ),
    out_prefix = "Phyto_RUEN_CD",
    main_ng = "Phytoplankton RUE per N (ng)",
    main_g  = "Phytoplankton RUE per N",
    cols_g = c('royalblue3','#278B9AFF','#1D2645FF','#44A57CFF','#819A7AFF','#58A449FF','#CEC917FF','violetred3','salmon2','gold3','black')
  ),
  
  Zoo_SR = list(
    file = "BRM_zoo_SR_CD_P_260313.RDATA",
    groups = list(
      Abiotic = c("polyTemperature2rawEQTRUE1","polyTemperature2rawEQTRUE2",
                  "polySalinity2rawEQTRUE1","polySalinity2rawEQTRUE2",
                  "polypH2rawEQTRUE1","polypH2rawEQTRUE2",
                  "polyZooplankton_MaxDepth2rawEQTRUE1", "polyZooplankton_MaxDepth2rawEQTRUE2"),
      Phyto_div = c("genus_richness", "FEve", "FDis"),
      Cyano_dom = c("cyano_prop"),
      DOY = c("polyDOY2rawEQTRUE1","polyDOY2rawEQTRUE2",
              "polydaylight_hours2rawEQTRUE1", "polydaylight_hours2rawEQTRUE2"),
      Site = "Site_Id",
      Year = "Year",
      Residuals = "residuals"
    ),
    out_prefix = "Zoo_SR_CD",
    main_ng = "Zooplankton Species Richness (ng)",
    main_g  = "Zooplankton Species Richness",
    cols_g = c('royalblue3','#278B9AFF','#44A57CFF','violetred3','salmon2','gold3','black')
  ),
  
  Zoo_FEve = list(
    file = "BRM_zoo_FE_CD_P_260313.RDATA",
    groups = list(
      Abiotic = c("polyTemperature2rawEQTRUE1","polyTemperature2rawEQTRUE2",
                  "polySalinity2rawEQTRUE1","polySalinity2rawEQTRUE2",
                  "polypH2rawEQTRUE1","polypH2rawEQTRUE2",
                  "polyZooplankton_MaxDepth2rawEQTRUE1", "polyZooplankton_MaxDepth2rawEQTRUE2"),
      Phyto_div = c("genus_richness", "FEve", "FDis"),
      Cyano_dom = c("cyano_prop"),
      DOY = c("polyDOY2rawEQTRUE1","polyDOY2rawEQTRUE2",
              "polydaylight_hours2rawEQTRUE1", "polydaylight_hours2rawEQTRUE2"),
      Site = "Site_Id",
      Year = "Year",
      Residuals = "residuals"
    ),
    out_prefix = "Zoo_FE_CD",
    main_ng = "Zooplankton Functional Evenness (ng)",
    main_g  = "Zooplankton Functional Evenness",
    cols_g = c('royalblue3','#278B9AFF','#44A57CFF','violetred3','salmon2','gold3','black')
  ),
  
  Zoo_FDis = list(
    file = "BRM_zoo_FD_CD_P_260313.RDATA",
    groups = list(
      Abiotic = c("polyTemperature2rawEQTRUE1","polyTemperature2rawEQTRUE2",
                  "polySalinity2rawEQTRUE1","polySalinity2rawEQTRUE2",
                  "polypH2rawEQTRUE1","polypH2rawEQTRUE2",
                  "polyZooplankton_MaxDepth2rawEQTRUE1", "polyZooplankton_MaxDepth2rawEQTRUE2"),
      Phyto_div = c("genus_richness", "FEve", "FDis"),
      Cyano_dom = c("cyano_prop"),
      DOY = c("polyDOY2rawEQTRUE1","polyDOY2rawEQTRUE2",
              "polydaylight_hours2rawEQTRUE1", "polydaylight_hours2rawEQTRUE2"),
      Site = "Site_Id",
      Year = "Year",
      Residuals = "residuals"
    ),
    out_prefix = "Zoo_FD_CD",
    main_ng = "Zooplankton Functional Dispersion (ng)",
    main_g  = "Zooplankton Functional Dispersion",
    cols_g = c('royalblue3','#278B9AFF','#44A57CFF','violetred3','salmon2','gold3','black')
  ),
  
  Zoo_RUE = list(
    file = "BRM_zoo_RUE_CD_PI_260313.RDATA",
    groups = list(
      Abiotic = c("Temperature","ITemperatureE2",
                  "Salinity","ISalinityE2",
                  "pH","IpHE2",
                  "MaxDepth","IMaxDepthE2"),
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
      DOY = c("DOY","IDOYE2",
              "daylight_hours", "Idaylight_hoursE2"),
      Site = "Site_Id",
      Year = "Year",
      Residuals = "residuals"
    ),
    out_prefix = "Zoo_RUE_CD",
    main_ng = "Zooplankton RUE (ng)",
    main_g  = "Zooplankton RUE",
    cols_g = c('royalblue3','#278B9AFF','#1D2645FF','#44A57CFF','#819A7AFF','#58A449FF','#CEC917FF','violetred3','salmon2','gold3','black')
  )
)

#..............................................................................
#run everything#
setwd(base_dir)

for (nm in names(configs)) {
  cfg <- configs[[nm]]
  
  message("Running VP for: ", nm)
  
  load(cfg$file)   # loads object: model
  
  run_vp(
    model = model,
    groups = cfg$groups,
    out_prefix = cfg$out_prefix,
    out_dir = out_dir,
    main_ng = cfg$main_ng,
    main_g  = cfg$main_g,
    cols_g  = cfg$cols_g
  )
}


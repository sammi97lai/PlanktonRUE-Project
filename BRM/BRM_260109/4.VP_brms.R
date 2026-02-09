#
#               General script to plot results for brms results (variance partition)
#                         08.12.2025
#                   Hurdle models  - with vegetation information & Salinity
#                         Emy Guilbault
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
#library(brms)

## plotting
library("ggplot2")
library(vioplot)
library(ggpubr)


##...........................................................................................................

##.........................................................
setwd("~/GitHub/Phytoplankton-project/phyto_RUE/BRM_260109")

## get the model
##.....................................................................................................
##.....................................................................................................

##  models
model <- readRDS("BRM_phyto_GD_PI_260109.rds")

# 1- we extract the fixed and random effect coefficients:
# Fixed part 
fixeff_mod = fixef(model, summary = FALSE) ## mcmc / coef ab coef pa
dim(fixeff_mod) # 8000   15

# random part
ranef_mod = ranef(model, summary = FALSE) ##list of random effec for which mcmcm / groups of rd effects / ab pa
#length(ranef_mod) # 3
#dim(ranef_mod[[1]]) # 8000  8   1 
#dim(ranef_mod[[2]]) # 8000  25   1 
#dim(ranef_mod[[2]]) # 8000  22   2 

## Since we have a hurdle model there are two separate parts occurrence (hu_ part) and abundance (by default)
# prep beta - coef for fixed effects
beta = fixeff_mod
colnames(beta)

# prep lambda - coef for random effects
# for the combine randomf effects we need intercept and slope coefficients side by side as column
# rd3 = aperm(ranef_mod[[3]], c(1,3,2))
# rd3_mat = matrix(rd3, nrow = dim(ranef_mod[[3]])[1], byrow = FALSE)
# colnames(rd3_mat) = rep(dimnames(ranef_mod[[3]])[[2]], each = 2)

rd3_mat <- ranef_mod[[3]][,,1]
colnames(rd3_mat) <- dimnames(ranef_mod[[3]])[[2]]

# lambda = cbind(ranef_mod[[2]][,,1], rd3_mat,  
#                ranef_mod[[1]][,,1]) # just intercept part..?
lambda = cbind(
  ranef_mod[[1]][,,1],  # Site_Id
  ranef_mod[[2]][,,1],  # Year
  rd3_mat               # MaxDepth
)
colnames(lambda) # 

# select the names of the covariates
## there are the same for the two parts of the hurdle model so this is fine
cov_names =  colnames(fixeff_mod)

## define some covariates to group
# groups <- list(
#   abiotic = c("Temperature", "Salinity", "pH", "log_N_P_ratio"), ## 
#   phyto_div = c("genus_richness", "FEve", "FDis", "cyano_prop"),  ## 
#   DOY = c('DOY'),
#   site = c('Site_Id'),
#   year = c('Year'),
#   maxdepth = c('MaxDepth')#,
#   #residuals = 'residuals'
# ) #group variables into categories#
groups <- list(
  abiotic = grep("^poly(Temperature|Salinity|pH|log_N_P_ratio)", colnames(X_Beta), value = TRUE),
  phyto_div = c("genus_richness", "FEve", "FDis", "cyano_prop"),
  DOY = grep("DOY", colnames(X_Beta), value = TRUE),
  site = "Site_Id",
  year = "Year",
  maxdepth = "MaxDepth"
)

# define some needed parameters
linear_terms <- unname(unlist(groups)) #flattened vector of all variable names
n_lt <- length(linear_terms)
n_group <- length(groups)


# 2-  we prepare the data matrix with all environemntal and random covariates values
# this cannot be accessed directly within brms so we build it based on the model formula and the data fed into the model
form_mod = model$formula

## envrionmental variable matrix
foo <- lFormula(eval(form_mod), model$data)
X_Beta <- as.matrix(foo$X)
colnames(X_Beta)[1] =  paste('(Intercept)')

## random variable matrix
Eta <- t(as.matrix(foo$reTrms$Zt))
colnames(Eta)  
## to fitEta## to fit the existing function we used in rstan we need to rename this to match
## we still have the same random effects but if different random effect both the names here (area2, year) and in the function need to change
# prepare the column names to match the function
# Lsite = length(unique(model$data$Site_Id))
# Lyear = length(unique(model$data$Year)) ## we need it twice because of the intercept and slope
# Lmaxd = length(unique(model$data$MaxDepth))

# colnames(Eta)[1:Lsite] = paste('b[(Intercept) Site_Id', colnames(Eta)[1:Lsite], sep=' ')
# colnames(Eta)[(Lsite+1):(Lsite+Lyear)] =
#   paste('b[(Intercept) Year', colnames(Eta)[(Lsite+1):(Lsite+Lyear)], sep=' ')
# # colnames(Eta)[(Lsite+1):(Lsite+2*Lyear)] = paste('b[(Intercept) Year', colnames(Eta)[(Lsite+1):(Lsite+2*Lyear)], sep=' ')
# # colnames(Eta)[(Lsite+2*Lyear+1):(Lsite+2*Lyear+Lmaxd)] = paste('b[(Intercept) MaxDepth', colnames(Eta)[(Lsite+2*Lyear+1):(Lsite+2*Lyear+Lmaxd)], sep=' ')
# colnames(Eta)[(Lsite+Lyear+1):(Lsite+Lyear+Lmaxd)] = paste('b[(Intercept) MaxDepth', colnames(Eta)[(Lsite+Lyear+1):(Lsite+Lyear+Lmaxd)], sep=' ')
site_idx <- grep("Site", colnames(Eta))
year_idx <- grep("Year", colnames(Eta))
maxd_idx <- grep("MaxDepth", colnames(Eta))

colnames(Eta)[site_idx] <- paste('b[(Intercept) Site_Id', colnames(Eta)[site_idx], sep=' ')
colnames(Eta)[year_idx] <- paste('b[(Intercept) Year', colnames(Eta)[year_idx], sep=' ')
colnames(Eta)[maxd_idx] <- paste('b[(Intercept) MaxDepth', colnames(Eta)[maxd_idx], sep=' ')

### prepare the parameters and matrix to be used in Torsti's function

X <- cbind(X_Beta, Eta) ## combine environmental part and random effects
y <- foo$fr[1] #extract y#


params <- as.matrix(cbind(beta, lambda)) #converts to matrix
colnames(params) = colnames(X)
dimnames(params) = list('iterations' = 1:dim(params)[1], 
                        'parameters' = colnames(params))
n_iter <- nrow(params)


## VP functions adapted from Torsti .............................................

## Functions .........................

#builds matrix A of linear terms bu combining effects from covariates, random effects (ID and year), and residuals
#represent specific contributions
construct_lin_term_mat <- function(params, X, y) {
  
  ## Linear terms for the growth rate and density dependence
  X_cov <- X[, !grepl('(Intercept)', colnames(X), fixed = TRUE)]
  A_cov <- as.matrix(X_cov %*% Matrix::Diagonal(x = params[colnames(X_cov)]))
  colnames(A_cov) <- colnames(X_cov)
  
  ## Linear term for the ID random effect
  X_iid <- X[, grep('b[(Intercept) Site_Id', colnames(X), fixed = TRUE)]
  A_iid <- X_iid %*% params[colnames(X_iid)]
  colnames(A_iid) <- c('Site_Id')
  
  ## Linear term for the yearly random effect
  X_iid2 <- X[, grep('b[(Intercept) Year', colnames(X), fixed = TRUE)]
  A_iid2 <- X_iid2 %*% params[colnames(X_iid2)]
  colnames(A_iid2) <- c('Year')
  
  ## Linear term for the yearly random effect
  X_iid3 <- X[, grep('b[(Intercept) MaxDepth', colnames(X), fixed = TRUE)]
  A_iid3 <- X_iid3 %*% params[colnames(X_iid3)]
  colnames(A_iid3) <- c('MaxDepth')
  
  A <- cbind(A_cov, A_iid, A_iid2, A_iid3)
  
  ## Linear term for residuals (omiting the intercept above does not affect variance)
  #A <- cbind(A, residuals = y - rowSums(A))
  
  return(as.matrix(A))
}


## Covariance matrix
#calculate covariance matrix of linear terms
#measures relationship between linear terms
lin_term_cov <- function(params, X, y) {
  A <- construct_lin_term_mat(params, X, y)
  A <- A[,linear_terms]
  K <- cov(A)
  dimnames(K) <- list(colnames(A), colnames(A))
  return(K)
}


vp_plot <- function(x, col = 'grey80', nm, main, las =1) {
  vp <- as.data.frame(aperm(x, c('iterations', 'variance')))
  colnames(vp) <- nm
  
  if (max(vp) <= 1 & min(vp) >= 0) {
    ylim <- 0:1
  } else {
    ylim <- NULL
  }
  
  vioplot(
    vp,
    col = col,
    las=las,
    drawRect = TRUE,
    pchMed = NA,
    ylim = ylim,
    main = main
  )
}


## VP functions adapted from Torsti .............................................


### Variance partition $V_x$
#quantify how much variance each term explains relative to total variance

# A # No grouping
K <- apply(params, 'iterations', lin_term_cov, X, y)
dim(K) <- c(n_lt, n_lt, n_iter)
dimnames(K) <- list(rows = linear_terms, cols = linear_terms, iterations = NULL)

P_x <- apply(K, 'iterations', function(x) {diag(x) / sum(diag(x))})
dimnames(P_x) <- list(variance = linear_terms, iterations = NULL)

vp_plot(P_x, nm = linear_terms, main = 'no grouping', las = 2)



# B # Variance partitioning for groups of linear terms
### Variance partition $V_b$

B <- sapply(groups, function(x) { as.integer(linear_terms %in% x) })
dimnames(B) <- list(linear_terms = linear_terms, groups = colnames(B))

K_B <- apply(K, 'iterations', function(x, B) { t(B) %*% x %*% B }, B)
dim(K_B) <- c(n_group, n_group, n_iter)
dimnames(K_B) <- list(rows = names(groups), cols = names(groups), iterations = NULL)

P <- apply(K_B, 'iterations', function(x) {diag(x) / sum(diag(x))})
dimnames(P) <- list(variance = names(groups), iterations = NULL)

vp_plot(P, col = c('turquoise4', 'purple4', 'orange2', 'beige', 'red4', 'gray30'),
        nm = names(groups), main = 'grouping', las=2)


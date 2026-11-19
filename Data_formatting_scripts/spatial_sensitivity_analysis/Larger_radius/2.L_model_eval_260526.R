################################################################################
###################### Model Evaluation and Results ############################
################################################################################

library(dplyr)
library(tidyr)
library(brms)
library(ggplot2)
library(posterior)
library(bayesplot)
library(patchwork)
library(loo)
library(tidybayes)
library(car)
library(corrplot)

rm(list = ls())

setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/spatial_sensitivity_analysis/Larger_radius/RDATA")
load('plankton_RUE_data_scaled_260526.RDATA')

setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/spatial_sensitivity_analysis/Larger_radius/RDATA")

################# Getting results for all models ################################

# ===========================================================
# for phytoplankton div models
# Automatic output folder
date_tag <- format(Sys.Date(), "%Y%m%d")
outdir <- file.path("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/spatial_sensitivity_analysis/Larger_radius/Results",
                    paste0("BRM_phytodiv_results_", date_tag))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat("\n Results will be saved to:", outdir, "\n")

# Load model list
source("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/spatial_sensitivity_analysis/Larger_radius/2.1.phytodiv_model_list.R")

# Check it loaded correctly
length(model_list)
names(model_list[[1]])
# ===========================================================

# ===========================================================
# for phytoplankton RUE models
# Automatic output folder
date_tag <- format(Sys.Date(), "%Y%m%d")
outdir <- file.path("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/spatial_sensitivity_analysis/Larger_radius/Results",
                    paste0("BRM_phytoRUE_results_", date_tag))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat("\n Results will be saved to:", outdir, "\n")

# Load model list
source("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/spatial_sensitivity_analysis/Larger_radius/2.2.phytoRUE_model_list.R")

# Check it loaded correctly
length(model_list)
names(model_list[[1]])
# ===========================================================

# ===========================================================
# for zooplankton models
# Automatic output folder
date_tag <- format(Sys.Date(), "%Y%m%d")
outdir <- file.path("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/spatial_sensitivity_analysis/Larger_radius/Results",
                    paste0("BRM_zoo_results_", date_tag))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat("\n Results will be saved to:", outdir, "\n")

# Load model list
source("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/spatial_sensitivity_analysis/Larger_radius/2.3.zoo_model_list.R")

# Check it loaded correctly
length(model_list)
names(model_list[[1]])
# ===========================================================

# ===========================================================
# for cyanobacteria dominance models
# Automatic output folder
date_tag <- format(Sys.Date(), "%Y%m%d")
outdir <- file.path("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/spatial_sensitivity_analysis/Larger_radius/Results",
                    paste0("BRM_cyanodom_results_", date_tag))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat("\n Results will be saved to:", outdir, "\n")

# Load model list
source("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/spatial_sensitivity_analysis/Larger_radius/2.4.cyanodom_model_list.R")

# Check it loaded correctly
length(model_list)
names(model_list[[1]])
# ===========================================================

# =============================================================================
# Function: Check beta trend from Emy's code
# =============================================================================
check_beta_trend <-function(LpredDeriv){
  
  pn<-dim(LpredDeriv)[2]
  response = matrix(NA,nrow = nrow(LpredDeriv),ncol = 1)
  
  for (i in 1:nrow(LpredDeriv)){
    if (LpredDeriv[i,1]>0 & LpredDeriv[i,pn]>0){   
      # derivative is positive at both ends -> clear case 1
      response[i] = 1 
    } else if (LpredDeriv[i,1]<0 & LpredDeriv[i,pn]<0) {
      # derivative is negative at both ends -> clear case 2
      response[i] = 2
    } else if (LpredDeriv[i,1]>0 & LpredDeriv[i,pn]<0) {
      # derivative changes from postitive to negative -> have to check in more detail
      
      if (sum(LpredDeriv[i,]>0)/pn<0.2){
        # derivative is positive over less than 20% of the x-gradient extent -> classify as case 2
        response[i] = 2
        
      } else if (sum(LpredDeriv[i,]>0)/pn>0.8){
        # derivative is positive over more than 80% of the x-gradient extent -> classify as case 1
        response[i] = 1
        
      } else {
        # derivative is postive or negative only over at most 60% of the x-gradient extent -> classify as case 3
        response[i] = 3
      }
      
    } else if (LpredDeriv[i,1]<0 & LpredDeriv[i,pn]>0) {
      # derivative changes from negative to positive -> have to check in more detail
      
      if (sum(LpredDeriv[i,]<0)/pn<0.2){
        # derivative is negative over less than 20% of the x-gradient extent -> classify as case 1
        response[i] = 1
        
      } else if (sum(LpredDeriv[i,]<0)/pn>0.8){
        # derivative is negative over more than 80% of the x-gradient extent -> classify as case 2
        response[i] = 2
        
      } else {response[i] = 4 
      # derivative is positive or negative only over at most 60% of the x-gradient extent 
      
      }
    }
    
  }
  # calculate the probabilities for the cases
  Pr1 = sum(response==1,na.rm=T)/length(response) #positive; derivative is positive over more than 80%
  Pr2 = sum(response==2,na.rm=T)/length(response) #negative, derivative is negative over more than 80%
  Pr3 = sum(response==3,na.rm=T)/length(response) #change; derivative is positive or negative only over at most 60% 
  Pr4 = sum(response==4,na.rm=T)/length(response) #change; Ushape 
  
  return(c(Pr1,Pr2,Pr3,Pr4))}

# =============================================================================
# Function: compute derivative + classify
# =============================================================================
classify_predictor <- function(model, data, var, n = 200) {
  
  # Create a grid over the predictor range
  xgrid <- seq(min(data[[var]], na.rm = TRUE),
               max(data[[var]], na.rm = TRUE),
               length.out = n)
  
  # Extract posterior samples
  post <- as_draws_df(model)
  
  # Identify linear & quadratic coefficient names
  coef_names <- colnames(post)
  
  # linear term
  b1_name <- coef_names[
    coef_names == paste0("b_", var) |
      grepl(paste0("^b_poly", var, ".*1$"), coef_names)
  ]
  
  # quadratic term
  b2_name <- coef_names[
    grepl(paste0("^b_poly", var, ".*2$"), coef_names) |
      grepl(paste0("^b_I.*", var, ".*2"), coef_names)
  ]
  
  if (length(b1_name) != 1 | length(b2_name) != 1) {
    stop(paste("Could not find coefficients for", var))
  }
  
  b1 <- post[[b1_name]]
  b2 <- post[[b2_name]]
  
  # Derivative for quadratic poly(raw=TRUE)
  # d/dx = b1 + 2*b2*x
  LpredDeriv <- outer(seq_along(b1), xgrid,
                      function(i, x) b1[i] + 2 * b2[i] * x)
  
  # classify shapes
  res <- check_beta_trend(LpredDeriv)
  names(res) <- c("Pr_positive", "Pr_negative", "Pr_hump", "Pr_Ushape")
  
  return(res)
}

# ===========================================================
# Function: Full model analysis
# ===========================================================
analyze_model <- function(model_file, name, resp = NULL, predictors = NULL) {
  
  cat("\n========================================\n")
  cat("Running model:", name, "\n")
  cat("File:", model_file, "\n")
  cat("========================================\n")
  
  if (!file.exists(model_file)) {
    message("Model file not found: ", model_file)
    return(invisible(NULL))
  }
  
  # --- Load model ---
  e <- new.env()
  load(model_file, envir = e)
  
  if (!"model" %in% ls(e)) {
    stop("No object named 'model' found in ", model_file)
  }
  
  m <- e$model
  
  # --- Summary sink ---
  summary_path <- file.path(outdir, paste0("summary_", name, ".txt"))
  sink(summary_path)
  
  cat("========================================\n")
  cat(" Model Summary -", name, "\n")
  cat("========================================\n\n")
  print(summary(m))
  sink()
  
  # --- Posterior predictive check ---
  ppc_types <- c(
    "dens_overlay",
    "ecdf_overlay",
    "intervals",
    "scatter_avg",
    "loo_pit_overlay"
  )
  
  for (ppc_type in ppc_types) {
    
    cat("Generating PPC:", ppc_type, "for model:", name, "\n")
    
    p_ppc <- tryCatch(
      {
        pp_check(m, type = ppc_type) +
          ggtitle(paste("PPC:", ppc_type, "-", name))
      },
      error = function(e) {
        message("Skipping PPC type '", ppc_type, "': ", e$message)
        return(NULL)
      }
    )
    
    if (!is.null(p_ppc)) {
      ggsave(
        file.path(outdir, paste0("ppcheck_", ppc_type, "_", name, ".png")),
        plot = p_ppc,
        width = 6,
        height = 4,
        dpi = 300
      )
    }
  }
  
  # --- LOO ---
  cat("Calculating LOO and ELPD for model:", name, "\n")
  loo_res <- loo(m, moment_match = TRUE)
  
  # if (any(loo_res$diagnostics$pareto_k > 0.7, na.rm = TRUE)) {
  #   message("High Pareto-k detected → running reloo() for ", name)
  #   loo_res <- loo(
  #     m,
  #     reloo = TRUE,
  #     cores = 1
  #   )
  # }
  
  sink(summary_path, append = TRUE)
  cat("\n\n========================================\n")
  cat(" LOO Cross-Validation Results\n")
  cat("========================================\n\n")
  print(loo_res)
  sink()
  
  # --- Derivative classification (poly terms only) ---
  if (!is.null(predictors) && length(predictors) > 0) {
    coef_names <- rownames(fixef(m))
    post <- as_draws_df(m)
    poly_predictors <- predictors[sapply(predictors, function(var) {
      any(
        grepl(paste0("poly", var), colnames(post)) |      # poly(x,2)
          grepl(paste0("^b_", var, "$"), colnames(post)) |  # linear term
          grepl(paste0("I.*", var, ".*2"), colnames(post))  # I(x^2)
      )
    })]
    
    if (length(poly_predictors) == 0) {
      cat("\nNo poly() predictors found → skipping derivative trends.\n")
      sink(summary_path, append = TRUE)
      cat("\nDerivative-based Predictor Classification\n")
      cat("No poly() predictors. No derivative trends computed.\n")
      sink()
    } else {
      cat("\nCalculating derivative-based trends for poly() predictors...\n")
      
      derivative_results <- lapply(poly_predictors, function(var) {
        tryCatch(
          classify_predictor(model = m, data = m$data, var = var),
          error = function(e) {
            message("Skipping derivative classification for ", var, ": ", e$message)
            return(c(NA, NA, NA, NA))
          }
        )
      })
      
      derivative_df <- as.data.frame(do.call(rbind, derivative_results))
      derivative_df$predictor <- poly_predictors
      derivative_df <- derivative_df[, c("predictor", "Pr_positive",
                                         "Pr_negative", "Pr_hump", "Pr_Ushape")]
      
      # Remove predictors with all NA to avoid plotting issues
      derivative_df <- derivative_df[!apply(is.na(derivative_df[,2:5]), 1, all), ]
      
      # Append to summary
      sink(summary_path, append = TRUE)
      cat("\nDerivative-based Predictor Classification\n")
      if (nrow(derivative_df) > 0) {
        cat("Predictor | Pr_positive | Pr_negative | Pr_hump | Pr_Ushape\n")
        cat("-----------------------------------------------------------\n")
        for (i in 1:nrow(derivative_df)) {
          cat(sprintf("%-15s | %10.3f | %10.3f | %7.3f | %9.3f\n",
                      derivative_df$predictor[i],
                      derivative_df$Pr_positive[i],
                      derivative_df$Pr_negative[i],
                      derivative_df$Pr_hump[i],
                      derivative_df$Pr_Ushape[i]))
        }
      } else {
        cat("No valid derivative classifications were computed.\n")
      }
      sink()
      
      print(derivative_df)
    }
  }
  
  # --- MCMC plot ---
  p_mcmc <- mcmc_plot(m, prob = 0.8, prob_outer = 0.95) +
    ggtitle(paste("MCMC Coefficient Plot -", name))
  ggsave(file.path(outdir, paste0("mcmc_coef_", name, ".png")),
         plot = p_mcmc, width = 6, height = 4, dpi = 300)
  
  # R-hat
  rhats <- brms::rhat(m)
  fixed_rhats <- rhats[grep("^b_", names(rhats))]
  
  if (length(fixed_rhats) > 0) {
    p_rhat <- mcmc_rhat(fixed_rhats) +
      yaxis_text(hjust = 1) +
      ggtitle(paste("MCMC R-hat (Fixed Effects) -", name))
    
    ggsave(
      file.path(outdir, paste0("mcmc_rhat_fixed_", name, ".png")),
      plot = p_rhat,
      width = 6,
      height = 4,
      dpi = 300
    )
  } else {
    message("No fixed effects R-hat values found for ", name)
  }
  
  # Effective sample size
  neff_ratios <- bayesplot::neff_ratio(m)
  
  fixed_neff <- neff_ratios[grep("^b_", names(neff_ratios))]
  
  if (length(fixed_neff) > 0) {
    p_neff <- mcmc_neff(fixed_neff) +
      yaxis_text(hjust = 1) +
      ggtitle(paste("MCMC Effective Sample Size (Fixed Effects) -", name))
    
    ggsave(
      file.path(outdir, paste0("mcmc_neff_fixed_", name, ".png")),
      plot = p_neff,
      width = 6,
      height = 4,
      dpi = 300
    )
  } else {
    message("No fixed-effect Neff values found for ", name)
  }
  
  # Autocorrelation (lags = 9)
  p_acf <- mcmc_plot(m, type = "acf", lags = 9) +
    ggtitle(paste("MCMC Autocorrelation (lags = 9) -", name))
  
  ggsave(
    file.path(outdir, paste0("mcmc_acf_", name, ".png")),
    plot = p_acf,
    width = 6,
    height = 4,
    dpi = 300
  )
  
  
  # --- Define predictor color groups ---
  abiotic_vars <- c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered", "MaxDepth", 
                    "Zooplankton_MaxDepth", "DOY", "daylight_hours")
  phyto_vars <- c("genus_richness", "phyto_evenness", "phyto_FEve_PD", "phyto_FDis_PD", "phyto_FEve_ZD", "phyto_FDis_ZD", "phyto_FEve_PRUE", 
                  "phyto_FDis_PRUE", "phyto_FEve_ZRUE", "phyto_FDis_ZRUE", "cyano_prop",
                  "genus_richness:Temperature", "genus_richness:Salinity", "genus_richness:TotalNitrogen_unfiltered",
                  "phyto_evenness:Temperature", "phyto_evenness:Salinity", "phyto_evenness:TotalNitrogen_unfiltered",
                  "phyto_FEve_PRUE:Temperature", "phyto_FEve_PRUE:Salinity",  "phyto_FEve_PRUE:TotalNitrogen_unfiltered", 
                  "phyto_FDis_PRUE:Temperature", "phyto_FDis_PRUE:Salinity",  "phyto_FDis_PRUE:TotalNitrogen_unfiltered",
                  "phyto_FEve_ZRUE:Temperature", "phyto_FEve_ZRUE:Salinity",  "phyto_FEve_ZRUE:TotalNitrogen_unfiltered", 
                  "phyto_FDis_ZRUE:Temperature", "phyto_FDis_ZRUE:Salinity",  "phyto_FDis_ZRUE:TotalNitrogen_unfiltered",
                  "cyano_prop:Temperature", "cyano_prop:Salinity", "cyano_prop:TotalNitrogen_unfiltered")
  zoo_vars <- c("species_richness", "zoo_evenness", "zoo_FEve", "zoo_FDis",
                "species_richness:Temperature", "species_richness:Salinity", "species_richness:TotalNitrogen_unfiltered",
                "zoo_evenness:Temperature", "zoo_evenness:Salinity", "zoo_evenness:TotalNitrogen_unfiltered",
                "zoo_FEve:Temperature", "zoo_FEve:Salinity", "zoo_FEve:TotalNitrogen_unfiltered",
                "zoo_FDis:Temperature", "zoo_FDis:Salinity", "zoo_FDis:TotalNitrogen_unfiltered")
  
  group_colors <- c(
    "Abiotic" = "#1f78b4",               
    "Phytoplankton diversity" = "#33a02c", 
    "Zooplankton diversity" = "darkorchid1"    
  )
  
  get_group_color <- function(var) {
    if (var %in% abiotic_vars) return(group_colors["Abiotic"])
    if (var %in% phyto_vars) return(group_colors["Phytoplankton diversity"])
    if (var %in% zoo_vars) return(group_colors["Zooplankton diversity"])
    return("grey50")
  }
  
  
  cat("Completed:", name, "\n\n")
  invisible(NULL)
}

# ===========================================================
# Run all models
# ===========================================================
for (mdl in model_list) {
  tryCatch(
    analyze_model(
      model_file = mdl$file,
      name = mdl$name,
      resp = mdl$resp,
      predictors = mdl$predictors
    ),
    error = function(e) message("Error in ", mdl$name, ": ", e$message)
  )
}

cat("\n All model analyses completed.\nResults saved in:", outdir, "\n")


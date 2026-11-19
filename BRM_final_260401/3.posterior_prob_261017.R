
################################################################################
############  Calculating posterior probabilities  ############################
################################################################################

library(dplyr)
library(brms)
library(posterior)
library(writexl)

rm(list = ls())

# Directories
model_dir <- path.expand(
  "~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA"
)

results_dir <- path.expand(
  "~/GitHub/PlanktonRUE-Project/Results"
)


# ------------------------------------------------------------------------------
# Function to load a specific RDATA file and calculate probabilities
# ------------------------------------------------------------------------------

calculate_posterior_probabilities <- function(file) {
  
  model_name <- tools::file_path_sans_ext(basename(file))
  
  message("Processing: ", model_name)
  
  # Load model into an isolated environment
  env <- new.env()
  load(file, envir = env)
  
  # Find brms model objects
  objects <- ls(env)
  
  model_objects <- objects[
    vapply(
      objects,
      function(nm) inherits(env[[nm]], "brmsfit"),
      logical(1)
    )
  ]
  
  if (length(model_objects) == 0) {
    stop("No fitted brms model found in: ", file)
  }
  
  if (length(model_objects) > 1) {
    stop(
      "Multiple brms models found in ", file,
      ". Please specify which model object to use."
    )
  }
  
  model <- env[[model_objects[1]]]
  
  # Extract posterior draws
  draws <- posterior::as_draws_df(model)
  
  # Fixed-effect coefficients
  coef_names <- grep("^b_", names(draws), value = TRUE)
  
  if (length(coef_names) == 0) {
    stop("No fixed-effect coefficients found in: ", file)
  }
  
  # Calculate probabilities and credible intervals
  results <- lapply(coef_names, function(coef) {
    
    x <- draws[[coef]]
    
    prob_positive <- mean(x > 0)
    prob_negative <- mean(x < 0)
    
    data.frame(
      Model = model_name,
      Parameter = sub("^b_", "", coef),
      Mean = mean(x),
      Median = median(x),
      CI_lower = unname(quantile(x, 0.025)),
      CI_upper = unname(quantile(x, 0.975)),
      Prob_positive = prob_positive,
      Prob_negative = prob_negative,
      Prob_direction = max(prob_positive, prob_negative)
    )
  })
  
  bind_rows(results) %>%
    mutate(
      Probability_label = case_when(
        Prob_direction >= 0.95 ~ ">95%",
        Prob_direction >= 0.90 ~ ">90%",
        Prob_direction >= 0.85 ~ ">85%",
        TRUE ~ "<85%"
      ),
      CI_excludes_zero = CI_lower > 0 | CI_upper < 0
    )
}


# ------------------------------------------------------------------------------
# Specify ONLY the models you want to analyse
# ------------------------------------------------------------------------------

model_files <- c(
  #phyto div models
  "BRM_phyto_GD_P_260916.RDATA",
  "BRM_phyto_SE_P_260916.RDATA",
  "BRM_phyto_FEve_P_260916.RDATA",
  "BRM_phyto_FDis_P_260916.RDATA",
  "BRM_cyanodom_P_261009_z.RDATA",
  #phyto RUE models
  "BRM_phytoRUEP_CD_PI_261013_All.RDATA",
  "BRM_phytoRUEN_CD_PI_261013_All.RDATA",
  #zoo div models
  "BRM_zoo_SR_CD_P_261013_All.RDATA",
  "BRM_zoo_SE_CD_P_261013_All.RDATA",
  "BRM_zoo_FE_CD_P_261013_All.RDATA",
  "BRM_zoo_FD_CD_P_261013_All.RDATA",
  #zoo RUE models
  "BRM_zoo_RUE_CD_PI_261013_All.RDATA"
  # Add or remove filenames as needed
)

# Build full paths
model_paths <- file.path(model_dir, model_files)

# Check that all specified files exist
missing_files <- model_paths[!file.exists(model_paths)]

if (length(missing_files) > 0) {
  stop(
    "The following model files were not found:\n",
    paste(missing_files, collapse = "\n")
  )
}


# ------------------------------------------------------------------------------
# Run analysis on the specified models only
# ------------------------------------------------------------------------------

all_results <- bind_rows(
  lapply(model_paths, calculate_posterior_probabilities)
)


# ------------------------------------------------------------------------------
# View and export results
# ------------------------------------------------------------------------------

print(as.data.frame(all_results), row.names = FALSE)

# Export as Excel
output_file <- file.path(
  results_dir,
  "All_models_posterior_probabilities.xlsx"
)

write_xlsx(
  all_results,
  path = output_file
)

message("Excel file saved to: ", output_file)

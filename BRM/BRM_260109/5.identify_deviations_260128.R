################################################################################
############### Identify points that deviate from model ########################
################################################################################
library(brms)
library(loo)

# ---- Load the model ----
model_file <- "~/GitHub/Phytoplankton-project/phyto_RUE/BRM_260109/BRM_cyanodom_L_260109_Z.RDATA"
load(model_file)  # assumes the object inside is named 'model'
m <- model

# ---- Specify response variable ----
resp <- "log_RUE_A_P"  # adjust to your model's response

# ---- Compute LOO ----
loo_res <- loo(m)

# ---- Only use rows actually used in the model ----
used_rows <- !is.na(m$data[[resp]])
data_used <- m$data[used_rows, ]

# ---- Posterior predictions and fitted values ----
yrep <- posterior_predict(m)          # draws x observations
mu_hat <- fitted(m, summary = FALSE)
mu_mean <- colMeans(mu_hat)
mu_sd   <- apply(mu_hat, 2, sd)

# ---- Standardised residuals ----
resid_std <- (data_used[[resp]] - mu_mean) / mu_sd

# ---- Posterior predictive tail probability ----
pp_tail <- apply(yrep, 2, function(x) mean(x < data_used[[resp]]))

# ---- Pareto-k ----
pareto_k <- loo_res$diagnostics$pareto_k
if(length(pareto_k) != nrow(data_used)){
  warning("Pareto-k length mismatch: filling with NA")
  pareto_k <- rep(NA, nrow(data_used))
}

# ---- Build deviation table ----
deviation_table <- data.frame(
  Site_Id           = data_used$Site_Id,      # site ID
  Latitude          = data_used$Latitude,     # latitude
  Longitude         = data_used$Longitude,    # longitude
  Date              = data_used$Date,
  observed_response = data_used[[resp]],
  resid_std         = resid_std,
  pp_tail           = pp_tail,
  pareto_k          = pareto_k,
  deviates          = abs(resid_std) > 3 |
    pp_tail < 0.025 | pp_tail > 0.975 |
    (!is.na(pareto_k) & pareto_k > 0.7)
)

# ---- Optional: keep only deviating points ----
deviation_table <- deviation_table[deviation_table$deviates, ]

# ---- Save table ----
write.csv(deviation_table, "deviating_points_cyanodom_L_z.csv", row.names = FALSE)

# ---- View table ----
head(deviation_table)
#################################################################################
####################### Conditional effect plots ###############################
################################################################################

library(dplyr)
library(tidyr)
library(brms)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(ragg)

rm(list = ls())

setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA")

#plotting ce function ----------------------------------------------
plot_ce_simple <- function(var, model, resp = NULL,
                           xlab = var,
                           ylab = resp,
                           title = var,
                           legend_title = NULL) {
  
  ce <- conditional_effects(model, effects = var, resp = resp)
  p <- plot(ce, plot = FALSE)[[1]]
  
  is_interaction <- grepl(":", var)
  
  rawdat <- model$data
  
  if (!is.null(resp)) {
    rawdat <- rawdat[!is.na(rawdat[[resp]]), ]
  }
  
  set.seed(123)
  rawdat <- rawdat[sample(nrow(rawdat),
                          min(5000, nrow(rawdat))), ]
  
  group_color <- get_group_color(var)
  
  # Add raw data points for main effects only
  if (!is_interaction && var %in% names(rawdat)) {
    p <- p +
      geom_point(
        data = rawdat,
        aes(x = .data[[var]], y = .data[[resp]]),
        inherit.aes = FALSE,
        color = "black",
        alpha = 0.3
      )
  }
  
  # Legend title
  if (!is.null(legend_title)) {
    p <- p +
      labs(colour = legend_title,
           fill = legend_title)
  }
  
  # Labels and theme
  p <- p +
    labs(
      title = gsub(":", " × ", title),
      x = xlab,
      y = ylab
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(color = group_color,
                                face = "bold"),
      axis.title = element_text(color = group_color)
    )
  
  return(p)
}

abiotic_vars <- c("Temperature", "Salinity", "pH", "TotalNitrogen_unfiltered", "MaxDepth", 
                  "Zooplankton_MaxDepth", "DOY", "daylight_hours")
phyto_vars <- c("genus_richness", "FEve", "FDis", "cyano_prop",
                "genus_richness:Temperature", "genus_richness:Salinity", "genus_richness:TotalNitrogen_unfiltered",
                "FEve:Temperature", "FEve:Salinity",  "FEve:TotalNitrogen_unfiltered", 
                "FDis:Temperature", "FDis:Salinity",  "FDis:TotalNitrogen_unfiltered",
                "cyano_prop:Temperature", "cyano_prop:Salinity", "cyano_prop:TotalNitrogen_unfiltered")
zoo_vars <- c("species_richness", "zoo_FEve", "zoo_FDis",
              "species_richness:Temperature", "species_richness:Salinity", "species_richness:TotalNitrogen_unfiltered",
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

#Phyto RUE per P #--------------------------------------------------
setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA")
load("BRM_phytoRUEP_CD_PI_260408.RDATA")

p1 <- plot_ce_simple("Temperature", model, "log_RUE_B_P", xlab = "Temperature", 
                     ylab = "RUE per P", title = "Temperature")
p2 <- plot_ce_simple("Salinity", model, "log_RUE_B_P", xlab = "Salinity", ylab = "RUE per P", title = "Salinity")
p3 <- plot_ce_simple("pH", model, "log_RUE_B_P",xlab = "pH", ylab = "RUE per P", title = "pH")
p4 <- plot_ce_simple("TotalNitrogen_unfiltered", model, "log_RUE_B_P",
                     xlab = "Nitrogen concentrations", ylab = "RUE per P", title = "Nitrogen concentrations")
p5 <- plot_ce_simple("genus_richness", model, "log_RUE_B_P", 
                     xlab = "Phytoplankton genus richness", ylab = "RUE per P", 
                     title = "Phytoplankton genus richness")
p6 <- plot_ce_simple("FEve", model, "log_RUE_B_P", xlab = "Phytoplankton FEve", 
                     ylab = "RUE per P", title = "Phytoplankton FEve")
p7 <- plot_ce_simple("FDis", model, "log_RUE_B_P", xlab = "Phytoplankton FDis", 
                     ylab = "RUE per P", title = "Phytoplankton FDis")
p8 <- plot_ce_simple("cyano_prop", model, "log_RUE_B_P", xlab = "Cyanobacteria Dominance", 
                     ylab = "RUE per P", title = "Cyanobacteria Dominance")
p9 <- plot_ce_simple("species_richness", model, "log_RUE_B_P", xlab = "Zooplankton species richness", 
                     ylab = "RUE per P", title = "Zooplankton species richness")
p10 <- plot_ce_simple("zoo_FEve", model, "log_RUE_B_P", xlab = "Zooplankton FEve", 
                      ylab = "RUE per P", title = "Zooplankton FEve")
p11 <- plot_ce_simple("zoo_FDis", model, "log_RUE_B_P", xlab = "Zooplankton FDis", 
                      ylab = "RUE per P", title = "Zooplankton FDis")
p12 <- plot_ce_simple("FEve:TotalNitrogen_unfiltered", model, "log_RUE_B_P", xlab = "Phytoplankton FEve", 
                      ylab = "RUE per P", title = "Phytoplankton FEve:Nitrogen concentrations", 
                      legend_title = "Nitrogen concentrations")
p13 <- plot_ce_simple("FDis:Salinity", model, "log_RUE_B_P", xlab = "Phytoplankton FDis", 
                      ylab = "RUE per P", title = "Phytoplankton FDis:Salinity", 
                      legend_title = "Salinity")
p14 <- plot_ce_simple("cyano_prop:Salinity", model, "log_RUE_B_P", xlab = "Cyanobacteria Dominance", 
                      ylab = "RUE per P", title = "Cyanobacteria Dominance:Salinity", 
                      legend_title = "Salinity")
p15 <- plot_ce_simple("cyano_prop:TotalNitrogen_unfiltered", model, "log_RUE_B_P", xlab = "Cyanobacteria Dominance", 
                      ylab = "RUE per P", title = "Cyanobacteria Dominance:Nitrogen concentrations", 
                      legend_title = "Nitrogen concentrations")
p16 <- plot_ce_simple("zoo_FDis:TotalNitrogen_unfiltered", model, "log_RUE_B_P", xlab = "Zooplankton FDis", 
                      ylab = "RUE per P", title = "Zooplankton FDis:Nitrogen concentrations", 
                      legend_title = "Nitrogen concentrations")
p17 <- plot_ce_simple("MaxDepth", model, "log_RUE_B_P", xlab = "Depth", 
                      ylab = "RUE per P", title = "Depth")
p18 <- plot_ce_simple("DOY", model, "log_RUE_B_P", xlab = "DOY", 
                      ylab = "RUE per P", title = "DOY")
p19 <- plot_ce_simple("daylight_hours", model, "log_RUE_B_P", xlab = "Daylight Hours", 
                      ylab = "RUE per P", title = "Daylight Hours")

combined_plot <- wrap_plots(
  p1, p2, p3, p4, p5,
  p6, p7, p8, p9, p10,
  p11,
  
  p12 + plot_layout(widths = 2),
  p13 + plot_layout(widths = 2),
  p14 + plot_layout(widths = 2),
  p15 + plot_layout(widths = 2),
  p16 + plot_layout(widths = 2),
  
  p17, p18, p19,
  
  ncol = 5
)

combined_plot

setwd("~/GitHub/PlanktonRUE-Project/Figures/CE_plots")
ggsave("CE_RUEP.png", 
       plot = combined_plot, width = 25, height = 11, dpi = 300, device = ragg::agg_png)

#Phyto RUE per N #--------------------------------------------------
setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA")
load("BRM_phytoRUEN_CD_PI_260408.RDATA")

p1 <- plot_ce_simple("Temperature", model, "log_RUE_B_N", xlab = "Temperature", 
                     ylab = "RUE per N", title = "Temperature")
p2 <- plot_ce_simple("Salinity", model, "log_RUE_B_N", xlab = "Salinity", 
                     ylab = "RUE per N", title = "Salinity")
p3 <- plot_ce_simple("pH", model, "log_RUE_B_N", xlab = "pH", 
                     ylab = "RUE per N", title = "pH")
p4 <- plot_ce_simple("TotalNitrogen_unfiltered", model, "log_RUE_B_N", xlab = "Nitrogen concentrations", 
                     ylab = "RUE per N", title = "Nitrogen concentrations")
p5 <- plot_ce_simple("genus_richness", model, "log_RUE_B_N", xlab = "Phytoplankton genus richness", 
                     ylab = "RUE per N", title = "Phytoplankton genus richness")
p6 <- plot_ce_simple("FEve", model, "log_RUE_B_N", xlab = "Phytoplankton FEve", 
                     ylab = "RUE per N", title = "Phytoplankton FEve")
p7 <- plot_ce_simple("FDis", model, "log_RUE_B_N", xlab = "Phytoplankton FDis", 
                     ylab = "RUE per N", title = "Phytoplankton FDis")
p8 <- plot_ce_simple("cyano_prop", model, "log_RUE_B_N", xlab = "Cyanobacteria Dominance", 
                     ylab = "RUE per N", title = "Cyanobacteria Dominance")
p9 <- plot_ce_simple("species_richness", model, "log_RUE_B_N", xlab = "Zooplankton species richness", 
                     ylab = "RUE per N", title = "Zooplankton species richness")
p10 <- plot_ce_simple("zoo_FEve", model, "log_RUE_B_N", xlab = "Zooplankton FEve", 
                      ylab = "RUE per N", title = "Zooplankton FEve")
p11 <- plot_ce_simple("zoo_FDis", model, "log_RUE_B_N", xlab = "Zooplankton FDis", 
                      ylab = "RUE per N", title = "Zooplankton FDis")
p12 <- plot_ce_simple("FEve:TotalNitrogen_unfiltered", model, "log_RUE_B_N", xlab = "Phytoplankton FEve", 
                      ylab = "RUE per N", title = "phytoplankton FEve:Nitrogen concentrations",
                      legend_title = "Nitrogen concentrations")
p13 <- plot_ce_simple("FDis:Salinity", model, "log_RUE_B_N", xlab = "Phytoplankton FDis", 
                      ylab = "RUE per N", title = "phytoplankton FDis:Salinity",
                      legend_title = "Salinity")
p14 <- plot_ce_simple("FDis:TotalNitrogen_unfiltered", model, "log_RUE_B_N", xlab = "Phytoplankton FDis", 
                      ylab = "RUE per N", title = "phytoplankton FDis:Nitrogen concentrations",
                      legend_title = "Nitrogen concentrations")
p15 <- plot_ce_simple("cyano_prop:Salinity", model, "log_RUE_B_N", xlab = "Cyanobacteria Dominance", 
                      ylab = "RUE per N", title = "Cyanobacteria Dominance:Salinity",
                      legend_title = "Salinity")
p16 <- plot_ce_simple("cyano_prop:TotalNitrogen_unfiltered", model, "log_RUE_B_N", xlab = "Cyanobacteria Dominance", 
                      ylab = "RUE per N", title = "Cyanobacteria Dominance:Nitrogen concentrations",
                      legend_title = "Nitrogen concentrations")
p17 <- plot_ce_simple("zoo_FDis:TotalNitrogen_unfiltered", model, "log_RUE_B_N", xlab = "Zooplankton FDis", 
                      ylab = "RUE per N", title = "Zooplankton FDis:Nitrogen concentrations",
                      legend_title = "Nitrogen concentrations")
p18 <- plot_ce_simple("MaxDepth", model, "log_RUE_B_N", xlab = "Depth", 
                      ylab = "RUE per N", title = "Depth")
p19 <- plot_ce_simple("DOY", model, "log_RUE_B_N", xlab = "DOY", 
                      ylab = "RUE per N", title = "DOY")
p20 <- plot_ce_simple("daylight_hours", model, "log_RUE_B_N", xlab = "Daylight Hours", 
                      ylab = "RUE per N", title = "Daylight Hours")

combined_plot <- wrap_plots(
  p1, p2, p3, p4, p5,
  p6, p7, p8, p9, p10,
  p11,
  
  p12 + plot_layout(widths = 2),
  p13 + plot_layout(widths = 2),
  p14 + plot_layout(widths = 2),
  p15 + plot_layout(widths = 2),
  p16 + plot_layout(widths = 2),
  p17 + plot_layout(widths = 2),
  
  p18, p19, p20,
  
  ncol = 5
)

combined_plot

setwd("~/GitHub/PlanktonRUE-Project/Figures/CE_plots")
ggsave("CE_RUEN.png", 
       plot = combined_plot, width = 25, height = 11, dpi = 300, device = ragg::agg_png)

#Zoo RUE #--------------------------------------------------
setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA")
load("BRM_zoo_RUE_CD_PI_260427.RDATA")

p1 <- plot_ce_simple("Temperature", model, "log_RUE_zooplankton", xlab = "Temperature", 
                     ylab = "Zooplankton RUE", title = "Temperature")
p2 <- plot_ce_simple("Salinity", model, "log_RUE_zooplankton", xlab = "Salinity", 
                     ylab = "Zooplankton RUE", title = "Salinity")
p3 <- plot_ce_simple("pH", model, "log_RUE_zooplankton", xlab = "pH", 
                     ylab = "Zooplankton RUE", title = "pH")
p4 <- plot_ce_simple("TotalNitrogen_unfiltered", model, "log_RUE_zooplankton", xlab = "Nitrogen concentrations", 
                     ylab = "Zooplankton RUE", title = "Nitrogen concentrations")
p5 <- plot_ce_simple("genus_richness", model, "log_RUE_zooplankton", xlab = "Phytoplankton genus richness", 
                     ylab = "Zooplankton RUE", title = "Phytoplankton genu richness")
p6 <- plot_ce_simple("FEve", model, "log_RUE_zooplankton", xlab = "Phytoplankton FEve", 
                     ylab = "Zooplankton RUE", title = "Phytoplankton FEve")
p7 <- plot_ce_simple("FDis", model, "log_RUE_zooplankton", xlab = "Phytoplankton FDis", 
                     ylab = "Zooplankton RUE", title = "Phytoplankton FDis")
p8 <- plot_ce_simple("cyano_prop", model, "log_RUE_zooplankton", xlab = "Cyanobacteria Dominance", 
                     ylab = "Zooplankton RUE", title = "Cyanobacteria Dominance")
p9 <- plot_ce_simple("species_richness", model, "log_RUE_zooplankton", xlab = "Zooplankton species richness", 
                     ylab = "Zooplankton RUE", title = "Zooplankton species richness")
p10 <- plot_ce_simple("zoo_FEve", model, "log_RUE_zooplankton", xlab = "Zooplankton FEve", 
                      ylab = "Zooplankton RUE", title = "Zooplankton Feve")
p11 <- plot_ce_simple("zoo_FDis", model, "log_RUE_zooplankton", xlab = "Zooplankton FDis", 
                      ylab = "Zooplankton RUE", title = "Zooplankton FDis")
p12 <- plot_ce_simple("genus_richness:Temperature", model, "log_RUE_zooplankton", xlab = "Phytoplankton genus richness", 
                      ylab = "Zooplankton RUE", title = "Phytoplankton genus richness:Temperature",
                      legend_title = "Temperature")
p13 <- plot_ce_simple("genus_richness:Salinity", model, "log_RUE_zooplankton", xlab = "Phytoplankton genus richness", 
                      ylab = "Zooplankton RUE", title = "Phytoplankton genus richness:Salinity",
                      legend_title = "Salinity")
p14 <- plot_ce_simple("FDis:Salinity", model, "log_RUE_zooplankton", xlab = "Phytoplankton FDis", 
                      ylab = "Zooplankton RUE", title = "Phytoplankton FDis:Salinity",
                      legend_title = "Salinity")
p15 <- plot_ce_simple("cyano_prop:Temperature", model, "log_RUE_zooplankton", xlab = "Cyanobacteria Dominance", 
                      ylab = "Zooplankton RUE", title = "Cyanobacteria Dominance:Temperature",
                      legend_title = "Temperature")
p16 <- plot_ce_simple("Zooplankton_MaxDepth", model, "log_RUE_zooplankton", xlab = "Depth", 
                      ylab = "Zooplankton RUE", title = "Depth")
p17 <- plot_ce_simple("DOY", model, "log_RUE_zooplankton", xlab = "DOY", 
                      ylab = "Zooplankton RUE", title = "DOY")
p18 <- plot_ce_simple("daylight_hours", model, "log_RUE_zooplankton", xlab = "Daylight Hours", 
                      ylab = "Zooplankton RUE", title = "Daylight Hours")

combined_plot <- wrap_plots(
  p1, p2, p3, p4, p5,
  p6, p7, p8, p9, p10,
  p11,
  
  p12 + plot_layout(widths = 2),
  p13 + plot_layout(widths = 2),
  p14 + plot_layout(widths = 2),
  p15 + plot_layout(widths = 2),
  
  p16, p17, p18,
  
  ncol = 5
)

combined_plot

setwd("~/GitHub/PlanktonRUE-Project/Figures/CE_plots")
ggsave("CE_zoo_RUE.png", 
       plot = combined_plot, width = 25, height = 11, dpi = 300, device = ragg::agg_png)

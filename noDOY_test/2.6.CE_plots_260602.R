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

setwd("~/GitHub/PlanktonRUE-Project/noDOY_test/RDATA")

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
                  "Zooplankton_MaxDepth", "daylight_hours")
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

#Phyto genus richness #--------------------------------------------------
setwd("~/GitHub/PlanktonRUE-Project/noDOY_test/RDATA")
load("BRM_phyto_GD_P_260602.RDATA")

p1 <- plot_ce_simple("Temperature", model, "genus_richness", xlab = "Temperature", 
                     ylab = "Phytoplankton genus richness", title = "Temperature")
p2 <- plot_ce_simple("Salinity", model, "genus_richness", xlab = "Salinity", 
                     ylab = "Phytoplankton genus richness", title = "Salinity")
p3 <- plot_ce_simple("pH", model, "genus_richness",xlab = "pH", 
                     ylab = "Phytoplankton genus richness", title = "pH")
p4 <- plot_ce_simple("TotalNitrogen_unfiltered", model, "genus_richness",
                     xlab = "Nitrogen concentrations", 
                     ylab = "Phytoplankton genus richness", title = "Nitrogen concentrations")
p5 <- plot_ce_simple("MaxDepth", model, "genus_richness", xlab = "Depth", 
                      ylab = "Phytoplankton genus richness", title = "Depth")
p7 <- plot_ce_simple("daylight_hours", model, "genus_richness", xlab = "Daylight Hours", 
                      ylab = "Phytoplankton genus richness", title = "Daylight Hours")

combined_plot <- wrap_plots(
  p1, p2, p3, p4, p5,
  p7,
  
  ncol = 4
)

combined_plot

setwd("~/GitHub/PlanktonRUE-Project/noDOY_test/Results")
ggsave("CE_Phyto_GD.png", 
       plot = combined_plot, width = 20, height = 8, dpi = 300, device = ragg::agg_png)

#Phyto FEve #-----------------------------------------------------------------
setwd("~/GitHub/PlanktonRUE-Project/noDOY_test/RDATA")
load("BRM_phyto_FEve_P_260602.RDATA")

p1 <- plot_ce_simple("Temperature", model, "FEve", xlab = "Temperature", 
                     ylab = "Phytoplankton FEve", title = "Temperature")
p2 <- plot_ce_simple("Salinity", model, "FEve", xlab = "Salinity", 
                     ylab = "Phytoplankton FEve", title = "Salinity")
p3 <- plot_ce_simple("pH", model, "FEve",xlab = "pH", 
                     ylab = "Phytoplankton FEve", title = "pH")
p4 <- plot_ce_simple("TotalNitrogen_unfiltered", model, "FEve",
                     xlab = "Nitrogen concentrations", 
                     ylab = "Phytoplankton FEve", title = "Nitrogen concentrations")
p5 <- plot_ce_simple("MaxDepth", model, "FEve", xlab = "Depth", 
                     ylab = "Phytoplankton FEve", title = "Depth")
p7 <- plot_ce_simple("daylight_hours", model, "FEve", xlab = "Daylight Hours", 
                     ylab = "Phytoplankton FEve", title = "Daylight Hours")

combined_plot <- wrap_plots(
  p1, p2, p3, p4, p5,
  p7,
  
  ncol = 4
)

combined_plot

setwd("~/GitHub/PlanktonRUE-Project/noDOY_test/Results")
ggsave("CE_Phyto_FEve.png", 
       plot = combined_plot, width = 20, height = 8, dpi = 300, device = ragg::agg_png)

#Phyto FDis #-----------------------------------------------------------------
setwd("~/GitHub/PlanktonRUE-Project/noDOY_test/RDATA")
load("BRM_phyto_FDis_P_260602.RDATA")

p1 <- plot_ce_simple("Temperature", model, "FDis", xlab = "Temperature", 
                     ylab = "Phytoplankton FDis", title = "Temperature")
p2 <- plot_ce_simple("Salinity", model, "FDis", xlab = "Salinity", 
                     ylab = "Phytoplankton FDis", title = "Salinity")
p3 <- plot_ce_simple("pH", model, "FDis",xlab = "pH", 
                     ylab = "Phytoplankton FDis", title = "pH")
p4 <- plot_ce_simple("TotalNitrogen_unfiltered", model, "FDis",
                     xlab = "Nitrogen concentrations", 
                     ylab = "Phytoplankton FDis", title = "Nitrogen concentrations")
p5 <- plot_ce_simple("MaxDepth", model, "FDis", xlab = "Depth", 
                     ylab = "Phytoplankton FDis", title = "Depth")
p7 <- plot_ce_simple("daylight_hours", model, "FDis", xlab = "Daylight Hours", 
                     ylab = "Phytoplankton FDis", title = "Daylight Hours")

combined_plot <- wrap_plots(
  p1, p2, p3, p4, p5,
  p7,
  
  ncol = 4
)

combined_plot

setwd("~/GitHub/PlanktonRUE-Project/noDOY_test/Results")
ggsave("CE_Phyto_FDis.png", 
       plot = combined_plot, width = 20, height = 8, dpi = 300, device = ragg::agg_png)

#Cyanobacteria Dominance #------------------------------------------------------
setwd("~/GitHub/PlanktonRUE-Project/noDOY_test/RDATA")
load("BRM_cyanodom_P_260602_Z.RDATA")

p1 <- plot_ce_simple("Temperature", model, "cyano_prop", xlab = "Temperature", 
                     ylab = "Cyanobacteria Dominance", title = "Temperature")
p2 <- plot_ce_simple("Salinity", model, "cyano_prop", xlab = "Salinity", 
                     ylab = "Cyanobacteria Dominance", title = "Salinity")
p3 <- plot_ce_simple("pH", model, "cyano_prop",xlab = "pH", 
                     ylab = "Cyanobacteria Dominance", title = "pH")
p4 <- plot_ce_simple("TotalNitrogen_unfiltered", model, "cyano_prop",
                     xlab = "Nitrogen concentrations", 
                     ylab = "Cyanobacteria Dominance", title = "Nitrogen concentrations")
p5 <- plot_ce_simple("MaxDepth", model, "cyano_prop", xlab = "Depth", 
                     ylab = "Cyanobacteria Dominance", title = "Depth")
p7 <- plot_ce_simple("daylight_hours", model, "cyano_prop", xlab = "Daylight Hours", 
                     ylab = "Cyanobacteria Dominance", title = "Daylight Hours")

combined_plot <- wrap_plots(
  p1, p2, p3, p4, p5,
  p7,
  
  ncol = 4
)

combined_plot

setwd("~/GitHub/PlanktonRUE-Project/noDOY_test/Results")
ggsave("CE_Phyto_CD.png", 
       plot = combined_plot, width = 20, height = 8, dpi = 300, device = ragg::agg_png)

#Zooplankton Species richness #------------------------------------------------------
setwd("~/GitHub/PlanktonRUE-Project/noDOY_test/RDATA")
load("BRM_zoo_SR_CD_P_260602.RDATA")

p1 <- plot_ce_simple("Temperature", model, "species_richness", xlab = "Temperature", 
                     ylab = "Zooplankton Species Richness", title = "Temperature")
p2 <- plot_ce_simple("Salinity", model, "species_richness", xlab = "Salinity", 
                     ylab = "Zooplankton Species Richness", title = "Salinity")
p3 <- plot_ce_simple("pH", model, "species_richness",xlab = "pH", 
                     ylab = "Zooplankton Species Richness", title = "pH")
p4 <- plot_ce_simple("TotalNitrogen_unfiltered", model, "species_richness",
                     xlab = "Nitrogen concentrations", 
                     ylab = "Zooplankton Species Richness", title = "Nitrogen concentrations")
p5 <- plot_ce_simple("genus_richness", model, "species_richness", 
                     xlab = "Phytoplankton genus richness", ylab = "Zooplankton Species Richness", 
                     title = "Phytoplankton genus richness")
p6 <- plot_ce_simple("FEve", model, "species_richness", xlab = "Phytoplankton FEve", 
                     ylab = "Zooplankton Species Richness", title = "Phytoplankton FEve")
p7 <- plot_ce_simple("FDis", model, "species_richness", xlab = "Phytoplankton FDis", 
                     ylab = "Zooplankton Species Richness", title = "Phytoplankton FDis")
p8 <- plot_ce_simple("cyano_prop", model, "species_richness", xlab = "Cyanobacteria Dominance", 
                     ylab = "Zooplankton Species Richness", title = "Cyanobacteria Dominance")
p9 <- plot_ce_simple("Zooplankton_MaxDepth", model, "species_richness", xlab = "Depth", 
                     ylab = "Zooplankton Species Richness", title = "Depth")
p11 <- plot_ce_simple("daylight_hours", model, "species_richness", xlab = "Daylight Hours", 
                     ylab = "Zooplankton Species Richness", title = "Daylight Hours")

combined_plot <- wrap_plots(
  p1, p2, p3, p4, p5,
  p6, p7, p8, p9, p11,
  
  ncol = 5
)

combined_plot

setwd("~/GitHub/PlanktonRUE-Project/noDOY_test/Results")
ggsave("CE_Zoo_SD.png", 
       plot = combined_plot, width = 20, height = 10, dpi = 300, device = ragg::agg_png)

#Zooplankton FEve #-------------------------------------------------------------
setwd("~/GitHub/PlanktonRUE-Project/noDOY_test/RDATA")
load("BRM_zoo_FE_CD_P_260602.RDATA")

p1 <- plot_ce_simple("Temperature", model, "zoo_FEve", xlab = "Temperature", 
                     ylab = "Zooplankton FEve", title = "Temperature")
p2 <- plot_ce_simple("Salinity", model, "zoo_FEve", xlab = "Salinity", 
                     ylab = "Zooplankton FEve", title = "Salinity")
p3 <- plot_ce_simple("pH", model, "zoo_FEve",xlab = "pH", 
                     ylab = "Zooplankton FEve", title = "pH")
p4 <- plot_ce_simple("TotalNitrogen_unfiltered", model, "zoo_FEve",
                     xlab = "Nitrogen concentrations", 
                     ylab = "Zooplankton FEve", title = "Nitrogen concentrations")
p5 <- plot_ce_simple("genus_richness", model, "zoo_FEve", 
                     xlab = "Phytoplankton genus richness", ylab = "Zooplankton FEve", 
                     title = "Phytoplankton genus richness")
p6 <- plot_ce_simple("FEve", model, "zoo_FEve", xlab = "Phytoplankton FEve", 
                     ylab = "Zooplankton FEve", title = "Phytoplankton FEve")
p7 <- plot_ce_simple("FDis", model, "zoo_FEve", xlab = "Phytoplankton FDis", 
                     ylab = "Zooplankton FEve", title = "Phytoplankton FDis")
p8 <- plot_ce_simple("cyano_prop", model, "zoo_FEve", xlab = "Cyanobacteria Dominance", 
                     ylab = "Zooplankton FEve", title = "Cyanobacteria Dominance")
p9 <- plot_ce_simple("Zooplankton_MaxDepth", model, "zoo_FEve", xlab = "Depth", 
                     ylab = "Zooplankton FEve", title = "Depth")
p11 <- plot_ce_simple("daylight_hours", model, "zoo_FEve", xlab = "Daylight Hours", 
                      ylab = "Zooplankton FEve", title = "Daylight Hours")

combined_plot <- wrap_plots(
  p1, p2, p3, p4, p5,
  p6, p7, p8, p9, p11,
  
  ncol = 5
)

combined_plot

setwd("~/GitHub/PlanktonRUE-Project/noDOY_test/Results")
ggsave("CE_Zoo_FEve.png", 
       plot = combined_plot, width = 20, height = 10, dpi = 300, device = ragg::agg_png)

#Zooplankton FDis #-------------------------------------------------------------
setwd("~/GitHub/PlanktonRUE-Project/noDOY_test/RDATA")
load("BRM_zoo_FD_CD_P_260602.RDATA")

p1 <- plot_ce_simple("Temperature", model, "zoo_FDis", xlab = "Temperature", 
                     ylab = "Zooplankton FDis", title = "Temperature")
p2 <- plot_ce_simple("Salinity", model, "zoo_FDis", xlab = "Salinity", 
                     ylab = "Zooplankton FDis", title = "Salinity")
p3 <- plot_ce_simple("pH", model, "zoo_FDis",xlab = "pH", 
                     ylab = "Zooplankton FDis", title = "pH")
p4 <- plot_ce_simple("TotalNitrogen_unfiltered", model, "zoo_FDis",
                     xlab = "Nitrogen concentrations", 
                     ylab = "Zooplankton FDis", title = "Nitrogen concentrations")
p5 <- plot_ce_simple("genus_richness", model, "zoo_FDis", 
                     xlab = "Phytoplankton genus richness", ylab = "Zooplankton FDis", 
                     title = "Phytoplankton genus richness")
p6 <- plot_ce_simple("FEve", model, "zoo_FDis", xlab = "Phytoplankton FEve", 
                     ylab = "Zooplankton FDis", title = "Phytoplankton FEve")
p7 <- plot_ce_simple("FDis", model, "zoo_FDis", xlab = "Phytoplankton FDis", 
                     ylab = "Zooplankton FDis", title = "Phytoplankton FDis")
p8 <- plot_ce_simple("cyano_prop", model, "zoo_FDis", xlab = "Cyanobacteria Dominance", 
                     ylab = "Zooplankton FDis", title = "Cyanobacteria Dominance")
p9 <- plot_ce_simple("Zooplankton_MaxDepth", model, "zoo_FDis", xlab = "Depth", 
                     ylab = "Zooplankton FDis", title = "Depth")
p11 <- plot_ce_simple("daylight_hours", model, "zoo_FDis", xlab = "Daylight Hours", 
                      ylab = "Zooplankton FDis", title = "Daylight Hours")

combined_plot <- wrap_plots(
  p1, p2, p3, p4, p5,
  p6, p7, p8, p9, p11,
  
  ncol = 5
)

combined_plot

setwd("~/GitHub/PlanktonRUE-Project/noDOY_test/Results")
ggsave("CE_Zoo_FDis.png", 
       plot = combined_plot, width = 20, height = 10, dpi = 300, device = ragg::agg_png)

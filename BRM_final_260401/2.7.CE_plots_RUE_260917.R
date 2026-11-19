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
phyto_vars <- c("genus_richness", "phyto_evenness", "phyto_FEve_ZD", "phyto_FDis_ZD", "phyto_FEve_PRUE", 
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

#Phyto RUE per P #--------------------------------------------------
setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA")
load("BRM_phytoRUEP_CD_PI_261013_All.RDATA")

p1 <- plot_ce_simple("Temperature", model, "log_RUE_B_P", xlab = "Temperature", 
                     ylab = "RUE per P", title = "Temperature")
p2 <- plot_ce_simple("Salinity", model, "log_RUE_B_P", xlab = "Salinity", ylab = "RUE per P", title = "Salinity")
p3 <- plot_ce_simple("pH", model, "log_RUE_B_P",xlab = "pH", ylab = "RUE per P", title = "pH")
p4 <- plot_ce_simple("TotalNitrogen_unfiltered", model, "log_RUE_B_P",
                     xlab = "Nitrogen concentrations", ylab = "RUE per P", title = "Nitrogen concentrations")
p5 <- plot_ce_simple("genus_richness", model, "log_RUE_B_P", 
                     xlab = "Phytoplankton Genus Richness", ylab = "RUE per P", 
                     title = "Phytoplankton Genus Richness")
p6 <- plot_ce_simple("phyto_evenness", model, "log_RUE_B_P", 
                     xlab = "Phytoplankton Evenness", ylab = "RUE per P", 
                     title = "Phytoplankton Evenness")
p7 <- plot_ce_simple("phyto_FEve_PRUE", model, "log_RUE_B_P", xlab = "Phytoplankton FEve", 
                     ylab = "RUE per P", title = "Phytoplankton FEve")
p8 <- plot_ce_simple("phyto_FDis_PRUE", model, "log_RUE_B_P", xlab = "Phytoplankton FDis", 
                     ylab = "RUE per P", title = "Phytoplankton FDis")
p9 <- plot_ce_simple("cyano_prop", model, "log_RUE_B_P", xlab = "Cyanobacteria Dominance", 
                     ylab = "RUE per P", title = "Cyanobacteria Dominance")
p10 <- plot_ce_simple("species_richness", model, "log_RUE_B_P", xlab = "Zooplankton Species Richness", 
                     ylab = "RUE per P", title = "Zooplankton Species Richness")
p11 <- plot_ce_simple("zoo_evenness", model, "log_RUE_B_P", xlab = "Zooplankton Evenness", 
                     ylab = "RUE per P", title = "Zooplankton Evenness")
p12 <- plot_ce_simple("zoo_FEve", model, "log_RUE_B_P", xlab = "Zooplankton FEve", 
                      ylab = "RUE per P", title = "Zooplankton FEve")
p13 <- plot_ce_simple("zoo_FDis", model, "log_RUE_B_P", xlab = "Zooplankton FDis", 
                      ylab = "RUE per P", title = "Zooplankton FDis")
p14 <- plot_ce_simple("genus_richness:Temperature", model, "log_RUE_B_P", xlab = "Phytoplankton Genus Richness", 
                      ylab = "RUE per P", title = "Phytoplankton Genus Richness:Temperature", 
                      legend_title = "Temperature")
p15 <- plot_ce_simple("genus_richness:Salinity", model, "log_RUE_B_P", xlab = "Phytoplankton Genus Richness", 
                      ylab = "RUE per P", title = "Phytoplankton Genus Richness:Salinity", 
                      legend_title = "Salinity")
p16 <- plot_ce_simple("phyto_evenness:Temperature", model, "log_RUE_B_P", xlab = "Phytoplankton Evenness", 
                      ylab = "RUE per P", title = "Phytoplankton Evenness:Temperature", 
                      legend_title = "Temperature")
p17 <- plot_ce_simple("phyto_evenness:Salinity", model, "log_RUE_B_P", xlab = "Phytoplankton Evenness", 
                      ylab = "RUE per P", title = "Phytoplankton Evenness:Salinity", 
                      legend_title = "Salinity")
p18 <- plot_ce_simple("phyto_evenness:TotalNitrogen_unfiltered", model, "log_RUE_B_P", xlab = "Phytoplankton Evenness", 
                      ylab = "RUE per P", title = "Phytoplankton Evenness:Nitrogen concentrations", 
                      legend_title = "Nitrogen concentrations")
p19 <- plot_ce_simple("phyto_FEve_PRUE:Temperature", model, "log_RUE_B_P", xlab = "Phytoplankton FEve", 
                      ylab = "RUE per P", title = "Phytoplankton FEve:Temperature", 
                      legend_title = "Temperature")
p20 <- plot_ce_simple("phyto_FEve_PRUE:Salinity", model, "log_RUE_B_P", xlab = "Phytoplankton FEve", 
                      ylab = "RUE per P", title = "Phytoplankton FEve:Salinity", 
                      legend_title = "Salinity")
p21 <- plot_ce_simple("phyto_FEve_PRUE:TotalNitrogen_unfiltered", model, "log_RUE_B_P", xlab = "Phytoplankton FEve", 
                      ylab = "RUE per P", title = "Phytoplankton FEve:Nitrogen concentrations", 
                      legend_title = "Nitrogen concentrations")
p22 <- plot_ce_simple("phyto_FDis_PRUE:Temperature", model, "log_RUE_B_P", xlab = "Phytoplankton FDis", 
                      ylab = "RUE per P", title = "Phytoplankton FDis:Temperature", 
                      legend_title = "Temperature")
p23 <- plot_ce_simple("phyto_FDis_PRUE:Salinity", model, "log_RUE_B_P", xlab = "Phytoplankton FDis", 
                      ylab = "RUE per P", title = "Phytoplankton FDis:Salinity", 
                      legend_title = "Salinity")
p24 <- plot_ce_simple("phyto_FDis_PRUE:TotalNitrogen_unfiltered", model, "log_RUE_B_P", xlab = "Phytoplankton FDis", 
                      ylab = "RUE per P", title = "Phytoplankton FDis:Nitrogen concentrations", 
                      legend_title = "Nitrogen concentrations")
p25 <- plot_ce_simple("cyano_prop:Temperature", model, "log_RUE_B_P", xlab = "Cyanobacteria Dominance", 
                      ylab = "RUE per P", title = "Cyanobacteria Dominance:Temperature", 
                      legend_title = "Temperature")
p26 <- plot_ce_simple("species_richness:Temperature", model, "log_RUE_B_P", xlab = "Zooplankton Species Richness", 
                      ylab = "RUE per P", title = "Zooplankton Species Richness:Temperature", 
                      legend_title = "Temperature")
p27 <- plot_ce_simple("zoo_evenness:Temperature", model, "log_RUE_B_P", xlab = "Zooplankton Evenness", 
                      ylab = "RUE per P", title = "Zooplankton Evenness:Temperature", 
                      legend_title = "Temperature")
p28 <- plot_ce_simple("zoo_evenness:TotalNitrogen_unfiltered", model, "log_RUE_B_P", xlab = "Zooplankton Evenness", 
                      ylab = "RUE per P", title = "Zooplankton Evenness:Nitrogen concentrations", 
                      legend_title = "Nitrogen concentrations")
p29 <- plot_ce_simple("zoo_FEve:Temperature", model, "log_RUE_B_P", xlab = "Zooplankton FEve", 
                      ylab = "RUE per P", title = "Zooplankton FEve:Temperature", 
                      legend_title = "Temperature")
p30 <- plot_ce_simple("zoo_FEve:Salinity", model, "log_RUE_B_P", xlab = "Zooplankton FEve", 
                      ylab = "RUE per P", title = "Zooplankton FEve:Salinity", 
                      legend_title = "Salinity")
p31 <- plot_ce_simple("zoo_FDis:Temperature", model, "log_RUE_B_P", xlab = "Zooplankton FDis", 
                      ylab = "RUE per P", title = "Zooplankton FDis:Temperature", 
                      legend_title = "Temperature")
p32 <- plot_ce_simple("zoo_FDis:TotalNitrogen_unfiltered", model, "log_RUE_B_P", xlab = "Zooplankton FDis", 
                      ylab = "RUE per P", title = "Zooplankton FDis:Nitrogen concentrations", 
                      legend_title = "Nitrogen concentrations")
p33 <- plot_ce_simple("MaxDepth", model, "log_RUE_B_P", xlab = "Depth", 
                      ylab = "RUE per P", title = "Depth")
p34 <- plot_ce_simple("daylight_hours", model, "log_RUE_B_P", xlab = "Daylight Hours", 
                      ylab = "RUE per P", title = "Daylight Hours")

combined_plot <- wrap_plots(
  p1, p2, p3, p4, p5, p6, 
  p7, p8, p9, p10, p11, p12, p13, 
  
  p14 + plot_layout(widths = 2), p15 + plot_layout(widths = 2), 
  p16 + plot_layout(widths = 2), p17 + plot_layout(widths = 2), 
  p18 + plot_layout(widths = 2), p19 + plot_layout(widths = 2), 
  p20 + plot_layout(widths = 2), p21 + plot_layout(widths = 2), 
  p22 + plot_layout(widths = 2), p23 + plot_layout(widths = 2), 
  p24 + plot_layout(widths = 2), p25 + plot_layout(widths = 2), 
  p26 + plot_layout(widths = 2), p27 + plot_layout(widths = 2), 
  p28 + plot_layout(widths = 2), p29 + plot_layout(widths = 2), 
  p30 + plot_layout(widths = 2), p31 + plot_layout(widths = 2), 
  p32 + plot_layout(widths = 2),
  
  p33, p34,
  
  ncol = 6
)

combined_plot

setwd("~/GitHub/PlanktonRUE-Project/Figures/CE_plots")
ggsave("CE_RUEP.png", 
       plot = combined_plot, width = 30, height = 20, dpi = 300, device = ragg::agg_png)

#Phyto RUE per N #--------------------------------------------------
setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA")
load("BRM_phytoRUEN_CD_PI_261013_All.RDATA")

p1 <- plot_ce_simple("Temperature", model, "log_RUE_B_N", xlab = "Temperature", 
                     ylab = "RUE per N", title = "Temperature")
p2 <- plot_ce_simple("Salinity", model, "log_RUE_B_N", xlab = "Salinity", 
                     ylab = "RUE per N", title = "Salinity")
p3 <- plot_ce_simple("pH", model, "log_RUE_B_N", xlab = "pH", 
                     ylab = "RUE per N", title = "pH")
p4 <- plot_ce_simple("TotalNitrogen_unfiltered", model, "log_RUE_B_N", xlab = "Nitrogen concentrations", 
                     ylab = "RUE per N", title = "Nitrogen concentrations")
p5 <- plot_ce_simple("genus_richness", model, "log_RUE_B_N", 
                     xlab = "Phytoplankton Genus Richness", ylab = "RUE per N", 
                     title = "Phytoplankton Genus Richness")
p6 <- plot_ce_simple("phyto_evenness", model, "log_RUE_B_N", 
                     xlab = "Phytoplankton Evenness", ylab = "RUE per N", 
                     title = "Phytoplankton Evenness")
p7 <- plot_ce_simple("phyto_FEve_PRUE", model, "log_RUE_B_N", xlab = "Phytoplankton FEve", 
                     ylab = "RUE per N", title = "Phytoplankton FEve")
p8 <- plot_ce_simple("phyto_FDis_PRUE", model, "log_RUE_B_N", xlab = "Phytoplankton FDis", 
                     ylab = "RUE per N", title = "Phytoplankton FDis")
p9 <- plot_ce_simple("cyano_prop", model, "log_RUE_B_N", xlab = "Cyanobacteria Dominance", 
                     ylab = "RUE per N", title = "Cyanobacteria Dominance")
p10 <- plot_ce_simple("species_richness", model, "log_RUE_B_N", xlab = "Zooplankton Species Richness", 
                      ylab = "RUE per N", title = "Zooplankton Species Richness")
p11 <- plot_ce_simple("zoo_evenness", model, "log_RUE_B_N", xlab = "Zooplankton Evenness", 
                      ylab = "RUE per N", title = "Zooplankton Evenness")
p12 <- plot_ce_simple("zoo_FEve", model, "log_RUE_B_N", xlab = "Zooplankton FEve", 
                      ylab = "RUE per N", title = "Zooplankton FEve")
p13 <- plot_ce_simple("zoo_FDis", model, "log_RUE_B_N", xlab = "Zooplankton FDis", 
                      ylab = "RUE per N", title = "Zooplankton FDis")
p14 <- plot_ce_simple("genus_richness:Temperature", model, "log_RUE_B_N", xlab = "Phytoplankton Genus Richness", 
                      ylab = "RUE per N", title = "Phytoplankton Genus Richness:Temperature", 
                      legend_title = "Temperature")
p15 <- plot_ce_simple("genus_richness:Salinity", model, "log_RUE_B_N", xlab = "Phytoplankton Genus Richness", 
                      ylab = "RUE per N", title = "Phytoplankton Genus Richness:Salinity", 
                      legend_title = "Salinity")
p16 <- plot_ce_simple("genus_richness:TotalNitrogen_unfiltered", model, "log_RUE_B_N", xlab = "Phytoplankton Genus Richness", 
                      ylab = "RUE per N", title = "Phytoplankton Genus Richness:Nitrogen concentrations", 
                      legend_title = "Nitrogen concentrations")
p17 <- plot_ce_simple("phyto_evenness:Salinity", model, "log_RUE_B_N", xlab = "Phytoplankton Evenness", 
                      ylab = "RUE per N", title = "Phytoplankton Evenness:Salinity", 
                      legend_title = "Salinity")
p18 <- plot_ce_simple("phyto_evenness:TotalNitrogen_unfiltered", model, "log_RUE_B_N", xlab = "Phytoplankton Evenness", 
                      ylab = "RUE per N", title = "Phytoplankton Evenness:Nitrogen concentrations", 
                      legend_title = "Nitrogen concentrations")
p19 <- plot_ce_simple("phyto_FEve_PRUE:Salinity", model, "log_RUE_B_N", xlab = "Phytoplankton FEve", 
                      ylab = "RUE per N", title = "Phytoplankton FEve:Salinity", 
                      legend_title = "Salinity")
p20 <- plot_ce_simple("phyto_FEve_PRUE:TotalNitrogen_unfiltered", model, "log_RUE_B_N", xlab = "Phytoplankton FEve", 
                      ylab = "RUE per N", title = "Phytoplankton FEve:Nitrogen concentrations", 
                      legend_title = "Nitrogen concentrations")
p21 <- plot_ce_simple("phyto_FDis_PRUE:Salinity", model, "log_RUE_B_N", xlab = "Phytoplankton FDis", 
                      ylab = "RUE per N", title = "Phytoplankton FDis:Salinity", 
                      legend_title = "Salinity")
p22 <- plot_ce_simple("phyto_FDis_PRUE:TotalNitrogen_unfiltered", model, "log_RUE_B_N", xlab = "Phytoplankton FDis", 
                      ylab = "RUE per N", title = "Phytoplankton FDis:Nitrogen concentrations", 
                      legend_title = "Nitrogen concentrations")
p23 <- plot_ce_simple("cyano_prop:Temperature", model, "log_RUE_B_N", xlab = "Cyanobacteria Dominance", 
                      ylab = "RUE per N", title = "Cyanobacteria Dominance:Temperature", 
                      legend_title = "Temperature")
p24 <- plot_ce_simple("cyano_prop:Salinity", model, "log_RUE_B_N", xlab = "Cyanobacteria Dominance", 
                      ylab = "RUE per N", title = "Cyanobacteria Dominance:Salinity", 
                      legend_title = "Salinity")
p25 <- plot_ce_simple("species_richness:Temperature", model, "log_RUE_B_N", xlab = "Zooplankton Species Richness", 
                      ylab = "RUE per N", title = "Zooplankton Species Richness:Temperature", 
                      legend_title = "Temperature")
p26 <- plot_ce_simple("zoo_evenness:Temperature", model, "log_RUE_B_N", xlab = "Zooplankton Evenness", 
                      ylab = "RUE per N", title = "Zooplankton Evenness:Temperature", 
                      legend_title = "Temperature")
p27 <- plot_ce_simple("zoo_evenness:TotalNitrogen_unfiltered", model, "log_RUE_B_N", xlab = "Zooplankton Evenness", 
                      ylab = "RUE per N", title = "Zooplankton Evenness:Nitrogen concentrations", 
                      legend_title = "Nitrogen concentrations")
p28 <- plot_ce_simple("zoo_FEve:Salinity", model, "log_RUE_B_N", xlab = "Zooplankton FEve", 
                      ylab = "RUE per N", title = "Zooplankton FEve:Salinity", 
                      legend_title = "Salinity")
p29 <- plot_ce_simple("zoo_FDis:Temperature", model, "log_RUE_B_N", xlab = "Zooplankton FDis", 
                      ylab = "RUE per N", title = "Zooplankton FDis:Temperature", 
                      legend_title = "Temperature")
p30 <- plot_ce_simple("zoo_FDis:TotalNitrogen_unfiltered", model, "log_RUE_B_N", xlab = "Zooplankton FDis", 
                      ylab = "RUE per N", title = "Zooplankton FDis:Nitrogen concentrations", 
                      legend_title = "Nitrogen concentrations")
p31 <- plot_ce_simple("MaxDepth", model, "log_RUE_B_N", xlab = "Depth", 
                      ylab = "RUE per N", title = "Depth")
p32 <- plot_ce_simple("daylight_hours", model, "log_RUE_B_N", xlab = "Daylight Hours", 
                      ylab = "RUE per N", title = "Daylight Hours")

combined_plot <- wrap_plots(
  p1, p2, p3, p4, p5,
  p6, p7, p8, p9, p10,
  p11, p12, p13,
  
  p14 + plot_layout(widths = 2), p15 + plot_layout(widths = 2), 
  p16 + plot_layout(widths = 2), p17 + plot_layout(widths = 2), 
  p18 + plot_layout(widths = 2), p19 + plot_layout(widths = 2), 
  p20 + plot_layout(widths = 2), p21 + plot_layout(widths = 2), 
  p22 + plot_layout(widths = 2), p23 + plot_layout(widths = 2), 
  p24 + plot_layout(widths = 2), p25 + plot_layout(widths = 2), 
  p26 + plot_layout(widths = 2), p27 + plot_layout(widths = 2), 
  p28 + plot_layout(widths = 2), p29 + plot_layout(widths = 2), 
  p30 + plot_layout(widths = 2), 
  
  p31, p32,
  
  ncol = 6
)

combined_plot

setwd("~/GitHub/PlanktonRUE-Project/Figures/CE_plots")
ggsave("CE_RUEN.png", 
       plot = combined_plot, width = 30, height = 20, dpi = 300, device = ragg::agg_png)

#Zoo RUE #--------------------------------------------------
setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA")
load("BRM_zoo_RUE_CD_PI_261013_All.RDATA")

p1 <- plot_ce_simple("Temperature", model, "log_RUE_zooplankton", xlab = "Temperature", 
                     ylab = "Zooplankton RUE", title = "Temperature")
p2 <- plot_ce_simple("Salinity", model, "log_RUE_zooplankton", xlab = "Salinity", 
                     ylab = "Zooplankton RUE", title = "Salinity")
p3 <- plot_ce_simple("pH", model, "log_RUE_zooplankton", xlab = "pH", 
                     ylab = "Zooplankton RUE", title = "pH")
p4 <- plot_ce_simple("TotalNitrogen_unfiltered", model, "log_RUE_zooplankton", xlab = "Nitrogen concentrations", 
                     ylab = "Zooplankton RUE", title = "Nitrogen concentrations")
p5 <- plot_ce_simple("genus_richness", model, "log_RUE_zooplankton", xlab = "Phytoplankton Genus Richness", 
                     ylab = "Zooplankton RUE", title = "Phytoplankton Genus Richness")
p6 <- plot_ce_simple("phyto_evenness", model, "log_RUE_zooplankton", 
                     xlab = "Phytoplankton Evenness", ylab = "Zooplankton RUE", 
                     title = "Phytoplankton Evenness")
p7 <- plot_ce_simple("phyto_FEve_ZRUE", model, "log_RUE_zooplankton", xlab = "Phytoplankton FEve", 
                     ylab = "Zooplankton RUE", title = "Phytoplankton FEve")
p8 <- plot_ce_simple("phyto_FDis_ZRUE", model, "log_RUE_zooplankton", xlab = "Phytoplankton FDis", 
                     ylab = "Zooplankton RUE", title = "Phytoplankton FDis")
p9 <- plot_ce_simple("cyano_prop", model, "log_RUE_zooplankton", xlab = "Cyanobacteria Dominance", 
                     ylab = "Zooplankton RUE", title = "Cyanobacteria Dominance")
p10 <- plot_ce_simple("species_richness", model, "log_RUE_zooplankton", xlab = "Zooplankton Species Richness", 
                     ylab = "Zooplankton RUE", title = "Zooplankton Species Richness")
p11 <- plot_ce_simple("zoo_evenness", model, "log_RUE_zooplankton", xlab = "Zooplankton Evenness", 
                      ylab = "Zooplankton RUE", title = "Zooplankton Evenness")
p12 <- plot_ce_simple("zoo_FEve", model, "log_RUE_zooplankton", xlab = "Zooplankton FEve", 
                      ylab = "Zooplankton RUE", title = "Zooplankton Feve")
p13 <- plot_ce_simple("zoo_FDis", model, "log_RUE_zooplankton", xlab = "Zooplankton FDis", 
                      ylab = "Zooplankton RUE", title = "Zooplankton FDis")
p14 <- plot_ce_simple("genus_richness:Temperature", model, "log_RUE_zooplankton", xlab = "Phytoplankton Genus Richness", 
                      ylab = "Zooplankton RUE", title = "Phytoplankton Genus Richness:Temperature",
                      legend_title = "Temperature")
p15 <- plot_ce_simple("genus_richness:Salinity", model, "log_RUE_zooplankton", xlab = "Phytoplankton Genus Richness", 
                      ylab = "Zooplankton RUE", title = "Phytoplankton Genus Richness:Salinity",
                      legend_title = "Salinity")
p16 <- plot_ce_simple("genus_richness:TotalNitrogen_unfiltered", model, "log_RUE_zooplankton", xlab = "Phytoplankton Genus Richness", 
                      ylab = "Zooplankton RUE", title = "Phytoplankton Genus Richness:Nitrogen concentrations",
                      legend_title = "Nitrogen concentrations")
p17 <- plot_ce_simple("phyto_evenness:Salinity", model, "log_RUE_zooplankton", xlab = "Phytoplankton Evenness", 
                      ylab = "Zooplankton RUE", title = "Phytoplankton Evenness:Salinity",
                      legend_title = "Salinity")
p18 <- plot_ce_simple("phyto_FEve_ZRUE:Salinity", model, "log_RUE_zooplankton", xlab = "Phytoplankton FEve", 
                      ylab = "Zooplankton RUE", title = "Phytoplankton FEve:Salinity",
                      legend_title = "Salinity")
p19 <- plot_ce_simple("phyto_FEve_ZRUE:TotalNitrogen_unfiltered", model, "log_RUE_zooplankton", xlab = "Phytoplankton FEve", 
                      ylab = "Zooplankton RUE", title = "Phytoplankton FEve:Nitrogen concentrations",
                      legend_title = "Nitrogen concentrations")
p20 <- plot_ce_simple("phyto_FDis_ZRUE:Temperature", model, "log_RUE_zooplankton", xlab = "Phytoplankton FDis", 
                      ylab = "Zooplankton RUE", title = "Phytoplankton FDis:Temperature",
                      legend_title = "Temperature")
p21 <- plot_ce_simple("cyano_prop:TotalNitrogen_unfiltered", model, "log_RUE_zooplankton", xlab = "Cyanobacteria Dominance", 
                      ylab = "Zooplankton RUE", title = "Cyanobacteria Dominance:Nitrogen concentrations",
                      legend_title = "Nitrogen concentrations")
p22 <- plot_ce_simple("species_richness:TotalNitrogen_unfiltered", model, "log_RUE_zooplankton", xlab = "Zooplankton Species Richness", 
                      ylab = "Zooplankton RUE", title = "Zooplankton Species Richness:Nitrogen concentrations",
                      legend_title = "Nitrogen concentrations")
p23 <- plot_ce_simple("zoo_evenness:Temperature", model, "log_RUE_zooplankton", xlab = "Zooplankton Evenness", 
                      ylab = "Zooplankton RUE", title = "Zooplankton Evenness:Temperature",
                      legend_title = "Temperature")
p24 <- plot_ce_simple("zoo_evenness:TotalNitrogen_unfiltered", model, "log_RUE_zooplankton", xlab = "Zooplankton Evenness", 
                      ylab = "Zooplankton RUE", title = "Zooplankton Evenness:Nitrogen concentrations",
                      legend_title = "Nitrogen concentrations")
p25 <- plot_ce_simple("zoo_FEve:Temperature", model, "log_RUE_zooplankton", xlab = "Zooplankton FEve", 
                      ylab = "Zooplankton RUE", title = "Zooplankton FEve:Temperature",
                      legend_title = "Temperature")
p26 <- plot_ce_simple("zoo_FEve:Salinity", model, "log_RUE_zooplankton", xlab = "Zooplankton FEve", 
                      ylab = "Zooplankton RUE", title = "Zooplankton FEve:Salinity",
                      legend_title = "Salinity")
p27 <- plot_ce_simple("zoo_FDis:Temperature", model, "log_RUE_zooplankton", xlab = "Zooplankton FDis", 
                      ylab = "Zooplankton RUE", title = "Zooplankton FDis:Temperature",
                      legend_title = "Temperature")
p28 <- plot_ce_simple("zoo_FDis:TotalNitrogen_unfiltered", model, "log_RUE_zooplankton", xlab = "Zooplankton FDis", 
                      ylab = "Zooplankton RUE", title = "Zooplankton FDis:Nitrogen concentrations",
                      legend_title = "Nitrogen concentrations")
p29 <- plot_ce_simple("Zooplankton_MaxDepth", model, "log_RUE_zooplankton", xlab = "Depth", 
                      ylab = "Zooplankton RUE", title = "Depth")
p30 <- plot_ce_simple("daylight_hours", model, "log_RUE_zooplankton", xlab = "Daylight Hours", 
                      ylab = "Zooplankton RUE", title = "Daylight Hours")

combined_plot <- wrap_plots(
  p1, p2, p3, p4, p5,
  p6, p7, p8, p9, p10,
  p11, p12, p13,
  
  p14 + plot_layout(widths = 2), p15 + plot_layout(widths = 2), 
  p16 + plot_layout(widths = 2), p17 + plot_layout(widths = 2), 
  p18 + plot_layout(widths = 2), p19 + plot_layout(widths = 2), 
  p20 + plot_layout(widths = 2), p21 + plot_layout(widths = 2), 
  p22 + plot_layout(widths = 2), p23 + plot_layout(widths = 2), 
  p24 + plot_layout(widths = 2), p25 + plot_layout(widths = 2), 
  p26 + plot_layout(widths = 2), p27 + plot_layout(widths = 2), 
  p28 + plot_layout(widths = 2), 
  
  p29, p30,
  
  ncol = 6
)

combined_plot

setwd("~/GitHub/PlanktonRUE-Project/Figures/CE_plots")
ggsave("CE_zoo_RUE.png", 
       plot = combined_plot, width = 30, height = 15, dpi = 300, device = ragg::agg_png)

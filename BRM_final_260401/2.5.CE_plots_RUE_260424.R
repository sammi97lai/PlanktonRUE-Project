#################################################################################
####################### Conditional effect plots ###############################
################################################################################

library(dplyr)
library(tidyr)
library(brms)
library(ggplot2)
library(ggpubr)
library(patchwork)

rm(list = ls())

setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA")

#plotting ce function ----------------------------------------------
plot_ce_simple <- function(var, model, resp = NULL) {
  
  ce <- conditional_effects(model, effects = var, resp = resp)
  p <- plot(ce, plot = FALSE)[[1]]
  
  is_interaction <- grepl(":", var)
  
  rawdat <- model$data
  rawdat <- rawdat[!is.na(rawdat[[resp]]), ]
  rawdat <- rawdat[sample(nrow(rawdat), min(5000, nrow(rawdat))), ]
  
  group_color <- get_group_color(var)
  
  p +
    {if(!is_interaction && var %in% names(rawdat)) 
      geom_point(
        data = rawdat,
        aes(x = .data[[var]], y = .data[[resp]]),
        inherit.aes = FALSE,
        color = "black",
        alpha = 0.3
      ) else NULL} +
    ggtitle(var) +
    theme_minimal() +
    theme(
      plot.title = element_text(color = group_color, face = "bold"),
      axis.title = element_text(color = group_color)
    )
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

p1 <- plot_ce_simple("Temperature", model, "log_RUE_B_P")
p2 <- plot_ce_simple("Salinity", model, "log_RUE_B_P")
p3 <- plot_ce_simple("pH", model, "log_RUE_B_P")
p4 <- plot_ce_simple("TotalNitrogen_unfiltered", model, "log_RUE_B_P")
p5 <- plot_ce_simple("genus_richness", model, "log_RUE_B_P")
p6 <- plot_ce_simple("FEve", model, "log_RUE_B_P")
p7 <- plot_ce_simple("FDis", model, "log_RUE_B_P")
p8 <- plot_ce_simple("cyano_prop", model, "log_RUE_B_P")
p9 <- plot_ce_simple("species_richness", model, "log_RUE_B_P")
p10 <- plot_ce_simple("zoo_FEve", model, "log_RUE_B_P")
p11 <- plot_ce_simple("zoo_FDis", model, "log_RUE_B_P")
p12 <- plot_ce_simple("FEve:TotalNitrogen_unfiltered", model, "log_RUE_B_P")
p13 <- plot_ce_simple("FDis:Salinity", model, "log_RUE_B_P")
p14 <- plot_ce_simple("cyano_prop:Salinity", model, "log_RUE_B_P")
p15 <- plot_ce_simple("cyano_prop:TotalNitrogen_unfiltered", model, "log_RUE_B_P")
p16 <- plot_ce_simple("zoo_FDis:TotalNitrogen_unfiltered", model, "log_RUE_B_P")
p17 <- plot_ce_simple("MaxDepth", model, "log_RUE_B_P")
p18 <- plot_ce_simple("DOY", model, "log_RUE_B_P")
p19 <- plot_ce_simple("daylight_hours", model, "log_RUE_B_P")

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

setwd("~/GitHub/PlanktonRUE-Project/Figures")
ggsave("CE_RUEP.png", 
       plot = combined_plot, width = 18, height = 8, dpi = 300)

#Phyto RUE per N #--------------------------------------------------
setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA")
load("BRM_phytoRUEN_CD_PI_260408.RDATA")

p1 <- plot_ce_simple("Temperature", model, "log_RUE_B_N")
p2 <- plot_ce_simple("Salinity", model, "log_RUE_B_N")
p3 <- plot_ce_simple("pH", model, "log_RUE_B_N")
p4 <- plot_ce_simple("TotalNitrogen_unfiltered", model, "log_RUE_B_N")
p5 <- plot_ce_simple("genus_richness", model, "log_RUE_B_N")
p6 <- plot_ce_simple("FEve", model, "log_RUE_B_N")
p7 <- plot_ce_simple("FDis", model, "log_RUE_B_N")
p8 <- plot_ce_simple("cyano_prop", model, "log_RUE_B_N")
p9 <- plot_ce_simple("species_richness", model, "log_RUE_B_N")
p10 <- plot_ce_simple("zoo_FEve", model, "log_RUE_B_N")
p11 <- plot_ce_simple("zoo_FDis", model, "log_RUE_B_N")
p12 <- plot_ce_simple("FEve:TotalNitrogen_unfiltered", model, "log_RUE_B_N")
p13 <- plot_ce_simple("FDis:Salinity", model, "log_RUE_B_N")
p14 <- plot_ce_simple("FDis:TotalNitrogen_unfiltered", model, "log_RUE_B_N")
p15 <- plot_ce_simple("cyano_prop:Salinity", model, "log_RUE_B_N")
p16 <- plot_ce_simple("cyano_prop:TotalNitrogen_unfiltered", model, "log_RUE_B_N")
p17 <- plot_ce_simple("zoo_FDis:TotalNitrogen_unfiltered", model, "log_RUE_B_N")
p18 <- plot_ce_simple("MaxDepth", model, "log_RUE_B_N")
p19 <- plot_ce_simple("DOY", model, "log_RUE_B_N")
p20 <- plot_ce_simple("daylight_hours", model, "log_RUE_B_N")

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

setwd("~/GitHub/PlanktonRUE-Project/Figures")
ggsave("CE_RUEN.png", 
       plot = combined_plot, width = 18, height = 8, dpi = 300)

#Zoo RUE #--------------------------------------------------
setwd("~/GitHub/PlanktonRUE-Project/BRM_final_260401/RDATA")
load("BRM_zoo_RUE_CD_PI_260427.RDATA")

p1 <- plot_ce_simple("Temperature", model, "log_RUE_zooplankton")
p2 <- plot_ce_simple("Salinity", model, "log_RUE_zooplankton")
p3 <- plot_ce_simple("pH", model, "log_RUE_zooplankton")
p4 <- plot_ce_simple("TotalNitrogen_unfiltered", model, "log_RUE_zooplankton")
p5 <- plot_ce_simple("genus_richness", model, "log_RUE_zooplankton")
p6 <- plot_ce_simple("FEve", model, "log_RUE_zooplankton")
p7 <- plot_ce_simple("FDis", model, "log_RUE_zooplankton")
p8 <- plot_ce_simple("cyano_prop", model, "log_RUE_zooplankton")
p9 <- plot_ce_simple("species_richness", model, "log_RUE_zooplankton")
p10 <- plot_ce_simple("zoo_FEve", model, "log_RUE_zooplankton")
p11 <- plot_ce_simple("zoo_FDis", model, "log_RUE_zooplankton")
p12 <- plot_ce_simple("genus_richness:Temperature", model, "log_RUE_zooplankton")
p13 <- plot_ce_simple("genus_richness:Salinity", model, "log_RUE_zooplankton")
p14 <- plot_ce_simple("FDis:Salinity", model, "log_RUE_zooplankton")
p15 <- plot_ce_simple("cyano_prop:Temperature", model, "log_RUE_zooplankton")
p16 <- plot_ce_simple("Zooplankton_MaxDepth", model, "log_RUE_zooplankton")
p17 <- plot_ce_simple("DOY", model, "log_RUE_zooplankton")
p18 <- plot_ce_simple("daylight_hours", model, "log_RUE_zooplankton")

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

setwd("~/GitHub/PlanktonRUE-Project/Figures")
ggsave("CE_zoo_RUE.png", 
       plot = combined_plot, width = 18, height = 8, dpi = 300)

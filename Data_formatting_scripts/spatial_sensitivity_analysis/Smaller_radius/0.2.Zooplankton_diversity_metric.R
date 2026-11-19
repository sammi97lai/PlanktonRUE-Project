library(dplyr)
library(readxl)
library(tidyr)
library(vegan)
library(FD)
library(lubridate)
library(stringr)

setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/RDATA")
load('Zooplankton_REC_clean3_260307.RDATA')

setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/spatial_sensitivity_analysis/Smaller_radius/RDATA")
load('plankton_RUE_data2_260526.RDATA')

####################### Zooplankton Hill Diversity #############################

#clean Zooplankton dataset
Zooplankton_biomass_data <- Zooplankton_REC_clean3 %>%
  select(SampleEvent, Date, latitude, longitude, 
         Zooplankton_MinDepth, Zooplankton_MaxDepth, Species, WetWeight) %>%
  pivot_wider(
    names_from = Species,
    values_from = WetWeight,
    values_fill = 0
  ) # Fills missing combinations with 0

biomass_data <- Zooplankton_biomass_data[-c(1:6)] #keep only abundance data, remove env factors#

Zooplankton_biomass_data2 <- Zooplankton_biomass_data #make new data set

# Hill numbers
Zooplankton_biomass_data2$shannon_zoo <- diversity(biomass_data, index = "shannon")  # q = 1 (Shannon diversity)
Zooplankton_biomass_data2$zoo_Hill2 <- 1 / diversity(biomass_data, index = "invsimpson")  # q = 2 (Simpson diversity)

#Species richness#
# Calculate species richness per sample
Zooplankton_biomass_data2$species_richness <- apply(biomass_data, 1, function(x) sum(x > 0))

# Shannon/Pielou's evenness
Zooplankton_biomass_data2$zoo_evenness <- ifelse(
  Zooplankton_biomass_data2$species_richness > 1,
  Zooplankton_biomass_data2$shannon_zoo / log(Zooplankton_biomass_data2$species_richness),
  NA
)

save(Zooplankton_biomass_data2, file = 'Zooplankton_biomass_data2_260526.RDATA')

zoo_diversity<-Zooplankton_biomass_data2%>%
  select(c(SampleEvent, shannon_zoo, zoo_Hill2, species_richness, zoo_evenness))

#adding species richness data to zooplankton dataset#
plankton_RUE_data3 <- left_join(plankton_RUE_data2, zoo_diversity, by = "SampleEvent")
save(plankton_RUE_data3, file = 'plankton_RUE_data3_260526.RDATA')

################################################################################
######### Functional diversity - functional evenness and disparity #############
#Functional Diversity Metric Calculations----------------------------------------
rm(list = ls())
setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/spatial_sensitivity_analysis/Smaller_radius/RDATA")
load("plankton_RUE_data3_260526.RDATA")
setwd("~/GitHub/PlanktonRUE-Project/Data/Zooplankton/Zooplankton_functional_traits")
Traits <- read_excel("Zooplankton_traits_clean_290925.xlsx")

#turn species name to dataframe and reorder
selected_col_names <- colnames(plankton_RUE_data3)[116:min(155, ncol(plankton_RUE_data3))]
species_table <- data.frame(Species = selected_col_names)
species_table$Species <- sub("_z$", "", species_table$Species)
Traits2 <- Traits %>% 
  filter(Species %in% species_table$Species)%>%
  mutate(Species = factor(Species, levels = species_table$Species)) %>%
  arrange(Species)  # Arrange rows based on the ordered Species

biomass_sub <-plankton_RUE_data3 %>% 
  select(116:155)
colnames(biomass_sub) <- sub("_z$", "", colnames(biomass_sub))
biomass <- biomass_sub %>%
  mutate(Site = paste("Site_", row_number())) %>%
  relocate(Site)
rownames(biomass) <- biomass$Site           
biomass$Site <- NULL                           
biomass <- as.data.frame(biomass)

setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/spatial_sensitivity_analysis/Smaller_radius/RDATA")
write.csv(biomass, "biomass.csv",row.names = FALSE)

# Order species labels alphabetically
zoo_traits <- Traits2[order(Traits2$Species), ]

#substitute NAs in meanBL with values from BLMax
zoo_traits <- zoo_traits %>%
  mutate(meanBodyLength_mm = ifelse(is.na(meanBodyLength_mm),
                                    BodyLengthMax_mm,
                                    meanBodyLength_mm))

#log transform body length traits#
hist(zoo_traits$meanBodyLength_mm)
zoo_traits$log_meanBL <- log(zoo_traits$meanBodyLength_mm) 
hist(zoo_traits$log_meanBL)
hist(zoo_traits$BodyLengthMax_mm)
zoo_traits$log_BLMax <- log(zoo_traits$BodyLengthMax_mm) 
hist(zoo_traits$log_BLMax)

zoo_traits <- zoo_traits %>% select(-c(meanBodyLength_mm, BodyLengthMax_mm, 
                                       Genus, SpeciesEpithet, References, Notes))
zoo_traits <- zoo_traits %>% select(-c(ReproductionMode, log_BLMax)) #too many NAs#
zoo_traits$FeedingMode <- as.factor(zoo_traits$FeedingMode)
# zoo_traits$ReproductionMode <- as.factor(zoo_traits$ReproductionMode)
zoo_traits$TrophicGroup <- as.factor(zoo_traits$TrophicGroup)

colSums(is.na(zoo_traits))

zoo_traits <- as.data.frame(zoo_traits)              # drop tibble class
rownames(zoo_traits) <- zoo_traits$Species           # assign species as row names
zoo_traits$Species <- NULL                           # remove the Species column

##calculating functional diversity##
dist.traits <- gowdis(zoo_traits) #with gower distance#
dbFD_trait <- dbFD(dist.traits, biomass, corr = "cailliez",
                   calc.FDiv = TRUE, #calc.FGR = TRUE,
                   calc.FRic = TRUE, calc.CWM = FALSE)
save(dbFD_trait, file = 'FD_traits_biomass_260526.RDATA')

setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/spatial_sensitivity_analysis/Smaller_radius/RDATA")
load('FD_traits_biomass_260526.RDATA')
par(mfrow=c(2,2))
hist(dbFD_trait$FEve)
hist(dbFD_trait$FDiv) 
hist(dbFD_trait$FDis)
hist(dbFD_trait$RaoQ)

## combine all info and plot
#siteInfo = filtered_carbon_env_data2[, c(1:7)]

#Divphyto = data.frame(siteInfo, FEve = dbFD_trait$FEve, 
#FDis = dbFD_trait$FDis, RaoQ = dbFD_trait$RaoQ)
Divphyto = data.frame(zoo_FEve = dbFD_trait$FEve,zoo_FDis = dbFD_trait$FDis)

summary(Divphyto)
rownames(Divphyto) <- NULL

setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/spatial_sensitivity_analysis/Smaller_radius/RDATA")
save(Divphyto, file = 'FD_zoo_traits_biomass_260526.RDATA')

#combine with biomass dataset#
setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/spatial_sensitivity_analysis/Smaller_radius/RDATA")
load('FD_zoo_traits_biomass_260526.RDATA')
setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/spatial_sensitivity_analysis/Smaller_radius/RDATA")
load("plankton_RUE_data3_260526.RDATA")
plankton_RUE_data4 <- cbind(plankton_RUE_data3,Divphyto)
save(plankton_RUE_data4, file = 'plankton_RUE_data4_260526.RDATA')

#clean the dataset or remove species info, and remove NAs #.................
load("plankton_RUE_data4_260526.RDATA")
# plankton_RUE_data5 <- plankton_RUE_data4 %>% select(-(c(16:84, 113:152))) #remove species biomass info

colSums(is.na(plankton_RUE_data4)) > 0
sum(is.na(plankton_RUE_data4$zoo_FEve))
plankton_RUE_data4[is.na(plankton_RUE_data4$zoo_FEve), ]
#only two species of zooplankton appear in these two rows (Bosmina and Limnocalanus macrurus), found offshore, aug 2003

#remove rows where zoo FEve are found
plankton_RUE_data5 <- plankton_RUE_data4[!is.na(plankton_RUE_data4$zoo_FEve), ]
save(plankton_RUE_data5, file = 'plankton_RUE_data5_260526.RDATA')

#.............................................................................
#### Plotting Results ###################################################
library(ggplot2)
library(ggpubr)
library(patchwork)
#FD#
FE <- ggplot(data = plankton_RUE_data4, aes(x=longitude, y=latitude,
                                                   color = zoo_FEve)) +
  geom_point(size = 3) +
  scale_color_viridis_c(option="turbo", direction=-1) +
  theme_bw() + ggtitle('Zooplankton Functional Evenness') +
  theme(legend.position="bottom")
FD <- ggplot(data = plankton_RUE_data4, aes(x=longitude, y=latitude,
                                                   color = zoo_FDis)) +
  geom_point(size = 3) +
  scale_color_viridis_c(option="turbo", direction=-1) +
  theme_bw() + ggtitle('Zooplankton Functional Dispersion') +
  theme(legend.position="bottom")

plot = ggarrange(FE, FD, nrow = 1, ncol = 2)

plot

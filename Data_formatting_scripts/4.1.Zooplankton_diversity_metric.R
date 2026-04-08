library(dplyr)
library(readxl)
library(tidyr)
library(vegan)
library(FD)
library(lubridate)
library(stringr)

setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/RDATA")
load('plankton_RUE_data2_260307.RDATA')
load('Zooplankton_REC_clean3_260307.RDATA')

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
Zooplankton_biomass_data2$zoo_Hill1 <- exp(diversity(biomass_data, index = "shannon"))  # q = 1 (Shannon diversity)
Zooplankton_biomass_data2$zoo_Hill2 <- 1 / diversity(biomass_data, index = "invsimpson")  # q = 2 (Simpson diversity)

#Species richness#
# Calculate species richness per sample
Zooplankton_biomass_data2$species_richness <- apply(biomass_data, 1, function(x) sum(x > 0))
Zooplankton_biomass_data2$log_SR <- log(Zooplankton_biomass_data2$species_richness)

save(Zooplankton_biomass_data2, file = 'Zooplankton_biomass_data2_260307.RDATA')

zoo_diversity<-Zooplankton_biomass_data2%>%select(c(SampleEvent, zoo_Hill1, zoo_Hill2, species_richness, log_SR))

#adding species richness data to zooplankton dataset#
plankton_RUE_data3 <- left_join(plankton_RUE_data2, zoo_diversity, by = "SampleEvent")
save(plankton_RUE_data3, file = 'plankton_RUE_data3_260307.RDATA')

################################################################################
######### Functional diversity - functional evenness and disparity #############
rm(list = ls())
setwd("~/GitHub/PlanktonRUE-Project/Data/Zooplankton/Zooplankton_functional_traits")
Species <- read_excel("Zooplankton_species.xlsx")

#extracting traits##------------------------------------------------------------
#extract traits from pata and hunt#....................................................
setwd("~/GitHub/PlanktonRUE-Project/Data/Zooplankton/Zooplankton_functional_traits/Pelagic-Ecosystems-Zooplankton_trait_database-5e78979")
s.format <- read.csv("data_input/trait_dataset_standard_format_20230628.csv")[-1,]
trait.directory <- read.csv("data_input/trait_directory_20230628.csv") %>% 
  distinct(traitID, .keep_all = TRUE)
# taxonomy table
taxonomy <- read.csv("data_input/taxonomy_table_20230628.csv") 
# stage table
lifestagelist <- read.csv("data_input/lifestage_directory_20230628.csv") %>% 
  select(-c(majorgroup, notes))
#traits dataset
traits.lvl2 <- read.csv("data_input/Trait_dataset_level2/trait_dataset_level2-2023-09-14.csv")
# Change assocTemperature column class to allow row binds
traits.lvl2$assocTemperature <- as.character(traits.lvl2$assocTemperature)
traits.lvl2$verbatimTemperature <- as.character(traits.lvl2$verbatimTemperature)

# Select only the column with species names for matching
species.list <- Species %>% 
  select(Species) %>% 
  # Join with taxonomy to get taxonID
  left_join(taxonomy, by = c("Species" = "verbatimScientificName")) %>% 
  filter(!is.na(taxonID)) %>% 
  distinct()

#filter trait database
trait.list <- c("bodyLengthMax", "reproductionMode",
                "trophicGroup", "feedingMode")
trait.subset <- traits.lvl2 %>% 
  filter(traitName %in% trait.list) %>% 
  filter(taxonID %in% species.list$taxonID) %>% 
  left_join(species.list, by = "taxonID") %>% 
  arrange(Species, traitName)

#identify missing traits#
# Find species with missing bodyLength
missing.length <- species.list %>% 
  filter(!(taxonID %in% filter(trait.subset, traitName == "bodyLengthMax")$taxonID))

#generalize missing traits at genus or family level
traits.generalized <- s.format  # empty template
trait <- "bodyLengthMax"
gen.level <- "genus"
for (taxon in missing.length$scientificName){
  gen <- getGroupLevelValue(taxon, trait, gen.level = gen.level, 
                            trait.df = traits.lvl2, 
                            taxonomy.df = taxonomy)
  if(is.data.frame(gen)) {
    gen <- gen %>% mutate(traitValue = as.character(traitValue))
    traits.generalized <- traits.generalized %>% bind_rows(gen)
  }
}
# Add generalized traits to your dataset
traits.lvl2.ed <- bind_rows(trait.subset, traits.generalized)

#finalize dataset
# Final subset of traits for your zooplankton
trait.subset <- traits.lvl2.ed %>% 
  filter(traitName %in% trait.list,
         taxonID %in% species.list$taxonID)
# Check coverage
trait.subset %>% 
  distinct(traitName, taxonID) %>% 
  group_by(traitName) %>% 
  summarise(Nrecords = n(), .groups = "drop") %>% 
  mutate(Perc.sp = Nrecords / nrow(species.list) * 100)

traits.wide <- trait.subset %>%
  select(Species, traitName, traitValue) %>%
  group_by(Species, traitName) %>%
  mutate(row = row_number()) %>%     # distinguish duplicates
  ungroup() %>%
  pivot_wider(
    names_from = traitName,
    values_from = traitValue
  ) %>%
  select(-row)   # remove the helper column if not needed

#find max of bodylengthmax for each genus#
genus_means <- traits.lvl2 %>%
  filter(traitName == "bodyLengthMax") %>%
  mutate(traitValue = as.numeric(traitValue)) %>%   # ensure numeric
  group_by(genus) %>%                             # use the genus column from your taxonomy
  summarise(max_bodyLength = max(traitValue, na.rm = TRUE),
            n_species = n())

#find most common trait for each genus#
trait.list <- c("reproductionMode", "trophicGroup", "feedingMode")
dominant_traits <- traits.lvl2 %>%
  filter(traitName %in% trait.list) %>%              # only categorical traits of interest
  group_by(genus, traitName, traitValue) %>%         # group by genus + trait + value
  summarise(n = n(), .groups = "drop") %>%
  group_by(genus, traitName) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%         # keep most frequent traitValue per trait
  ungroup()

#find most common trait for each family#
dominant_traits <- traits.lvl2 %>%
  filter(traitName %in% trait.list) %>%              # only categorical traits of interest
  group_by(family, traitName, traitValue) %>%         # group by genus + trait + value
  summarise(n = n(), .groups = "drop") %>%
  group_by(family, traitName) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%         # keep most frequent traitValue per trait
  ungroup()

#find most common trait for each order#
dominant_traits <- traits.lvl2 %>%
  filter(traitName %in% trait.list) %>%              # only categorical traits of interest
  group_by(order, traitName, traitValue) %>%         # group by genus + trait + value
  summarise(n = n(), .groups = "drop") %>%
  group_by(order, traitName) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%         # keep most frequent traitValue per trait
  ungroup()

#find most common trait for each phylum#
dominant_traits <- traits.lvl2 %>%
  filter(traitName %in% trait.list) %>%              # only categorical traits of interest
  group_by(phylum, traitName, traitValue) %>%         # group by genus + trait + value
  summarise(n = n(), .groups = "drop") %>%
  group_by(phylum, traitName) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%         # keep most frequent traitValue per trait
  ungroup()

#find most common trait for each major group#
dominant_traits <- traits.lvl2 %>%
  filter(traitName %in% trait.list) %>%              # only categorical traits of interest
  group_by(majorgroup, traitName, traitValue) %>%         # group by genus + trait + value
  summarise(n = n(), .groups = "drop") %>%
  group_by(majorgroup, traitName) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%         # keep most frequent traitValue per trait
  ungroup()

#find max of body length for each phylum
phylum_means <- traits.lvl2 %>%
  filter(traitName == "bodyLengthMax") %>%
  mutate(traitValue = as.numeric(traitValue)) %>%   # ensure numeric
  group_by(phylum) %>%                            
  summarise(max_bodyLength = max(traitValue, na.rm = TRUE),
            n_species = n())

#find max of body length for each order
order_means <- traits.lvl2 %>%
  filter(traitName == "bodyLengthMax") %>%
  mutate(traitValue = as.numeric(traitValue)) %>%   # ensure numeric
  group_by(order) %>%                             
  summarise(max_bodyLength = max(traitValue, na.rm = TRUE),
            n_species = n())

#find max of body length for each family
family_means <- traits.lvl2 %>%
  filter(traitName == "bodyLengthMax") %>%
  mutate(traitValue = as.numeric(traitValue)) %>%   # ensure numeric
  group_by(family) %>%                             
  summarise(max_bodyLength = max(traitValue, na.rm = TRUE),
            n_species = n())
#...............................................................................

#extract traits from Beisner 2021
setwd("~/GitHub/PlanktonRUE-Project/Data/Zooplankton/Zooplankton_functional_traits/Beisner_etal_2021")
Beisner_traits<- read.csv("Functionaltraitsmatrix.csv",sep=";")
Beisner_summary <- Beisner_traits %>%
  group_by(Species.name, feeding.type, trophic.group, habitat) %>%
  summarise(mean_body_length = mean(mean.length_mm, na.rm = TRUE),
            n = n(),   # optional: number of records averaged
            .groups = "drop")
Beisner_genus_summary <- Beisner_traits %>%
  mutate(Genus = word(Species.name, 1, sep = fixed("."))) %>%   # extract genus before the dot
  group_by(Genus) %>%
  summarise(mean_body_length = mean(mean.length_mm, na.rm = TRUE),
            n = n(),  # how many records per genus
            .groups = "drop")

#extract traits from Hebert etal 2016
setwd("~/GitHub/PlanktonRUE-Project/Data/Zooplankton/Zooplankton_functional_traits/Hebert_etal_2016_zooplankton_traits/ecy1337-sup-0001-data")
Hebert_traits<- read.csv("zooplankton_traits.csv",sep=";", dec=",")
Hebert_summary <- Hebert_traits %>%
  mutate(
    Body.length = ifelse(
      is.na(Body.length) & !is.na(Min.bl) & !is.na(Max.bl),
      (Min.bl + Max.bl) / 2,
      Body.length
    )
  ) %>%
  group_by(Genus, Species) %>%
  summarise(
    mean_bodylength = mean(Body.length, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )
#genus level
Hebert_genus_summary <- Hebert_traits %>%
  mutate(
    Body.length = ifelse(
      is.na(Body.length) & !is.na(Min.bl) & !is.na(Max.bl),
      (Min.bl + Max.bl) / 2,
      Body.length
    )
  ) %>%
  group_by(Genus) %>%
  summarise(
    mean_bodylength = mean(Body.length, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )
#group level
Hebert_group_summary <- Hebert_traits %>%
  mutate(
    Body.length = ifelse(
      is.na(Body.length) & !is.na(Min.bl) & !is.na(Max.bl),
      (Min.bl + Max.bl) / 2,
      Body.length
    )
  ) %>%
  group_by(Group) %>%
  summarise(
    mean_bodylength = mean(Body.length, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

Hebert_traits %>%
  filter(!is.na(Trophic.group)) %>%     # drop NA trophic groups
  group_by(Group) %>%
  count(Trophic.group, sort = TRUE) %>%
  slice_max(n, n = 1, with_ties = FALSE)

#.............................................................................

#Functional Diversity Metric Calculations----------------------------------------
rm(list = ls())
setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/RDATA")
load("plankton_RUE_data3_260307.RDATA")
setwd("~/GitHub/PlanktonRUE-Project/Data/Zooplankton/Zooplankton_functional_traits")
Traits <- read_excel("Zooplankton_traits_clean_290925.xlsx")

#turn species name to dataframe and reorder
selected_col_names <- colnames(plankton_RUE_data3)[113:min(152, ncol(plankton_RUE_data3))]
species_table <- data.frame(Species = selected_col_names)
species_table$Species <- sub("_z$", "", species_table$Species)
Traits2 <- Traits %>% 
  filter(Species %in% species_table$Species)%>%
  mutate(Species = factor(Species, levels = species_table$Species)) %>%
  arrange(Species)  # Arrange rows based on the ordered Species

biomass_sub <-plankton_RUE_data3 %>% 
  select(113:152)
colnames(biomass_sub) <- sub("_z$", "", colnames(biomass_sub))
biomass <- biomass_sub %>%
  mutate(Site = paste("Site_", row_number())) %>%
  relocate(Site)
rownames(biomass) <- biomass$Site           
biomass$Site <- NULL                           
biomass <- as.data.frame(biomass)

setwd("~/GitHub/PlanktonRUE-Project/Data/Zooplankton/Zooplankton_functional_traits/RDATA_and_csv")
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
save(dbFD_trait, file = 'FD_traits_biomass_260307.RDATA')

setwd("~/GitHub/PlanktonRUE-Project/Data/Zooplankton/Zooplankton_functional_traits/RDATA_and_csv")
load('FD_traits_biomass_260307.RDATA')
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

setwd("~/GitHub/PlanktonRUE-Project/Data/Zooplankton/Zooplankton_functional_traits/RDATA_and_csv")
save(Divphyto, file = 'FD_zoo_traits_biomass_260307.RDATA')

#combine with biomass dataset#
setwd("~/GitHub/PlanktonRUE-Project/Data/Zooplankton/Zooplankton_functional_traits/RDATA_and_csv")
load('FD_zoo_traits_biomass_260307.RDATA')
setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/RDATA")
load("plankton_RUE_data3_260307.RDATA")
plankton_RUE_data4 <- cbind(plankton_RUE_data3,Divphyto)
save(plankton_RUE_data4, file = 'plankton_RUE_data4_260307.RDATA')

#clean the dataset or remove species info, and remove NAs #.................
load("plankton_RUE_data4_260307.RDATA")
# plankton_RUE_data5 <- plankton_RUE_data4 %>% select(-(c(16:84, 113:152))) #remove species biomass info

colSums(is.na(plankton_RUE_data4)) > 0
sum(is.na(plankton_RUE_data4$zoo_FEve))
plankton_RUE_data4[is.na(plankton_RUE_data4$zoo_FEve), ]
#only two species of zooplankton appear in these two rows (Bosmina and Limnocalanus macrurus), found offshore, aug 2003

#remove rows where zoo FEve are found
plankton_RUE_data5 <- plankton_RUE_data4[!is.na(plankton_RUE_data4$zoo_FEve), ]
save(plankton_RUE_data5, file = 'plankton_RUE_data5_260309.RDATA')

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

setwd("~/GitHub/PlanktonRUE-Project/Figures")
png("Functional_Diversity_zoo.png", width = 8, height = 5, units = "in", res = 500)
plot
dev.off()

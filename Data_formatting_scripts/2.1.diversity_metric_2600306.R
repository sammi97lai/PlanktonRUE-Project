######################### Diversity metric calculations #####################################

library(dplyr)
library(readxl)
library(tidyr)
library(vegan)
library(FD)

setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/RDATA")

############## Hill Diversity, Genus richness and Rarity-weighted Richness #################----------------------------------

load('filtered_biomass_env_data.RDATA')
# load('non_filtered_biomass_data.RDATA') #include rare species

bio_data <- filtered_biomass_env_data[-c(1:17)] #keep only abundance data, remove env factors#

filtered_biomass_env_data2 <- filtered_biomass_env_data #make new data set

# Hill numbers
filtered_biomass_env_data2$Hill1 <- exp(diversity(bio_data, index = "shannon"))  # q = 1 (Shannon diversity)
filtered_biomass_env_data2$Hill2 <- 1 / diversity(bio_data, index = "invsimpson")  # q = 2 (Simpson diversity)

#Genus richness#
# Calculate genus richness per sample
filtered_biomass_env_data2$genus_richness <- apply(bio_data, 1, function(x) sum(x > 0))
filtered_biomass_env_data2$log_G <- log(filtered_biomass_env_data2$genus_richness)

#Rarity-weighted richness#
pa_matrix <- ifelse(bio_data > 0, 1, 0) #covert to presence-absence data
species_freq <- colSums(pa_matrix) #calculate species frequency (number of sites species is present)
rarity_weights <- 1/species_freq
weighted_matrix <- sweep(pa_matrix, 2, rarity_weights, `*`)
rwr <- rowSums(weighted_matrix)
filtered_biomass_env_data2$rwr <- rwr

save(filtered_biomass_env_data2, file = 'filtered_biomass_env_data2.RDATA')

phyto_diversity<-filtered_biomass_env_data2%>%select(c(Hill1,Hill2,log_G,genus_richness,rwr))


#####################################################################################################
###############Functional diversity - functional evenness and disparity ###########################------------------
rm(list = ls())
setwd("~/GitHub/PlanktonRUE-Project/Data/Phytoplankton/Phytoplankton_functional_traits")

trait_klais <- read_excel("Klais_et_al/phytoplankton_trait_compilation.xlsx", sheet = "taxonomy & traits")
trait_Lehtinen <- read_excel("Supplement Table S1_Trait data table_Lehtinen et al.xlsx", sheet = "Trait Data Table")
Edwards_et_al_2015 <- read_excel("Edwards_et_al_2015.xlsx",sheet = "Table1")
trait_Thomas <- read.csv("Thomas_et_al/traits_derived_20140606.csv", comment.char="#")

#selecting certain traits#
trait_Edwards<- Edwards_et_al_2015 %>% 
  select (species, vmax_amm, qmax_amm, vmax_nit, 
          qmax_nit, , vmax_p, qmax_p)
trait_Thomas_sub<-trait_Thomas %>% 
  select (phylum, class, order, family, genus, species, tmin, tmax, emp.max.growth)

#substituting ones without genus with species name or putting a genus column in#
trait_Lehtinen2 <- trait_Lehtinen %>%
  mutate(Genus = coalesce(Genus, Species))
anyNA(trait_Lehtinen2$Genus) #check NAs# #FALSE#
trait_klais2 <- trait_klais %>%
  mutate(Genus = coalesce(Genus, valid_name))
anyNA(trait_klais2$Genus) #check NAs# #FALSE#
trait_Thomas_sub <- trait_Thomas_sub %>%
  filter(genus != "" & !is.na(genus)) #removing blanks in genus#
trait_Edwards2 <- trait_Edwards %>% 
  separate(species, into = c("Genus", "species"), sep = " ", extra = "merge", fill = "right")
anyNA(trait_Edwards2$Genus) #check NAs# #FALSE#

#data organization#
names(trait_klais2)[names(trait_klais2) == "nfix"] <- "Nfixing"
names(trait_klais2)[names(trait_klais2) == "buoyancy"] <- "Buoyancy"
names(trait_klais2)[names(trait_klais2) == "motility"] <- "Motility"
names(trait_Thomas_sub)[names(trait_Thomas_sub) == "genus"] <- "Genus"

#select columns#
trait_Lehtinen2_filter<- trait_Lehtinen2 %>% select(-c(Class, Order, Species, SizeDescription_CountingUnit))
trait_klais2_filter<- trait_klais2 %>% 
  select(Genus, Nfixing, Buoyancy, silica, Motility, autotrophy, heterotrophy, fyco, chlc, chlb, chain)
trait_Edwards2_filter<- trait_Edwards2 %>% 
  select(Genus, vmax_amm, vmax_nit, vmax_p, qmax_amm, 
         qmax_nit, qmax_p)
trait_Thomas_filter<- trait_Thomas_sub %>%
  select(Genus, tmin, tmax, emp.max.growth)

#because of similar genuses, had to keep binary values that was majority#
majority <- function(x) {
  if (length(x) == 0) return(NA)  # Return NA if no valid values
  table_x <- table(na.omit(x))  # Create a frequency table excluding NAs
  if (length(table_x) == 0) return(NA)  # If only NAs, return NA
  names(table_x)[which.max(table_x)]  # Return the most frequent value
}

trait_klais2_summary <- trait_klais2_filter %>%
  group_by(Genus) %>%
  summarise(across(c(Nfixing, Buoyancy, silica, Motility, autotrophy, heterotrophy, fyco, chlc, chlb, chain), 
                   majority, .names = "{.col}"))
print(trait_klais2_summary)

trait_Lehtinen2_summary <- trait_Lehtinen2_filter %>%
  group_by(Genus) %>%
  summarise(across(c(CountingUnit, Nfixing, Buoyancy, Motility, Mixotrophy, Harmfulness, SizeGroup), 
                   majority, .names = "{.col}"),
            across(c(`Biovolume_µm-3`,`ESD(Equivalent Shperical Diameter)`), mean, .names = "{.col}"))
print(trait_Lehtinen2_summary)

#prop.table(table(trait_klais2_filter$Genus, trait_klais2_filter$autotrophy), margin = 1) #datachecking#
#most of the functions are similar to majority of 1 that most of the species within the same genus has the ability#

# Custom function to convert non-numeric values to NA and replace NaN with NA
convert_to_na <- function(x) {
  result <- as.numeric(as.character(x))  # First convert to character, then to numeric
  result[is.nan(result)] <- NA  # Replace NaN with NA
  return(result)
}

trait_Edwards2_summary <- trait_Edwards2_filter %>%
  mutate(across(-Genus, as.numeric))%>%
  group_by(Genus)%>% 
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE), .names = "{.col}"))%>%
  mutate(across(everything(), ~ifelse(is.nan(.), NA, .)))  # Replace any remaining NaN with NA after summarizing
print(trait_Edwards2_summary)

trait_Thomas_summary <- trait_Thomas_filter %>%
  group_by(Genus) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))
print(trait_Thomas_summary)

#merge datasets#
phyto_traits <- full_join(trait_klais2_summary, trait_Lehtinen2_summary, by = "Genus", suffix=c("_K","_L"))%>%
  mutate(across(ends_with("_K"), 
                ~ coalesce(.x, get(sub("_K", "_L", cur_column()))), 
                .names = "{.col}")) %>%
  rename_with(~ sub("_K", "", .), ends_with("_K")) %>%
  select(-ends_with("_L"))  # Drop unnecessary _df2 columns
phyto_traits <- full_join(phyto_traits, trait_Edwards2_summary, by = "Genus")
phyto_traits <- full_join(phyto_traits, trait_Thomas_summary, by = "Genus")

#sort genera together like dataset#
#remove certain Genus
filtered_phyto_traits <- phyto_traits %>% 
  filter(!Genus %in% c("Katablepharis", "Leucocryptos", "Telonema",
                       "Choanoflagellatea", "Diaphanoeca", "Cysts")) %>%
  select(-heterotrophy) ##removing heterotrophs and cysts#

#group Genus and relabel and average traits#....................................
filtered_phyto_traits <- filtered_phyto_traits %>%
  mutate(Genus_grouped = case_when(
    Genus %in% c("Monoraphidium", "Koliella") ~ 
      "Monoraphidium complex", #inconsistencies in ID#
    Genus %in% c("Oscillatoriales", "Pseudanabaena", "Planktolyngbya") ~ 
      "Oscillatoriales/Synechococcales",#group all non-nitrogen-fixing cyanobacteria#
    Genus %in% c("Paraphysomonas", "Ochromonadales cf.", "Chromulina", 
                 "Chrysolykos", "Chrysidiastrum", "Chrysosphaerella",
                 "Chrysamoeba", "Monochrysis", "Phaeaster", "Epipyxis", "Bitrichia") ~ 
      "Ochromonadales", #everything under order Ochromonadales#
    Genus %in% c("Plagioselmis", "Cryptomonas", "Teleaulax","Hemiselmis",
                 "Rhodomonas", "Chroomonas", "Cryptomonadales cf.") ~ 
      "Cryptomonadales", #genus level ID not always possible#
    Genus %in% c("Flagellates oval", "Flagellates round", "Unicell", "Monad",
                 "Flagellates oval cf.", "Flagellates (oval)", "Flagellates (sphere)",
                 "Flagellate biflagella", "Flagellate uniflagella") ~ 
      "Unidentified flagellates unicell and monads", #usually unidentifiable#
    Genus %in% c("Anabaena", "Dolichospermum") ~ 
      "Anabaena/Dolichsopsermum", #name change#
    Genus %in% c("Chroococcales", "Cyanodictyon", "Merismopedia", "Lemmermanniella",
                 "Aphanocapsa", "Romeria",
                 "Synechococcus", "Chroococcales cf.") ~ 
      "Chroococcales/Synechococcales", #genus level ID not always possible"
    Genus %in% c("Centrales", "Cyclotella", "Eupodiscales") ~ "Centrales spp.", 
    Genus %in% c("Eutreptiella", "Trachelomonas", "Lepocinclis",
                 "Eutreptia", "Euglena", "Phacus", "Monomorphina", "Strombomonas",
                 "Colacium", "Euglenales cf.", "Cryptoglena") ~ 
      "Euglenales", #everything under order Euglenales#
    Genus %in% c("Dinophyceae cysts", "Dinophyceae") ~ "Unidentified Dinophyceae", 
    TRUE ~ Genus  # Keep others unchanged
  ))

# Select only numeric columns for averaging
numeric_cols <- c("Biovolume_µm-3", "ESD(Equivalent Shperical Diameter)", "vmax_amm", "vmax_nit",
                  "vmax_p", "qmax_amm", "qmax_nit", "qmax_p", "tmin", "tmax", "emp.max.growth")

# Group by your renamed group and calculate means
grouped_averages <- filtered_phyto_traits %>%
  group_by(Genus_grouped) %>%
  summarise(across(all_of(numeric_cols), ~ mean(.x, na.rm = TRUE))) %>%
  mutate(across(where(is.numeric), ~replace(.x, is.nan(.x), NA)))#replace NaN with NA

get_mode <- function(x) {
  x <- x[!is.na(x)]  # remove NAs
  if(length(x) == 0) return(NA)  # if all were NA, return NA
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
} #get most frequent value in group for categorical values#

categorical_cols <- c("Nfixing", "Buoyancy", "silica", "Motility", "autotrophy", "fyco",
                      "chlc", "chlb", "chain", "CountingUnit", "Mixotrophy", "Harmfulness", "SizeGroup")

categorical_summary <- filtered_phyto_traits %>%
  group_by(Genus_grouped) %>%
  summarise(across(all_of(categorical_cols), get_mode))

final_phyto_traits <- left_join(categorical_summary, grouped_averages, by = "Genus_grouped")
colnames(final_phyto_traits)[colnames(final_phyto_traits) == "Genus_grouped"] <- "Genus" #rename to Genus

#match with biomass genus data#..............................................
setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/RDATA")
load("filtered_biomass_env_data2.RDATA")
selected_col_names <- colnames(filtered_biomass_env_data2)[18:min(91, ncol(filtered_biomass_env_data2))]
genus_table <- data.frame(Genus = selected_col_names)
final_phyto_traits2 <- final_phyto_traits %>%
  filter(Genus %in% genus_table$Genus)%>%
  mutate(Genus = factor(Genus, levels = genus_table$Genus)) %>%
  arrange(Genus)  # Arrange rows based on the ordered Genus

setwd("~/GitHub/PlanktonRUE-Project/Data/Phytoplankton/Phytoplankton_functional_traits/RDATA_and_csv")
save(final_phyto_traits2, file = 'final_phyto_traits_biomass.RDATA')

#calculating functional diversity#..........................................
rm(list = ls())
setwd("~/GitHub/PlanktonRUE-Project/Data/Phytoplankton/Phytoplankton_functional_traits/RDATA_and_csv")
load('final_phyto_traits_biomass.RDATA')

setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/RDATA")
load("filtered_biomass_env_data2.RDATA")

setwd("~/GitHub/PlanktonRUE-Project/Data/Phytoplankton/Phytoplankton_functional_traits/RDATA_and_csv")

biomass_sub <-filtered_biomass_env_data2 %>%
  select(-c(Site_Id, Sampling_Id, SettledSample_Id, Date, Time,
            latitude, longitude, MinDepth, MaxDepth,
            `Chlorophyll-a`, Temperature, pH, Salinity, TotalPhosphorus_unfiltered,
            TotalNitrogen_unfiltered, Oxygen_soluble, Hill1, Hill2,
            log_G, genus_richness, rwr, Region))
biomass <- biomass_sub %>%
  mutate(Site = paste("Site_", row_number())) %>%
  relocate(Site)
biomass <- as.data.frame(biomass)
write.csv(biomass, "biomass.csv",row.names = FALSE)

# Order species labels alphabetically
traits <- final_phyto_traits2[order(final_phyto_traits2$Genus), ]

#log transform biovolume#
hist(traits$`Biovolume_µm-3`)
traits$log_Biovolume <- log(traits$`Biovolume_µm-3`) 
hist(traits$log_Biovolume)
traits <- traits %>% select(-c(`Biovolume_µm-3`))

#make sure trait data is in right format
# str(traits)

write.csv(traits, "Traits_260306.csv",row.names = FALSE)
colSums(is.na(traits))

#filter out traits#
new_traits <- traits %>%
  select(Genus,Nfixing, Buoyancy, Motility, chain, Mixotrophy, SizeGroup, log_Biovolume,
         vmax_nit, vmax_p,
         silica, fyco, chlc, chlb)

write.csv(new_traits, "new_Traits_260306.csv",row.names = FALSE)

#create PCA of functional traits#
rm(list = ls())
setwd("~/GitHub/PlanktonRUE-Project/Data/Phytoplankton/Phytoplankton_functional_traits/RDATA_and_csv")
Traits <- read.csv("new_Traits_260306.csv", header = T, row.names = 1)

Traits_clean <- Traits #replace NA with column mean (just for PCA)
Traits_clean[] <- lapply(Traits_clean, function(x) {
  if (is.numeric(x)) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
  }
  return(x)
})
trait_scaled <- scale(Traits_clean)

pca_result <- prcomp(trait_scaled, center = TRUE, scale. = TRUE)

library(factoextra)
# Biplot
fviz_pca_biplot(pca_result,
                repel = TRUE, # avoids label overlap
                label = "var", 
                col.var = "black", # trait arrows
                col.ind = "steelblue",
                pointshape = 21,
                addEllipses = FALSE) # points

#Traits_clean["Coscinodiscus", ] #large Diatom, very high in Biovolume

##calculating functional diversity##
setwd("~/GitHub/PlanktonRUE-Project/Data/Phytoplankton/Phytoplankton_functional_traits/RDATA_and_csv")
Biomass <- read.csv("biomass.csv", header = T, row.names = 1, check.names=FALSE)
Traits <- read.csv("new_Traits_260306.csv", header = T, row.names = 1)

# Check which rows sum to zero
zero_rows <- which(rowSums(Biomass, na.rm = TRUE) == 0)
zero_rows #1555

dist.traits <- gowdis(Traits) #with gower distance#

dbFD_trait <- dbFD(dist.traits, Biomass, corr = "cailliez",
                   calc.FDiv = TRUE, #calc.FGR = TRUE,
                   calc.FRic = TRUE, calc.CWM = FALSE)
save(dbFD_trait, file = 'FD_traits_biomass_260306.RDATA')

setwd("~/GitHub/PlanktonRUE-Project/Data/Phytoplankton/Phytoplankton_functional_traits/RDATA_and_csv")
load('FD_traits_biomass_260306.RDATA')
par(mfrow=c(2,2))
hist(dbFD_trait$FEve)
hist(dbFD_trait$FDiv) 
hist(dbFD_trait$FDis)
hist(dbFD_trait$RaoQ)

## combine all info and plot
#siteInfo = filtered_carbon_env_data2[, c(1:7)]

#Divphyto = data.frame(siteInfo, FEve = dbFD_trait$FEve, 
#FDis = dbFD_trait$FDis, RaoQ = dbFD_trait$RaoQ)
Divphyto = data.frame(FEve = dbFD_trait$FEve,FDis = dbFD_trait$FDis)

summary(Divphyto)
rownames(Divphyto) <- NULL


#combine with rest of phyto dataset#
setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/RDATA")
load('filtered_biomass_env_data2.RDATA')
filtered_biomass_env_data3 <- cbind(filtered_biomass_env_data2,Divphyto)
save(filtered_biomass_env_data3, file = 'filtered_biomass_env_data3_260306.RDATA')

#...............................................................................
####### Functional Diversity- Selected traits for Predation #############
setwd("~/GitHub/PlanktonRUE-Project/Data/Phytoplankton/Phytoplankton_functional_traits/RDATA_and_csv")
Biomass <- read.csv("biomass.csv", header = T, row.names = 1, check.names=FALSE)
Traits <- read.csv("new_Traits_260306.csv", header = T, row.names = 1)
colnames(Traits)

#selected traits#
selected_trait <- Traits %>%
  select(Buoyancy, silica, Motility, log_Biovolume, chain)
# str(selected_trait)
# summary(selected_trait)
selected_trait$Buoyancy <- as.numeric(selected_trait$Buoyancy)
selected_trait$silica <- as.numeric(selected_trait$silica)
selected_trait$Motility <- as.numeric(selected_trait$Motility)
selected_trait$chain <- as.numeric(selected_trait$chain)

dist.selectedtraits <- gowdis(selected_trait) #with gower distance#

#calculate functional diversity metrics#
#for biomass
dbFD_trait2 <- dbFD(dist.selectedtraits, Biomass, corr = "cailliez",
                   calc.FDiv = TRUE, #calc.FGR = TRUE,
                   calc.FRic = TRUE, calc.CWM = FALSE)
save(dbFD_trait2, file = 'FD_traits_selected_biomass_260306.RDATA')

setwd("~/GitHub/PlanktonRUE-Project/Data/Phytoplankton/Phytoplankton_functional_traits/RDATA_and_csv")
load('FD_traits_selected_biomass_260306.RDATA')
par(mfrow=c(2,2))
hist(dbFD_trait2$FEve)
hist(dbFD_trait2$FDiv) 
hist(dbFD_trait2$FDis)
hist(dbFD_trait2$RaoQ)

Divphyto2 = data.frame(FEve_selected = dbFD_trait2$FEve,FDis_selected = dbFD_trait2$FDis)

summary(Divphyto2)
rownames(Divphyto2) <- NULL

#combine with biomass dataset#
setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/RDATA")
load('filtered_biomass_env_data3_260306.RDATA')
filtered_biomass_env_data4 <- cbind(filtered_biomass_env_data3,Divphyto2)
save(filtered_biomass_env_data4, file = 'filtered_biomass_env_data4_260306.RDATA')

#...............................................................................
#################### Cyanobacteria Dominance ###################################
setwd("~/GitHub/PlanktonRUE-Project/Data/Phytoplankton/Phytoplankton_functional_traits/RDATA_and_csv")
Biomass <- read.csv("biomass.csv", header = T, row.names = 1, check.names=FALSE)

#identify cyanobacteria species#
Cyano <- c("Anabaena/Dolichsopsermum", "Aphanizomenon", "Chroococcales/Synechococcales",
           "Nodularia", "Snowella", "Woronichinia", "Eucapsis", "Aphanothece",
           "Planktothrix")

Cyano_Dom <- Biomass %>%
  rowwise() %>%
  mutate(
    total = sum(c_across(everything())),               # total across all genus in the row
    cyano_sum = sum(c_across(all_of(Cyano))),        # sum of cyanobacteria in the row
    cyano_prop = (cyano_sum / total)
  ) %>%
  ungroup() %>%
  select(cyano_prop)  # optional: remove helper columns

setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/RDATA")
load('filtered_biomass_env_data4_260306.RDATA')
filtered_biomass_env_data5 <- cbind(filtered_biomass_env_data4,Cyano_Dom)
save(filtered_biomass_env_data5, file = 'filtered_biomass_env_data5_260306.RDATA')

#...............................................................................
####################plotting results############################################
library(ggplot2)
library(ggpubr)
library(patchwork)
setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/RDATA")
load('filtered_biomass_env_data5_260306.RDATA')

#FD#
FE <- ggplot(data = filtered_biomass_env_data5, aes(x=longitude, y=latitude, 
                                                   color = FEve)) + 
  geom_point(size = 1.2) +
  scale_color_viridis_c(option="turbo", direction=-1) +
  theme_bw() + ggtitle('Functional Evenness') +
  theme(legend.position="bottom")
FD <- ggplot(data = filtered_biomass_env_data5, aes(x=longitude, y=latitude,
                                                   color = FDis)) + 
  geom_point(size = 1.2) +
  scale_color_viridis_c(option="turbo", direction=-1) +
  theme_bw() + ggtitle('Functional Dispersion') +
  theme(legend.position="bottom")

#CWM#
Cyano <- ggplot(data = filtered_biomass_env_data5, aes(x=longitude, y=latitude,
                                                      color = cyano_prop)) + 
  geom_point(size = 1.2) +
  scale_color_viridis_c(option="turbo",direction=-1) +
  theme_bw() + ggtitle('Cyanobacteria Dominance') +
  theme(legend.position="bottom")

#FD#
FE_selected <- ggplot(data = filtered_biomass_env_data5, aes(x=longitude, y=latitude, 
                                                    color = FEve_selected)) + 
  geom_point(size = 1.2) +
  scale_color_viridis_c(option="turbo", direction=-1) +
  theme_bw() + ggtitle('Functional Evenness: selected') +
  theme(legend.position="bottom")
FD_selected <- ggplot(data = filtered_biomass_env_data5, aes(x=longitude, y=latitude,
                                                    color = FDis_selected)) + 
  geom_point(size = 1.2) +
  scale_color_viridis_c(option="turbo", direction=-1) +
  theme_bw() + ggtitle('Functional Dispersion: selected') +
  theme(legend.position="bottom")

plot = ggarrange(FE, FE_selected, FD_selected, FD, Cyano, nrow = 2, ncol = 3)

setwd("~/GitHub/PlanktonRUE-Project/Figures")
png("Functional_Diversity_phyto.png", width = 12, height = 8, units = "in", res = 500)
plot
dev.off()

#with patchwork
(FE | FE_selected | Cyano) /
  (FD | FD_selected )


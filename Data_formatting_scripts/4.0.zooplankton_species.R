########################## Zooplankton Species #################################
library(dplyr)
library(readxl)
library(tidyr)
library(lubridate)
library(writexl)

setwd("~/GitHub/Phytoplankton-project/Data/Zooplankton/Zooplankton_data_Syke")
Zooplankton_REC_2 <- read.csv("ZooplanktonMarine_V3EXT_11Sep2025.csv")
setwd("~/GitHub/Phytoplankton-project/Data/Zooplankton/Zooplankton_functional_traits")

#get species list#..............................................................
#clean dataset
Zooplankton_REC_clean <- Zooplankton_REC_2 %>%
  filter(Year >= 2000 & Year <= 2024) %>%
  mutate(Date = make_date(Year, Month, Day)) %>%
  rename(latitude = Latitude,    
         longitude = Longitude,   
         Zooplankton_MinDepth = MinDepth,
         Zooplankton_MaxDepth = MaxDepth,
         WetWeight = WetWeight_mg.m3
  ) %>%
  filter(!is.na(latitude) & !is.na(longitude))%>%
  mutate(
    EnvironmentType = recode(
      EnvironmentType,
      "Coast, outer archipelago" = "Coastal",
      "Coast, inner archipelago" = "Coastal",
      "Open sea" = "Open water"
    ))

# list unique species names
Species <- Zooplankton_REC_clean %>%
  select(Genus, SpeciesEpithet, Species) %>%
  distinct() %>%
  arrange(Genus, SpeciesEpithet)

write_xlsx(Species, "Zooplankton_species.xlsx")

#remove not "true" zooplankton "free living"
remove_taxa <- c(
  # sessile/benthic invertebrates
  "Amphibalanus", "Argulus", "Balanus", "Bylgides", "Einhornia", "Gammarus", 
  "Harmothoe", "Marenzelleria", "Monoporeia", "Mytilus", "Neomysis", 
  "Palaemon", "Pontoporeia", "Pygospio", "Rhithropanopeus", 
  "Saduria", "Spionidae",
  # sessile ciliates
  "Acineta", "Cothurnia", "Vorticella", "Zoothamnium",
  # higher-level placeholders (not real zooplankton)
  "Actinopteri", "Animalia", "Bivalvia", "Bryozoa", 
  "Chironomidae", "Ciliophora", "Cirripedia", 
  "Decapoda", "Gastropoda",
  "Harpacticoida", "Nematoda", "Nemertea", "Oligochaeta",
  "Polychaeta", "Sessilida"
)

Zooplankton_REC_clean2 <- Zooplankton_REC_clean %>%
  filter(!Genus %in% remove_taxa,
         !Species %in% remove_taxa) 

#check species list again#
# Species <- Zooplankton_REC_clean2 %>%
#   select(Genus, SpeciesEpithet, Species) %>%
#   distinct() %>%
#   arrange(Genus, SpeciesEpithet)

write_xlsx(Species, "Zooplankton_species.xlsx")

setwd("~/GitHub/Phytoplankton-project/phyto_RUE/dataformatting")
save(Zooplankton_REC_clean2, file = 'Zooplankton_REC_clean2_230925.RDATA')

#keep only crustacean zooplankton#
load('Zooplankton_REC_clean2_230925.RDATA')
remove_taxa2 <- c(
  #rotifers
  "Asplanchna", "Brachionus", "Collotheca", "Conochilus", "Euchlanis", "Filinia", 
  "Kellicottia", "Keratella", "Lecane", "Notholca", "Polyarthra", "Synchaeta", 
  "Trichocerca", "Trichotria", "Rotifera",
  #ciliates
  "Askenasia", "Didinium", "Helicostomella", "Leprotintinnus", "Strombidium", "Tintinnopsis",
  #protozoa
  "Radiosperma", "Protozoa",
  #missing traits
  "Cercopagis", "Pleopis",
  #remove species that are 0s
  "Hyperia", "Paracalanus", "Bythotrephes", "Caligus", "Oikopleura", "Mertensia",
  "Podonidae", "Cnidaria", "Mnemiopsis", "Eudiaptomus"
)

Zooplankton_REC_clean3 <- Zooplankton_REC_clean2 %>%
  filter(!Genus %in% remove_taxa2,
         !Species %in% remove_taxa2) 

Species <- Zooplankton_REC_clean3 %>%
  select(Genus, SpeciesEpithet, Species) %>%
  distinct() %>%
  arrange(Genus, SpeciesEpithet)

setwd("~/GitHub/Phytoplankton-project/phyto_RUE/dataformatting")
save(Zooplankton_REC_clean3, file = 'Zooplankton_REC_clean3_260925.RDATA')

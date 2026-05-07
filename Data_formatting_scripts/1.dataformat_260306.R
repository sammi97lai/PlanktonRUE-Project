library(dplyr)
library(readxl)
library(tidyr)
library(vegan)
library(openxlsx) # Save Excel file
library(lubridate)
library(ggplot2)
library(sf)         # For spatial features
library(ggspatial)  # For better map styling
library(viridis)    # For color scale

setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts")

################# Data Organisation and formatting ###############################

##Offshore##---------------
##phytoplankton###
Finland_offshore_phyto_since2000 <- 
  read_excel("~/GitHub/PlanktonRUE-Project/Data/Phytoplankton/Phytoplankton_data_Syke/Finland_offshore_phyto_since2000.xlsx")
#water quality#
Finland_offshore_wq_since2000_from_Vesla_database <- 
  read_excel("~/GitHub/PlanktonRUE-Project/Data/Water quality/Finland_offshore_wq_since2000_from Vesla database.xlsx")
Finland_offshore_wq_since2000_from_Vesla_database$Time <- format(as.POSIXct(Finland_offshore_wq_since2000_from_Vesla_database$Time), format = "%H:%M")
Finland_offshore_wq_since2000_from_Vesla_database<-Finland_offshore_wq_since2000_from_Vesla_database%>%
  rename(`Coordinates_ EUREF_WGS84_lat`=`Coordinates_WGS84 lat`,
         `Coordinates_ EUREF_WGS84_long`=`Coordinates_WGS84 long`)

#merge wq dataset with phyto#
Offshore_phyto<-semi_join(Finland_offshore_phyto_since2000,Finland_offshore_wq_since2000_from_Vesla_database,
                          by=c("Date","Time","Coordinates_ EUREF_WGS84_lat","Coordinates_ EUREF_WGS84_long")) 
Offshore_wq<-semi_join(Finland_offshore_wq_since2000_from_Vesla_database,Finland_offshore_phyto_since2000,
                       by=c("Date","Time","Coordinates_ EUREF_WGS84_lat","Coordinates_ EUREF_WGS84_long")) 
merged_offshore <- Offshore_wq %>%
  left_join(Offshore_phyto, by = c("Station","Site_Id","Date","Time","Coordinates_ EUREF_WGS84_lat","Coordinates_ EUREF_WGS84_long"),
            relationship = "many-to-many")
merged_offshore <- merged_offshore %>% ungroup()%>% mutate(Region = "Offshore")
merged_offshore$`Coordinates_ EUREF_WGS84_lat` <- as.numeric(merged_offshore$`Coordinates_ EUREF_WGS84_lat`) #change lat and lon to numeric#
merged_offshore$`Coordinates_ EUREF_WGS84_long` <- as.numeric(merged_offshore$`Coordinates_ EUREF_WGS84_long`)
n_distinct(merged_offshore$Site_Id)

##Outer coastal##------------------
#water quality#
# Load the data file
Finland_outer_coastal_wq_since20001 <- 
  read_excel("~/GitHub/PlanktonRUE-Project/Data/Water quality/Finland_outer coastal_wq_since2000.xlsx",sheet = "Outer coastal wq part 1")
Finland_outer_coastal_wq_since20002 <- 
  read_excel("~/GitHub/PlanktonRUE-Project/Data/Water quality/Finland_outer_coastal_wq_since20002.xlsx")
Finland_outer_coastal_wq_since20003 <- 
  read_excel("~/GitHub/PlanktonRUE-Project/Data/Water quality/Finland_outer_coastal_wq_since20003.xlsx")

#combine water quality datasets#
Outer_coastal_wq <- bind_rows(Finland_outer_coastal_wq_since20001, Finland_outer_coastal_wq_since20002)
Outer_coastal_wq <- bind_rows(Outer_coastal_wq, Finland_outer_coastal_wq_since20003)
Outer_coastal_wq$Date <- as.Date(Outer_coastal_wq$Date, format = "%d.%m.%Y")

#phytoplankton#
Finland_outer_coastal_phyto_since2000 <- 
  read_excel("~/GitHub/PlanktonRUE-Project/Data/Phytoplankton/Phytoplankton_data_Syke/Finland_outer coastal_phyto_since2000.xlsx")
Finland_outer_coastal_phyto_since2000$Date <- as.POSIXct(Finland_outer_coastal_phyto_since2000$Date, format = "%d.%m.%Y", tz = "UTC")
Finland_outer_coastal_phyto_since2000 <- Finland_outer_coastal_phyto_since2000 %>%
  separate(`Coordinates_ EUREF-FIN/WGS85`, into = c("Coordinates_ EUREF_WGS84_lat", "Coordinates_ EUREF_WGS84_long"), sep = "-")

#merge datasets#
Outer_coastal_phyto<-semi_join(Finland_outer_coastal_phyto_since2000,Outer_coastal_wq,
                               by=c("Site_Id","Sampling_Id","SettledSample_Id","Date","Time")) 
OUTCOAST_wq_filtered<-semi_join(Outer_coastal_wq,Finland_outer_coastal_phyto_since2000,
                                by=c("Site_Id","Sampling_Id","SettledSample_Id","Date","Time"))
#merge
merged_OUTCOAST <- OUTCOAST_wq_filtered %>%
  left_join(Outer_coastal_phyto, by = c("Site_Id","Sampling_Id","SettledSample_Id","Date","Time"))
merged_OUTCOAST <- merged_OUTCOAST %>% ungroup()%>% mutate(Region = "Outer coastal")
merged_OUTCOAST$`Coordinates_ EUREF_WGS84_lat` <- as.numeric(merged_OUTCOAST$`Coordinates_ EUREF_WGS84_lat`) #change lat and lon to numeric#
merged_OUTCOAST$`Coordinates_ EUREF_WGS84_long` <- as.numeric(merged_OUTCOAST$`Coordinates_ EUREF_WGS84_long`)
n_distinct(merged_OUTCOAST$Site_Id) #find number of sites#

##Inner coastal##-------------------------
#water quality#
# Load the data file
Finland_inner_coastal_wq_since20001 <- 
  read_excel("~/GitHub/PlanktonRUE-Project/Data/Water quality/Finland_Inner coastal_wq_since2000.xlsx",sheet = "Inner coastal wq part 1")
Finland_inner_coastal_wq_since20002 <- 
  read_excel("~/GitHub/PlanktonRUE-Project/Data/Water quality/Finland_Inner coastal_wq_since2000.xlsx",sheet = "Inner coastal wq part 2")
Finland_inner_coastal_wq_since20003 <- 
  read_excel("~/GitHub/PlanktonRUE-Project/Data/Water quality/Finland_Inner coastal_wq_since2000.xlsx", sheet = "Inner coastal wq part 3")
Finland_inner_coastal_wq_since20004 <- 
  read_excel("~/GitHub/PlanktonRUE-Project/Data/Water quality/Finland_Inner coastal_wq_since2000.xlsx", sheet = "Inner coastal wq part 4")
Finland_inner_coastal_wq_since20005 <- 
  read_excel("~/GitHub/PlanktonRUE-Project/Data/Water quality/Finland_Inner coastal_wq_since2000.xlsx", sheet = "Inner coastal wq part 5")
Finland_inner_coastal_wq_since20006 <- 
  read_excel("~/GitHub/PlanktonRUE-Project/Data/Water quality/Finland_Inner coastal_wq_since2000.xlsx", sheet = "Inner coastal wq part 6")
# Load the mapping table
mapping <- read_excel("~/GitHub/PlanktonRUE-Project/Data/Water quality/Finland_Inner coastal_wq_since2000_mapping.xlsx")
# Convert mapping table to named vector
name_mapping <- setNames(mapping$Translation, mapping$OriginalHeadline)
# Replace column names
colnames(Finland_inner_coastal_wq_since20001) <- 
  recode(colnames(Finland_inner_coastal_wq_since20001), !!!name_mapping)
colnames(Finland_inner_coastal_wq_since20002) <- 
  recode(colnames(Finland_inner_coastal_wq_since20002), !!!name_mapping)
colnames(Finland_inner_coastal_wq_since20003) <- 
  recode(colnames(Finland_inner_coastal_wq_since20003), !!!name_mapping)
colnames(Finland_inner_coastal_wq_since20004) <- 
  recode(colnames(Finland_inner_coastal_wq_since20004), !!!name_mapping)
colnames(Finland_inner_coastal_wq_since20005) <- 
  recode(colnames(Finland_inner_coastal_wq_since20005), !!!name_mapping)
colnames(Finland_inner_coastal_wq_since20006) <- 
  recode(colnames(Finland_inner_coastal_wq_since20006), !!!name_mapping)
#combine water quality datasets#
Inner_coastal_wq <- bind_rows(Finland_inner_coastal_wq_since20001, Finland_inner_coastal_wq_since20002)
Inner_coastal_wq <- bind_rows(Inner_coastal_wq, Finland_inner_coastal_wq_since20003)
Inner_coastal_wq <- bind_rows(Inner_coastal_wq, Finland_inner_coastal_wq_since20004)
Inner_coastal_wq <- bind_rows(Inner_coastal_wq, Finland_inner_coastal_wq_since20005)
Inner_coastal_wq <- bind_rows(Inner_coastal_wq, Finland_inner_coastal_wq_since20006)
Inner_coastal_wq$Date <- as.Date(Inner_coastal_wq$Date, format = "%d.%m.%Y")

#phytoplankton#
Finland_inner_coastal_phyto_since2000 <- 
  read_excel("~/GitHub/PlanktonRUE-Project/Data/Phytoplankton/Phytoplankton_data_Syke/Finland_inner coastal_phyto_since2000.xlsx")
Finland_inner_coastal_phyto_since2000$Date <- as.Date(Finland_inner_coastal_phyto_since2000$Date, format = "%d.%m.%Y")
Finland_inner_coastal_phyto_since2000 <- Finland_inner_coastal_phyto_since2000 %>%
  separate(`Coordinates_ EUREF-FIN/WGS85`, into = c("Coordinates_ EUREF_WGS84_lat", "Coordinates_ EUREF_WGS84_long"), sep = "-")

#merge datasets#
Inner_coastal_phyto<-semi_join(Finland_inner_coastal_phyto_since2000,Inner_coastal_wq,
                               by=c("Site_Id","Sampling_Id","SettledSample_Id","Date","Time")) 
Inner_wq_filtered<-semi_join(Inner_coastal_wq,Finland_inner_coastal_phyto_since2000,
                             by=c("Site_Id","Sampling_Id","SettledSample_Id","Date","Time"))
#merge
merged_inner_coastal <- Inner_wq_filtered %>%
  left_join(Inner_coastal_phyto, by = c("Site_Id","Sampling_Id","SettledSample_Id","Date","Time"))
merged_inner_coastal <- merged_inner_coastal %>% mutate(Region = "Inner coastal")
merged_inner_coastal$`Coordinates_ EUREF_WGS84_lat` <- as.numeric(merged_inner_coastal$`Coordinates_ EUREF_WGS84_lat`) #change lat and lon to numeric#
merged_inner_coastal$`Coordinates_ EUREF_WGS84_long` <- as.numeric(merged_inner_coastal$`Coordinates_ EUREF_WGS84_long`)
n_distinct(merged_inner_coastal$Site_Id) #find number of sites#

########################### merge datasets  #######################################
#merge all the phytoplankton datasets together#
Phytoplankton_data <- bind_rows(merged_offshore, merged_OUTCOAST)
Phytoplankton_data <- bind_rows(Phytoplankton_data, merged_inner_coastal)

setwd("~/GitHub/PlanktonRUE-Project/Data_formatting_scripts/RDATA")
save(Phytoplankton_data, file = 'Phytoplankton_data_format.RDATA')
load('Phytoplankton_data_format.RDATA')
sum(is.na(Phytoplankton_data$AmmoniumNitrogen_unfiltered))
#data selection and formatting#
Phytoplankton_data2<- Phytoplankton_data%>%
  select(Site_Id, Sampling_Id, SettledSample_Id, Station, Date, Time,`Coordinates_ EUREF_WGS84_lat`,
         `Coordinates_ EUREF_WGS84_long`, StationDepth, Region,
         Temperature, Salinity, pH, Oxygen_soluble, TotalPhosphorus_unfiltered, TotalNitrogen_unfiltered, 
         `Chlorophyll-a`, Class, Order, Genus, Species, Trophy,
         Biovolume_um3, UnitsPerLiter, Biomass_ugPerLiter, Carbon_ugPerLiter, SamplingDepth) #subsetting data#
Phytoplankton_data2<-Phytoplankton_data2%>%
  rename(`latitude`=`Coordinates_ EUREF_WGS84_lat`,
         `longitude`=`Coordinates_ EUREF_WGS84_long`)
Phytoplankton_data2 <- Phytoplankton_data2[Phytoplankton_data2$SamplingDepth != "Syvyysväli", ] # Drop rows 
#Phytoplankton_data2$Depth <- as.numeric(sub(".*-", "", Phytoplankton_data2$SamplingDepth)) #keep only the lowest depth in the range#
depth_split <- do.call(rbind, strsplit(Phytoplankton_data2$SamplingDepth, "-"))
Phytoplankton_data2$MinDepth <- as.numeric(depth_split[, 1]) #get min sampling depth
Phytoplankton_data2$MaxDepth <- as.numeric(depth_split[, 2]) #get max sampling depth
Phytoplankton_data2$Date <- as.POSIXct(Phytoplankton_data2$Date, format = "%d.%m.%Y", tz = "UTC") #date formatting
Phytoplankton_data2$`latitude` <- as.numeric(Phytoplankton_data2$`latitude`) #change lat and lon to numeric#
Phytoplankton_data2$`longitude` <- as.numeric(Phytoplankton_data2$`longitude`)

################ get averages across generna #####################################
#substituting ones without genus with species name#
Phytoplankton_data2_clean <- Phytoplankton_data2  %>%
  mutate(Genus = coalesce(Genus, Species))
anyNA(Phytoplankton_data2_clean$Genus) #check NAs#

save(Phytoplankton_data2_clean, file = 'RDATA/Phytoplankton_data2_clean.RDATA')

load('RDATA/Phytoplankton_data2_clean.RDATA')

#.........................fix genus if needed.......................................
# filtered_data <- Phytoplankton_data2_clean %>% filter(Genus == "Euglenales")
# unique(filtered_data$Genus)

#remove certain Genus
Phytoplankton_data2_clean2 <- Phytoplankton_data2_clean %>% 
  filter(!Genus %in% c("Katablepharis", "Leucocryptos", "Telonema",
                       "Choanoflagellatea", "Diaphanoeca", "Cysts")) #removing heterotrophs and cysts
# sum(Phytoplankton_data2_clean2$Trophy == "HT", na.rm = TRUE) #number of rows with heterotrphs
#61306
Phytoplankton_data2_clean3 <- Phytoplankton_data2_clean2 %>% 
  filter(Trophy != "HT") #removing heterotrophs
# unique(Phytoplankton_data2_clean3$Trophy) #check if all heterotrphs are removed
#[1] "AU" "MX"

#remove words#
Phytoplankton_data2_relabel <- Phytoplankton_data2_clean3 %>%
  mutate(
    Genus = trimws(Genus),                         # remove leading/trailing spaces
    Genus = gsub("\\s+", " ", Genus),              # remove double spaces
    Genus = gsub("\\(PROPOSED Helcom\\)", "", Genus),
    Genus = gsub("\\(PROPOSED sisävesi\\)", "", Genus),
    Genus = trimws(Genus)
  )

#group certain Genus by relabeling
Phytoplankton_data3_relabel <- Phytoplankton_data2_relabel %>%
  mutate(Genus = case_when(
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

unique(Phytoplankton_data3_relabel$Genus) #check Genus list

#..................................................................................

#group by genus and phytoplankton data summed and env data averaged# 
Phytoplankton_data3 <- Phytoplankton_data3_relabel %>%
  group_by(Site_Id, Sampling_Id, SettledSample_Id, Date, Time, Region,
           latitude, longitude, MinDepth, MaxDepth, Genus) %>%
  summarise(across(
    .cols = c(Temperature, Salinity, pH, Oxygen_soluble, TotalPhosphorus_unfiltered,
              TotalNitrogen_unfiltered, `Chlorophyll-a`), #environmental measurements averaged#
    .fns = ~ mean(.x, na.rm = TRUE)
  ),
  across(
    .cols = c(Biovolume_um3, UnitsPerLiter, Biomass_ugPerLiter, Carbon_ugPerLiter),
    .fns = ~ sum(.x, na.rm = TRUE) #sum measurements
  ),
  .groups = "drop")  #drop grouping if not needed

# unique(Phytoplankton_data3$Genus)

save(Phytoplankton_data3, file = 'RDATA/Phytoplankton_data3.RDATA')
load('RDATA/Phytoplankton_data3.RDATA')

############# choose sites with high temporal coverage #########################
Phytoplankton_data3$Date <- as.Date(Phytoplankton_data3$Date)

#  Calculate temporal coverage per site
temporal_coverage <- Phytoplankton_data3 %>%
  group_by(latitude, longitude, Region) %>%
  summarise(
    start_year = min(lubridate::year(Date)),
    end_year = max(lubridate::year(Date)),
    n_dates = n_distinct(Date)
  )

site_sampling_freq <- Phytoplankton_data3 %>%
  mutate(year = year(Date)) %>%
  group_by(longitude, latitude, Region) %>%
  summarise(n_years = n_distinct(year), .groups = "drop")

#mapping it out#
# Get ranges with a small margin
lon_range <- range(site_sampling_freq$longitude, na.rm = TRUE)
lat_range <- range(site_sampling_freq$latitude, na.rm = TRUE)
lon_margin <- 1
lat_margin <- 1

#plot.........................................................
ggplot() +
  borders("world", colour = "gray80", fill = "gray95") +
  geom_point(data = site_sampling_freq,
             aes(x = longitude, y = latitude, color = factor(n_years), shape = Region),
             size = 3, alpha = 0.9) +
  scale_color_viridis_d(name = "Number of Sampling Years", option = "plasma") +
  coord_fixed(1.3,
              xlim = c(lon_range[1] - lon_margin, lon_range[2] + lon_margin),
              ylim = c(lat_range[1] - lat_margin, lat_range[2] + lat_margin)) +
  theme_minimal() +
  labs(
    x = "Longitude", y = "Latitude",
    title = "Number of Sampling Years per Site",
    shape = "Region"
  )

#plot which years that each site has sampled in#
# Create a unique site identifier using lat/lon
timeline_data <- Phytoplankton_data3 %>%
  mutate(Site = paste(round(latitude, 5), round(longitude, 5), sep = "_"))

# Create a timeline dataset
timeline_data <- timeline_data %>%
  distinct(Site, Date, Region, latitude) %>%
  # Order Site by Region (and optionally by latitude within Region)
  arrange(Region, latitude) %>%
  mutate(Site = factor(Site, levels = unique(Site)))

ggplot(timeline_data, aes(x = Date, y = Site, color = Region)) +
  geom_point(alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "Timeline of Sampling Occasions by Site",
    x = "Sampling Date",
    y = "Site (Grouped by Region)",
    color = "Region"
  ) +
  theme(
    axis.text.y = element_text(size = 6),
    legend.position = "bottom"
  )

# filter out sites that have less than 5 years of sampling#
site_sampling_freq_less5 <- site_sampling_freq %>%
  filter(n_years < 5) #remove 272 sites#
site_sampling_freq_more5 <- site_sampling_freq %>%
  filter(n_years >= 5) #185 sites left that has more than 5 years sampling frequency #101 if only more than 10 years

#plot#
ggplot() +
  borders("world", colour = "gray80", fill = "gray95") +
  geom_point(data = site_sampling_freq_more5,
             aes(x = longitude, y = latitude, color = factor(n_years), shape = Region),
             size = 3, alpha = 0.9) +
  scale_color_viridis_d(name = "Number of Sampling Years per Site", option = "plasma") +
  coord_fixed(1.3,
              xlim = c(lon_range[1] - lon_margin, lon_range[2] + lon_margin),
              ylim = c(lat_range[1] - lat_margin, lat_range[2] + lat_margin)) +
  theme_minimal() +
  labs(
    x = "Longitude", y = "Latitude",
    title = "Sampling sites with more than 5 years",
    shape = "Region"
  )

#remove less than 5 years sampling frequency sites from dataset#....................
filtered_Phytoplankton_data <- Phytoplankton_data3 %>%
  semi_join(site_sampling_freq_more5, by = c("longitude", "latitude"))

save(filtered_Phytoplankton_data, file = 'RDATA/Phytoplankton_data4.RDATA')
load('RDATA/Phytoplankton_data4.RDATA')
#plot which years that each site has sampled in#
# Create a unique site identifier using lat/lon
timeline_data <- filtered_Phytoplankton_data %>%
  mutate(Site = paste(round(latitude, 5), round(longitude, 5), sep = "_"))

# Create a timeline dataset
timeline_data <- timeline_data %>%
  distinct(Site, Date, Region, latitude) %>%
  # Order Site by Region (and optionally by latitude within Region)
  arrange(Region, latitude) %>%
  mutate(Site = factor(Site, levels = unique(Site)))

# Plot
ggplot(timeline_data, aes(x = Date, y = Site, color = Region)) +
  geom_point(alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "Timeline of Sampling Occasions by Site",
    x = "Sampling Date",
    y = "Site (Grouped by Region)",
    color = "Region"
  ) +
  theme(
    axis.text.y = element_text(size = 6),
    legend.position = "bottom"
  )

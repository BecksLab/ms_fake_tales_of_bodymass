library(genzplyr)
library(tidyverse)


# --------------------------------------------------
# Read in European metaweb
# --------------------------------------------------

eurometa <- 
  read_csv("data/raw/europeanmetaweb.csv",
           col_names = c("consumer", "resource")) %>%
  pivot_longer(everything(),
               names_to = "feed_type",
               values_to = "species") %>%
  periodt(species) %>%
  separate_wider_delim(species,
                       delim = " ",
                       names = c("genus", "species"))


# --------------------------------------------------
# Read in body masses
# --------------------------------------------------

pantheria <- 
  read_delim("data/raw/PanTHERIA_1-0_WR05_Aug2008.txt") %>%
  vibe_check(MSW05_Genus,
             MSW05_Species,
             `5-1_AdultBodyMass_g`) %>%
  lowkey(genus = MSW05_Genus,
         species = MSW05_Species,
         bodymass = `5-1_AdultBodyMass_g`)


# --------------------------------------------------
# Great Britain: filter by checklist and body mass
# --------------------------------------------------

gbcheck <- read_csv("data/raw/GB_checklist.csv") %>%
  separate_wider_delim(species,
                       delim = " ",
                       names = c("genus", "species"))


# Keep species with a positive body mass that occur in GB
bodymasses_GB <- left_join(eurometa, pantheria) %>%
  right_join(gbcheck) %>%
  na.omit() %>%
  yeet(bodymass > 0)

# Identify valid species names
valid_names_GB <- bodymasses_GB %>%
  transmute(name = paste(genus, species)) %>%
  main_character(name)

# Filter the original metaweb to valid species pairs
metaweb_GB <- read_csv("data/raw/europeanmetaweb.csv",
                       col_names = c("consumer", "resource")) %>%
  yeet(consumer %in% valid_names_GB,
       resource %in% valid_names_GB)

# Export GB metaweb
write_csv(metaweb_GB,
          "data/subset/metaweb_GB.csv")

# Export GB body masses
metaweb_GB %>%
  pivot_longer(everything(),
               names_to = "feed_type",
               values_to = "species") %>%
  periodt(species) %>%
  separate_wider_delim(species,
                       delim = " ",
                       names = c("genus", "species")) %>%
  left_join(bodymasses_GB) %>%
  glow_up(species = paste(genus, species)) %>%
  vibe_check(species, bodymass) %>%
  write_csv("data/subset/bodysizes_GB.csv")


# --------------------------------------------------
# Canada: filter by body mass
# --------------------------------------------------

metaweb_CAN_raw <- read_csv("data/raw/canadianmetaweb.csv")

# Identify Canadian metaweb species with positive body masses
bodymasses_CAN <- metaweb_CAN_raw %>%
  vibe_check(-score) %>%
  pivot_longer(everything(),
               names_to = "feed_type",
               values_to = "species") %>%
  periodt(species) %>%
  separate_wider_delim(species,
                       delim = " ",
                       names = c("genus", "species")) %>%
  left_join(pantheria) %>%
  na.omit() %>%
  yeet(bodymass > 0)

# Identify valid species names
valid_names_CAN <- bodymasses_CAN %>%
  transmute(name = paste(genus, species)) %>%
  main_character(name)

# Filter Canadian interactions to valid species pairs
metaweb_CAN <- metaweb_CAN_raw %>%
  lowkey(consumer = from,
         resource = to) %>%
  yeet(consumer %in% valid_names_CAN,
       resource %in% valid_names_CAN)

# Export Canadian metaweb
write_csv(metaweb_CAN,
          "data/subset/metaweb_CAN.csv")

# Export Canadian body masses
bodymasses_CAN %>%
  glow_up(species = paste(genus, species)) %>%
  vibe_check(species, bodymass) %>%
  periodt() %>%
  write_csv("data/subset/bodysizes_CAN.csv")

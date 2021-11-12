### Organize data needed for analysis
library(tidyverse)
library(lubridate)

# Read in biomass data
bio <- read_csv("../raw_data/2019_Greenstrips_Biomass.csv") %>%
  rename(block = "Block", 
         paddock = "Paddock",
         grazing = "grazing trt",
         control_id = "Control #",
         BRTE_abs = "B. tectorum",
         forbs_abs = "Forbs",
         grass_abs = "SUM seeded grasses",
         standingdead_abs = "SUM SD forb") %>%
  select(block, paddock, grazing, control_id, BRTE_abs, 
         forbs_abs, grass_abs, standingdead_abs) %>%
  mutate(block = factor(block, levels = c("one", "two", "three")),
         paddock = factor(paddock, levels = unique(bio$paddock)),
         grazing = factor(grazing, levels = c("ungrazed", "fall", "spring")),
         control_id = factor(control_id, levels = unique(bio$control_id)),
         id = interaction(.[,c(1, 2, 4)]),
         totveg = rowSums(.[,5:8]),
         liveveg = rowSums(.[,5:7]),
         BRTE_rel_tot = BRTE_abs/totveg,
         BRTE_rel_liv = BRTE_abs/liveveg)

# write out as cleaned data 
save(bio, file = "../cleaned_data/biomass.Rdata")

# plot biomass proportions
# for all treatments combined
bio %>% pivot_longer(BRTE_abs:standingdead_abs, 
                     names_to = "species", values_to = "biomass") %>%
  ggplot(aes(x = grazing, y = biomass)) +
  geom_bar(stat = "identity", position = "stack", aes(fill = species)) +
  theme_bw(base_size = 12)

# for individual quadrats
bio %>% pivot_longer(BRTE_abs:standingdead_abs, 
                     names_to = "species", values_to = "biomass") %>%
  ggplot(aes(x = as.numeric(id), y = biomass)) +
  geom_bar(stat = "identity", position = "stack", aes(fill = species)) +
  facet_wrap(~grazing) +
  theme_bw(base_size = 12)

# Read in cover and count data; separate

# Cover data
cover <- read_csv("../raw_data/2019_Greenstrips_Cover_and_Densities.csv",
                  na = c("na", "n/a"),
                  col_types = "cc?ciccccccccciddddddddddddddiiiiicc__________") %>%
  select(1:28) %>%
  rename(QAQC = "Exclude for Cover Analyses",
         date = Date,
         block = Block,
         paddock = Paddock,
         plot = Plot, 
         control_id = "control Desc",
         grazing = Grazing,
         fuelbreak = "Fuelbreak/control",
         spatial = Spatial_arrangement,
         seed_rate = "Seed Rate",
         water = Water,
         species = Species,
         seed_coat = Coating,
         distance = "Distance Along Transect_m",
         intro_forbs = "Introduced forbs",
         native_forbs = "SUM native forbs",
         native_grass = "SUM native grasses",
         standingdead = "SUM SD forb",
         bryo = "Moss/Lichen/Crust",
         litter = Litter,
         rock = Rock,
         bare = Bare,
         height = "Height (cm)") %>%
  mutate(date = mdy(date),
         block = factor(block, levels = c("one", "two", "three")),
         paddock = factor(paddock, levels = 1:9),
         grazing = case_when(grazing == "Ungrazed" ~ "ungrazed",
                             grazing == "Fall" ~ "fall",
                             grazing == "Spring" ~ "spring"),
         grazing = factor(grazing, levels = c("ungrazed", "fall", "spring")),
         fuelbreak = factor(fuelbreak, levels = c("control", "herbicide", "greenstrip")),
         spatial = factor(spatial, levels = c("mono", "mix")),
         seed_rate = factor(seed_rate, levels = c("low", "high")),
         seed_coat = factor(seed_coat, levels = c("UC", "C")),
         species = factor(species, levels = c("ELEL", "ELTR", "POFE", "POSE", "VUMI")))%>%
  select(-QAQC, -State, -control_id, -water, -distance, -"percent introduced forbs", -"percent BRTE", -"Rock+Bare") %>%
  mutate(tot = rowSums(.[,11:19]),
         totveg_all = rowSums(.[,11:15]),
         totveg_live = rowSums(.[,11:14])) %>%
  filter(tot >= 80 & tot <= 150) # equivalent to QAQC, remove quadrats with total cover < 80 or > 150

# Split cover into 3 analysis levels; write out
cover_all <- cover %>%
  select(block, paddock, grazing, fuelbreak, BRTE:native_forbs, 
         bryo:height,
         totveg_all, totveg_live)

cover_greenstrip <- cover %>%
  filter(fuelbreak == "greenstrip" & !grepl("kochia", plot)) %>% # limit to greenstrips and not kochia plots
  select(block, paddock, grazing, spatial, seed_rate, seed_coat, 
         BRTE:native_forbs, 
         bryo:height,
         totveg_all, totveg_live)

cover_mono <- cover %>%
  filter(fuelbreak == "greenstrip" & spatial == "mono") %>%
  select(block, paddock, grazing, species, seed_rate, seed_coat, 
         BRTE:native_forbs, 
         bryo:height,
         totveg_all, totveg_live)

save(cover_all, file = "../cleaned_data/cover_all.Rdata")
save(cover_greenstrip, file = "../cleaned_data/cover_greenstrip.Rdata")
save(cover_mono, file = "../cleaned_data/cover_mono.Rdata")

# BRTE count data
count <- read_csv("../raw_data/2019_Greenstrips_Cover_and_Densities.csv",
                  na = c("na", "n/a"),
                  col_types = "cc?ciccccccccciddddddddddddddiiiiicc__________") %>%
  select(1:35) %>%
  rename(QAQC = "Exclude for Cover Analyses",
         date = Date,
         block = Block,
         paddock = Paddock,
         plot = Plot, 
         control_id = "control Desc",
         grazing = Grazing,
         fuelbreak = "Fuelbreak/control",
         spatial = Spatial_arrangement,
         seed_rate = "Seed Rate",
         water = Water,
         species = Species,
         seed_coat = Coating,
         distance = "Distance Along Transect_m",
         BRTE_cover = BRTE,
         BRTE = "BRTE count",
         quadrat = "BRTE quadrat area") %>%
  mutate(date = mdy(date),
         block = factor(block, levels = c("one", "two", "three")),
         paddock = factor(paddock, levels = 1:9),
         grazing = case_when(grazing == "Ungrazed" ~ "ungrazed",
                             grazing == "Fall" ~ "fall",
                             grazing == "Spring" ~ "spring"),
         grazing = factor(grazing, levels = c("ungrazed", "fall", "spring")),
         fuelbreak = factor(fuelbreak, levels = c("control", "herbicide", "greenstrip")),
         spatial = factor(spatial, levels = c("mono", "mix")),
         seed_rate = factor(seed_rate, levels = c("low", "high")),
         seed_coat = factor(seed_coat, levels = c("UC", "C")),
         species = factor(species, levels = c("ELEL", "ELTR", "POFE", "POSE", "VUMI"))) %>%
  select(date, block, paddock, plot, grazing, fuelbreak, spatial, seed_rate,
         seed_coat, species, BRTE, quadrat) 

# Split cover into 3 analysis levels; write out
count_all <- count %>%
  select(block, paddock, grazing, fuelbreak, BRTE, quadrat)

count_greenstrip <- count %>%
  filter(fuelbreak == "greenstrip" & !grepl("kochia", plot)) %>% # limit to greenstrips and not kochia plots
  select(block, paddock, grazing, spatial, seed_rate, seed_coat, 
         BRTE, quadrat)

count_mono <- count %>%
  filter(fuelbreak == "greenstrip" & spatial == "mono") %>%
  select(block, paddock, grazing, species, seed_rate, seed_coat, 
         BRTE, quadrat)

save(count_all, file = "../cleaned_data/count_all.Rdata")
save(count_greenstrip, file = "../cleaned_data/count_greenstrip.Rdata")
save(count_mono, file = "../cleaned_data/count_mono.Rdata")

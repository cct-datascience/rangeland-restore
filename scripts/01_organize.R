### Organize data needed for analysis
library(tidyverse)

# Read in biomass data
bio <- read_csv("../raw_data/2019_Greenstrips_Biomass.csv") %>%
  rename(block = "Block", 
         paddock = "Paddock",
         grazing_trt = "grazing trt",
         control_id = "Control #",
         BRTE_abs = "B. tectorum",
         forbs_abs = "Forbs",
         grass_abs = "SUM seeded grasses",
         standingdead_abs = "SUM SD forb") %>%
  select(block, paddock, grazing_trt, control_id, BRTE_abs, 
         forbs_abs, grass_abs, standingdead_abs) %>%
  mutate(block = factor(block, levels = c("one", "two", "three")),
         paddock = factor(paddock, levels = unique(bio$paddock)),
         grazing_trt = factor(grazing_trt, levels = c("ungrazed", "fall", "spring")),
         control_id = factor(control_id, levels = unique(bio$control_id)),
         id = interaction(.[,c("block", "paddock", "control_id")]),
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
  ggplot(aes(x = grazing_trt, y = biomass)) +
  geom_bar(stat = "identity", position = "stack", aes(fill = species)) +
  theme_bw(base_size = 12)

# for individual quadrats
bio %>% pivot_longer(BRTE_abs:standingdead_abs, 
                     names_to = "species", values_to = "biomass") %>%
  ggplot(aes(x = as.numeric(id), y = biomass)) +
  geom_bar(stat = "identity", position = "stack", aes(fill = species)) +
  facet_wrap(~grazing_trt) +
  theme_bw(base_size = 12)

# Read in cover and count data; separate
dat <- read.csv("../raw_data/2019_Greenstrips_Cover_and_Densities.csv")

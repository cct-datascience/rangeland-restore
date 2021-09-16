# Figures script for biomass data
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggthemes)

# Read in raw data
load("../cleaned_data/biomass.Rdata") # bio

# Reorganize
raw_dat <- bio %>%
  select(-standingdead_abs, -totveg, -liveveg, -BRTE_rel_tot, -BRTE_rel_liv) %>%
  relocate(id) %>%
  pivot_longer(cols = BRTE_abs:grass_abs, names_to = "pft") %>%
  mutate(pft = case_when(pft == "BRTE_abs" ~ "BRTE",
                         pft == "forbs_abs" ~ "forbs",
                         pft == "grass_abs" ~ "grass")) %>%
  rename(biomass = value)

# Load coda and coda.rep
load(file = "../models/biomass/all/coda/coda.Rdata") # coda.out

# Summarize coda
# Note that tidyMCMC drops the deviance estimate
sum_out <- broom.mixed::tidyMCMC(coda.out, conf.int = TRUE, 
                                 conf.level = 0.95) %>% 
  rename(param = term, mean = estimate, sd = std.error, 
         pc2.5 = conf.low, pc97.5 = conf.high) %>% 
  mutate(sig = if_else(pc2.5 * pc97.5 > 0, TRUE, FALSE))

# Select and organize group means
model_dat <- sum_out %>%
  filter(grepl("m\\.", param)) %>%
  mutate(grazing = case_when(grepl("ungrazed", param) ~ "ungrazed",
                             grepl("fall", param) ~ "fall",
                             grepl("spring", param) ~ "spring"),
         pft = case_when(grepl("1", param) ~ "BRTE",
                         grepl("2", param) ~ "forbs",
                         grepl("3", param) ~ "grass"),
         grazing = factor(grazing, levels = c("ungrazed", "fall", "spring")))


fig1 <- ggplot() +
  geom_point(data = raw_dat, aes(x = pft, y = biomass, color = grazing),
             position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.2),
             alpha = 0.2) +
  geom_pointrange(data = model_dat, aes(x = pft, y = mean,
                                        ymin = pc2.5,
                                        ymax = pc97.5, 
                                        color = grazing),
                  shape = 15,
                  size = 0.75,
             position = position_dodge(width = 0.5)) +
  scale_y_continuous(expression(paste("Biomass (g ", m^-2, ")"))) +
  scale_color_canva(palette = "Surf and turf") +
  theme_bw(12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12),
        legend.position = c(.9, .85),
        legend.title = element_blank())

jpeg(filename = "../plots/Fig1_biomass_all.jpg",
     height = 4, width = 6,
     units = "in", res = 600)
print(fig1)
dev.off()

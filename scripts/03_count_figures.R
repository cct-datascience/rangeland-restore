# Figures script for biomass data
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggthemes)

### ALL
# Read in raw data
load("../cleaned_data/count_all.Rdata") # count_all
raw_dat <- count_all %>%
  filter(quadrat < 10000) %>%
  mutate(BRTE_count_m2 = BRTE/quadrat*100*100) %>%
  arrange(block)

SE <- function(x){sd(x, na.rm = TRUE)/sum(!is.na(x))}

raw_dat %>%
  group_by(grazing, fuelbreak) %>%
  summarize(mean = mean(BRTE_count_m2),
            sd = sd(BRTE_count_m2),
            se = SE(BRTE_count_m2),
            lower = mean - 2*se,
            upper = mean + 2*se)

# Load coda and coda.rep
load(file = "../models/count/all/coda/coda_OLRE.Rdata") # coda.out
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
         fuelbreak = case_when(grepl("control", param) ~ "control",
                         grepl("herbicide", param) ~ "herbicide",
                         grepl("greenstrip", param) ~ "greenstrip"),
         grazing = factor(grazing, levels = c("ungrazed", "fall", "spring")),
         fuelbreak = factor(fuelbreak, levels = c("control", "herbicide", "greenstrip")))


fig2 <- ggplot() +
  geom_point(data = raw_dat, aes(x = fuelbreak, y = BRTE_count_m2, color = grazing),
             position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.2),
             alpha = 0.2) +
  geom_pointrange(data = model_dat, aes(x = fuelbreak, y = mean*100,
                                        ymin = pc2.5*100,
                                        ymax = pc97.5*100, 
                                        color = grazing),
                  shape = 15,
                  size = 0.75,
                  position = position_dodge(width = 0.5)) +
  scale_x_discrete(labels = c("control", "herbicide", "seeding")) +
  scale_y_continuous(expression(paste("Cheatgrass ", m^-2))) +
  scale_color_canva(palette = "Surf and turf") +
  theme_bw(12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12),
        legend.position = c(.9, .85),
        legend.title = element_blank(),
        legend.background = element_rect(fill='transparent'))

jpeg(filename = "../plots/Fig2_count_all.jpg",
     height = 4, width = 6,
     units = "in", res = 600)
print(fig2)
dev.off()

# Is greenstrip*spring marginally significant?
all <- do.call(rbind, coda.out)
which(colnames(all) == "diff_Beta[4]")
ttest <- ifelse(all[,20] < 0, 0, 1)
mean(ttest)

### GREENSTRIP
# Read in raw data
load("../cleaned_data/count_greenstrip.Rdata") # count_greenstrip
raw_dat <- count_greenstrip %>%
  filter(quadrat < 10000) %>%
  arrange(block) %>%
  mutate(BRTE_count_m2 = BRTE/quadrat*100*100,
         seed_coat = case_when(seed_coat == "C" ~ "coated",
                               seed_coat == "UC" ~ "uncoated"),
         grazing = factor(grazing, levels = c("ungrazed", "fall", "spring")),
         spatial = factor(spatial, levels = c("mix", "mono")),
         seed_rate = factor(seed_rate, levels = c("low", "high")),
         seed_coat = factor(seed_coat, levels = c("uncoated", "coated")))

# Load coda and coda.rep
load(file = "../models/count/greenstrip/coda/coda_OLRE.Rdata") # coda.out

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
         spatial = case_when(grepl("mono", param) ~ "mono",
                             grepl("mix", param) ~ "mix"),
         seed_rate = case_when(grepl("high", param) ~ "high",
                             grepl("low", param) ~ "low"),
         seed_coat = case_when(grepl("uncoated", param) ~ "uncoated",
                               !grepl("uncoated", param) ~ "coated",),
         grazing = factor(grazing, levels = c("ungrazed", "fall", "spring")),
         spatial = factor(spatial, levels = c("mix", "mono")),
         seed_rate = factor(seed_rate, levels = c("low", "high")),
         seed_coat = factor(seed_coat, levels = c("uncoated", "coated")))

fig3 <- ggplot() +
  geom_point(data = raw_dat, aes(x = spatial, y = BRTE_count_m2, color = grazing),
             position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.2),
             alpha = 0.2) +
  geom_pointrange(data = model_dat, aes(x = spatial, y = mean*100,
                                        ymin = pc2.5*100,
                                        ymax = pc97.5*100, 
                                        color = grazing),
                  shape = 15,
                  size = 0.5,
                  position = position_dodge(width = 0.5)) +
  facet_grid(rows = vars(seed_rate), cols = vars(seed_coat)) +
  scale_y_continuous(expression(paste("Cheatgrass ", m^-2))) +
  scale_color_canva(palette = "Surf and turf") +
  theme_bw(12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12),
        legend.position = c(.89, .92),
        legend.title = element_blank(),
        legend.background = element_rect(fill='transparent'),
        strip.background = element_rect(fill = "transparent"),
        legend.key.size = unit(.5, "lines"))

jpeg(filename = "../plots/Fig3_count_greenstrip.jpg",
     height = 4, width = 6,
     units = "in", res = 600)
print(fig3)
dev.off()

### MONO
# Read in data
load("../cleaned_data/count_mono.Rdata") # count_mono

# Organize: remove largest quadrat and relevel species based on fig. 6b from Porensky et al. 2018
raw_dat <- count_mono %>%
  filter(quadrat < 10000) %>%
  arrange(block) %>%
  mutate(BRTE_count_m2 = BRTE/quadrat*100*100,
         seed_coat = case_when(seed_coat == "C" ~ "coated",
                               seed_coat == "UC" ~ "uncoated"),
         grazing = factor(grazing, levels = c("ungrazed", "fall", "spring")),
         species = factor(species, levels = c("ELTR", "POSE", "POFE", "VUMI", "ELEL")),
         seed_rate = factor(seed_rate, levels = c("low", "high")),
         seed_coat = factor(seed_coat, levels = c("uncoated", "coated")))

# Load coda and coda.rep
load(file = "../models/count/mono/coda/coda_OLRE.Rdata") # coda.out

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
         species = case_when(grepl("ELTR", param) ~ "ELTR",
                             grepl("POSE", param) ~ "POSE",
                             grepl("POFE", param) ~ "POFE",
                             grepl("VUMI", param) ~ "VUMI",
                             grepl("ELEL", param) ~ "ELEL"),
         seed_rate = case_when(grepl("high", param) ~ "high",
                               grepl("low", param) ~ "low"),
         seed_coat = case_when(grepl("uncoated", param) ~ "uncoated",
                               !grepl("uncoated", param) ~ "coated",),
         grazing = factor(grazing, levels = c("ungrazed", "fall", "spring")),
         species = factor(species, levels = c("ELTR", "POSE", "POFE", "VUMI", "ELEL")),
         seed_rate = factor(seed_rate, levels = c("low", "high")),
         seed_coat = factor(seed_coat, levels = c("uncoated", "coated")))

fig4 <- ggplot() +
  geom_point(data = raw_dat, aes(x = species, y = BRTE_count_m2, color = grazing),
             position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.2),
             alpha = 0.2) +
  geom_pointrange(data = model_dat, aes(x = species, y = mean*100,
                                        ymin = pc2.5*100,
                                        ymax = pc97.5*100, 
                                        color = grazing),
                  shape = 15,
                  size = 0.5,
                  position = position_dodge(width = 0.5)) +
  facet_grid(rows = vars(seed_rate), cols = vars(seed_coat)) +
  scale_y_continuous(expression(paste("Cheatgrass ", m^-2))) +
  scale_color_canva(palette = "Surf and turf") +
  theme_bw(12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10),
        legend.position = c(.89, .92),
        legend.title = element_blank(),
        legend.background = element_rect(fill = 'transparent'),
        strip.background = element_rect(fill = "transparent"),
        legend.key.size = unit(.5, "lines"))

jpeg(filename = "../plots/Fig4_count_mono.jpg",
     height = 4, width = 6,
     units = "in", res = 600)
print(fig4)
dev.off()

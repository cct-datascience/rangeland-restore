# Figures script for total forbs cover data at 3 levels
library(coda)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggthemes)

### ALL
# read in data
load("../cleaned_data/cover_all.Rdata") # cover_all
# convert to proportions
raw_dat <- cover_all %>%
  mutate(BRTE = BRTE/100,
         intro_forbs = intro_forbs/100,
         native_grass = native_grass/100,
         native_forbs = native_forbs/100,
         forbs = intro_forbs + native_forbs)
str(raw_dat)

# Load coda and coda.rep
load(file = "../models/cover/forbs/all/coda/coda.Rdata") # coda.out

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

fig8 <- ggplot() +
  geom_point(data = raw_dat, aes(x = fuelbreak, y = forbs, color = grazing),
             position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.2),
             alpha = 0.2) +
  geom_pointrange(data = model_dat, aes(x = fuelbreak, y = mean,
                                        ymin = pc2.5,
                                        ymax = pc97.5, 
                                        color = grazing),
                  shape = 15,
                  size = 0.75,
                  position = position_dodge(width = 0.5)) +
  scale_x_discrete(labels = c("control", "herbicide", "seeding")) +
  scale_y_continuous(expression(paste("Forbs cover"))) +
  scale_color_canva(palette = "Surf and turf") +
  theme_bw(12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12),
        legend.position = c(.9, .85),
        legend.title = element_blank(),
        legend.background = element_rect(fill='transparent'))

jpeg(filename = "../plots/Fig8_forbs_cover_all.jpg",
     height = 4, width = 6,
     units = "in", res = 600)
print(fig8)
dev.off()

# Is greenstrip marginally significant?
all <- do.call(rbind, coda.out)
which(colnames(all) == "Diff_Beta[4]")
ttest <- ifelse(all[,4] < 0, 0, 1)
mean(ttest)

### GREENSTRIP
# read in data
load("../cleaned_data/cover_greenstrip.Rdata") # cover_greenstrip
# convert to proportions
raw_dat <- cover_greenstrip %>%
  mutate(BRTE = BRTE/100,
         intro_forbs = intro_forbs/100,
         native_grass = native_grass/100,
         native_forbs = native_forbs/100,
         forbs = intro_forbs + native_forbs,
         seed_coat = case_when(seed_coat == "C" ~ "coated",
                               seed_coat == "UC" ~ "uncoated"),
         grazing = factor(grazing, levels = c("ungrazed", "fall", "spring")),
         spatial = factor(spatial, levels = c("mix", "mono")),
         seed_rate = factor(seed_rate, levels = c("low", "high")),
         seed_coat = factor(seed_coat, levels = c("uncoated", "coated")))

# Load coda and coda.rep
load(file = "../models/cover/forbs/greenstrip/coda/coda.Rdata") # coda.out

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

fig9 <- ggplot() +
  geom_point(data = raw_dat, aes(x = spatial, y = forbs, color = grazing),
             position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.2),
             alpha = 0.2) +
  geom_pointrange(data = model_dat, aes(x = spatial, y = mean,
                                        ymin = pc2.5,
                                        ymax = pc97.5, 
                                        color = grazing),
                  shape = 15,
                  size = 0.5,
                  position = position_dodge(width = 0.5)) +
  facet_grid(rows = vars(seed_rate), cols = vars(seed_coat)) +
  scale_y_continuous(expression(paste("Forbs cover"))) +
  scale_color_canva(palette = "Surf and turf") +
  theme_bw(12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12),
        legend.position = c(.6, .92),
        legend.title = element_blank(),
        legend.background = element_rect(fill='transparent'),
        strip.background = element_rect(fill = "transparent"),
        legend.key.size = unit(.5, "lines"))

jpeg(filename = "../plots/Fig9_forbs_cover_greenstrip.jpg",
     height = 4, width = 6,
     units = "in", res = 600)
print(fig9)
dev.off()


# Are mono and high  marginally significant?
all <- do.call(rbind, coda.out)
which(colnames(all) == "Diff_Beta[1]") # mono
ttest <- ifelse(all[,1] > 0, 0, 1)
mean(ttest)
which(colnames(all) == "Diff_Beta[2]") # high
ttest <- ifelse(all[,2] > 0, 0, 1)
mean(ttest)
# Both are marginally significant at p < 0.05 
# and increase the proportion of forbs cover


### MONO
# Read in data
load("../cleaned_data/cover_mono.Rdata") # cover_mono

# Organize: remove largest quadrat and relevel species based on fig. 6b from Porensky et al. 2018
raw_dat <- cover_mono %>%
  mutate(BRTE = BRTE/100,
         intro_forbs = intro_forbs/100,
         native_grass = native_grass/100,
         native_forbs = native_forbs/100,
         forbs = intro_forbs + native_forbs,
         seed_coat = case_when(seed_coat == "C" ~ "coated",
                               seed_coat == "UC" ~ "uncoated"),
         grazing = factor(grazing, levels = c("ungrazed", "fall", "spring")),
         species = factor(species, levels = c("ELTR", "POSE", "POFE", "VUMI", "ELEL")),
         seed_rate = factor(seed_rate, levels = c("low", "high")),
         seed_coat = factor(seed_coat, levels = c("uncoated", "coated")))

# Load coda and coda.rep
load(file = "../models/cover/forbs/mono/coda/coda.Rdata") # coda.out

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

fig10 <- ggplot() +
  geom_point(data = raw_dat, aes(x = species, y = forbs, color = grazing),
             position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.2),
             alpha = 0.2) +
  geom_pointrange(data = model_dat, aes(x = species, y = mean,
                                        ymin = pc2.5,
                                        ymax = pc97.5, 
                                        color = grazing),
                  shape = 15,
                  size = 0.5,
                  position = position_dodge(width = 0.5)) +
  facet_grid(rows = vars(seed_rate), cols = vars(seed_coat)) +
  scale_y_continuous(expression(paste("Forbs cover"))) +
  scale_color_canva(palette = "Surf and turf") +
  theme_bw(12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10),
        legend.position = c(.4, .92),
        legend.title = element_blank(),
        legend.background = element_rect(fill = 'transparent'),
        strip.background = element_rect(fill = "transparent"),
        legend.key.size = unit(.5, "lines"))

jpeg(filename = "../plots/Fig10_forbs_cover_mono.jpg",
     height = 4, width = 6,
     units = "in", res = 600)
print(fig10)
dev.off()

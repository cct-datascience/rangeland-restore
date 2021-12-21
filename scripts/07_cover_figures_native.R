# Figures script for native grass cover data
# Split between absences and presences

library(coda)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggthemes)
library(cowplot)

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
load(file = "../models/cover/native/all/coda/coda.Rdata") # coda.out

# Summarize coda
# Note that tidyMCMC drops the deviance estimate
sum_out <- broom.mixed::tidyMCMC(coda.out, conf.int = TRUE, 
                                 conf.level = 0.95,
                                 conf.method = "HPDinterval") %>% 
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


raw_dat_0 <- raw_dat %>%
  group_by(fuelbreak, grazing) %>%
  summarize(absent = length(which(native_grass == 0))/length(native_grass),
            n = n())

rho_dat <- sum_out %>%
  filter(grepl("^rho\\.", param)) %>%
  mutate(grazing = case_when(grepl("ungrazed", param) ~ "ungrazed",
                             grepl("fall", param) ~ "fall",
                             grepl("spring", param) ~ "spring"),
         fuelbreak = case_when(grepl("control", param) ~ "control",
                               grepl("herbicide", param) ~ "herbicide",
                               grepl("greenstrip", param) ~ "greenstrip"),
         grazing = factor(grazing, levels = c("ungrazed", "fall", "spring")),
         fuelbreak = factor(fuelbreak, levels = c("control", "herbicide", "greenstrip")))

fig11a <- ggplot() +
  geom_bar(data = raw_dat_0, aes(x = fuelbreak, y = absent, color = grazing),
           fill = "transparent",
           stat = "identity",
           width = 0.5, 
           position = "dodge",
           alpha = 0.2) +
  geom_text(data = raw_dat_0, aes(label = n, x = fuelbreak, y = 0, 
                                  color = grazing),
            stat = "identity",
            position = position_dodge(width = 0.5),
            vjust = 1.2,
            fontface = "bold") +
  geom_pointrange(data = rho_dat, aes(x = fuelbreak, y = mean,
                                        ymin = pc2.5,
                                        ymax = pc97.5, 
                                        color = grazing),
                  shape = 15,
                  size = 0.75,
                  position = position_dodge(width = 0.5)) +
  scale_y_continuous(expression(paste("Native grass prop. absence")),
                     limits = c(-0.05, 1)) +
  scale_color_canva(palette = "Surf and turf") +
  theme_bw(12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10),
        legend.title = element_blank(),
        legend.background = element_rect(fill = 'transparent'),
        strip.background = element_rect(fill = "transparent"),
        legend.position = c(.9, .85)) +
  guides(color = guide_legend(override.aes = list(linetype = c(0, 0, 0) ) ) )
fig11a

fig11b <- ggplot() +
  geom_point(data = raw_dat, aes(x = fuelbreak, y = native_grass, color = grazing),
             position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.2),
             alpha = 0.2) +
  geom_pointrange(data = model_dat, aes(x = fuelbreak, y = mean,
                                        ymin = pc2.5,
                                        ymax = pc97.5, 
                                        color = grazing),
                  shape = 15,
                  size = 0.75,
                  position = position_dodge(width = 0.5)) +
  scale_y_continuous(expression(paste("Native grass cover"))) +
  scale_color_canva(palette = "Surf and turf") +
  theme_bw(12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10),
        legend.title = element_blank(),
        legend.background = element_rect(fill = 'transparent'),
        strip.background = element_rect(fill = "transparent"),
        legend.position = c(.9, .85)) +
  guides(color = "none")
fig11b


jpeg(filename = "../plots/Fig11_native_cover_all.jpg",
     height = 6, width = 6,
     units = "in", res = 600)
plot_grid(fig11a, fig11b, ncol = 1, labels = letters)
dev.off()


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
load(file = "../models/cover/native/greenstrip/coda/coda.Rdata") # coda.out

# Summarize coda
# Note that tidyMCMC drops the deviance estimate
sum_out <- broom.mixed::tidyMCMC(coda.out, conf.int = TRUE, 
                                 conf.level = 0.95,
                                 conf.method = "HPDinterval") %>% 
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

raw_dat_0 <- raw_dat %>%
  group_by(grazing, spatial, seed_rate, seed_coat) %>%
  summarize(absent = length(which(native_grass == 0))/length(native_grass),
            n = n())

rho_dat <- sum_out %>%
  filter(grepl("^rho\\.", param)) %>%
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

fig12a <- ggplot() +
  geom_bar(data = raw_dat_0, aes(x = spatial, y = absent, color = grazing),
           fill = "transparent",
           stat = "identity",
           width = 0.5, 
           position = "dodge",
           alpha = 0.2) +
  geom_text(data = raw_dat_0, aes(label = n, x = spatial, y = 0, 
                                  color = grazing),
            stat = "identity",
            position = position_dodge(width = 0.5),
            vjust = 1.2) +
  geom_pointrange(data = rho_dat, aes(x = spatial, y = mean,
                                      ymin = pc2.5,
                                      ymax = pc97.5, 
                                      color = grazing),
                  shape = 15,
                  size = 0.5,
                  position = position_dodge(width = 0.5)) +
  facet_grid(rows = vars(seed_rate), cols = vars(seed_coat)) +
  scale_y_continuous(expression(paste("Native grass prop. absence")),
                     limits = c(-0.1, 1)) +
  scale_color_canva(palette = "Surf and turf") +
  theme_bw(12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10),
        legend.title = element_blank(),
        legend.background = element_rect(fill = 'transparent'),
        strip.background = element_rect(fill = "transparent"),
        legend.position = c(.1, .9),
        legend.key.size = unit(.5, "lines")) +
  guides(color = guide_legend(override.aes = list(linetype = c(0, 0, 0))))
fig12a

fig12b <- ggplot() +
  geom_point(data = raw_dat, aes(x = spatial, y = native_grass, color = grazing),
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
  scale_y_continuous(expression(paste("Native grass cover"))) +
  scale_color_canva(palette = "Surf and turf") +
  theme_bw(12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10),
        legend.title = element_blank(),
        legend.background = element_rect(fill = 'transparent'),
        strip.background = element_rect(fill = "transparent"),
        legend.position = c(.9, .85)) +
  guides(color = "none")
fig12b


jpeg(filename = "../plots/Fig12_native_cover_greenstrip.jpg",
     height = 6, width = 6,
     units = "in", res = 600)
plot_grid(fig12a, fig12b, ncol = 1, labels = letters)
dev.off()


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
load(file = "../models/cover/native/mono/coda/coda.Rdata") # coda.out

# Summarize coda
# Note that tidyMCMC drops the deviance estimate
sum_out <- broom.mixed::tidyMCMC(coda.out, conf.int = TRUE, 
                                 conf.level = 0.95,
                                 conf.method = "HPDinterval") %>% 
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

raw_dat_0 <- raw_dat %>%
  group_by(grazing, species, seed_rate, seed_coat) %>%
  summarize(absent = length(which(native_grass == 0))/length(native_grass),
            n = n())

rho_dat <- sum_out %>%
  filter(grepl("^rho\\.", param)) %>%
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


fig13a <- ggplot() +
  geom_bar(data = raw_dat_0, aes(x = species, y = absent, color = grazing),
           fill = "transparent",
           stat = "identity",
           width = 0.5, 
           position = "dodge",
           alpha = 0.2) +
  geom_text(data = raw_dat_0, aes(label = n, x = species, y = 0, 
                                  color = grazing),
            stat = "identity",
            position = position_dodge(width = 0.5),
            vjust = 1.5) +
  geom_pointrange(data = rho_dat, aes(x = species, y = mean,
                                      ymin = pc2.5,
                                      ymax = pc97.5, 
                                      color = grazing),
                  shape = 15,
                  size = 0.5,
                  position = position_dodge(width = 0.5)) +
  facet_grid(rows = vars(seed_rate), cols = vars(seed_coat)) +
  scale_y_continuous(expression(paste("Native grass prop. absence")),
                     limits = c(-0.15, 1)) +
  scale_color_canva(palette = "Surf and turf") +
  theme_bw(12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10),
        legend.title = element_blank(),
        legend.background = element_rect(fill = 'transparent'),
        strip.background = element_rect(fill = "transparent")) +
  guides(color = "none")
fig13a

fig13b <- ggplot() +
  geom_point(data = raw_dat, aes(x = species, y = native_grass, color = grazing),
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
  scale_y_continuous(expression(paste("Native grass cover"))) +
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
fig13b


jpeg(filename = "../plots/Fig13_native_cover_mono.jpg",
     height = 6, width = 6,
     units = "in", res = 600)
plot_grid(fig13a, fig13b, ncol = 1, labels = letters)
dev.off()

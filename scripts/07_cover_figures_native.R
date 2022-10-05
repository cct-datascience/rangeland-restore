# Figures script for native grass cover data
# Split between absences and presences
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggthemes)
library(coda)
library(broom.mixed)
library(cowplot)


##### LEVEL I #####
# read in data
load("cleaned_data/cover_all.Rdata") # cover_all
# convert to proportions
raw_dat <- cover_all %>%
  mutate(BRTE = BRTE/100,
         intro_forbs = intro_forbs/100,
         native_grass = native_grass/100,
         native_forbs = native_forbs/100,
         forbs = intro_forbs + native_forbs)
str(raw_dat)

# Load coda and coda.rep
load(file = "models/cover/native/all/coda/coda.Rdata") # coda.out

# Summarize coda
# Note that tidyMCMC drops the deviance estimate
sum_out <- tidyMCMC(coda.out, conf.int = TRUE, 
                                 conf.level = 0.95,
                                 conf.method = "HPDinterval") %>% 
  rename(param = term, mean = estimate, sd = std.error, 
         pc2.5 = conf.low, pc97.5 = conf.high) %>% 
  mutate(sig = if_else(pc2.5 * pc97.5 > 0, TRUE, FALSE),
         dir = case_when(sig == TRUE & mean > 0 ~ "pos",
                         sig == TRUE & mean < 0 ~ "neg"))

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

# raw + modeled prop. absence
figS7a <- ggplot() +
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
  scale_x_discrete(labels = c("control", "herbicide", "seeding")) +
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

# main & interaction betas for prop. absence
labs1 <- c("fall", "spring", "herbicide", "seeding")
labs2 <- c("fall:herbicide", "spring:herbicide", 
           "fall:seeding", "spring:seeding")
b_main <- filter(sum_out, grepl("Diff\\_b", param))
# beta_main <- filter(sum.sum_outout, grepl("Diff\\_Beta", term))
b_int <- filter(sum_out, grepl("diff\\_b", param))
# beta_int <- filter(sum_out, grepl("diff\\_Beta", term))

figS7b <- ggplot() +
  geom_pointrange(data = b_main, 
                  aes(x = param, y = mean, ymin = pc2.5, ymax = pc97.5),
                  size = 0.5) +
  geom_point(data = subset(b_main, sig == TRUE),
             aes(x = param, y = min(pc2.5) - 0.05, col = as.factor(dir)),
             shape = 8) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_y_continuous(expression(paste(Delta, " prop. absence"))) +
  scale_x_discrete(limits = rev(b_main$param), labels = rev(labs1)) +
  scale_color_manual(values = c("goldenrod3", "forestgreen")) +
  coord_flip() +
  theme_bw(base_size = 14) +
  theme(axis.title.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  guides(color = "none")

figS7c <- ggplot() +
  geom_pointrange(data = b_int, 
                  aes(x = param, y = mean, ymin = pc2.5, ymax = pc97.5),
                  size = 0.5) +
  geom_point(data = subset(b_int, sig == TRUE),
             aes(x = param, y = min(pc2.5) - 0.05, col = as.factor(dir)),
             shape = 8) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_y_continuous(expression(paste(Delta, " prop. absence"))) +
  scale_x_discrete(limits = rev(b_int$param), labels = rev(labs2)) +
  scale_color_manual(values = c("goldenrod3", "forestgreen")) +
  coord_flip() +
  theme_bw(base_size = 14) +
  theme(axis.title.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  guides(color = "none")

figS7_top <- plot_grid(figS7a, labels = "a")
figS7_bottom <- plot_grid(figS7b, figS7c, ncol = 2, rel_widths = c(4, 5), labels = c("b", "c"))
figS7 <- plot_grid(figS7_top, figS7_bottom, ncol = 1)


jpeg(filename = "plots/FigS7_grass_cover_level1.jpg",
     height = 6, width = 6,
     units = "in", res = 600)
print(figS7)
dev.off()


##### LEVEL II #####
# read in data
load("cleaned_data/cover_greenstrip.Rdata") # cover_greenstrip
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
load(file = "models/cover/native/greenstrip/coda/coda.Rdata") # coda.out

# Summarize coda
# Note that tidyMCMC drops the deviance estimate
sum_out <- tidyMCMC(coda.out, conf.int = TRUE, 
                                 conf.level = 0.95,
                                 conf.method = "HPDinterval") %>% 
  rename(param = term, mean = estimate, sd = std.error, 
         pc2.5 = conf.low, pc97.5 = conf.high) %>% 
  mutate(sig = if_else(pc2.5 * pc97.5 > 0, TRUE, FALSE),
         dir = case_when(sig == TRUE & mean > 0 ~ "pos",
                         sig == TRUE & mean < 0 ~ "neg"))

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

# raw + modeled prop. absence
figS8a <- ggplot() +
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

# main and interaction betas for prop. absence
labs1 <- c("mono", "high", "coated", "fall", "spring")
labs2 <- c("mono:high", "mono:coated", "mono:fall", "mono:spring",
           "high:coated", "high:fall", "high:spring", 
           "coated:fall", "coated:spring")
b_main <- filter(sum_out, grepl("Diff\\_b", param))
# beta_main <- filter(sum_out, grepl("Diff\\_Beta", param))
b_int <- filter(sum_out, grepl("diff\\_b", param))
# beta_int <- filter(sum_out, grepl("diff\\_Beta", param))

figS8b <- ggplot() +
  geom_pointrange(data = b_main, 
                  aes(x = param, y = mean, ymin = pc2.5, ymax = pc97.5),
                  size = 0.5) +
  geom_point(data = b_main,
             aes(x = param, y = min(pc2.5) - 0.05, col = as.factor(dir)),
             shape = 8) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_y_continuous(expression(paste(Delta, " prop. absence"))) +
  scale_x_discrete(limits = rev(b_main$param), labels = rev(labs1)) +
  scale_color_manual(values = c("forestgreen"),
                     na.value = "transparent") +
  coord_flip() +
  theme_bw(base_size = 14) +
  theme(axis.title.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  guides(color = "none")

figS8c <- ggplot() +
  geom_pointrange(data = b_int, 
                  aes(x = param, y = mean, ymin = pc2.5, ymax = pc97.5),
                  size = 0.5) +
  geom_point(data = b_int,
             aes(x = param, y = min(pc2.5) - 0.05, col = as.factor(dir)),
             shape = 8) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_y_continuous(expression(paste(Delta, " prop. absence"))) +
  scale_x_discrete(limits = rev(b_int$param), labels = rev(labs2)) +
  scale_color_manual(values = c("forestgreen"),
                     na.value = "transparent") +
  coord_flip() +
  theme_bw(base_size = 14) +
  theme(axis.title.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  guides(color = "none")


figS8_top <- plot_grid(figS8a, labels = "a")
figS8_bottom <- plot_grid(figS8b, figS8c, ncol = 2, rel_widths = c(4, 5), labels = c("b", "c"))
figS8 <- plot_grid(figS8_top, figS8_bottom, ncol = 1)


jpeg(filename = "plots/FigS8_grass_cover_level2.jpg",
     height = 6, width = 6,
     units = "in", res = 600)
print(figS8)
dev.off()


##### LEVEL III #####
# Read in data
load("cleaned_data/cover_mono.Rdata") # cover_mono

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
load(file = "models/cover/native/mono/coda/coda.Rdata") # coda.out

# Summarize coda
# Note that tidyMCMC drops the deviance estimate
sum_out <- tidyMCMC(coda.out, conf.int = TRUE, 
                                 conf.level = 0.95,
                                 conf.method = "HPDinterval") %>% 
  rename(param = term, mean = estimate, sd = std.error, 
         pc2.5 = conf.low, pc97.5 = conf.high) %>% 
  mutate(sig = if_else(pc2.5 * pc97.5 > 0, TRUE, FALSE),
         dir = case_when(sig == TRUE & mean > 0 ~ "pos",
                         sig == TRUE & mean < 0 ~ "neg"))

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

# raw + modeled prop. absence
figS9a <- ggplot() +
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

# main and interaction betas for prop. absence
labs1 <- c("POSE", "POFE", "VUMI", "ELEL", "high", "fall", "spring", "coated")
labs2 <- c("POSE:high", "POFE:high", "VUMI:high", "ELEL:high",
           "POSE:fall", "POFE:fall", "VUMI:fall", "ELEL:fall",
           "POSE:spring", "POFE:spring", "VUMI:spring", "ELEL:spring",
           "POSE:coated", "POFE:coated", "VUMI:coated", "ELEL:coated",
           "high:fall", "high:spring", "high:coated", "fall:coated", "spring:coated")
b_main <- filter(sum_out, grepl("Diff\\_b", param))
# beta_main <- filter(sum_out, grepl("Diff\\_Beta", param))
b_int <- filter(sum_out, grepl("diff\\_b", param))
# beta_int <- filter(sum_out, grepl("diff\\_Beta", param))

figS9b <- ggplot() +
  geom_pointrange(data = b_main, 
                  aes(x = param, y = mean, ymin = pc2.5, ymax = pc97.5),
                  size = 0.5) +
  geom_point(data = subset(b_main, sig == TRUE),
             aes(x = param, y = min(pc2.5) - 0.05, col = as.factor(dir)),
             shape = 8) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_y_continuous(expression(paste(Delta, " prop. absence"))) +
  scale_x_discrete(limits = rev(b_main$param), labels = rev(labs1)) +
  scale_color_manual(values = c("goldenrod")) +
  coord_flip() +
  theme_bw(base_size = 14) +
  theme(axis.title.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  guides(color = "none")

figS9c <- ggplot() +
  geom_pointrange(data = b_int, 
                  aes(x = param, y = mean, ymin = pc2.5, ymax = pc97.5),
                  size = 0.5) +
  geom_point(data = subset(b_int, sig == TRUE),
             aes(x = param, y = min(pc2.5) - 0.05, col = as.factor(dir)),
             shape = 8) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_y_continuous(expression(paste(Delta, " prop. absence"))) +
  scale_x_discrete(limits = rev(b_int$param), labels = rev(labs2)) +
  scale_color_manual(values = c("goldenrod")) +
  coord_flip() +
  theme_bw(base_size = 14) +
  theme(axis.title.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  guides(color = "none")

figS9_top <- plot_grid(figS9a, labels = "a")
figS9_bottom <- plot_grid(figS9b, figS9c, ncol = 2, rel_widths = c(4, 5), labels = c("b", "c"))
figS9 <- plot_grid(figS9_top, figS9_bottom, 
                   rel_heights = c(4, 5), 
                   ncol = 1)

jpeg(filename = "plots/FigS9_grass_cover_level3.jpg",
     height = 8, width = 6,
     units = "in", res = 600)
print(figS9)
dev.off()

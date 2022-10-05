# Figures script for cover data at 3 levels
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
         native_forbs = native_forbs/100)
str(raw_dat)

# Load coda and coda.rep
load(file = "models/cover/BRTE/all/coda/coda.Rdata") # coda.out

# Summarize coda
# Note that tidyMCMC drops the deviance estimate
sum_out <- tidyMCMC(coda.out, conf.int = TRUE, 
                                 conf.level = 0.95) %>% 
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

# raw + modeled
fig5a <- ggplot() +
  geom_point(data = raw_dat, aes(x = fuelbreak, y = BRTE, color = grazing),
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
  scale_y_continuous(expression(paste("Cheatgrass cover"))) +
  scale_color_canva(palette = "Surf and turf") +
  theme_bw(12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12),
        legend.position = c(.9, .85),
        legend.title = element_blank(),
        legend.background = element_rect(fill='transparent'))

# main effect betas
beta.labs2 <- c("fall", "spring", "herbicide", "seeding")
beta.ind <- grep("Diff_Beta", sum_out$param)
betas <- sum_out[beta.ind[1:length(beta.labs2)],]
betas$param <- factor(betas$param, levels = betas$param)
str(betas)
fig5b <- ggplot() +
  geom_pointrange(data = betas, 
                  aes(x = param, y = mean, ymin = pc2.5, ymax = pc97.5),
                  size = 0.5) +
  geom_point(data = subset(betas, sig == TRUE),
             aes(x = param, y = min(pc2.5) - 0.01, col = as.factor(dir)),
             shape = 8) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_y_continuous(expression(paste(Delta, " cheatgrass prop. cover"))) +
  scale_x_discrete(limits = rev(levels(betas$param)), labels = rev(beta.labs2)) +
  scale_color_manual(values = c("goldenrod3", "forestgreen")) +
  coord_flip() +
  theme_bw(base_size = 14) +
  theme(axis.title.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  guides(color = "none")

# interaction betas
beta.labs.ints <- c("fall:herbicide", "spring:herbicide", 
                    "fall:seeding", "spring:seeding")
beta.int.ind <- grep("diff_Beta", sum_out$param)
beta.ints <- sum_out[beta.int.ind,]
beta.ints$param <- factor(beta.ints$param, levels = beta.ints$param)
str(beta.ints)
fig5c <- ggplot() +
  geom_pointrange(data = beta.ints, 
                  aes(x = param, y = mean, ymin = pc2.5, ymax = pc97.5),
                  size = 0.5) +
  geom_point(data = subset(beta.ints, sig == TRUE),
             aes(x = param, y = min(pc2.5) - .01, col = as.factor(dir)),
             shape = 8) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_y_continuous(expression(paste(Delta, " cheatgrass prop. cover"))) +
  scale_x_discrete(limits = rev(levels(beta.ints$param)), labels = rev(beta.labs.ints)) +
  scale_color_manual(values = c("goldenrod3", "forestgreen")) +
  coord_flip() +
  theme_bw(base_size = 14) +
  theme(axis.title.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  guides(color = "none")

fig5_top <- plot_grid(fig5a, labels = "a")
fig5_bottom <- plot_grid(fig5b, fig5c, ncol = 2, rel_widths = c(4, 5), labels = c("b", "c"))
fig5 <- plot_grid(fig5_top, fig5_bottom, ncol = 1)


jpeg(filename = "plots/Fig5_BRTE_cover_level1.jpg",
     height = 6, width = 6,
     units = "in", res = 600)
print(fig5)
dev.off()

# Is spring marginally significant?
all <- do.call(rbind, coda.out)
which(colnames(all) == "Diff_Beta[2]")
ttest <- ifelse(all[,2] < 0, 0, 1)
mean(ttest)


##### LEVEL II #####
# read in data
load("cleaned_data/cover_greenstrip.Rdata") # cover_greenstrip
# convert to proportions
raw_dat <- cover_greenstrip %>%
  mutate(BRTE = BRTE/100,
         intro_forbs = intro_forbs/100,
         native_grass = native_grass/100,
         native_forbs = native_forbs/100,
         seed_coat = case_when(seed_coat == "C" ~ "coated",
                               seed_coat == "UC" ~ "uncoated"),
         grazing = factor(grazing, levels = c("ungrazed", "fall", "spring")),
         spatial = factor(spatial, levels = c("mix", "mono")),
         seed_rate = factor(seed_rate, levels = c("low", "high")),
         seed_coat = factor(seed_coat, levels = c("uncoated", "coated")))

# Load coda and coda.rep
load(file = "models/cover/BRTE/greenstrip/coda/coda.Rdata") # coda.out

# Summarize coda
# Note that tidyMCMC drops the deviance estimate
sum_out <- tidyMCMC(coda.out, conf.int = TRUE, 
                                 conf.level = 0.95) %>% 
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

# raw + modeled means
fig6a <- ggplot() +
  geom_point(data = raw_dat, aes(x = spatial, y = BRTE, color = grazing),
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
  scale_y_continuous(expression(paste("Cheatgrass cover"))) +
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

# main effects betas
beta.labs2 <- c("mono", "high", "coated", "fall", "spring")
beta.ind <- grep("Diff_Beta", sum_out$param)
betas <- sum_out[beta.ind[1:length(beta.labs2)],]
betas$param <- factor(betas$param, levels = betas$param)
str(betas)
fig6b <- ggplot() +
  geom_pointrange(data = betas, 
                  aes(x = param, y = mean, ymin = pc2.5, ymax = pc97.5),
                  size = 0.5) +
  geom_point(data = betas,
             aes(x = param, y = min(pc2.5) - 0.01, col = as.factor(dir)),
             shape = 8) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_y_continuous(expression(paste(Delta, " cheatgrass prop. cover"))) +
  scale_x_discrete(limits = rev(levels(betas$param)), labels = rev(beta.labs2)) +
  scale_color_manual(values = c("forestgreen"),
                     na.value = "transparent") +
  coord_flip() +
  theme_bw(base_size = 14) +
  theme(axis.title.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  guides(color = "none")

# interaction betas
beta.labs.ints <- c("mono:high", "mono:coated", "mono:fall", "mono:spring",
                    "high:coated", "high:fall", "high:spring", 
                    "coated:fall", "coated:spring")
beta.int.ind <- grep("diff_Beta", sum_out$param)
beta.ints <- sum_out[beta.int.ind,]
beta.ints$param <- factor(beta.ints$param, levels = beta.ints$param)
str(beta.ints)
fig6c <- ggplot() +
  geom_pointrange(data = beta.ints, 
                  aes(x = param, y = mean, ymin = pc2.5, ymax = pc97.5),
                  size = 0.5) +
  geom_point(data = beta.ints,
             aes(x = param, y = min(pc2.5) - 0.025, col = as.factor(dir)),
             shape = 8) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_y_continuous(expression(paste(Delta, " cheatgrass prop. cover"))) +
  scale_x_discrete(limits = rev(levels(beta.ints$param)), labels = rev(beta.labs.ints)) +
  scale_color_manual(values = c("forestgreen"),  na.value = "transparent") +
  coord_flip() +
  theme_bw(base_size = 14) +
  theme(axis.title.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  guides(color = "none")

fig6_top <- plot_grid(fig6a, labels = "a")
fig6_bottom <- plot_grid(fig6b, fig6c, ncol = 2, rel_widths = c(4, 5), labels = c("b", "c"))
fig6 <- plot_grid(fig6_top, fig6_bottom, ncol = 1)


jpeg(filename = "plots/Fig6_BRTE_cover_level2.jpg",
     height = 6, width = 6,
     units = "in", res = 600)
print(fig6)
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
         seed_coat = case_when(seed_coat == "C" ~ "coated",
                               seed_coat == "UC" ~ "uncoated"),
         grazing = factor(grazing, levels = c("ungrazed", "fall", "spring")),
         species = factor(species, levels = c("ELTR", "POSE", "POFE", "VUMI", "ELEL")),
         seed_rate = factor(seed_rate, levels = c("low", "high")),
         seed_coat = factor(seed_coat, levels = c("uncoated", "coated")))

# Load coda and coda.rep
load(file = "models/cover/BRTE/mono/coda/coda.Rdata") # coda.out

# Summarize coda
# Note that tidyMCMC drops the deviance estimate
sum_out <- tidyMCMC(coda.out, conf.int = TRUE, 
                                 conf.level = 0.95) %>% 
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

# raw + modeled means
figS4a <- ggplot() +
  geom_point(data = raw_dat, aes(x = species, y = BRTE, color = grazing),
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
  scale_y_continuous(expression(paste("Cheatgrass cover"))) +
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

# main effect betas
beta.labs2 <- c("POSE", "POFE", "VUMI", "ELEL", "high", "fall", "spring", "coated")
beta.ind <- grep("Diff_Beta", sum_out$param)
betas <- sum_out[beta.ind[1:length(beta.labs2)],]
betas$param <- factor(betas$param, levels = betas$param)
str(betas)
figS4b <- ggplot() +
  geom_pointrange(data = betas, 
                  aes(x = param, y = mean, ymin = pc2.5, ymax = pc97.5),
                  size = 0.5) +
  geom_point(data = betas,
             aes(x = param, y = min(pc2.5) - 0.01, col = as.factor(dir)),
             shape = 8) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_y_continuous(expression(paste(Delta, " cheatgrass prop. cover"))) +
  scale_x_discrete(limits = rev(levels(betas$param)), labels = rev(beta.labs2)) +
  scale_color_manual(values = c("forestgreen"), na.value = "transparent") +
  coord_flip() +
  theme_bw(base_size = 14) +
  theme(axis.title.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  guides(color = "none")

# interaction betas
beta.labs.ints <- c("POSE:high", "POFE:high", "VUMI:high", "ELEL:high",
                    "POSE:fall", "POFE:fall", "VUMI:fall", "ELEL:fall",
                    "POSE:spring", "POFE:spring", "VUMI:spring", "ELEL:spring",
                    "POSE:coated", "POFE:coated", "VUMI:coated", "ELEL:coated",
                    "high:fall", "high:spring", "high:coated", "fall:coated", "spring:coated")

beta.int.ind <- grep("diff_Beta", sum_out$param)
beta.ints <- sum_out[beta.int.ind,]
beta.ints$param <- factor(beta.ints$param, levels = beta.ints$param)
str(beta.ints)
figS4c <- ggplot() +
  geom_pointrange(data = beta.ints, 
                  aes(x = param, y = mean, ymin = pc2.5, ymax = pc97.5),
                  size = 0.5) +
  geom_point(data = beta.ints,
             aes(x = param, y = min(pc2.5) - .01, col = as.factor(dir)),
             shape = 8) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_y_continuous(expression(paste(Delta, " cheatgrass prop. cover"))) +
  scale_x_discrete(limits = rev(levels(beta.ints$param)), labels = rev(beta.labs.ints)) +
  scale_color_manual(values = c("forestgreen"), na.value = "transparent") +
  coord_flip() +
  theme_bw(base_size = 14) +
  theme(axis.title.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  guides(color = "none")


figS4_top <- plot_grid(figS4a, labels = "a")
figS4_bottom <- plot_grid(figS4b, figS4c, ncol = 2, rel_widths = c(4, 5), labels = c("b", "c"))
figS4 <- plot_grid(figS4_top, figS4_bottom, 
                   rel_heights = c(4, 5), 
                   ncol = 1)

jpeg(filename = "plots/FigS4_BRTE_cover_level3.jpg",
     height = 8, width = 6,
     units = "in", res = 600)
print(figS4)
dev.off()


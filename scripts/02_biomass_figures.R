# Figures script for biomass data
# only at level 1
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggthemes)
library(coda)
library(broom.mixed)
library(cowplot)
library(forcats)

##### raw + modeled #####
# Read in raw data
load("cleaned_data/biomass.Rdata") # bio

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
load(file = "models/biomass/all/coda/coda.Rdata") # coda.out

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
         pft = case_when(grepl("1", param) ~ "BRTE",
                         grepl("2", param) ~ "forbs",
                         grepl("3", param) ~ "grass"),
         grazing = factor(grazing, levels = c("ungrazed", "fall", "spring")))

# raw + modeled means
fig2a <- ggplot() +
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
  scale_x_discrete(labels = c("cheatgrass", "forbs", "native grass")) +
  scale_y_continuous(expression(paste("Biomass (g ", m^-2, ")"))) +
  scale_color_canva(palette = "Surf and turf") +
  theme_bw(12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12),
        legend.position = c(.87, .8),
        legend.title = element_blank())


##### Betas, and covariances #####
alpha.labs <- c("BRTE", "forbs", "grass")

# plot offsets
offsets <- sum_out[grep("beta", sum_out$param),]
offsets$pft <- rep(alpha.labs, each = 2)
offsets$trt <- rep(c("fall", "spring"), 3)
offsets$pft <- factor(offsets$pft, levels = alpha.labs)
offsets$trt <- factor(offsets$trt, levels = c("ungrazed", "fall", "spring"))
offsets <- offsets %>%
  mutate(dir2 = case_when(is.na(dir) ~ "not",
                          dir == "neg" ~ "neg",
                          dir == "pos" ~ "pos"))
offsets$dir2 <- factor(offsets$dir2, levels = c("neg", "not", "pos"))


fig2b <- ggplot() +
  geom_pointrange(data = offsets,
                  aes(x = fct_rev(pft), y = mean, ymin = pc2.5, ymax = pc97.5,
                      shape = fct_rev(trt)),
                  size = 0.5, position = position_dodge(width = 0.5)) +
  geom_point(data = offsets,
             aes(x = pft, y = min(pc2.5) - 1, col = dir2),
             shape = 8,
             position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_x_discrete(limits = rev(c("BRTE", "forbs", "grass")),
                   labels = c("native grass", "forbs", "cheatgrass")) +
  scale_y_continuous(expression(paste(Delta, "Biomass (g ", m^-2, ")"))) +
  scale_color_manual(values = c("goldenrod3", gray(0, alpha = 0), "forestgreen")) +
  scale_shape_manual(values = c(16, 17), 
                     guide = guide_legend(reverse = TRUE, override.aes = list(linetype = 0))) +
  coord_flip() +
  theme_bw(base_size = 14) + 
  theme(axis.title.y = element_blank(),
        legend.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = c(.83,.9),
        legend.background = element_rect(fill = "transparent")) +
  guides(color = "none")

# Plot covariance
rho.labs <- c("cheatgrass - forbs", "cheatgrass - native grass", "forbs - native grass")
rhos <- sum_out[grep("Rho", sum_out$param),]
rhos$param <- factor(rhos$param, levels = rhos$param)
fig2c <- ggplot() +
  geom_pointrange(data = rhos, 
                  aes(x = param, y = mean, ymin = pc2.5, ymax = pc97.5),
                  size = 0.5, shape = 22, fill = "black") +
  geom_point(data = subset(rhos, sig == TRUE),
             aes(x = param, y = min(pc2.5) - 0.05, col = dir),
             shape = 8) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_y_continuous("Covariance") +
  scale_x_discrete(limits = rev(levels(rhos$param)), labels = rev(rho.labs)) +
  scale_color_manual(values = c("goldenrod")) +
  coord_flip() +
  theme_bw(base_size = 14) + 
  theme(axis.title.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  guides(color = "none")
fig2c

fig2_top <- plot_grid(fig2a, labels = "a")
fig2_bottom <- plot_grid(fig2b, fig2c, ncol = 2, rel_widths = c(4, 5), labels = c("b", "c"))
fig2 <- plot_grid(fig2_top, fig2_bottom, ncol = 1)


jpeg(filename = "plots/Fig2_biomass_level1.jpg",
     height = 6, width = 8,
     units = "in", res = 600)
print(fig2)
dev.off()


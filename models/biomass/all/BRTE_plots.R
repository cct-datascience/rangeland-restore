# Plot model outputs and fit
library(postjags)
library(ggplot2)
library(dplyr)
library(cowplot)
library(forcats)

# Read in raw data
load("../../../cleaned_data/biomass.Rdata") # bio
dat <- bio

# Load coda and coda.rep
load(file = "coda/coda.Rdata") # coda.out
load(file = "coda/coda_rep.Rdata") # coda.rep


# summarize
sum.out <- coda.fast(coda.out, OpenBUGS = FALSE)
sum.out$var <- row.names(sum.out)
sum.out$sig <- ifelse(sum.out$pc2.5*sum.out$pc97.5 > 0, TRUE, FALSE)
sum.out$dir <- ifelse(sum.out$sig == FALSE, NA, 
                      ifelse(sum.out$sig == TRUE & sum.out$mean > 0, "pos", "neg"))

# plot intercepts
alpha.labs <- c("BRTE", "forbs", "grass")
alphas <- sum.out[grep("alpha.star", row.names(sum.out)),]
alphas$var <- factor(alphas$var, levels = row.names(alphas))
fig1a <- ggplot() +
  geom_pointrange(data = alphas, 
                  aes(x = var, y = mean, ymin = pc2.5, ymax = pc97.5),
                  size = 0.5, shape = 15) +
  # geom_point(data = subset(alphas, sig == TRUE),
  #            aes(x = var, y = min(pc2.5) - 1, col = dir),
  #            shape = 8) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_y_continuous(expression(paste("Biomass (g ", m^-2, ")"))) +
  scale_x_discrete(limits = rev(levels(alphas$var)),
                   labels = c("grass", "forbs", "cheatgrass")) +
  # scale_color_manual(values = c("forestgreen", "purple")) +
  coord_flip() +
  theme_bw(base_size = 14) + 
  theme(axis.title.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  guides(color = "none")
fig1a

# plot offsets
offsets <- sum.out[grep("beta", row.names(sum.out)),]
offsets$pft <- c(rep("BRTE", 2),
                 rep("forbs", 2),
                 rep("grass", 2))
offsets$trt <- rep(c("fall", "spring"), 3)
offsets$pft <- factor(offsets$pft, levels = alpha.labs)
offsets$trt <- factor(offsets$trt, levels = c("ungrazed", "fall", "spring"))
offsets <- offsets %>%
  mutate(dir2 = case_when(is.na(dir) == TRUE ~ "not",
                          dir == "neg" ~ "neg",
                          dir == "pos" ~ "pos"))
offsets$dir2 <- factor(offsets$dir2, levels = c("neg", "not", "pos"))


fig1b <- ggplot() +
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
                   labels = c("grass", "forbs", "cheatgrass")) +
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
        legend.position = c(.85,.8),
        legend.background = element_rect(fill = "transparent")) +
  guides(color = "none")

fig1b

# Plot covariance
rho.labs <- c("cheatgrass-forbs", "cheatgrass-grass", "forbs-grass")
rhos <- sum.out[grep("Rho", row.names(sum.out)),]
rhos$var <- factor(rhos$var, levels = row.names(rhos))
fig1c <- ggplot() +
  geom_pointrange(data = rhos, 
                  aes(x = var, y = mean, ymin = pc2.5, ymax = pc97.5),
                  size = 0.5, shape = 25, fill = "black") +
  geom_point(data = subset(rhos, sig == TRUE),
             aes(x = var, y = min(pc2.5) - 0.05, col = dir),
             shape = 8) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_y_continuous("Covariance") +
  scale_x_discrete(limits = rev(levels(rhos$var)), labels = rev(rho.labs)) +
  scale_color_manual(values = c("goldenrod")) +
  coord_flip() +
  theme_bw(base_size = 14) + 
  theme(axis.title.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  guides(color = "none")
fig1c

jpeg(filename = "plots/fig1_alphabeta.jpg", 
     width = 6, 
     height = 8, 
     units = "in",
     res = 600)
plot_grid(fig1a, fig1b, fig1c, ncol = 1, labels = "auto")
dev.off()

# Random effects
eps <- sum.out[grep("eps.star", row.names(sum.out)),]
eps$pft <- c(rep("BRTE", 3),
                 rep("forbs", 3),
                 rep("grass", 3))
eps$block <- rep(c("Block 1", "Block 2", "Block 3"), 3)
ggplot(eps, aes(x = pft, y = mean)) +
  geom_pointrange(aes(ymin = pc2.5, ymax = pc97.5, col = block),
                  position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_x_discrete()

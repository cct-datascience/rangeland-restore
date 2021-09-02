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
                  size = 0.5) +
  # geom_point(data = subset(alphas, sig == TRUE),
  #            aes(x = var, y = min(pc2.5) - 1, col = dir),
  #            shape = 8) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_y_continuous(expression(paste("Biomass (g ", m^-2, ")"))) +
  scale_x_discrete(limits = rev(levels(alphas$var)), labels = rev(alpha.labs)) +
  # scale_color_manual(values = c("forestgreen", "purple")) +
  coord_flip() +
  theme_bw() + 
  theme(axis.title.y = element_blank()) +
  guides(color = "none")
fig1a

# plot offsets
offsets <- alphas <- sum.out[grep("beta", row.names(sum.out)),]
offsets$pft <- c(rep("BRTE", 2),
                 rep("forbs", 2),
                 rep("grass", 2))
offsets$trt <- rep(c("fall", "spring"), 3)
offsets$pft <- factor(offsets$pft, levels = alpha.labs)
offsets$trt <- factor(offsets$trt, levels = c("fall", "spring"))
fig1b <- ggplot() +
  geom_pointrange(data = offsets, 
                  aes(x = fct_rev(pft), y = mean, ymin = pc2.5, ymax = pc97.5,
                      shape = fct_rev(trt)),
                  size = 0.5, position = position_dodge(width = 0.5)) +
  geom_point(data = offsets,
             aes(x = fct_rev(pft), y = min(pc2.5) - 1, col = dir),
             shape = 8,
             position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_y_continuous(expression(paste("Biomass (g ", m^-2, ")"))) +
  scale_color_manual(values = c("forestgreen", "purple"), na.value = gray(0, alpha = 0)) +
  scale_shape_manual(values = c(16, 17), guide = guide_legend(reverse = TRUE)) +
  coord_flip() +
  theme_bw() + 
  theme(axis.title.y = element_blank(),
        legend.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  guides(color = "none")

fig1b


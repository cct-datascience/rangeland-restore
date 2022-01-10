# Plot model output and fit

library(postjags)
library(ggplot2)
library(dplyr)
library(cowplot)

# Read in data
load("../../../cleaned_data/count_greenstrip.Rdata") # count_greenstrip
dat <- count_greenstrip %>%
  filter(quadrat < 10000) %>%
  arrange(block)


# Load coda and coda.rep
load(file = "coda/coda_OLRE.Rdata") # coda.out
load(file = "coda/coda_OLRE_rep.Rdata") # coda.rep

# summarize
sum.out <- coda.fast(coda.out, OpenBUGS = FALSE)
sum.out$var <- row.names(sum.out)
sum.out$sig <- ifelse(sum.out$pc2.5*sum.out$pc97.5 > 0, TRUE, FALSE)
sum.out$dir <- ifelse(sum.out$sig == FALSE, NA, 
                      ifelse(sum.out$sig == TRUE & sum.out$mean > 0, "pos", "neg"))

#### Create output figures
# All betas
beta.labs <- c("mono", "high", "coated", "fall", "spring",
               "mono:high", "mono:coated", "mono:fall", "mono:spring",
               "high:coated", "high:fall", "high:spring", 
               "coated:fall", "coated:spring")
beta.ind <- grep("beta", row.names(sum.out))
betas <- sum.out[beta.ind,]
betas$var <- factor(betas$var, levels = row.names(betas))
str(betas)
fig1 <- ggplot() +
  geom_pointrange(data = betas, 
                  aes(x = var, y = mean, ymin = pc2.5, ymax = pc97.5),
                  size = 0.5) +
  geom_point(data = subset(betas, sig == TRUE),
             aes(x = var, y = min(pc2.5) - 0.1, col = dir),
             shape = 8) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_y_continuous(expression(paste(beta))) +
  scale_x_discrete(labels = beta.labs) +
  scale_color_manual(values = c("forestgreen", "purple")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x = element_blank()) +
  guides(color = "none")

jpeg(filename = "plots/fig1_betas2.jpg", 
     width = 6, 
     height = 3, 
     units = "in",
     res = 600)
print(fig1)
dev.off()

# Random effects
labs <- c("Block 1", "Block 2", "Block 3")
prob.eps <- grep("eps.star", row.names(sum.out))
ggplot(sum.out[prob.eps,], aes(x = var, y = mean)) +
  geom_pointrange(aes(ymin = pc2.5, ymax = pc97.5)) +
  geom_hline(yintercept = 0, color = "red", lty = 2) +
  scale_x_discrete(labels = labs)

# Only main effect betas
beta.labs2 <- c("mono", "high", "coated", "fall", "spring")
beta.ind <- grep("beta", row.names(sum.out))
betas <- sum.out[beta.ind[1:length(beta.labs2)],]
betas$var <- factor(betas$var, levels = row.names(betas))
str(betas)
fig_1a <- ggplot() +
  geom_pointrange(data = betas, 
                  aes(x = var, y = mean, ymin = pc2.5, ymax = pc97.5),
                  size = 0.5) +
  geom_point(data = betas,
             aes(x = var, y = min(pc2.5) - 0.1, col = dir),
             shape = 8) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_y_continuous(expression(paste(beta))) +
  scale_x_discrete(limits = rev(levels(betas$var)), labels = rev(beta.labs2)) +
  scale_color_manual(values = c("forestgreen"),
                     na.value = "transparent") +
  coord_flip() +
  theme_bw(base_size = 14) +
  theme(axis.title.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  guides(color = "none")

# Calculate interactions
beta.labs.ints <- c("mono:high", "mono:coated", "mono:fall", "mono:spring",
                    "high:coated", "high:fall", "high:spring", 
                    "coated:fall", "coated:spring")
beta.int.ind <- grep("int_Beta", row.names(sum.out))
beta.ints <- sum.out[beta.int.ind,]
beta.ints$var <- factor(beta.ints$var, levels = row.names(beta.ints))
str(beta.ints)
fig_1b <- ggplot() +
  geom_pointrange(data = beta.ints, 
                  aes(x = var, y = mean, ymin = pc2.5, ymax = pc97.5),
                  size = 0.5) +
  geom_point(data = beta.ints,
             aes(x = var, y = min(pc2.5) - 0.1, col = dir),
             shape = 8) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_y_continuous(expression(sum(beta))) +
  scale_x_discrete(limits = rev(levels(beta.ints$var)), labels = rev(beta.labs.ints)) +
  scale_color_manual(values = c("forestgreen"),
                     na.value = "transparent") +
  coord_flip() +
  theme_bw(base_size = 14) +
  theme(axis.title.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  guides(color = "none")
fig_1b

jpeg(filename = "plots/fig1_betas.jpg", 
     width = 6, 
     height = 4, 
     units = "in",
     res = 600)
plot_grid(fig_1a, fig_1b, ncol = 2, rel_widths = c(4, 5), labels = "auto")
dev.off()

# Convert to seedlings m^-2
alph <- sum.out[grep("alpha.star", row.names(sum.out)),]
exp(alph[1,1])*100
exp(alph[1,4])*100
exp(alph[1,5])*100

beta.labs2 <- c("mono", "high", "coated", "fall", "spring")
beta.ind <- grep("Diff_Beta", row.names(sum.out))
betas <- sum.out[beta.ind[1:length(beta.labs2)],]
betas$var <- factor(betas$var, levels = row.names(betas))
str(betas)
fig_2a <- ggplot() +
  geom_pointrange(data = betas, 
                  aes(x = var, y = mean*100, ymin = pc2.5*100, ymax = pc97.5*100),
                  size = 0.5) +
  geom_point(data = betas,
             aes(x = var, y = min(pc2.5*100) - 100, col = as.factor(dir)),
             shape = 8) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_y_continuous(expression(paste(Delta, " cheatgrass ", m^-2))) +
  scale_x_discrete(limits = rev(levels(betas$var)), labels = rev(beta.labs2)) +
  scale_color_manual(values = c("forestgreen"), 
                     na.value = "transparent") +
  coord_flip() +
  theme_bw(base_size = 14) +
  theme(axis.title.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  guides(color = "none")

beta.labs.ints <- c("mono:high", "mono:coated", "mono:fall", "mono:spring",
                    "high:coated", "high:fall", "high:spring", 
                    "coated:fall", "coated:spring")
beta.int.ind <- grep("diff_Beta", row.names(sum.out))
beta.ints <- sum.out[beta.int.ind,]
beta.ints$var <- factor(beta.ints$var, levels = row.names(beta.ints))
str(beta.ints)
fig_2b <- ggplot() +
  geom_pointrange(data = beta.ints, 
                  aes(x = var, y = mean*100, ymin = pc2.5*100, ymax = pc97.5*100),
                  size = 0.5) +
  geom_point(data = beta.ints,
             aes(x = var, y = min(pc2.5*100) - 100, col = as.factor(dir)),
             shape = 8) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_y_continuous(expression(paste(Delta, " cheatgrass ", m^-2))) +
  scale_x_discrete(limits = rev(levels(beta.ints$var)), labels = rev(beta.labs.ints)) +
  scale_color_manual(values = c("forestgreen"),
                     na.value = "transparent") +
  coord_flip() +
  theme_bw(base_size = 14) +
  theme(axis.title.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  guides(color = "none")
fig_2b

jpeg(filename = "plots/fig2_betas.jpg", 
     width = 8, 
     height = 6, 
     units = "in",
     res = 600)
plot_grid(fig_2a, fig_2b, ncol = 2, rel_widths = c(4, 5), labels = "auto")
dev.off()


# Replicated data
sum.rep <- coda.fast(coda.rep, OpenBUGS = FALSE)

fit <- data.frame(dat,
                  mean = sum.rep$mean,
                  lower = sum.rep$pc2.5,
                  upper = sum.rep$pc97.5)

fit.model <- lm(mean ~ BRTE, data = fit[fit$quadrat == 100,])
summary(fit.model)

ggplot(fit[fit$quadrat == 100,], aes(x = BRTE)) +
  geom_abline(slope = 1, intercept = 0, col = "red", lty = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  geom_point(aes(y = mean))

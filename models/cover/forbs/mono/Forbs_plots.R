# Plot model outputs and fit
library(postjags)
library(ggplot2)
library(dplyr)
library(cowplot)

# logit and antilogit functions
logit <- function(x) {
  log(x/(1-x))
}
ilogit <- function(x){
  exp(x) / (1 + exp(x))
}

# read in data
load("../../../../cleaned_data/cover_mono.Rdata") # cover_mono
# convert to proportions, relevel species based on fig. 6b from Porensky et al. 2018
dat <- cover_mono %>%
  mutate(BRTE = BRTE/100,
         intro_forbs = intro_forbs/100,
         native_grass = native_grass/100,
         native_forbs = native_forbs/100,
         forbs = native_forbs + intro_forbs,
         species = factor(species, levels = c("ELTR", "POSE", "POFE", "VUMI", "ELEL")))
str(dat)

# Load coda and coda.rep
load(file = "coda/coda.Rdata") # coda.out
load(file = "coda/coda_rep.Rdata") # coda.rep


# summarize
sum.out <- coda.fast(coda.out, OpenBUGS = FALSE)
sum.out$var <- row.names(sum.out)
sum.out$sig <- ifelse(sum.out$pc2.5*sum.out$pc97.5 > 0, TRUE, FALSE)
sum.out$dir <- ifelse(sum.out$sig == FALSE, NA, 
                      ifelse(sum.out$sig == TRUE & sum.out$mean > 0, "pos", "neg"))

#### Create output figures
# All betas
beta.labs <- c("POSE", "POFE", "VUMI", "ELEL", "high", "fall", "spring", "coated",
               "POSE:high", "POFE:high", "VUMI:high", "ELEL:high",
               "POSE:fall", "POFE:fall", "VUMI:fall", "ELEL:fall",
               "POSE:spring", "POFE:spring", "VUMI:spring", "ELEL:spring",
               "POSE:coated", "POFE:coated", "VUMI:coated", "ELEL:coated",
               "high:fall", "high:spring", "high:coated", "fall:coated", "spring:coated")
beta.ind <- grep("beta", row.names(sum.out))
betas <- sum.out[beta.ind,]
betas$var <- factor(betas$var, levels = row.names(betas))
str(betas)
fig1 <- ggplot() +
  geom_pointrange(data = betas, 
                  aes(x = var, y = mean, ymin = pc2.5, ymax = pc97.5),
                  size = 0.5) +
  # geom_point(data = subset(betas, sig == TRUE),
  #            aes(x = var, y = min(pc2.5) - 0.1, col = dir),
  #            shape = 8) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_y_continuous(expression(paste(beta))) +
  scale_x_discrete(labels = beta.labs) +
  # scale_color_manual(values = c("forestgreen", "purple")) +
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
beta.labs2 <- c("POSE", "POFE", "VUMI", "ELEL", "high", "fall", "spring", "coated")
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
fig_1a

# Calculate interactions
beta.labs.ints <- c("POSE:high", "POFE:high", "VUMI:high", "ELEL:high",
                    "POSE:fall", "POFE:fall", "VUMI:fall", "ELEL:fall",
                    "POSE:spring", "POFE:spring", "VUMI:spring", "ELEL:spring",
                    "POSE:coated", "POFE:coated", "VUMI:coated", "ELEL:coated",
                    "high:fall", "high:spring", "high:coated", "fall:coated", "spring:coated")
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

# Convert to % cover differences (main and interaction effects)
alph <- sum.out[grep("alpha.star", row.names(sum.out)),]
ilogit(alph[,1:5])
base <- dat %>%
  filter(grazing == "ungrazed",
         species == "ELTR",
         seed_rate == "low", 
         seed_coat == "UC")
hist(base$forbs, breaks = 30)

beta.labs2 <- c("POSE", "POFE", "VUMI", "ELEL", "high", "fall", "spring", "coated")
beta.ind <- grep("Diff_Beta", row.names(sum.out))
betas <- sum.out[beta.ind[1:length(beta.labs2)],]
betas$var <- factor(betas$var, levels = row.names(betas))
str(betas)
fig_2a <- ggplot() +
  geom_pointrange(data = betas, 
                  aes(x = var, y = mean, ymin = pc2.5, ymax = pc97.5),
                  size = 0.5) +
  geom_point(data = betas,
             aes(x = var, y = min(pc2.5) - 0.01, col = as.factor(dir)),
             shape = 8) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_y_continuous(expression(paste(Delta, " forbs prop. cover"))) +
  scale_x_discrete(limits = rev(levels(betas$var)), labels = rev(beta.labs2)) +
  scale_color_manual(values = c("forestgreen"),
                     na.value = "transparent") +
  coord_flip() +
  theme_bw(base_size = 14) +
  theme(axis.title.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  guides(color = "none")
fig_2a

beta.labs.ints <- c("POSE:high", "POFE:high", "VUMI:high", "ELEL:high",
                    "POSE:fall", "POFE:fall", "VUMI:fall", "ELEL:fall",
                    "POSE:spring", "POFE:spring", "VUMI:spring", "ELEL:spring",
                    "POSE:coated", "POFE:coated", "VUMI:coated", "ELEL:coated",
                    "high:fall", "high:spring", "high:coated", "fall:coated", "spring:coated")
beta.int.ind <- grep("diff_Beta", row.names(sum.out))
beta.ints <- sum.out[beta.int.ind,]
beta.ints$var <- factor(beta.ints$var, levels = row.names(beta.ints))
str(beta.ints)
fig_2b <- ggplot() +
  geom_pointrange(data = beta.ints, 
                  aes(x = var, y = mean, ymin = pc2.5, ymax = pc97.5),
                  size = 0.5) +
  geom_point(data = beta.ints,
             aes(x = var, y = min(pc2.5) - .01, col = as.factor(dir)),
             shape = 8) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_y_continuous(expression(paste(Delta, " forbs prop. cover"))) +
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


# Replicated data summary and fit
sum.rep <- coda.fast(coda.rep, OpenBUGS = FALSE)
z.rep <- sum.rep[grep("y.0.rep", row.names(sum.rep)),]
cont.rep <- sum.rep[grep("y.c.rep", row.names(sum.rep)),]

#align
y.temp <- with(dat, ifelse(forbs == 1 | forbs == 0, forbs, NA))
y.0 <- ifelse(is.na(y.temp), 0, 1)
which.0 <- which(y.0 == 1)
which.cont <- which(y.0 == 0)


fit <- rbind.data.frame(cbind(dat[which.0, ],
                              mean = z.rep$mean[which.0],
                              median = z.rep$median[which.0],
                              lower = z.rep$pc2.5[which.0],
                              upper = z.rep$pc97.5[which.0]),
                        cbind(dat[which.cont, ],
                              mean = cont.rep$mean,
                              median = cont.rep$median,
                              lower = cont.rep$pc2.5,
                              upper = cont.rep$pc97.5)
)


fit.model <- lm(mean ~ forbs, data = fit)
summary(fit.model)

ggplot(fit, aes(x = forbs)) +
  geom_abline(slope = 1, intercept = 0, col = "red", lty = 2) +
  # geom_errorbar(aes(ymin = lower, ymax = upper)) +
  geom_point(aes(y = mean))

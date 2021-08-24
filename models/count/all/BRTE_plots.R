# Plot model outputs

# Read in raw data
load("../../../cleaned_data/count_all.Rdata") # count_all
dat <- count_all %>%
  filter(quadrat < 10000)

# Load coda and coda.re
load(file = "coda/coda.Rdata") # coda.out
load(file = "coda/coda_rep.Rdata") # coda.rep


# summarize
sum.out <- coda.fast(coda.out, OpenBUGS = FALSE)
sum.out$var <- row.names(sum.out)
sum.out$sig <- ifelse(sum.out$pc2.5*sum.out$pc97.5 > 0, TRUE, FALSE)
sum.out$dir <- ifelse(sum.out$sig == FALSE, NA, 
                      ifelse(sum.out$sig == TRUE & sum.out$mean > 0, "pos", "neg"))


beta.labs <- c("fall", "spring", "herbicide", "greenstrip", 
               "fall:herbicide", "spring:herbicide", "fall:greenstrip", "spring:greenstrip")
beta.ind <- grep("beta", row.names(sum.out))
betas <- sum.out[beta.ind,]
betas$var <- factor(betas$var, levels = row.names(betas))
str(betas)
ggplot() +
  geom_pointrange(data = betas, 
                  aes(x = var, y = mean, ymin = pc2.5, ymax = pc97.5)) +
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

# Random effects
# labs <- c()
# for(i in 1:9){
#   labs[i] <- paste0("Plot ", i)
# }
# prob.eps <- grep("eps.star", row.names(sum.out))
# df <- data.frame(sum.out[prob.eps,], block = rep(1:3, each = 3))
# ggplot(df, aes(x = var, y = mean)) +
#   geom_pointrange(aes(ymin = pc2.5, ymax = pc97.5, color = as.factor(block))) +
#   geom_hline(yintercept = 0, color = "red", lty = 2) +
#   scale_x_discrete(labels = labs)

labs <- c("Block 1", "Block 2", "Block 3")
prob.eps <- grep("eps.star", row.names(sum.out))
ggplot(sum.out[prob.eps,], aes(x = var, y = mean)) +
  geom_pointrange(aes(ymin = pc2.5, ymax = pc97.5)) +
  geom_hline(yintercept = 0, color = "red", lty = 2) +
  scale_x_discrete(labels = labs)

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
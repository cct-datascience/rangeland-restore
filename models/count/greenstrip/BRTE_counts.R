# Cheatgrass counts modeled as Poisson distribution

library(rjags)
load.module('dic')
load.module('glm')
library(mcmcplots)
library(postjags)
library(ggplot2)
library(dplyr)

# Read in data
load("../../../cleaned_data/count_greenstrip.Rdata") # count_greenstrip
dat <- count_greenstrip
str(dat)

# Quadrat sizes
dat <- dat %>%
  filter(quadrat < 10000)
table(dat$quadrat)

# Plot
dat %>%
  ggplot(aes(x = grazing, y = BRTE/quadrat)) +
  geom_jitter(aes(color = block))

dat %>%
  ggplot(aes(x = spatial, y = BRTE/quadrat)) +
  geom_jitter(aes(color = block)) +
  facet_wrap(~grazing, ncol = 3)

dat %>%
  ggplot(aes(x = seed_rate, y = BRTE/quadrat)) +
  geom_jitter(aes(color = block)) +
  facet_wrap(~grazing, ncol = 3)

dat %>%
  ggplot(aes(x = seed_coat, y = BRTE/quadrat)) +
  geom_jitter(aes(color = block)) +
  facet_wrap(~grazing, ncol = 3)


# model matrix
X <- model.matrix( ~ (spatial + seed_rate + seed_coat + grazing)^2, data = dat) 
colnames(X)

log(sd(tapply(dat$BRTE, dat$block, FUN = mean)))

# link paddocks to block
# link <- dat %>%
#   group_by(block) %>%
#   summarize(paddock = unique(paddock))

# Assemble model inputs
datlist <- list(counts = dat$BRTE,
                area = dat$quadrat/100, # convert to square decimeters
                N = nrow(dat),
                mix = as.numeric(X[,2]),
                high = as.numeric(X[,3]),
                coated = as.numeric(X[,4]),
                fall = as.numeric(X[,5]),
                spring = as.numeric(X[,6]),
                mix_high = as.numeric(X[,7]),
                mix_coated = as.numeric(X[,8]),
                mix_fall = as.numeric(X[,9]),
                mix_spring = as.numeric(X[,10]),
                high_coated = as.numeric(X[,11]),
                high_fall = as.numeric(X[,12]),
                high_spring = as.numeric(X[,13]),
                coated_fall = as.numeric(X[,14]),
                coated_spring = as.numeric(X[,15]),
                nL = ncol(X) - 1, # number of levels
                # pad = as.numeric(dat$paddock),
                # Np = length(unique(dat$paddock)),
                block = as.numeric(dat$block),
                Nb = length(unique(dat$block)),
                Ab = 5) # stand deviation among blocks
str(datlist)

# likely intercept value
base <- dat %>%
  filter(grazing == "ungrazed",
         spatial == "mono",
         seed_rate == "low", 
         seed_coat == "UC")
hist(base$BRTE, breaks = 30)
table(base$quadrat)
log(mean(base$BRTE))

# initials
inits <- function(){
  list(alpha = runif(1, 0, 5),
       beta = rnorm(ncol(X) - 1, 0, 10),
       tau.Eps = runif(1, 0, 3)
       )
}
initslist <- list(inits(), inits(), inits())


# Or, use previous starting values + set seed
load("inits/inits.Rdata")# saved.state, second element is inits
initslist <- list(append(saved.state[[2]][[1]], list(.RNG.name = array("base::Marsaglia-Multicarry"), .RNG.seed = array(13))),
                  append(saved.state[[2]][[2]], list(.RNG.name = array("base::Wichmann-Hill"), .RNG.seed = array(89))),
                  append(saved.state[[2]][[3]], list(.RNG.name = array("base::Super-Duper"), .RNG.seed = array(18))))

# model
jm <- jags.model(file = "BRTE_counts_Poisson.jags",
                 inits = initslist,
                 n.chains = 3,
                 data = datlist)
# update(jm, 10000)

# params to monitor
params <- c("deviance", "Dsum", # evaluate fit
            "alpha", "beta",  # parameters
            "tau.Eps", "sig.eps", # precision/variance terms
            "alpha.star", "eps.star") #identifiable intercept and random effects

coda.out <- coda.samples(jm, variable.names = params,
                         n.iter = 15000, thin = 5)

# plot chains
mcmcplot(coda.out, parms = c("deviance", "Dsum", "alpha.star", 
                             "beta", "eps.star",
                             "sig.eps"))

# dic samples
dic.out <- dic.samples(jm, n.iter = 5000)
dic.out

# convergence?
gel <- gelman.diag(coda.out, multivariate = FALSE)
gel

# If not converged, restart model from final iterations
# newinits <-  initfind(coda.out) 
# newinits[[1]]
# saved.state <- removevars(newinits, variables = c(1, 3, 5:6))
# saved.state[[1]]
# save(saved.state, file = "inits/inits.Rdata")

save(coda.out, file = "coda/coda.Rdata")

# summarize
sum.out <- coda.fast(coda.out, OpenBUGS = FALSE)
sum.out$var <- row.names(sum.out)
sum.out$sig <- ifelse(sum.out$pc2.5*sum.out$pc97.5 > 0, TRUE, FALSE)
sum.out$dir <- ifelse(sum.out$sig == FALSE, NA, 
                      ifelse(sum.out$sig == TRUE & sum.out$mean > 0, "pos", "neg"))


beta.labs <- c("mix", "high", "coated", "fall", "spring",
               "mix:high", "mix:coated", "mix:fall", "mix:spring",
               "high:coated", "high:fall", "high:spring", 
               "coated:fall", "coated:spring")
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

# Model fit
params <- c("counts.rep") #monitor replicated data
coda.rep <- coda.samples(jm, variable.names = params,
                         n.iter = 5000, thin = 1)
save(coda.rep, file = "coda/coda_rep.Rdata")

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

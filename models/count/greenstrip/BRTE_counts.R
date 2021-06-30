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
# dat2 <- dat %>%
#   filter(quadrat < 10000)
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
X <- model.matrix( ~ grazing + spatial + seed_rate + seed_coat, data = dat) 

# link paddocks to block
link <- dat %>%
  group_by(block) %>%
  summarize(paddock = unique(paddock))

# standard deviation among paddocks and blocks
sd(tapply(dat$BRTE, dat$paddock, FUN = mean))
sd(tapply(dat$BRTE, dat$block, FUN = mean))

# Assemble model inputs
datlist <- list(counts = dat$BRTE,
                area = dat$quadrat/100, # convert to square decimeters
                N = nrow(dat),
                fall = as.numeric(X[,2]),
                spring = as.numeric(X[,3]),
                mix = as.numeric(X[,4]),
                high = as.numeric(X[,5]),
                coated = as.numeric(X[,6]),
                nL = ncol(X) - 1, # number of levels
                pad = as.numeric(dat$paddock),
                Np = length(unique(dat$paddock)),
                block = as.numeric(link$block),
                Nb = length(unique(dat$block)),
                Ab = 10) # stand deviation among paddocks and blocks
str(datlist)

# initials
inits <- function(){
  list(alpha = runif(1, -3, 0),
       beta = runif(ncol(X) - 1, -3, 0),
       psi = runif(1, 0, 1),
       tau.Eps = runif(1, 0, 1),
       tau.Gam = runif(1, 0, 1))
}
initslist <- list(inits(), inits(), inits())

# model
jm <- jags.model(file = "BRTE_counts.jags",
                 inits = initslist,
                 n.chains = 3,
                 data = datlist)
update(jm, 10000)

# params to monitor
params <- c("deviance", "Dsum", # evaluate fit
            "alpha", "beta", "psi", # parameters
            "diff", # derived parameters
            "tau.Eps", "tau.Gam", "sig.eps", "sig.gam", # precision/variance terms
            "alpha.star", "eps.star", "gam.star",#identifiable intercept and random effects
            "eps", "gam", "eps.avg") 

coda.out <- coda.samples(jm, variable.names = params,
                         n.iter = 5000, thin = 1)

# plot chains
mcmcplot(coda.out, parms = c("deviance", "Dsum", "beta", "psi",
                             "alpha.star",  "eps.star", "gam.star",
                             "sig.eps", "sig.gam"))
mcmcplot(coda.out, parms = c("diff"))

# dic samples
dic.out <- dic.samples(jm, n.iter = 5000)
dic.out

# convergence?
gel <- gelman.diag(coda.out, multivariate = FALSE)
gel

# If not converged, restart model from final iterations
newinits <-  initfind(coda.out) 
newinits[[1]]
saved.state <- removevars(newinits, variables = c(1, 3, 5:8, 10:11))
its <- list(saved.state[[2]][[2]],
            saved.state[[2]][[2]],
            saved.state[[2]][[2]])

coda.out[[1]][5000,1]
coda.out[[2]][5000,1]
coda.out[[3]][5000,1]
coda.out[[1]][5000,8]
coda.out[[2]][5000,8]
coda.out[[3]][5000,8]
# mc <- coda.out[[3]]
# means <- round(colMeans(mc),2)

# do random effects add up correctly?
inds.eps <- grep("eps.star", colnames(mc))
head(rowSums(mc[,inds.eps[1:3]])) # paddocks within block 1
head(rowSums(mc[,inds.eps[4:6]])) # paddocks within block 2
head(rowSums(mc[,inds.eps[7:9]])) # paddocks within block 3
inds.gam <-grep("gam.star", colnames(mc))
head(rowSums(mc[,inds.gam])) # blocks

# plot trace plots
traplot(coda.out, parms = "alpha.star")
traplot(coda.out, parms = "eps")
traplot(coda.out, parms = "eps.star")
traplot(coda.out, parms = "gam")
traplot(coda.out, parms = "gam.star")
traplot(coda.out, parms = "diff")


# summarize
sum.out <- coda.fast(mcmc.list(mc), OpenBUGS = FALSE)
sum.out$var <- row.names(sum.out)




beta.ind <- grep("beta", row.names(sum.out))
ggplot(sum.out[beta.ind,], aes(x = var, y = mean)) +
  geom_pointrange(aes(ymin = pc2.5, ymax = pc97.5))

diff.ind <- grep("diff", row.names(sum.out))
ggplot(sum.out[diff.ind,], aes(x = var, y = mean)) +
  geom_pointrange(aes(ymin = pc2.5, ymax = pc97.5))

prob.ind <- grep("prob", row.names(sum.out))
ggplot(sum.out[prob.ind,], aes(x = var, y = mean)) +
  geom_pointrange(aes(ymin = pc2.5, ymax = pc97.5))

which(dat$quadrat > 100)

# Cheatgrass counts modeled as Poisson distribution

library(rjags)
load.module('dic')
library(mcmcplots)
library(postjags)
library(ggplot2)
library(dplyr)

# Read in data
load("../../../cleaned_data/count_all.Rdata") # count_all
dat <- count_all
str(dat)

# Remove 2 rows with largest quadrat size
dat2 <- dat %>%
  filter(quadrat < 10000)

# Plot
dat2 %>%
  ggplot(aes(x = fuelbreak, y = BRTE/quadrat)) +
  geom_jitter(aes(color = block)) +
  facet_wrap(~grazing, ncol = 3)

dat2 %>%
  filter(grazing == "fall") %>%
  ggplot(aes(x = fuelbreak, y = BRTE/quadrat)) +
  geom_jitter(aes(color = block))


# model matrix
X <- model.matrix( ~ grazing + fuelbreak, data = dat2) 

# link paddocks to block
link <- dat2 %>%
  group_by(block) %>%
  summarize(paddock = unique(paddock))

# standard deviation among paddocks and blocks
sd(tapply(dat2$BRTE, dat2$paddock, FUN = mean))
sd(tapply(dat2$BRTE, dat2$block, FUN = mean))

# Assemble model inputs
datlist <- list(counts = dat2$BRTE,
                area = dat2$quadrat/100, # convert to square decimeters
                N = nrow(dat2),
                fall = X[,2],
                spring = X[,3],
                herbicide = X[,4],
                greenstrip = X[,5],
                nL = ncol(X) - 1, # number of levels
                pad = as.numeric(dat2$paddock),
                Np = length(unique(dat2$paddock)),
                block = as.numeric(link$block),
                Nb = length(unique(dat2$block)),
                Ab = 10) # stand deviation among paddocks and blocks


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
                 inits = its,
                 n.chains = 3,
                 data = datlist)
update(jm, 10000)

# params to monitor
params <- c("deviance", "Dsum", # evaluate fit
            "alpha", "beta", "psi", # parameters
            "prob", "diff", # derived parameters
            "tau.Eps", "tau.Gam", "sig.eps", "sig.gam", # precision/variance terms
            "alpha.star", "eps.star", "gam.star") #identifiable intercept and random effects

coda.out <- coda.samples(jm, variable.names = params,
                         n.iter = 5000, thin = 1)

# plot chains
mcmcplot(coda.out, parms = c("deviance", "Dsum", "beta", "psi",
                             "alpha.star",  "eps.star", "gam.star",
                             "sig.eps", "sig.gam"))
mcmcplot(coda.out, parms = c("prob"))

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
save(saved.state, file = "inits/inits.Rdata")

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
grep("eps.star", colnames(mc))
head(rowSums(mc[,13:15])) # paddocks within block 1
head(rowSums(mc[,16:18])) # paddocks within block 2
head(rowSums(mc[,19:21])) # paddocks within block 3
grep("gam.star", colnames(mc))
head(rowSums(mc[,22:24])) # blocks

# plot trace plots
traplot(coda.out, parms = "alpha.star")
traplot(coda.out, parms = "eps.star")
traplot(coda.out, parms = "gam.star")
traplot(coda.out, parms = "prob")


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

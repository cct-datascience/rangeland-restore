# Cheatgrass counts modeled as Poisson distribution
# Ultimately using BRTE_counts_Poisson.jags file

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
dat <- dat %>%
  filter(quadrat < 10000)


# Plot
dat %>%
  ggplot(aes(x = fuelbreak, y = BRTE/quadrat)) +
  geom_jitter(aes(color = paddock, shape = block)) +
  facet_wrap(~grazing, ncol = 3)

dat %>%
  filter(grazing == "fall") %>%
  ggplot(aes(x = fuelbreak, y = BRTE/quadrat)) +
  geom_jitter(aes(color = block))


# model matrix
X <- model.matrix( ~ grazing * fuelbreak, data = dat) 
colnames(X)

# standard deviation among paddocks and blocks
log(sd(tapply(dat$BRTE, dat$paddock, FUN = mean)))
log(sd(tapply(dat$BRTE, dat$block, FUN = mean)))

# link paddocks to block
# link <- dat %>%
#   group_by(block) %>%
#   summarize(paddock = unique(paddock))

# Assemble model inputs
datlist <- list(counts = dat$BRTE,
                area = dat$quadrat/100, # convert to square decimeters
                N = nrow(dat),
                fall = X[,2],
                spring = X[,3],
                herbicide = X[,4],
                greenstrip = X[,5],
                fall_herbicide = X[,6],
                spring_herbicide = X[,7],
                fall_greenstrip = X[,8],
                spring_greenstrip = X[,9],
                nL = ncol(X) - 1, # number of offset levels
                # pad = as.numeric(dat$paddock),
                # Np = length(unique(dat$paddock)),
                block = as.numeric(dat$block),
                Nb = length(unique(dat$block)),
                Ab = 5) # stand deviation among paddocks and blocks

# likely intercept value
base <- dat %>%
  filter(grazing == "ungrazed",
         fuelbreak == "control")
hist(base$BRTE, breaks = 30)
table(base$quadrat)
log(median(base$BRTE))

# generate random initials
inits <- function(){
  list(alpha = rnorm(1, 0, 10),
       beta = rnorm(ncol(X) - 1, 0, 10),
       tau.Eps = runif(1, 0, 1))
}
initslist <- list(inits(), inits(), inits())

# Or, use previous starting values + set seed
load("inits/inits.Rdata")# saved.state, second element is inits
initslist <- list(append(saved.state[[2]][[1]], list(.RNG.name = array("base::Super-Duper"), .RNG.seed = array(13))),
                  append(saved.state[[2]][[2]], list(.RNG.name = array("base::Wichmann-Hill"), .RNG.seed = array(89))),
                  append(saved.state[[2]][[3]], list(.RNG.name = array("base::Mersenne-Twister"), .RNG.seed = array(18))))

# model
jm <- jags.model(file = "BRTE_counts_Poisson.jags",
                 inits = initslist,
                 n.chains = 3,
                 data = datlist)
# update(jm, 10000)

# params to monitor
params <- c("deviance", "Dsum", # evaluate fit
            "alpha", "beta", # parameters
            "tau.Eps", "sig.eps", # precision/variance terms
            "alpha.star", "eps.star")#identifiable intercept and random effects

coda.out <- coda.samples(jm, variable.names = params,
                         n.iter = 15000, thin = 5)

# plot chains
mcmcplot(coda.out, parms = c("deviance", "Dsum", "beta",
                             "alpha.star",  "eps.star", "sig.eps"))

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


# Model fit
params <- c("counts.rep") #monitor replicated data
coda.rep <- coda.samples(jm, variable.names = params,
                         n.iter = 5000, thin = 1)
save(coda.rep, file = "coda/coda_rep.Rdata")

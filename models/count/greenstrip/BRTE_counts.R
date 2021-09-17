# Cheatgrass counts modeled as zero-inflated Poisson distribution

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
jm <- jags.model(file = "BRTE_counts_ziPoisson.jags",
                 inits = initslist,
                 n.chains = 3,
                 data = datlist)
# update(jm, 10000)

# params to monitor
params <- c("deviance", "Dsum", # evaluate fit
            "alpha", "beta", "psi",  # parameters
            "tau.Eps", "sig.eps", # precision/variance terms
            "alpha.star", "eps.star", # identifiable intercept and random effects
            "int_Beta", "Diff_Beta", "diff_Beta", # monitored interaction terms
            "m.mono.low.uncoated.ungrazed", "m.mono.low.uncoated.fall", "m.mono.low.uncoated.spring", 
            "m.mono.low.coated.ungrazed", "m.mono.low.coated.fall", "m.mono.low.coated.spring",
            "m.mono.high.uncoated.ungrazed", "m.mono.high.uncoated.fall", "m.mono.high.uncoated.spring",
            "m.mono.high.coated.ungrazed", "m.mono.high.coated.fall", "m.mono.high.coated.spring",
            "m.mix.low.uncoated.ungrazed", "m.mix.low.uncoated.fall", "m.mix.low.uncoated.spring",
            "m.mix.low.coated.ungrazed", "m.mix.low.coated.fall", "m.mix.low.coated.spring",
            "m.mix.high.uncoated.ungrazed", "m.mix.high.uncoated.fall", "m.mix.high.uncoated.spring",
            "m.mix.high.coated.ungrazed", "m.mix.high.coated.fall", "m.mix.high.coated.spring") 

coda.out <- coda.samples(jm, variable.names = params,
                         n.iter = 15000, thin = 5)

# plot chains
mcmcplot(coda.out, parms = c("deviance", "Dsum", "alpha.star", 
                             "beta", "psi", "eps.star",
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
# saved.state <- removevars(newinits, variables = c(1:2, 4, 6:10))
# saved.state[[1]]
# save(saved.state, file = "inits/inits.Rdata")

save(coda.out, file = "coda/coda_zip.Rdata")

# Model fit
params <- c("counts.rep") #monitor replicated data
coda.rep <- coda.samples(jm, variable.names = params,
                         n.iter = 15000, thin = 5)
save(coda.rep, file = "coda/coda_zip_rep.Rdata")

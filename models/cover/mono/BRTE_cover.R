# BRTE cover modeled as zero-or-1 inflated beta
# cover proportions collected on 1x1 m quadrats

library(rjags)
load.module('dic')
library(mcmcplots)
library(postjags)
library(ggplot2)
library(dplyr)

# logit and antilogit functions
logit <- function(x) {
  log(x/(1-x))
}
ilogit <- function(x){
  exp(x) / (1 + exp(x))
}

# read in data
load("../../../cleaned_data/cover_mono.Rdata") # cover_mono
# convert to proportions
dat <- cover_mono %>%
  mutate(BRTE = BRTE/100,
         intro_forbs = intro_forbs/100,
         native_grass = native_grass/100,
         native_forbs = native_forbs/100)

# Relevel species based on fig. 6b from Porensky et al. 2018
dat$species <- factor(dat$species, levels = c("ELTR", "POSE", "POFE", "VUMI", "ELEL"))

str(dat)

# plot

ggplot(dat, aes(x = species, y = BRTE, color = seed_rate)) +
  geom_point() +
  facet_wrap(~grazing)

# model matrix
X <- model.matrix( ~ (species + seed_rate + grazing + seed_coat)^2, data = dat) 
colnames(X)

# split the data into discrete and continuous components
y.temp <- with(dat, ifelse(BRTE == 1 | BRTE == 0, BRTE, NA))
y.discrete <- ifelse(is.na(y.temp), 0, 1)

# group discrete response + predictors
y.d <- y.temp[!is.na(y.temp)]
x.d <- X[y.discrete == 1,]
n.discrete <- length(y.d)
block.d <- as.numeric(dat$block)[y.discrete == 1]

# group continuous response + predictors
which.cont <- which(y.discrete == 0)
y.c <- dat$BRTE[which.cont]
x.c <- X[which.cont,]
n.cont <- length(y.c)
block.c <- as.numeric(dat$block)[which.cont]

# Assemble model inputs
datlist <- list(N = nrow(dat),
                y.discrete = y.discrete,
                n.discrete = n.discrete,
                y.d = y.d, 
                POSE = x.d[,2],
                POFE = x.d[,3],
                VUMI = x.d[,4],
                ELEL = x.d[,5],
                high = x.d[,6],
                fall = x.d[,7], 
                spring = x.d[,8], 
                coated = x.d[,9],
                block.d = block.d,
                n.cont = n.cont,
                y.c = y.c,
                POSE2 = x.c[,2],
                POFE2 = x.c[,3],
                VUMI2 = x.c[,4],
                ELEL2 = x.c[,5],
                high2 = x.c[,6],
                fall2 = x.c[,7], 
                spring2 = x.c[,8], 
                coated2 = x.c[,9],
                block.c = block.c,
                Nb = length(unique(dat$block)),
                nL = ncol(X) - 1)

# likely intercept value
base <- dat %>%
  filter(grazing == "ungrazed",
         species == "ELTR",
         seed_rate == "low", 
         seed_coat == "UC")
hist(base$BRTE, breaks = 30)
logit(median(base$BRTE))

# generate random initials
inits <- function(){
  list(alpha = rnorm(1, 0, 10),
       beta = rnorm(ncol(X) - 1, 0, 10),
       tau = runif(1, 0, 1),
       tau.eps = runif(1, 0, 1),
       psi = runif(1, 0, 1))
}
initslist <- list(inits(), inits(), inits())

# Or, use previous starting values + set seed
load("inits/inits.Rdata")# saved.state, second element is inits
initslist <- list(append(saved.state[[2]][[1]], list(.RNG.name = array("base::Marsaglia-Multicarry"), .RNG.seed = array(14))),
                  append(saved.state[[2]][[2]], list(.RNG.name = array("base::Super-Duper"), .RNG.seed = array(89))),
                  append(saved.state[[2]][[3]], list(.RNG.name = array("base::Mersenne-Twister"), .RNG.seed = array(18))))

# model
jm <- jags.model(file = "BRTE_cover_zoib.jags",
                 inits = initslist,
                 n.chains = 3,
                 data = datlist)
# update(jm, 10000)

# params to monitor
params <- c("deviance", # evaluate fit
            "psi", "alpha", "beta", # parameters
            "tau", "sig", "tau.eps", "sig.eps", # precision/variance terms
            "alpha.star", "eps.star", # identifiable intercept and random effects
            "int_Beta", "Diff_Beta", "diff_Beta") # monitored main and two-way treatment effects

coda.out <- coda.samples(jm, variable.names = params,
                         n.iter = 15000, thin = 5)

# plot chains
mcmcplot(coda.out, parms = c("deviance", "psi", "beta",
                             "alpha.star",  "eps.star", 
                             "sig", "sig.eps"))
caterplot(coda.out, parms = "eps.star", reorder = FALSE)
caterplot(coda.out, parms = "beta", reorder = FALSE)
caterplot(coda.out, parms = "Diff_Beta", reorder = FALSE)
caterplot(coda.out, parms = "diff_Beta", reorder = FALSE)

# dic samples
dic.out <- dic.samples(jm, n.iter = 5000)
dic.out

# convergence?
gel <- gelman.diag(coda.out, multivariate = FALSE)
gel

# If not converged, restart model from final iterations
# newinits <-  initfind(coda.out)
# newinits[[1]]
# saved.state <- removevars(newinits, variables = c(1, 3, 5:7, 9:10))
# saved.state[[1]]
# save(saved.state, file = "inits/inits.Rdata")

save(coda.out, file = "coda/coda.Rdata")


# Model fit
params <- c("y.discrete.rep", "y.d.rep", "y.c.rep") #monitor replicated data
coda.rep <- coda.samples(jm, variable.names = params,
                         n.iter = 15000, thin = 5)
save(coda.rep, file = "coda/coda_rep.Rdata")


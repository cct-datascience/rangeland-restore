# Greenstrip forb cover modeled as 
# beta regression, no 0 or 1's present
# cover proportions collected on 1x1 m quadrats

library(rjags)
load.module('dic')
library(mcmcplots)
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
load("../../../../cleaned_data/cover_greenstrip.Rdata") # cover_greenstrip
# convert to proportions
dat <- cover_greenstrip %>%
  mutate(BRTE = BRTE/100,
         intro_forbs = intro_forbs/100,
         native_grass = native_grass/100,
         native_forbs = native_forbs/100,
         forbs = intro_forbs + native_forbs,
         spatial = factor(spatial, levels = c("mix", "mono")))
str(dat)

# check data for patterns
ggplot(dat, aes(x = spatial, y = forbs, color = block)) +
  geom_point() +
  facet_grid(rows = vars(grazing),
             cols = vars(block))

# model matrix
X <- model.matrix( ~ (spatial + seed_rate + seed_coat + grazing)^2, data = dat) 
colnames(X)

# split the data into discrete and continuous components
y.temp <- with(dat, ifelse(forbs == 1 | forbs == 0, 
                           forbs, NA))
y.discrete <- ifelse(is.na(y.temp), 0, 1)
sum(y.discrete)
# no discrete observations in this dataset
# only continuous response + predictors

# Assemble model inputs
datlist <- list(N = nrow(dat),
                y = dat$forbs,
                mono = X[,2],
                high = X[,3],
                coated = X[,4],
                fall = X[,5],
                spring = X[,6],
                mono_high = X[,7],
                mono_coated = X[,8],
                mono_fall = X[,9],
                mono_spring = X[,10],
                high_coated = X[,11],
                high_fall = X[,12],
                high_spring = X[,13],
                coated_fall = X[,14],
                coated_spring = X[,15],
                block = as.numeric(dat$block),
                Nb = length(unique(dat$block)),
                nL = ncol(X) - 1)

# likely intercept value
base <- dat %>%
  filter(grazing == "ungrazed",
         spatial == "mix",
         seed_rate == "low", 
         seed_coat == "UC")
logit(median(base$forbs))

# likely sig among blocks
sd(tapply(dat$forbs, dat$block, mean))

# generate random initials
inits <- function(){
  list(alpha = rnorm(1, 0, 10),
       beta = rnorm(ncol(X) - 1, 0, 10),
       tau = runif(1, 0, 1),
       tau.eps = runif(1, 0, 0.1))
}
initslist <- list(inits(), inits(), inits())

# Or, use previous starting values + set seed
load("inits/inits.Rdata")# saved.state, second element is inits
initslist <- list(append(saved.state[[2]][[1]], list(.RNG.name = array("base::Super-Duper"), .RNG.seed = array(13))),
                  append(saved.state[[2]][[2]], list(.RNG.name = array("base::Wichmann-Hill"), .RNG.seed = array(89))),
                  append(saved.state[[2]][[3]], list(.RNG.name = array("base::Mersenne-Twister"), .RNG.seed = array(18))))

# model
jm <- jags.model(file = "Forbs_cover_beta.jags",
                 inits = initslist,
                 n.chains = 3,
                 data = datlist)
# update(jm, 10000)

# params to monitor
params <- c("deviance", # evaluate fit
            "alpha", "beta", # parameters
            "tau", "sig", "tau.eps", "sig.eps", # precision/variance terms
            "alpha.star", "eps.star", # identifiable intercept and random effects
            "int_Beta", "Diff_Beta", "diff_Beta", # monitored main and two-way treatment effects
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
mcmcplot(coda.out, parms = c("deviance", "beta",
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
# saved.state <- removevars(newinits, variables = c(1, 3, 5:33))
# saved.state[[1]]
# save(saved.state, file = "inits/inits.Rdata")

save(coda.out, file = "coda/coda.Rdata")


# Model fit
params <- c("y.rep") #monitor replicated data
coda.rep <- coda.samples(jm, variable.names = params,
                         n.iter = 15000, thin = 5)
save(coda.rep, file = "coda/coda_rep.Rdata")

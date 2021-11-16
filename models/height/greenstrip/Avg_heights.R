# Average heights per plot
# modeled as a normal distribution

library(rjags)
load.module('dic')
library(mcmcplots)
library(postjags)
library(ggplot2)
library(dplyr)

# Read in data
load("../../../cleaned_data/cover_greenstrip.Rdata") # cover_greenstrip
dat <- cover_greenstrip%>%
  arrange(block) %>%
  mutate(spatial = factor(spatial, levels = c("mix", "mono")))
str(dat)


# Plot
hist(dat$height)
summary(dat$height)

dat %>%
  ggplot(aes(x = grazing, y = height)) +
  geom_jitter(aes(color = block)) +
  facet_grid(cols = vars(seed_rate),
             rows = vars(seed_coat))


# model matrix
X <- model.matrix( ~ (spatial + seed_rate + seed_coat + grazing)^2, data = dat) 
colnames(X)

# standard deviation among paddocks and blocks
log(sd(tapply(dat$BRTE, dat$block, FUN = mean)))

# Assemble model inputs
datlist <- list(height = dat$height,
                N = nrow(dat),
                mono = as.numeric(X[,2]),
                high = as.numeric(X[,3]),
                coated = as.numeric(X[,4]),
                fall = as.numeric(X[,5]),
                spring = as.numeric(X[,6]),
                mono_high = as.numeric(X[,7]),
                mono_coated = as.numeric(X[,8]),
                mono_fall = as.numeric(X[,9]),
                mono_spring = as.numeric(X[,10]),
                high_coated = as.numeric(X[,11]),
                high_fall = as.numeric(X[,12]),
                high_spring = as.numeric(X[,13]),
                coated_fall = as.numeric(X[,14]),
                coated_spring = as.numeric(X[,15]),
                nL = ncol(X) - 1, # number of offset levels
                block = as.numeric(dat$block),
                Nb = length(unique(dat$block)),
                Ab = 5) # stand deviation among paddocks and blocks

# likely intercept value
base <- dat %>%
  filter(grazing == "ungrazed",
         spatial == "mix",
         seed_rate == "low", 
         seed_coat == "UC")
hist(base$height, breaks = 30)
median(base$height)

# generate random initials
inits <- function(){
  list(alpha = rnorm(1, 0, 10),
       beta = rnorm(ncol(X) - 1, 0, 10),
       tau.Eps = runif(1, 0, 1),
       tau = runif(1, 0, 1))
}
initslist <- list(inits(), inits(), inits())

# Or, use previous starting values + set seed
load("inits/inits.Rdata")# saved.state, second element is inits
initslist <- list(append(saved.state[[2]][[1]], list(.RNG.name = array("base::Super-Duper"), .RNG.seed = array(13))),
                  append(saved.state[[2]][[2]], list(.RNG.name = array("base::Wichmann-Hill"), .RNG.seed = array(89))),
                  append(saved.state[[2]][[3]], list(.RNG.name = array("base::Mersenne-Twister"), .RNG.seed = array(18))))

# model
jm <- jags.model(file = "Avg_heights_norm.jags",
                 inits = initslist,
                 n.chains = 3,
                 data = datlist)
# update(jm, 10000)

# params to monitor
params <- c("deviance", "Dsum", # evaluate fit
            "alpha", "beta", # parameters
            "tau.Eps", "sig.eps", "tau", "sig", # precision/variance terms
            "alpha.star", "eps.star", # identifiable intercept and random effects
            "int_Beta", # monitored interaction effects
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
mcmcplot(coda.out, parms = c("deviance", "Dsum", "beta",
                             "alpha.star",  "eps.star", "sig.eps", 
                             "sig"))

caterplot(coda.out, parms = "beta", reorder = FALSE)
caterplot(coda.out, parms = "eps.star", reorder = FALSE)

# dic samples
dic.out <- dic.samples(jm, n.iter = 5000)
dic.out

# convergence?
gel <- gelman.diag(coda.out, multivariate = FALSE)
gel

# If not converged, restart model from final iterations
# newinits <-  initfind(coda.out)
# newinits[[1]]
# saved.state <- removevars(newinits, variables = c(1:2, 4, 6:34))
# saved.state[[1]]
# save(saved.state, file = "inits/inits.Rdata")

save(coda.out, file = "coda/coda.Rdata")


# Model fit
params <- c("height.rep") #monitor replicated data
coda.rep <- coda.samples(jm, variable.names = params,
                         n.iter = 15000, thin = 5)
save(coda.rep, file = "coda/coda_rep.Rdata")

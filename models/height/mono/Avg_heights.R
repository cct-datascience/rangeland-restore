# Average heights per plot
# modeled as a normal distribution

library(rjags)
load.module('dic')
library(mcmcplots)
library(postjags)
library(ggplot2)
library(dplyr)

# Read in data
load("../../../cleaned_data/cover_mono.Rdata") # cover_mono
dat <- cover_mono %>%
  mutate(species = factor(species, levels = c("ELTR", "POSE", "POFE", "VUMI", "ELEL")))

# Plot
hist(dat$height)
summary(dat$height)

dat %>%
  ggplot(aes(x = species, y = height)) +
  geom_jitter(aes(color = grazing)) +
  facet_grid(cols = vars(seed_rate),
             rows = vars(seed_coat))


# model matrix
X <- model.matrix( ~ (species + seed_rate + grazing + seed_coat)^2, data = dat) 
colnames(X)

# standard deviation among paddocks and blocks
log(sd(tapply(dat$height, dat$block, FUN = mean)))

# Assemble model inputs
datlist <- list(height = dat$height,
                N = nrow(dat),
                POSE = X[,2],
                POFE = X[,3],
                VUMI = X[,4],
                ELEL = X[,5],
                high = X[,6],
                fall = X[,7], 
                spring = X[,8], 
                coated = X[,9],
                nL = ncol(X) - 1, # number of offset levels
                block = as.numeric(dat$block),
                Nb = length(unique(dat$block)),
                Ab = 5) # stand deviation among paddocks and blocks

# likely intercept value
base <- dat %>%
  filter(grazing == "ungrazed",
         species == "ELTR",
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
            "int_Beta", # monitored effect combinations
            "m.ELEL.low.ungrazed.uncoated", "m.ELEL.low.ungrazed.coated", "m.ELEL.low.fall.uncoated", "m.ELEL.low.fall.coated", "m.ELEL.low.spring.uncoated", "m.ELEL.low.spring.coated", 
            "m.ELEL.high.ungrazed.uncoated", "m.ELEL.high.ungrazed.coated", "m.ELEL.high.fall.uncoated", "m.ELEL.high.fall.coated", "m.ELEL.high.spring.uncoated", "m.ELEL.high.spring.coated", 
            "m.VUMI.low.ungrazed.uncoated", "m.VUMI.low.ungrazed.coated", "m.VUMI.low.fall.uncoated", "m.VUMI.low.fall.coated", "m.VUMI.low.spring.uncoated", "m.VUMI.low.spring.coated", 
            "m.VUMI.high.ungrazed.uncoated", "m.VUMI.high.ungrazed.coated", "m.VUMI.high.fall.uncoated", "m.VUMI.high.fall.coated" , "m.VUMI.high.spring.uncoated", "m.VUMI.high.spring.coated", 
            "m.POFE.low.ungrazed.uncoated", "m.POFE.low.ungrazed.coated", "m.POFE.low.fall.uncoated", "m.POFE.low.fall.coated", "m.POFE.low.spring.uncoated", "m.POFE.low.spring.coated", 
            "m.POFE.high.ungrazed.uncoated", "m.POFE.high.ungrazed.coated", "m.POFE.high.fall.uncoated", "m.POFE.high.fall.coated", "m.POFE.high.spring.uncoated", "m.POFE.high.spring.coated", 
            "m.POSE.low.ungrazed.uncoated", "m.POSE.low.ungrazed.coated", "m.POSE.low.fall.uncoated", "m.POSE.low.fall.coated", "m.POSE.low.spring.uncoated", "m.POSE.low.spring.coated", 
            "m.POSE.high.ungrazed.uncoated", "m.POSE.high.ungrazed.coated", "m.POSE.high.fall.uncoated", "m.POSE.high.fall.coated", "m.POSE.high.spring.uncoated", "m.POSE.high.spring.coated", 
            "m.ELTR.low.ungrazed.uncoated", "m.ELTR.low.ungrazed.coated", "m.ELTR.low.fall.uncoated", "m.ELTR.low.fall.coated", "m.ELTR.low.spring.uncoated", "m.ELTR.low.spring.coated", 
            "m.ELTR.high.ungrazed.uncoated", "m.ELTR.high.ungrazed.coated", "m.ELTR.high.fall.uncoated", "m.ELTR.high.fall.coated", "m.ELTR.high.spring.uncoated", "m.ELTR.high.spring.coated"
            )
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
# saved.state <- removevars(newinits, variables = c(1, 3, 5:68))
# saved.state[[1]]
# save(saved.state, file = "inits/inits.Rdata")

save(coda.out, file = "coda/coda.Rdata")


# Model fit
params <- c("height.rep") #monitor replicated data
coda.rep <- coda.samples(jm, variable.names = params,
                         n.iter = 15000, thin = 5)
save(coda.rep, file = "coda/coda_rep.Rdata")

# Cheatgrass counts modeled as Poisson distribution
# ANOVA 4 factors, (5, 2, 3, 2 levels), two-way interactions
# only block RE, no paddocks nested within blocks

library(rjags)
load.module('dic')
load.module('glm')
library(mcmcplots)
library(postjags)
library(ggplot2)
library(dplyr)

# Read in data
load("../../../cleaned_data/count_mono.Rdata") # count_mono

# Organize: remove largest quadrat and relevel species based on fig. 6b from Porensky et al. 2018
dat <- count_mono %>%
  filter(quadrat < 10000) %>%
  mutate(species = factor(species, levels = c("ELTR", "POSE", "POFE", "VUMI", "ELEL"))) %>%
  arrange(block)
str(dat)

range(which(dat$block == "one"))
range(which(dat$block == "two"))
range(which(dat$block == "three"))


# Plot
dat %>%
  ggplot(aes(x = grazing, y = BRTE/quadrat)) +
  geom_jitter(aes(color = block))

dat %>%
  ggplot(aes(x = species, y = BRTE/quadrat)) +
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
X <- model.matrix( ~ (species + seed_rate + grazing + seed_coat)^2, data = dat) 
colnames(X)

# link paddocks to block
# link <- dat %>%
#   group_by(block) %>%
#   summarize(paddock = unique(paddock))

# standard deviation among paddocks and blocks
log(sd(tapply(dat$BRTE, dat$block, FUN = mean)))

# Assemble model inputs
datlist <- list(counts = dat$BRTE,
                area = dat$quadrat/100, # convert to square decimeters
                N = nrow(dat),
                POSE = as.numeric(X[,2]),
                POFE = as.numeric(X[,3]),
                VUMI = as.numeric(X[,4]),
                ELEL = as.numeric(X[,5]),
                high = as.numeric(X[,6]),
                fall = as.numeric(X[,7]), 
                spring = as.numeric(X[,8]), 
                coated = as.numeric(X[,9]),
                nL = ncol(X) - 1, # number of offset levels
                # pad = as.numeric(dat$paddock),
                # Np = length(unique(dat$paddock)),
                block = as.numeric(dat$block),
                Nb = length(unique(dat$block)),
                Ab = 5) # stand deviation among paddocks and blocks
str(datlist)

# likely intercept value
base <- dat %>%
  filter(grazing == "ungrazed",
         species == "ELTR",
         seed_rate == "low", 
         seed_coat == "UC")
hist(base$BRTE, breaks = 30)
summary(base$quadrat)
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
load("inits/inits_OLRE.Rdata")# saved.state, second element is inits
initslist <- list(append(saved.state[[2]][[1]], list(.RNG.name = array("base::Marsaglia-Multicarry"), .RNG.seed = array(13))),
                  append(saved.state[[2]][[2]], list(.RNG.name = array("base::Wichmann-Hill"), .RNG.seed = array(89))),
                  append(saved.state[[2]][[3]], list(.RNG.name = array("base::Mersenne-Twister"), .RNG.seed = array(18))))

# model
jm <- jags.model(file = "BRTE_counts_PoissonOLRE.jags",
                 inits = initslist,
                 n.chains = 3,
                 data = datlist)
# update(jm, 10000)

# params to monitor
params <- c("deviance", "Dsum", # evaluate fit
            "alpha", "beta", # parameters
            "tau.Eps", "sig.eps", # precision/variance terms for block RE
            "tau", "sig", # precision/variance terms for OLRE
            "alpha.star", "eps.star",  # identifiable intercept and block RE
            "int_Beta", "Diff_Beta", "diff_Beta", # calculated terms
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
                         n.iter = 150000, thin = 50)

# plot chains
mcmcplot(coda.out, parms = c("deviance", "Dsum","alpha.star", 
                             "beta", "eps.star", "sig.eps", "sig"))

traplot(coda.out, parms = "alpha.star")

# dic samples
dic.out <- dic.samples(jm, n.iter = 5000)
dic.out

# convergence?
gel <- gelman.diag(coda.out, multivariate = FALSE)
gel

# If not converged, restart model from final iterations
# newinits <-  initfind(coda.out)
# newinits[[1]]
# saved.state <- removevars(newinits, variables = c(1:2, 4, 6:70))
# saved.state[[1]]
# save(saved.state, file = "inits/inits_OLRE.Rdata")

save(coda.out, file = "coda/coda_OLRE.Rdata")

# Model fit
params <- c("counts.rep") #monitor replicated data
coda.rep <- coda.samples(jm, variable.names = params,
                         n.iter = 150000, thin = 50)
save(coda.rep, file = "coda/coda_OLRE_rep.Rdata")

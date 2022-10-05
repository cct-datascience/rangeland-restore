# BRTE cover modeled as zero inflated beta
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
load("cleaned_data/cover_mono.Rdata") # cover_mono
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
  facet_grid(cols = vars(grazing),
             rows = vars(seed_coat))

# No zeroes

# model matrix
X <- model.matrix( ~ (species + seed_rate + grazing + seed_coat)^2, data = dat) 
colnames(X)

# split the data into 0 and continuous components
y.temp <- with(dat, ifelse(BRTE == 0, 
                           BRTE, NA))
y.0 <- ifelse(is.na(y.temp), 0, 1)

# group continuous response + predictors
which.cont <- which(y.0 == 0)
y.c <- dat$BRTE[which.cont]
x.c <- X[which.cont,]
n.cont <- length(y.c)
block.c <- as.numeric(dat$block)[which.cont]

# Assemble model inputs
datlist <- list(N = nrow(dat),
                y.0 = y.0,
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

str(datlist)

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
       rho = runif(1, 0, 1))
}
initslist <- list(inits(), inits(), inits())

# Or, use previous starting values + set seed
load("models/cover/BRTE/mono/inits/inits.Rdata")# saved.state, second element is inits
initslist <- list(append(saved.state[[2]][[1]], list(.RNG.name = array("base::Marsaglia-Multicarry"), .RNG.seed = array(14))),
                  append(saved.state[[2]][[2]], list(.RNG.name = array("base::Super-Duper"), .RNG.seed = array(89))),
                  append(saved.state[[2]][[3]], list(.RNG.name = array("base::Mersenne-Twister"), .RNG.seed = array(18))))

# model
jm <- jags.model(file = "models/cover/BRTE/mono/BRTE_cover_zib.jags",
                 inits = initslist,
                 n.chains = 3,
                 data = datlist)
# update(jm, 10000)

# params to monitor
params <- c("deviance", # evaluate fit
            "rho", "alpha", "beta", # parameters
            "tau", "sig", "tau.eps", "sig.eps", # precision/variance terms
            "alpha.star", "eps.star", # identifiable intercept and random effects
            "int_Beta", "Diff_Beta", "diff_Beta", # monitored main and two-way treatment effects
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
mcmcplot(coda.out, parms = c("deviance", "rho", "beta",
                             "alpha.star",  "eps.star", 
                             "sig", "sig.eps"))
traplot(coda.out, parms = "rho")
caterplot(coda.out, parms = "eps.star", reorder = FALSE)
caterplot(coda.out, parms = "beta", reorder = FALSE)
caterplot(coda.out, parms = "Diff_Beta", reorder = FALSE)
caterplot(coda.out, parms = "diff_Beta", reorder = FALSE)
caterplot(coda.out, regex = "m\\.", reorder = FALSE)

# dic samples
dic.out <- dic.samples(jm, n.iter = 5000)
dic.out

# convergence?
gel <- gelman.diag(coda.out, multivariate = FALSE)
gel

# If not converged, restart model from final iterations
# newinits <-  initfind(coda.out)
# newinits[[1]]
# saved.state <- removevars(newinits, variables = c(1, 3, 5:67, 69:70))
# saved.state[[1]]
# save(saved.state, file = "models/cover/BRTE/mono/inits/inits.Rdata")

save(coda.out, file = "models/cover/BRTE/mono/coda/coda.Rdata")


# Model fit
params <- c("y.0.rep", "y.c.rep") #monitor replicated data
coda.rep <- coda.samples(jm, variable.names = params,
                         n.iter = 15000, thin = 5)
save(coda.rep, file = "models/cover/BRTE/mono/coda/coda_rep.Rdata")


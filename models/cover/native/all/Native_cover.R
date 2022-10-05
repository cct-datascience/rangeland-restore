# Native grass cover modeled as zero-or-1 inflated beta
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
load("cleaned_data/cover_all.Rdata") # cover_all
# convert to proportions
dat <- cover_all %>%
  mutate(BRTE = BRTE/100,
         intro_forbs = intro_forbs/100,
         native_grass = native_grass/100,
         native_forbs = native_forbs/100)
str(dat)

# model matrix
X <- model.matrix( ~ grazing * fuelbreak, data = dat) 
colnames(X)

# split the data into discrete and continuous components
y.temp <- with(dat, ifelse(native_grass == 1 | native_grass == 0, 
                           native_grass, NA))
y.0 <- ifelse(is.na(y.temp), 0, 1)

# group continuous response + predictors
which.cont <- which(y.0 == 0)
y.c <- dat$native_grass[which.cont]
x.c <- X[which.cont,]
n.cont <- length(y.c)
block.c <- as.numeric(dat$block)[which.cont]

# Assemble model inputs
datlist <- list(N = nrow(dat),
                y.0 = y.0,
                fall = X[,2],
                spring = X[,3],
                herbicide = X[,4],
                greenstrip = X[,5],
                fall_herbicide = X[,6],
                spring_herbicide = X[,7],
                fall_greenstrip = X[,8],
                spring_greenstrip = X[,9],
                block = as.numeric(dat$block),
                n.cont = n.cont,
                y.c = y.c,
                fall2 = x.c[,2],
                spring2 = x.c[,3],
                herbicide2 = x.c[,4],
                greenstrip2 = x.c[,5],
                fall_herbicide2 = x.c[,6],
                spring_herbicide2 = x.c[,7],
                fall_greenstrip2 = x.c[,8],
                spring_greenstrip2 = x.c[,9],
                block.c = block.c,
                Nb = length(unique(dat$block)),
                nL = ncol(X) - 1)

# likely intercept value
base <- dat %>%
  filter(grazing == "ungrazed",
         fuelbreak == "control")
logit(median(base$native_grass))

# generate random initials
inits <- function(){
  list(a = rnorm(1, 0, 10),
       b = rnorm(ncol(X) - 1, 0, 10),
       alpha = rnorm(1, 0, 10),
       beta = rnorm(ncol(X) - 1, 0, 10),
       tau = runif(1, 0, 1),
       tau.eps = runif(1, 0, 1),
       tau.eps.c = runif(1, 0, 1))
}
initslist <- list(inits(), inits(), inits())

# Or, use previous starting values + set seed
load("models/cover/native/all/inits/inits.Rdata")# saved.state, second element is inits
initslist <- list(append(saved.state[[2]][[1]], list(.RNG.name = array("base::Super-Duper"), .RNG.seed = array(13))),
                  append(saved.state[[2]][[2]], list(.RNG.name = array("base::Wichmann-Hill"), .RNG.seed = array(89))),
                  append(saved.state[[2]][[3]], list(.RNG.name = array("base::Mersenne-Twister"), .RNG.seed = array(18))))

# model
jm <- jags.model(file = "models/cover/native/all/Native_cover_zib.jags",
                 inits = initslist,
                 n.chains = 3,
                 data = datlist)
# update(jm, 10000)

# params to monitor
params <- c("deviance", # evaluate fit
            "a", "b", "alpha", "beta", # parameters
            "tau", "sig", # precision/variance terms
            "tau.eps", "sig.eps", "tau.eps.c", "sig.eps.c", 
            "a.star", "alpha.star", 
            "eps.star", "eps.star.c", # identifiable intercept and random effects
            "Diff_b", "diff_b", "Diff_Beta", "diff_Beta", # monitored main and two-way treatment effects
            "m.ungrazed.control", "m.ungrazed.herbicide", "m.ungrazed.greenstrip",
            "m.fall.control", "m.fall.herbicide", "m.fall.greenstrip",
            "m.spring.control", "m.spring.herbicide", "m.spring.greenstrip",
            "rho.ungrazed.control", "rho.ungrazed.herbicide", "rho.ungrazed.greenstrip",
            "rho.fall.control", "rho.fall.herbicide", "rho.fall.greenstrip",
            "rho.spring.control", "rho.spring.herbicide", "rho.spring.greenstrip") 
           

coda.out <- coda.samples(jm, variable.names = params,
                         n.iter = 15000, thin = 5)

# plot chains
mcmcplot(coda.out, parms = c("deviance", "a.star", "b", "beta",
                             "alpha.star",  "eps.star", "eps.star.c",
                             "sig", "sig.eps", "sig.eps.c"))
traplot(coda.out, parms = c("sig", "sig.eps", "sig.eps.c"))
caterplot(coda.out, parms = "eps.star", reorder = FALSE)
caterplot(coda.out, regex = "b\\[", reorder = FALSE)
caterplot(coda.out, parms = "beta", reorder = FALSE)
caterplot(coda.out, parms = "Diff_Beta", reorder = FALSE)
caterplot(coda.out, regex = "diff_b\\[", reorder = FALSE)
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
# saved.state <- removevars(newinits, variables = c(1, 3, 5, 8:31))
# saved.state[[1]]
# save(saved.state, file = "models/cover/native/all/inits/inits.Rdata")

save(coda.out, file = "models/cover/native/all/coda/coda.Rdata")


# Model fit
params <- c("y.0.rep", "y.c.rep") #monitor replicated data
coda.rep <- coda.samples(jm, variable.names = params,
                         n.iter = 15000, thin = 5)
save(coda.rep, file = "models/cover/native/all/coda/coda_rep.Rdata")

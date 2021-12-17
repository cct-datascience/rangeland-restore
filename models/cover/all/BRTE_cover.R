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
load("../../../cleaned_data/cover_all.Rdata") # cover_all
# convert to proportions
dat <- cover_all %>%
  mutate(BRTE = BRTE/100,
         intro_forbs = intro_forbs/100,
         native_grass = native_grass/100,
         native_forbs = native_forbs/100)
str(dat)


# plot
ggplot(dat, aes(x = grazing, y = BRTE, color = block)) +
  geom_point()

# Check distribution of zeros
dat %>%
  filter(BRTE == 0) %>%
ggplot(aes(x = fuelbreak, y = BRTE, color = block)) +
  geom_point() +
  facet_wrap(~grazing)

# model matrix
X <- model.matrix( ~ grazing * fuelbreak, data = dat) 
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
str(datlist)

# likely intercept value
base <- dat %>%
  filter(grazing == "ungrazed",
         fuelbreak == "control")
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
load("inits/inits.Rdata")# saved.state, second element is inits
initslist <- list(append(saved.state[[2]][[1]], list(.RNG.name = array("base::Super-Duper"), .RNG.seed = array(13))),
                  append(saved.state[[2]][[2]], list(.RNG.name = array("base::Wichmann-Hill"), .RNG.seed = array(89))),
                  append(saved.state[[2]][[3]], list(.RNG.name = array("base::Mersenne-Twister"), .RNG.seed = array(18))))

for(i in 1:3){
  initslist[[i]]$rho <- initslist[[i]]$psi
}

# model
jm <- jags.model(file = "BRTE_cover_zib.jags",
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
            "m.ungrazed.control", "m.ungrazed.herbicide", "m.ungrazed.greenstrip",
            "m.fall.control", "m.fall.herbicide", "m.fall.greenstrip",
            "m.spring.control", "m.spring.herbicide", "m.spring.greenstrip") 

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
caterplot(coda.out, regex = "m\\.")

# dic samples
dic.out <- dic.samples(jm, n.iter = 5000)
dic.out

# convergence?
gel <- gelman.diag(coda.out, multivariate = FALSE)
gel

# If not converged, restart model from final iterations
# newinits <-  initfind(coda.out)
# newinits[[1]]
# saved.state <- removevars(newinits, variables = c(1, 3, 5:16, 18:19))
# saved.state[[1]]
# save(saved.state, file = "inits/inits.Rdata")

save(coda.out, file = "coda/coda.Rdata")


# Model fit
params <- c("y.0.rep", "y.c.rep") #monitor replicated data
coda.rep <- coda.samples(jm, variable.names = params,
                         n.iter = 15000, thin = 5)
save(coda.rep, file = "coda/coda_rep.Rdata")


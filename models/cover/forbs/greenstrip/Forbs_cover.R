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
load("cleaned_data/cover_greenstrip.Rdata") # cover_greenstrip
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

summary(dat$forbs)
# No zeros

# model matrix
X <- model.matrix( ~ (spatial + seed_rate + seed_coat + grazing)^2, data = dat) 
colnames(X)

# split the data into 0 and continuous components
y.temp <- with(dat, ifelse(forbs == 0, 
                           forbs, NA))
y.0 <- ifelse(is.na(y.temp), 0, 1)

# group continuous response + predictors
which.cont <- which(y.0 == 0)
y.c <- dat$forbs[which.cont]
x.c <- X[which.cont,]
n.cont <- length(y.c)
block.c <- as.numeric(dat$block)[which.cont]

# Assemble model inputs
datlist <- list(N = nrow(dat),
                y.0 = y.0,
                n.cont = n.cont,
                y.c = y.c,
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
                block.c = block.c,
                Nb = length(unique(dat$block)),
                nL = ncol(X) - 1)
str(datlist)

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
       tau.eps = runif(1, 0, 0.1),
       rho = runif(1, 0, 1))
}
initslist <- list(inits(), inits(), inits())

# Or, use previous starting values + set seed
load("models/cover/forbs/greenstrip/inits/inits.Rdata")# saved.state, second element is inits
initslist <- list(append(saved.state[[2]][[1]], list(.RNG.name = array("base::Super-Duper"), .RNG.seed = array(13))),
                  append(saved.state[[2]][[2]], list(.RNG.name = array("base::Wichmann-Hill"), .RNG.seed = array(89))),
                  append(saved.state[[2]][[3]], list(.RNG.name = array("base::Mersenne-Twister"), .RNG.seed = array(18))))

# model
jm <- jags.model(file = "models/cover/forbs/greenstrip/Forbs_cover_zib.jags",
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
mcmcplot(coda.out, parms = c("deviance", "rho", "beta", 
                             "alpha.star",  "eps.star", 
                             "sig", "sig.eps"))
traplot(coda.out, parms = "rho")
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
# saved.state <- removevars(newinits, variables = c(1, 3, 5:31, 33:34))
# saved.state[[1]]
# save(saved.state, file = "models/cover/forbs/greenstrip/inits/inits.Rdata")

save(coda.out, file = "models/cover/forbs/greenstrip/coda/coda.Rdata")


# Model fit
params <- c("y.0.rep", "y.c.rep") #monitor replicated data
coda.rep <- coda.samples(jm, variable.names = params,
                         n.iter = 15000, thin = 5)
save(coda.rep, file = "models/cover/forbs/greenstrip/coda/coda_rep.Rdata")

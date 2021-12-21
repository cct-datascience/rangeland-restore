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
load("../../../../cleaned_data/cover_greenstrip.Rdata") # cover_greenstrip
# convert to proportions
dat <- cover_greenstrip %>%
  mutate(BRTE = BRTE/100,
         intro_forbs = intro_forbs/100,
         native_grass = native_grass/100,
         native_forbs = native_forbs/100,
         spatial = factor(spatial, levels = c("mix", "mono")))
str(dat)

#plot
dat %>%
  filter(native_grass != 0) %>%
  ggplot(aes(x = seed_rate, y = native_grass, color = seed_coat)) +
  geom_jitter(height = 0.1) +
  scale_y_continuous(limits = c(0, 0.5)) +
  facet_grid(cols = vars(grazing),
             rows = vars(spatial))

# Check distribution of zeros
dat %>%
  filter(native_grass == 0) %>%
  ggplot(aes(x = seed_rate, y = native_grass, color = seed_coat)) +
  geom_jitter(height = 0.1) +
  scale_y_continuous(limits = c(0, 0.5)) +
  facet_grid(cols = vars(grazing),
             rows = vars(spatial))

dat %>%
  filter(native_grass == 0) %>%
  count(spatial, grazing)

dat %>%
  group_by(spatial, seed_rate, seed_coat, grazing) %>%
  summarize(rho = length(which(native_grass == 0))/length(native_grass)) %>%
  ggplot(aes(x = seed_rate, y = rho, color = seed_coat)) +
  geom_point() +
  scale_y_continuous() +
  facet_grid(cols = vars(grazing),
             rows = vars(spatial))

dat %>%
  filter(native_grass != 0) %>%
  group_by(spatial, seed_rate, seed_coat, grazing) %>%
  summarize(mean_ng = mean(native_grass)) %>%
  ggplot(aes(x = seed_rate, y = mean_ng, color = seed_coat)) +
  geom_point() +
  scale_y_continuous() +
  facet_grid(cols = vars(grazing),
             rows = vars(spatial))

# model matrix
X <- model.matrix( ~ (spatial + seed_rate + seed_coat + grazing)^2, data = dat) 
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
                n.cont = n.cont,
                y.c = y.c,
                mono2 = x.c[,2],
                high2 = x.c[,3],
                coated2 = x.c[,4],
                fall2 = x.c[,5],
                spring2 = x.c[,6],
                mono_high2 = x.c[,7],
                mono_coated2 = x.c[,8],
                mono_fall2 = x.c[,9],
                mono_spring2 = x.c[,10],
                high_coated2 = x.c[,11],
                high_fall2 = x.c[,12],
                high_spring2 = x.c[,13],
                coated_fall2 = x.c[,14],
                coated_spring2 = x.c[,15],
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
hist(base$native_grass, breaks = 30)
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
load("inits/inits.Rdata")# saved.state, second element is inits
initslist <- list(append(saved.state[[2]][[1]], list(.RNG.name = array("base::Super-Duper"), .RNG.seed = array(41))),
                  append(saved.state[[2]][[2]], list(.RNG.name = array("base::Wichmann-Hill"), .RNG.seed = array(89))),
                  append(saved.state[[2]][[3]], list(.RNG.name = array("base::Mersenne-Twister"), .RNG.seed = array(11))))

# model
jm <- jags.model(file = "Native_cover_zib.jags",
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
            "m.mono.low.uncoated.ungrazed", "m.mono.low.uncoated.fall", "m.mono.low.uncoated.spring", 
            "m.mono.low.coated.ungrazed", "m.mono.low.coated.fall", "m.mono.low.coated.spring",
            "m.mono.high.uncoated.ungrazed", "m.mono.high.uncoated.fall", "m.mono.high.uncoated.spring",
            "m.mono.high.coated.ungrazed", "m.mono.high.coated.fall", "m.mono.high.coated.spring",
            "m.mix.low.uncoated.ungrazed", "m.mix.low.uncoated.fall", "m.mix.low.uncoated.spring",
            "m.mix.low.coated.ungrazed", "m.mix.low.coated.fall", "m.mix.low.coated.spring",
            "m.mix.high.uncoated.ungrazed", "m.mix.high.uncoated.fall", "m.mix.high.uncoated.spring",
            "m.mix.high.coated.ungrazed", "m.mix.high.coated.fall", "m.mix.high.coated.spring",
            "rho.mono.low.uncoated.ungrazed", "rho.mono.low.uncoated.fall", "rho.mono.low.uncoated.spring", 
            "rho.mono.low.coated.ungrazed", "rho.mono.low.coated.fall", "rho.mono.low.coated.spring",
            "rho.mono.high.uncoated.ungrazed", "rho.mono.high.uncoated.fall", "rho.mono.high.uncoated.spring",
            "rho.mono.high.coated.ungrazed", "rho.mono.high.coated.fall", "rho.mono.high.coated.spring",
            "rho.mix.low.uncoated.ungrazed", "rho.mix.low.uncoated.fall", "rho.mix.low.uncoated.spring",
            "rho.mix.low.coated.ungrazed", "rho.mix.low.coated.fall", "rho.mix.low.coated.spring",
            "rho.mix.high.uncoated.ungrazed", "rho.mix.high.uncoated.fall", "rho.mix.high.uncoated.spring",
            "rho.mix.high.coated.ungrazed", "rho.mix.high.coated.fall", "rho.mix.high.coated.spring") 
           

coda.out <- coda.samples(jm, variable.names = params,
                         n.iter = 15000, thin = 5)

# plot chains
mcmcplot(coda.out, parms = c("deviance", "a.star", "b", "beta",
                             "alpha.star",  "eps.star", "eps.star.c",
                             "sig", "sig.eps", "sig.eps.c"))
traplot(coda.out, parms = c("sig", "sig.eps", "sig.eps.c"))
traplot(coda.out, parms = c("a.star", "alpha.star"))
caterplot(coda.out, parms = "eps.star", reorder = FALSE)
caterplot(coda.out, parms = "eps.star.c", reorder = FALSE)
caterplot(coda.out, regex = "^b\\[", reorder = FALSE)
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
# saved.state <- removevars(newinits, variables = c(1:2, 4, 6, 9:63))
# saved.state[[1]]
# save(saved.state, file = "inits/inits.Rdata")

save(coda.out, file = "coda/coda.Rdata")


# Model fit
params <- c("y.0.rep", "y.c.rep") #monitor replicated data
coda.rep <- coda.samples(jm, variable.names = params,
                         n.iter = 15000, thin = 5)
save(coda.rep, file = "coda/coda_rep.Rdata")

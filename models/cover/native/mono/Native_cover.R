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
load("cleaned_data/cover_mono.Rdata") # cover_mono
# convert to proportions
dat <- cover_mono %>%
  mutate(BRTE = BRTE/100,
         intro_forbs = intro_forbs/100,
         native_grass = native_grass/100,
         native_forbs = native_forbs/100)

# Relevel species based on fig. 6b from Porensky et al. 2018
dat$species <- factor(dat$species, levels = c("ELTR", "POSE", "POFE", "VUMI", "ELEL"))


#plot
dat %>%
  filter(native_grass != 0) %>%
  ggplot(aes(x = seed_rate, y = native_grass, color = seed_coat)) +
  geom_jitter(height = 0.1) +
  scale_y_continuous(limits = c(0, 0.5)) +
  facet_grid(cols = vars(grazing),
             rows = vars(species))

# Check distribution of zeros
dat %>%
  filter(native_grass == 0) %>%
  ggplot(aes(x = seed_rate, y = native_grass, color = seed_coat)) +
  geom_jitter(height = 0.01, width = 0.2) +
  scale_y_continuous(limits = c(-0.1, 0.2)) +
  facet_grid(cols = vars(grazing),
             rows = vars(species))

dat %>%
  filter(native_grass == 0) %>%
  count(species, grazing)

dat %>%
  group_by(species, seed_rate, seed_coat, grazing) %>%
  summarize(rho = length(which(native_grass == 0))/length(native_grass)) %>%
  ggplot(aes(x = seed_rate, y = rho, color = seed_coat)) +
  geom_point() +
  scale_y_continuous() +
  facet_grid(cols = vars(grazing),
             rows = vars(species))

dat %>%
  filter(native_grass != 0) %>%
  group_by(species, seed_rate, seed_coat, grazing) %>%
  summarize(mean_ng = mean(native_grass)) %>%
  ggplot(aes(x = seed_rate, y = mean_ng, color = seed_coat)) +
  geom_point() +
  scale_y_continuous() +
  facet_grid(cols = vars(grazing),
             rows = vars(species))

# model matrix
X <- model.matrix( ~ (species + seed_rate + grazing + seed_coat)^2, data = dat) 
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
                POSE = X[,2],
                POFE = X[,3],
                VUMI = X[,4],
                ELEL = X[,5],
                high = X[,6],
                fall = X[,7], 
                spring = X[,8], 
                coated = X[,9],
                block = as.numeric(dat$block),
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
load("models/cover/native/mono/inits/inits.Rdata")# saved.state, second element is inits
initslist <- list(append(saved.state[[2]][[1]], list(.RNG.name = array("base::Super-Duper"), .RNG.seed = array(41))),
                  append(saved.state[[2]][[2]], list(.RNG.name = array("base::Wichmann-Hill"), .RNG.seed = array(89))),
                  append(saved.state[[2]][[3]], list(.RNG.name = array("base::Mersenne-Twister"), .RNG.seed = array(11))))

# model
jm <- jags.model(file = "models/cover/native/mono/Native_cover_zib.jags",
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
            "rho.ELEL.low.ungrazed.uncoated", "rho.ELEL.low.ungrazed.coated", "rho.ELEL.low.fall.uncoated", "rho.ELEL.low.fall.coated", "rho.ELEL.low.spring.uncoated", "rho.ELEL.low.spring.coated", 
            "rho.ELEL.high.ungrazed.uncoated", "rho.ELEL.high.ungrazed.coated", "rho.ELEL.high.fall.uncoated", "rho.ELEL.high.fall.coated", "rho.ELEL.high.spring.uncoated", "rho.ELEL.high.spring.coated", 
            "rho.VUMI.low.ungrazed.uncoated", "rho.VUMI.low.ungrazed.coated", "rho.VUMI.low.fall.uncoated", "rho.VUMI.low.fall.coated", "rho.VUMI.low.spring.uncoated", "rho.VUMI.low.spring.coated", 
            "rho.VUMI.high.ungrazed.uncoated", "rho.VUMI.high.ungrazed.coated", "rho.VUMI.high.fall.uncoated", "rho.VUMI.high.fall.coated" , "rho.VUMI.high.spring.uncoated", "rho.VUMI.high.spring.coated", 
            "rho.POFE.low.ungrazed.uncoated", "rho.POFE.low.ungrazed.coated", "rho.POFE.low.fall.uncoated", "rho.POFE.low.fall.coated", "rho.POFE.low.spring.uncoated", "rho.POFE.low.spring.coated", 
            "rho.POFE.high.ungrazed.uncoated", "rho.POFE.high.ungrazed.coated", "rho.POFE.high.fall.uncoated", "rho.POFE.high.fall.coated", "rho.POFE.high.spring.uncoated", "rho.POFE.high.spring.coated", 
            "rho.POSE.low.ungrazed.uncoated", "rho.POSE.low.ungrazed.coated", "rho.POSE.low.fall.uncoated", "rho.POSE.low.fall.coated", "rho.POSE.low.spring.uncoated", "rho.POSE.low.spring.coated", 
            "rho.POSE.high.ungrazed.uncoated", "rho.POSE.high.ungrazed.coated", "rho.POSE.high.fall.uncoated", "rho.POSE.high.fall.coated", "rho.POSE.high.spring.uncoated", "rho.POSE.high.spring.coated", 
            "rho.ELTR.low.ungrazed.uncoated", "rho.ELTR.low.ungrazed.coated", "rho.ELTR.low.fall.uncoated", "rho.ELTR.low.fall.coated", "rho.ELTR.low.spring.uncoated", "rho.ELTR.low.spring.coated", 
            "rho.ELTR.high.ungrazed.uncoated", "rho.ELTR.high.ungrazed.coated", "rho.ELTR.high.fall.uncoated", "rho.ELTR.high.fall.coated", "rho.ELTR.high.spring.uncoated", "rho.ELTR.high.spring.coated",
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
# saved.state <- removevars(newinits, variables = c(1:2, 4, 6, 9:135))
# saved.state[[1]]
# save(saved.state, file = "models/cover/native/mono/inits/inits.Rdata")

save(coda.out, file = "models/cover/native/mono/coda/coda.Rdata")


# Model fit
params <- c("y.0.rep", "y.c.rep") #monitor replicated data
coda.rep <- coda.samples(jm, variable.names = params,
                         n.iter = 15000, thin = 5)
save(coda.rep, file = "models/cover/native/mono/coda/coda_rep.Rdata")

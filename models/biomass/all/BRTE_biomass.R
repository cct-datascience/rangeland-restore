# Cheatgrass biomass, relative to total and total live


library(rjags)
load.module('dic')
load.module('glm')
library(mcmcplots)
library(postjags)
library(ggplot2)
library(dplyr)

# Read in data
load("../../../cleaned_data/biomass.Rdata") # bio
dat <- bio
str(dat)

# Plot
dat %>%
  ggplot(aes(x = grazing, y = BRTE_abs)) +
  geom_jitter(aes(color = paddock, shape = block))

dat %>%
  ggplot(aes(x = grazing, y = forbs_abs)) +
  geom_jitter(aes(color = paddock, shape = block))

dat %>%
  ggplot(aes(x = grazing, y = grass_abs)) +
  geom_jitter(aes(color = paddock, shape = block))

dat %>%
  ggplot() +
  geom_violin(aes(x = grazing, y = BRTE_abs, fill = "BRTE"), alpha = 0.5) +
  geom_violin(aes(x = grazing, y = forbs_abs, fill = "forbs"), alpha = 0.5) +
  geom_violin(aes(x = grazing, y = grass_abs, fill = "grass"), alpha = 0.5)

# model matrix
X <- model.matrix( ~ grazing, data = dat) 
colnames(X)

# standard deviation among paddocks and blocks
sd(tapply(dat$BRTE_abs, dat$block, FUN = mean))
sd(tapply(dat$forbs_abs, dat$block, FUN = mean))
sd(tapply(dat$grass_abs, dat$block, FUN = mean))

#create the R matrix, 3 by 43, 0's outside of diagonal
R <- diag(x = 1, 3, 3)

# Assemble model inputs
datlist <- list(biomass = dat[,5:7],
                N = nrow(dat),
                fall = X[,2],
                spring = X[,3],
                nL = ncol(X) - 1, # number of offset levels
                # pad = as.numeric(dat$paddock),
                # Np = length(unique(dat$paddock)),
                block = as.numeric(dat$block),
                Nb = length(unique(dat$block)),
                R = R,
                Ab = 5) # stand deviation among paddocks and blocks

# likely intercept value
base <- dat %>%
  filter(grazing == "ungrazed")
hist(base$BRTE_abs, breaks = 30)
colMeans(base[,5:7])

# Covariance of observed traits, for initials for Omega:
omega <- cov(dat[,5:7], use="complete.obs")
omega.gen <- function(x) {
  noise = rnorm(n = nrow(dat[,5:7])*ncol(dat[,5:7]), mean = 0, sd = 1)
  nois.mat = matrix(noise, ncol = ncol(dat[,5:7]))
  return(solve(cov(dat[,5:7] + nois.mat, use = "complete.obs")))
}

# Generate initials
inits <- function(){
  list(alpha = rnorm(3, 0, 10),
       beta = matrix(rnorm(3*(ncol(X) - 1), 0, 10), ncol = 3),
       tau.Eps = runif(3, 0, 1),
       omega = round(omega.gen(), 4))
}
initslist <- list(inits(), inits(), inits())


# Or, use previous starting values + set seed
load("inits/inits.Rdata")# saved.state, second element is inits
initslist <- list(append(saved.state[[2]][[1]], list(.RNG.name = array("base::Super-Duper"), .RNG.seed = array(13))),
                  append(saved.state[[2]][[2]], list(.RNG.name = array("base::Wichmann-Hill"), .RNG.seed = array(89))),
                  append(saved.state[[2]][[3]], list(.RNG.name = array("base::Mersenne-Twister"), .RNG.seed = array(18))))

# model
jm <- jags.model(file = "BRTE_biomass_mvnorm.jags",
                 inits = initslist,
                 n.chains = 3,
                 data = datlist)
# update(jm, 10000)

# params to monitor
params <- c("deviance", "Dsum", # evaluate fit
            "alpha", "beta", # parameters
            "tau.Eps", "sig.eps", # precision/variance terms
            "alpha.star", "eps.star", # identifiable intercept and random effects
            "omega", "Sig", "Rho") # monitored interaction effects

coda.out <- coda.samples(jm, variable.names = params,
                         n.iter = 15000, thin = 5)

# plot chains
mcmcplot(coda.out, parms = c("deviance", "Dsum", "beta",
                             "alpha.star",  "eps.star", "sig.eps",
                             "Sig", "Rho"))

traplot(coda.out, parms = "sig.eps")
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
# saved.state <- removevars(newinits, variables = c(1:3, 5, 7, 9))
# saved.state[[1]]
# save(saved.state, file = "inits/inits.Rdata")

save(coda.out, file = "coda/coda.Rdata")

# Model fit
params <- c("biomass.rep") #monitor replicated data
coda.rep <- coda.samples(jm, variable.names = params,
                         n.iter = 15000, thin = 5)
save(coda.rep, file = "coda/coda_rep.Rdata")

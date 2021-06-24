# Cheatgrass counts modeled as Poisson distribution

library(rjags)
load.module('dic')
library(mcmcplots)
library(postjags)
library(ggplot2)
library(data.table)

# Read in data
load("../../../cleaned_data/count_all.Rdata") # count_all
dat <- count_all
str(dat)
ggplot(dat) +
  geom_bar(aes(x = BRTE.count)) +
  geom_bar(aes(x = Native.grass.count..low.data.quality.), 
           fill = "blue", alpha = 0.2)

# model matrix
X <- model.matrix( ~ grazing + fuelbreak, data = dat) 

# Assemble model inputs
datlist <- list(counts = dat$BRTE,
                area = dat$quadrat,
                N = nrow(dat),
                spring = X[,2],
                fall = X[,3],
                herbicide = X[,4],
                greenstrip = X[,5],
                nL = ncol(X) - 1) # number of levels


# initials
inits <- function(){
  list(alpha = runif(1, -3, 0),
       beta = runif(ncol(X) - 1, -3, 0),
       psi = runif(1, 0, 1))
}
initslist <- list(inits(), inits(), inits())

# model
jm <- jags.model(file = "BRTE_counts.jags",
                 inits = initslist,
                 n.chains = 3,
                 data = datlist)
update(jm, 10000)
# params to monitor
params <- c("deviance", "Dsum", "alpha", "beta", "psi", "prob", "diff")
coda.out <- coda.samples(jm, variable.names = params,
                         n.iter = 5000, thin = 1)

# plot chains
mcmcplot(coda.out)

# convergence?
gel <- gelman.diag(coda.out, multivariate = FALSE)
gel

# summarize
sum.out <- coda.fast(coda.out, OpenBUGS = FALSE)
sum.out$var <- row.names(sum.out)

# plot
diff.ind <- grep("diff", row.names(sum.out))
ggplot(sum.out[diff.ind,], aes(x = var, y = mean)) +
  geom_pointrange(aes(ymin = pc2.5, ymax = pc97.5))

prob.ind <- grep("prob", row.names(sum.out))
ggplot(sum.out[prob.ind,], aes(x = var, y = mean)) +
  geom_pointrange(aes(ymin = pc2.5, ymax = pc97.5))

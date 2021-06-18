# Cheatgrass counts modeled as Poisson distribution

library(rjags)
load.module('dic')
library(mcmcplots)
library(postjags)
library(ggplot2)
library(data.table)

# Read in data
dat <- read.csv("../raw_data/2019_Greenstrips_Cover_and_Densities.csv")
colnames(dat)
ggplot(dat) +
  geom_bar(aes(x = BRTE.count)) +
  geom_bar(aes(x = Native.grass.count..low.data.quality.), 
           fill = "blue", alpha = 0.2)

# assign grazing levels
dat$grazing <- factor(dat$Grazing, 
                      levels = c("Ungrazed", "Fall", "Spring"))

datlist <- list(counts = dat$BRTE.count,
                N = nrow(dat),
                grazing = as.numeric(dat$grazing),
                Ng = length(unique(dat$grazing)))


# initials
inits <- function(){
  list(theta = runif(length(unique(dat$grazing)), 0, 1),
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
params <- c("deviance", "Dsum", "theta", "psi")
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
t.ind <- grep("theta", row.names(sum.out))
ggplot(sum.out[t.ind,], aes(x = var, y = mean)) +
  geom_pointrange(aes(ymin = pc2.5, ymax = pc97.5))
.
# Initial model control script for biomass of BRTE
library(rjags)
load.module('dic')
library(mcmcplots)
library(postjags)
library(ggplot2)
library(data.table)

# Read in data
dat <- read.csv("../raw_data/2019_greenstrips_Biomass.csv")
dat <- as.data.table(dat)
# exploratory plotting
plot(dat$B..tectorum)
points(dat$Forbs, col = "red")
hist(dat$B..tectorum)
hist(log(dat$B..tectorum))
boxplot(log(dat$B..tectorum)~dat$grazing.trt)

ggplot(dat, aes(x = grazing, y = log(B..tectorum), color = Paddock)) +
  geom_point() +
  facet_wrap(~Block)

# assign grazing levels
dat$grazing <- factor(dat$grazing.trt, 
                      levels = c("ungrazed", "fall", "spring"))
# link plot to block
dat$Block <- as.numeric(factor(dat$Block, levels = c("one", "two", "three")))
link <- dat[,.(block = unique(Block)), .(Paddock) ]
# make list of data
datlist <- list(y = log(dat$B..tectorum),
                N = nrow(dat),
                grazing = as.numeric(dat$grazing),
                Ng = length(unique(dat$grazing)),
                pad = as.numeric(factor(dat$Paddock)),
                Np = length(unique(dat$Paddock)),
                block = link$block,
                Nb = length(unique(dat$Block)))

# initials
inits <- function(){
  list(tau = runif(1, 0, 1),
       tau.eps = runif(1, 0, 1),
       tau.gam = runif(1, 0, 1),
       alpha = rnorm(datlist$Ng, 0, 10))
}
initslist <- list(inits(), inits(), inits())

# model
jm <- jags.model(file = "BRTE_biomass.jags",
                 inits = initslist,
                 n.chains = 3,
                 data = datlist)
update(jm, 10000)
# params to monitor
params <- c("deviance", "Dsum", "alpha", "tau", 
            "eps", "tau.eps", "gam", "tau.gam", "sigs",
            "diff")
coda.out <- coda.samples(jm, variable.names = params,
                      n.iter = 500000, thin = 10)

# plot chains
mcmcplot(coda.out)

# convergence?
gel <- gelman.diag(coda.out, multivariate = FALSE)
gel

# summarize
sum.out <- coda.fast(coda.out, OpenBUGS = FALSE)
sum.out$var <- row.names(sum.out)

# plot
a.ind <- grep("alpha", row.names(sum.out))
ggplot(sum.out[a.ind,], aes(x = var, y = mean)) +
  geom_pointrange(aes(ymin = pc2.5, ymax = pc97.5))

d.ind <- grep("diff", row.names(sum.out))
ggplot(sum.out[d.ind,], aes(x = var, y = mean)) +
  geom_pointrange(aes(ymin = pc2.5, ymax = pc97.5)) +
  geom_hline(yintercept = 0, col = "red")

# Calculate Bayesian pvalues for diff
diffs <- rbind(data.frame(coda.out[[1]][,d.ind]),
               data.frame(coda.out[[2]][,d.ind]),
               data.frame(coda.out[[3]][,d.ind]))

pval <- function(x){mean(ifelse(x > 0, 1, 0))}
1 - apply(diffs, 2, pval)

# replicated
params <- c("y.rep")
coda.rep <- coda.samples(jm, variable.names = params,
                         n.iter = 5000)
sum.rep <- coda.fast(coda.rep, OpenBUGS = FALSE)

pred <- data.frame(dat, 
                   pred = sum.rep$mean,
                   pred.lower = sum.rep$pc2.5,
                   pred.upper = sum.rep$pc97.5)

ggplot(pred, aes(x = log(B..tectorum), y = pred)) +
  geom_pointrange(aes(ymin = pred.lower, ymax = pred.upper)) +
  geom_abline(slope = 1, intercept = 0, color = "red")

m1 <- lm(pred ~ B..tectorum, data = pred)
summary(m1)

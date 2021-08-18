# Cheatgrass counts modeled as Poisson distribution
# ANOVA 4 factors, (5,2, 3, 2 levels)
# Remove grazing as a fixed affect, account for it in the nested RE

library(rjags)
load.module('dic')
load.module('glm')
library(mcmcplots)
library(postjags)
library(ggplot2)
library(dplyr)

# Read in data
load("../../../cleaned_data/count_mono.Rdata") # count_mono
dat <- count_mono
str(dat)

# Quadrat sizes
# dat2 <- dat %>%
#   filter(quadrat < 10000)
table(dat$quadrat)

# Plot
dat %>%
  ggplot(aes(x = grazing, y = BRTE/quadrat)) +
  geom_jitter(aes(color = block))

dat %>%
  ggplot(aes(x = species, y = BRTE/quadrat)) +
  geom_jitter(aes(color = block)) +
  facet_wrap(~grazing, ncol = 3)

dat %>%
  ggplot(aes(x = seed_rate, y = BRTE/quadrat)) +
  geom_jitter(aes(color = block)) +
  facet_wrap(~grazing, ncol = 3)

dat %>%
  ggplot(aes(x = seed_coat, y = BRTE/quadrat)) +
  geom_jitter(aes(color = block)) +
  facet_wrap(~grazing, ncol = 3)


# model matrix
X <- model.matrix( ~ (species + seed_rate + seed_coat), data = dat) 
dim(X)

# link paddocks to block
link <- dat %>%
  group_by(block) %>%
  summarize(paddock = unique(paddock))

# standard deviation among paddocks and blocks
sd(tapply(dat$BRTE, dat$paddock, FUN = mean))
sd(tapply(dat$BRTE, dat$block, FUN = mean))

# Assemble model inputs
datlist <- list(counts = dat$BRTE,
                area = dat$quadrat/100, # convert to square decimeters
                N = nrow(dat),
                ELTR = as.numeric(X[,2]),
                POFE = as.numeric(X[,3]),
                POSE = as.numeric(X[,4]),
                VUMI = as.numeric(X[,5]),
                high = as.numeric(X[,6]),
                coated = as.numeric(X[,7]),
                nL = ncol(X) - 1, # number of levels
                pad = as.numeric(dat$paddock),
                Np = length(unique(dat$paddock)),
                block = as.numeric(link$block),
                Nb = length(unique(dat$block)),
                Ab = 10) # stand deviation among paddocks and blocks
str(datlist)

# likely intercept value
base <- dat %>%
  filter(grazing == "ungrazed",
         spatial == "mono",
         seed_rate == "low", 
         seed_coat == "UC")
hist(base$BRTE, breaks = 30)
summary(base$quadrat)
mean(base$BRTE)

# initials
inits <- function(){
  list(alpha = rnorm(1, 0, 10),
       beta = rnorm(ncol(X) - 1, 0, 10),
       psi = runif(1, 0, 1),
       tau.Eps = runif(1, 0, 3),
       tau.Gam = runif(1, 0, 3)
       )
}
initslist <- list(inits(), inits(), inits())

# model
jm <- jags.model(file = "BRTE_counts_nograzing.jags",
                 inits = saved.state[[2]],
                 n.chains = 3,
                 data = datlist)
update(jm, 10000)

# params to monitor
params <- c("deviance", "Dsum", # evaluate fit
            "alpha", "beta", "psi", # parameters
            "diff", "prob", # derived parameters
            "tau.Eps", "tau.Gam", "sig.eps", "sig.gam", # precision/variance terms
            "alpha.star", "eps.star", "gam.star") #identifiable intercept and random effects
 
params <- c("deviance", "Dsum", # evaluate fit
            "alpha", "beta", #"psi", # parameters,
            "diff") # derived parameters)

coda.out <- coda.samples(jm, variable.names = params,
                         n.iter = 5000, thin = 1)

# plot chains
mcmcplot(coda.out, parms = c("deviance", "Dsum","psi","alpha.star", 
                             "beta", "eps.star", "gam.star",
                             "sig.eps", "sig.gam"))
mcmcplot(coda.out, parms = c("diff"))

# dic samples
dic.out <- dic.samples(jm, n.iter = 5000)
dic.out

# convergence?
gel <- gelman.diag(coda.out, multivariate = FALSE)
gel

# If not converged, restart model from final iterations
newinits <-  initfind(coda.out) 
newinits[[1]]
saved.state <- removevars(newinits, variables = c(1, 3, 5:8, 10:11))
saved.state[[1]]
save(saved.state, file = "inits/inits_nograzing.Rdata")

its <- list(saved.state[[2]][[3]],
            saved.state[[2]][[2]],
            saved.state[[2]][[2]])
grep("sig.gam", colnames(coda.out[[1]]))
coda.out[[1]][5000,1]
coda.out[[2]][5000,1]
coda.out[[3]][5000,1]
coda.out[[1]][5000,40]
coda.out[[2]][5000,40]
coda.out[[3]][5000,40]
coda.out[[1]][5000,39]
coda.out[[2]][5000,39]
coda.out[[3]][5000,39]
coda.out[[1]][5000,47]
coda.out[[2]][5000,47]
coda.out[[3]][5000,47]
# mc <- coda.out[[3]]
# means <- round(colMeans(mc),2)

# do random effects add up correctly?
inds.eps <- grep("eps.star", colnames(mc))
head(rowSums(mc[,inds.eps[1:3]])) # paddocks within block 1
head(rowSums(mc[,inds.eps[4:6]])) # paddocks within block 2
head(rowSums(mc[,inds.eps[7:9]])) # paddocks within block 3
inds.gam <-grep("gam.star", colnames(mc))
head(rowSums(mc[,inds.gam])) # blocks



# plot trace plots
traplot(coda.out, parms = "alpha.star")
traplot(coda.out, parms = "eps")
traplot(coda.out, parms = "eps.avg")
traplot(coda.out, parms = "eps.star")
traplot(coda.out, parms = "gam")
traplot(coda.out, parms = "gam.star")
traplot(coda.out, parms = "diff")
traplot(coda.out, parms = "tau.Gam")
traplot(coda.out, parms = "sig.gam")
traplot(coda.out, parms = "tau.Eps")
traplot(coda.out, parms = "sig.eps")


# summarize
sum.out <- coda.fast(coda.out, OpenBUGS = FALSE)
sum.out$var <- row.names(sum.out)


beta.ind <- grep("beta", row.names(sum.out))
ggplot(sum.out[beta.ind,], aes(x = var, y = mean)) +
  geom_pointrange(aes(ymin = pc2.5, ymax = pc97.5)) +
  geom_hline(yintercept = 0, color = "red", lty = 2)

diff.ind <- grep("diff", row.names(sum.out))
ggplot(sum.out[diff.ind,], aes(x = var, y = mean)) +
  geom_pointrange(aes(ymin = pc2.5, ymax = pc97.5)) +
  geom_hline(yintercept = 0, color = "red", lty = 2)

prob.ind <- grep("prob", row.names(sum.out))
ggplot(sum.out[prob.ind,], aes(x = var, y = mean)) +
  geom_pointrange(aes(ymin = pc2.5, ymax = pc97.5))


# Random effects
labs <- c()
for(i in 1:9){
  labs[i] <- paste0("Plot ", i)
}
prob.eps <- grep("eps.star", row.names(sum.out))
df <- data.frame(sum.out[prob.eps,], block = rep(1:3, each = 3))
ggplot(df, aes(x = var, y = mean)) +
  geom_pointrange(aes(ymin = pc2.5, ymax = pc97.5, color = as.factor(block))) +
  geom_hline(yintercept = 0, color = "red", lty = 2) +
  scale_x_discrete(labels = labs)

labs <- c("Block 1", "Block 2", "Block 3")
prob.gam <- grep("gam.star", row.names(sum.out))
ggplot(sum.out[prob.gam,], aes(x = var, y = mean)) +
  geom_pointrange(aes(ymin = pc2.5, ymax = pc97.5)) +
  geom_hline(yintercept = 0, color = "red", lty = 2) +
  scale_x_discrete(labels = labs)

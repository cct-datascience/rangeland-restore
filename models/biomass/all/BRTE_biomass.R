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

# Report summary (median) cover for
# BRTE, native grass, native forb, intro forb, litter, bareground
# Litter was low, did not analyze
# Native forb was low, lumped with intro forb

library(tidyverse)

load("../cleaned_data/cover_all.Rdata")
dat <- cover_all %>%
  mutate(BRTE = BRTE/100,
         intro_forbs = intro_forbs/100,
         native_grass = native_grass/100,
         native_forbs = native_forbs/100,
         bryo = bryo/100,
         litter = litter/100,
         rock = rock/100,
         bare = bare/100)

# Calculate medians
median(dat$BRTE)
median(dat$intro_forbs)
median(dat$native_grass)
median(dat$native_forbs)
median(dat$bryo)
median(dat$litter)
median(dat$rock)
median(dat$bare)

# Cover results:
# Of plant cover across all plots, BRTE had the highest median cover at 30% while 
# introduced forbs had the second highest median cover at 18%. Median cover of native 
# grasses and forbs were 1% and 0%, respectively. Litter and rock also had very low 
# median cover at 0.5%, while the median cover of bare ground was 47%. 
# Due to low cover of native forbs, native and introduced forbs were combined 
# for analysis. Litter and rock cover were too low for analysis. 
# 
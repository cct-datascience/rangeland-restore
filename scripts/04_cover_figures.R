# Figures script for cover data at 3 levels
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggthemes)

### ALL
# read in data
load("../cleaned_data/cover_all.Rdata") # cover_all
# convert to proportions
raw_dat <- cover_all %>%
  mutate(BRTE = BRTE/100,
         intro_forbs = intro_forbs/100,
         native_grass = native_grass/100,
         native_forbs = native_forbs/100)
str(raw_dat)

# Load coda and coda.rep
load(file = "../models/cover/all/coda/coda.Rdata") # coda.out

# Summarize coda
# Note that tidyMCMC drops the deviance estimate
sum_out <- broom.mixed::tidyMCMC(coda.out, conf.int = TRUE, 
                                 conf.level = 0.95) %>% 
  rename(param = term, mean = estimate, sd = std.error, 
         pc2.5 = conf.low, pc97.5 = conf.high) %>% 
  mutate(sig = if_else(pc2.5 * pc97.5 > 0, TRUE, FALSE))

# Select and organize group means
model_dat <- sum_out %>%
  filter(grepl("m\\.", param)) %>%
  mutate(grazing = case_when(grepl("ungrazed", param) ~ "ungrazed",
                             grepl("fall", param) ~ "fall",
                             grepl("spring", param) ~ "spring"),
         fuelbreak = case_when(grepl("control", param) ~ "control",
                         grepl("herbicide", param) ~ "herbicide",
                         grepl("greenstrip", param) ~ "greenstrip"),
         grazing = factor(grazing, levels = c("ungrazed", "fall", "spring")),
         fuelbreak = factor(fuelbreak, levels = c("control", "herbicide", "greenstrip")))

fig5 <- ggplot() +
  geom_point(data = raw_dat, aes(x = fuelbreak, y = BRTE, color = grazing),
             position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.2),
             alpha = 0.2) +
  geom_pointrange(data = model_dat, aes(x = fuelbreak, y = mean,
                                        ymin = pc2.5,
                                        ymax = pc97.5, 
                                        color = grazing),
                  shape = 15,
                  size = 0.75,
                  position = position_dodge(width = 0.5)) +
  scale_y_continuous(expression(paste("BRTE ", m^-2))) +
  scale_color_canva(palette = "Surf and turf") +
  theme_bw(12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12),
        legend.position = c(.9, .85),
        legend.title = element_blank(),
        legend.background = element_rect(fill='transparent'))

jpeg(filename = "../plots/Fig5_cover_all.jpg",
     height = 4, width = 6,
     units = "in", res = 600)
print(fig5)
dev.off()

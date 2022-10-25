# Use PRISM for 30 year normals and compare to local station
library(prism)
library(raster)
library(dplyr)
library(stringr)
library(ggplot2)
library(lubridate)
library(cowplot)

# Create vector of 2 points that bound TS Ranch (from manuscript methods section)
df <- data.frame(site = c("TS1","TS2"),
                 lon = c(-116.509, -116.554),
                 lat = c(40.843, 40.895)) %>%
  tibble::column_to_rownames(var = "site")

# custom function to extract month from label
get_sixth <- function(vec) {
  return(as.numeric(vec[6]))
}

get_fifth <- function(vec) {
  foo <- vec[5]
  return(as.numeric(substr(foo, 5, 6)))
}

get_fifth_year <- function(vec) {
  foo <- vec[5]
  return(as.numeric(substr(foo, 1, 4)))
}

##### PRISM PPT 2019 #####
# Resolution is a 4km
# Set directory for prism data
prism_set_dl_dir("raw_data/prism_2019_ppt")
get_prism_monthlys(type = "ppt",
                   years = 2019,
                   mon = 1:12,
                   keepZip = FALSE)
fn <- prism_archive_ls()
# Stack data frame
RS_2019_ppt <- pd_stack(fn)

# Extract from stack
out_2019_ppt <- extract(RS_2019_ppt, df) %>%
  as.data.frame() %>%
  mutate(site = c("TS1","TS2")) %>%
  relocate(site) %>%
  tidyr::pivot_longer(-site, 
                      names_to = "label",
                      values_to = "ppt") 

# add numeric month column
out_2019_ppt$month <- unlist(lapply(str_split(out_2019_ppt$label, "_"), FUN = get_fifth))
str(out_2019_ppt)

##### PRISM Tmax 2019 #####
# Set directory for prism data
prism_set_dl_dir("raw_data/prism_2019_tmax")
get_prism_monthlys(type = "tmax",
                  years = 2019,
                  mon = 1:12,
                  keepZip = FALSE)
fn <- prism_archive_ls()
# Stack data frame
RS_2019_tmax <- pd_stack(fn)

# Extract from stack
out_2019_tmax <- extract(RS_2019_tmax, df) %>%
  as.data.frame() %>%
  mutate(site = c("TS1","TS2")) %>%
  relocate(site) %>%
  tidyr::pivot_longer(-site, 
                      names_to = "label",
                      values_to = "tmax") 

# add numeric month column
out_2019_tmax$month <- unlist(lapply(str_split(out_2019_tmax$label, "_"), FUN = get_fifth))
str(out_2019_tmax)

##### PRISM Tmin 2019 #####
# Set directory for prism data
prism_set_dl_dir("raw_data/prism_2019_tmin")
get_prism_monthlys(type = "tmin",
                   years = 2019,
                   mon = 1:12,
                   keepZip = FALSE)
fn <- prism_archive_ls()
# Stack data frame
RS_2019_tmin <- pd_stack(fn)

# Extract from stack
out_2019_tmin <- extract(RS_2019_tmin, df) %>%
  as.data.frame() %>%
  mutate(site = c("TS1","TS2")) %>%
  relocate(site) %>%
  tidyr::pivot_longer(-site, 
                      names_to = "label",
                      values_to = "tmin") 

# add numeric month column
out_2019_tmin$month <- unlist(lapply(str_split(out_2019_tmin$label, "_"), FUN = get_fifth))
str(out_2019_tmin)

##### Combine 2019 and summarize #####
out_2019_all <- out_2019_ppt %>%
  left_join(select(out_2019_tmax, -label)) %>%
  left_join(select(out_2019_tmin, -label)) %>%
  select(-label) %>%
  relocate(month, .after = site) %>%
  group_by(month) %>%
  summarize(ppt = mean(ppt),
            tmax = mean(tmax),
            tmin = mean(tmin)) %>%
  ungroup() %>%
  mutate(month = month(month, label = TRUE))

save(out_2019_all, file = "cleaned_data/prism_2019.Rdata")

##### PRISM PPT NORMALS #####
# Set directory for prism data
prism_set_dl_dir("raw_data/prism_ppt")
get_prism_normals(type = "ppt",
                  resolution = "4km",
                  mon = 1:12,
                  annual = FALSE,
                  keepZip = FALSE)
fn <- prism_archive_ls()
# Stack data frame
RS_ppt <- pd_stack(fn)

# Extract from stack
out_ppt <- extract(RS_ppt, df) %>%
  as.data.frame() %>%
  mutate(site = c("TS1","TS2")) %>%
  relocate(site) %>%
  tidyr::pivot_longer(-site, 
                      names_to = "label",
                      values_to = "ppt") 

# add numeric month column
out_ppt$month <- unlist(lapply(str_split(out_ppt$label, "_"), FUN = get_sixth))
str(out_ppt)

##### PRISM Tmax NORMALS #####
# Set directory for prism data
prism_set_dl_dir("raw_data/prism_tmax")
get_prism_normals(type = "tmax",
                  resolution = "4km",
                  mon = 1:12,
                  annual = FALSE,
                  keepZip = FALSE)
fn <- prism_archive_ls()
# Stack data frame
RS_tmax <- pd_stack(fn)

# Extract from stack
out_tmax <- extract(RS_tmax, df) %>%
  as.data.frame() %>%
  mutate(site = c("TS1","TS2")) %>%
  relocate(site) %>%
  tidyr::pivot_longer(-site, 
                      names_to = "label",
                      values_to = "tmax") 

# add numeric month column
out_tmax$month <- unlist(lapply(str_split(out_tmax$label, "_"), FUN = get_sixth))
str(out_tmax)

##### PRISM Tmin NORMALS #####
# Set directory for prism data
prism_set_dl_dir("raw_data/prism_tmin")
get_prism_normals(type = "tmin",
                  resolution = "4km",
                  mon = 1:12,
                  annual = FALSE,
                  keepZip = FALSE)
fn <- prism_archive_ls()
# Stack data frame
RS_tmin <- pd_stack(fn)

# Extract from stack
out_tmin <- extract(RS_tmin, df) %>%
  as.data.frame() %>%
  mutate(site = c("TS1","TS2")) %>%
  relocate(site) %>%
  tidyr::pivot_longer(-site, 
                      names_to = "label",
                      values_to = "tmin") 

# add numeric month column
out_tmin$month <- unlist(lapply(str_split(out_tmin$label, "_"), FUN = get_sixth))
str(out_tmin)




##### Combine normals and summarize #####
out_all <- out_ppt %>%
  left_join(dplyr::select(out_tmax, -label)) %>%
  left_join(dplyr::select(out_tmin, -label)) %>%
  dplyr::select(-label) %>%
  relocate(month, .after = site) %>%
  group_by(month) %>%
  summarize(ppt = mean(ppt),
            tmax = mean(tmax),
            tmin = mean(tmin)) %>%
  ungroup() %>%
  mutate(month = month(month, label = TRUE))

save(out_all, file = "cleaned_data/prism_normals.Rdata")


##### Plotting all PRISM data #####
load("cleaned_data/prism_normals.Rdata") # out_all
load("cleaned_data/prism_2019.Rdata") # out_2019_all

# plot separately
fig1a <- ggplot(out_all, aes(x = month)) +
  geom_bar(aes(y = ppt), stat = "identity",
           fill = "gray90") +
  geom_point(aes(y = tmin, color = "Tmin")) +
  geom_path(aes(y = tmin, color = "Tmin", group = 1)) +
  geom_point(aes(y = tmax, color = "Tmax")) +
  geom_line(aes(y = tmax, color = "Tmax", group = 2)) +
  geom_text(aes(x = 0.5, y = 35, label = "1991-2020 Normals"),
            hjust = 0, 
            vjust = 0) +
  scale_y_continuous(name = "Temp (°C)| Precip (mm)") +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank())

fig1b <- ggplot(out_2019_all, aes(x = month)) +
  geom_bar(aes(y = ppt), stat = "identity",
           fill = "gray90") +
  geom_point(aes(y = tmin, color = "Tmin")) +
  geom_path(aes(y = tmin, color = "Tmin", group = 1)) +
  geom_point(aes(y = tmax, color = "Tmax")) +
  geom_line(aes(y = tmax, color = "Tmax", group = 2)) +
  geom_text(aes(x = 0.5, y = 80, label = "2019"),
            hjust = 0, 
            vjust = 0) +
  scale_y_continuous(name = "Temp (°C)| Precip (mm)") +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank())

# fig1 <- plot_grid(fig1a, fig1b, ncol = 1, align = TRUE)
# ggsave(filename = "plots/Fig1_weather.jpg",
#        plot = fig1,
#        width = 6,
#        height = 6, 
#        dpi = 600)


# plot together
out_all$type <- "1991-2020 Normals"
out_2019_all$type <- "2019"

out_df <- rbind.data.frame(out_all, out_2019_all)

fig1 <- out_df %>%
  ggplot(aes(x = month)) +
  geom_bar(aes(y = ppt), stat = "identity",
           fill = "gray90") +
  geom_point(aes(y = tmin, color = "Tmin")) +
  geom_path(aes(y = tmin, color = "Tmin", group = 1)) +
  geom_point(aes(y = tmax, color = "Tmax")) +
  geom_line(aes(y = tmax, color = "Tmax", group = 2)) +
  facet_wrap(~type, ncol = 1) +
  scale_y_continuous(name = "Temp (°C)| Precip (mm)") +
  scale_color_manual(values = c("chocolate3", "cornflowerblue"),
                     labels = c(~T[max], ~T[min])) +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.8, 0.9),
        legend.key = element_rect(color = NA),
        strip.background = element_blank())


ggsave(filename = "plots/Fig1_weather.jpg",
       plot = fig1,
       width = 4,
       height = 6, 
       dpi = 600)

#### PRISM water year ppt 2014-2019 #####
# Resolution is a 4km
# Set directory for prism data
prism_set_dl_dir("raw_data/prism_wateryear_ppt")
get_prism_monthlys(type = "ppt",
                   years = 2013,
                   mon = 10:12,
                   keepZip = FALSE)
                   
get_prism_monthlys(type = "ppt",
                   years = 2014:2018,
                   mon = 1:12,
                   keepZip = FALSE)

get_prism_monthlys(type = "ppt",
                   years = 2019,
                   mon = 1:9,
                   keepZip = FALSE)

fn <- prism_archive_ls()
# Stack data frame
RS_wateryear_ppt <- pd_stack(fn)

# Extract from stack
out_wateryear_ppt <- extract(RS_wateryear_ppt, df) %>%
  as.data.frame() %>%
  mutate(site = c("TS1","TS2")) %>%
  relocate(site) %>%
  tidyr::pivot_longer(-site, 
                      names_to = "label",
                      values_to = "ppt") 

# add numeric month, year, and water year
out_wateryear_ppt$month <- unlist(lapply(str_split(out_wateryear_ppt$label, "_"), FUN = get_fifth))
out_wateryear_ppt$year <- unlist(lapply(str_split(out_wateryear_ppt$label, "_"), FUN = get_fifth_year))
out_wateryear_ppt$wyear <- rep(rep(2014:2019, each = 12), 2)
str(out_wateryear_ppt)

# summarize by water year
out_wateryear_sum <- out_wateryear_ppt %>%
  group_by(site, wyear) %>%
  summarize(ann_precip = sum(ppt)) %>%
  ungroup() %>%
  group_by(wyear) %>%
  summarize(ann_precip = mean(ann_precip))

#### PRISM water year ppt normals #####
prism_set_dl_dir("raw_data/prism_ppt")
fn <- prism_archive_ls()
# Stack data frame
RS_ppt <- pd_stack(fn)

# Extract from stack
out_ppt <- extract(RS_ppt, df) %>%
  as.data.frame() %>%
  mutate(site = c("TS1","TS2")) %>%
  relocate(site) %>%
  tidyr::pivot_longer(-site, 
                      names_to = "label",
                      values_to = "ppt") 

# add numeric month column
out_ppt$month <- unlist(lapply(str_split(out_ppt$label, "_"), FUN = get_sixth))
out_ppt$normal <- "normal"
str(out_ppt)

# summarize by water year
out_ppt_sum <- out_ppt %>%
  group_by(site, normal) %>%
  summarize(ann_precip = sum(ppt)) %>%
  ungroup() %>%
  group_by(normal) %>%
  summarize(ann_precip = mean(ann_precip))

# Plot water year annual precip + normals
fig0 <- ggplot(data = out_wateryear_sum,
       aes(x = wyear, y = ann_precip)) +
  geom_col() +
  geom_hline(yintercept = out_ppt_sum$ann_precip,
             lty = 2) +
  scale_y_continuous("Water year precipitation (mm)") +
  scale_x_continuous(breaks = 2014:2019) +
  theme_bw(base_size = 12) +
  theme(axis.title.x = element_blank(),
        panel.grid = element_blank())

ggsave(filename = "plots/FigS0_precip.jpg",
       plot = fig0,
       width = 6,
       height = 3, 
       dpi = 600)

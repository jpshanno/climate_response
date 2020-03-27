library(tidyverse)
library(jpshanno)
library(lubridate)
library(tsibble)
library(broom)

dat <- 
  read_csv("Data/QAQC/Continuous_Sampling/well_levels.csv",
           col_types = cols(site = "c", 
                            sample_time = "T",
                            water_level_m = "d", 
                            .default = col_skip())) %>% 
  filter(site == "157",
         year(sample_time) == 2014) %>% 
  as_tsibble(key = site,
             index = sample_time) %>% 
  index_by(sample_date = as_date(sample_time)) %>% 
  summarize(water_level_m = mean(water_level_m)) %>% 
  head(-3)

plot(dat, type = "l")

# dat <- 
#   left_join(dat,
#             read_csv("Documents/Event_Response/data/precip_combined_external_sources.csv") %>% 
#               filter(site == "113") %>% 
#               mutate(precip_m = precip_mm / 1000) %>% 
#               select(sample_date,
#                      precip_m),
#             by = "sample_date")

# For 2012 in 156 this gives a good fit, and for 2014 at 157. This is probably
# the right way to do it

dat_ts <- ts(data = dat$water_level_m, 
             start = c(year(dat$sample_date[1]), 
                       yday(dat$sample_date[1])), 
             frequency = 365)

mod <- 
  lm(dat_ts ~ TSA::harmonic(dat_ts, m = 2))

mod <- 
  lm(dat$water_level_m ~ cbind(yday(dat$sample_date), 
                               sin(2*pi*(yday(min(dat$sample_date)):yday(max(dat$sample_date)))/365), 
                               cos(2*pi*(yday(min(dat$sample_date)):yday(max(dat$sample_date)))/365),
                               sin(pi*(yday(min(dat$sample_date)):yday(max(dat$sample_date)))/365), 
                               cos(pi*(yday(min(dat$sample_date)):yday(max(dat$sample_date)))/365)))

plot(water_level_m ~ sample_date, data = dat, type = 'l')
points(predict(mod)~dat$sample_date)

plot(predict(mod))
mod %>% 
  augment(.) %>% 
  ggplot(data = .,
         aes(x = 1:nrow(.),
             y = dat_ts)) +
  geom_line() +
  geom_line(aes(y = .fitted))  

dat_resid <- 
  ts(resid(mod),
     start = c(year(dat$sample_date[1]), 
               yday(dat$sample_date[1])), 
     frequency = 365)

plot(dat_resid)
acf(dat_resid)
pacf(dat_resid)
TSA::eacf(dat_resid)


# This approach did not work for 2012 at 156, but did work for 2014 at 157
# It probably doesn't work because it should be done on a time series. The help
# page from spectrum() states that it calculates the density over a range dependent
# on the frequency of x (which will equal 1 for a non-time series vectro)
# https://stats.stackexchange.com/a/61032
# ssp <- spectrum(dat$water_level_m)  
# per <- 1/ssp$freq[ssp$spec==max(ssp$spec)]
# per
# 
# X <-
#   (pi / per) * yday(dat$sample_date)
# 
# # https://stats.stackexchange.com/a/60504  
# mod <-
#   lm(water_level_m ~ sin(2 * X) + cos(2 * X) + sin(4 * X) + cos(4 * X) + sin(6 * X) + cos(6 * X),
#      data = dat)
# 
# augment(mod) %>% 
#   ggplot(data = .,
#          aes(x = 1:nrow(.),
#              y = water_level_m)) +
#   geom_line() +
#   geom_line(aes(y = .fitted))
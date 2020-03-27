library(jpshanno)
library(data.table)
library(ggplot2)
library(tidyverse)

raw_wl <- 
  fread("data/well_levels.csv", 
        select = c(site = "character", 
                   sample_time = "character", 
                   water_level_m = "numeric"))

raw_wl[, sample_date := as.IDate(sample_time)]

daily_wl <- 
  raw_wl[, .(water_level_cm = 100*mean(water_level_m)), by = .(site, sample_date)]
setkey(daily_wl, site, sample_date)

daily_wl <- 
  daily_wl %>%
  group_by(site, field_season = year(sample_date)) %>% 
  filter(sum(is.na(water_level_cm)) < 30) %>% 
  ungroup() %>% 
  mutate(treatment_period = set_treatment_period(sample_date,
                                                 study = ifelse(site %in% c("053", "113"),
                                                                 "pws",
                                                                 "eco")),
         treatment = case_when(site %in% c("113", "135", "157", "152", "111") ~ "Control",
                               site %in% c("053", "077", "156", "009", "139") ~ "Ash Cut",
                               TRUE ~ "Girdle"),
         treatment = factor(treatment, levels = c("Control", "Girdle", "Ash Cut")))

probs <- 
  seq(0.05, 0.95, 0.05)
  
exceedence <- 
  daily_wl %>% 
  group_by(site) %>% 
  summarize(quantiles = list(pivot_wider(enframe(quantile(water_level_cm, 
                                                          probs = probs, 
                                                          na.rm = TRUE))))) %>% 
  unnest(quantiles) %>% 
  set_names(c("site", paste0("exceedance_", sprintf("%0.2f", 1-probs))))


ggplot(data = filter(daily_wl, site == "119"), 
       aes(x = water_level_cm, 
           linetype = treatment_period,
           # group = as.factor(field_season),
           color = treatment)) + 
  stat_ecdf() + 
  scale_y_continuous("Probability of Exceedence", 
                     labels = function(x){paste0(as.integer(100*(1 - x)), "%")}) + 
  facet_wrap(~site, scales = "free_x") + 
  theme_minimal()

write_csv(exceedence,
          "output/tabular/wetland_water_level_exceedence_probabilities.csv")

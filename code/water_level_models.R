# To Do
# Some sites have very high delta wl (139 in October). Some things to still 
# check still
# Remove precip below the precip_threshold and see if you can get better model
# fits


library(jpshanno)
library(data.table)
library(ggplot2)
library(patchwork)
library(quantreg)
library(sirad)
library(sf)
library(gamm4)

raw_wl <- 
  fread("data/well_levels.csv", 
        select = c(site = "character", 
                   sample_time = "character", 
                   water_level_m = "numeric"))

raw_wl[, sample_date := as.IDate(sample_time)]

daily_wl <- 
  raw_wl[, 
         .(water_level_cm = 100*mean(water_level_m)), 
         keyby = .(site, sample_date)]

precip <- 
  fread("data/precipitation_daily.csv",
        select = c(site = "character",
                   sample_date = "character",
                   precip_cm = "numeric"))

precip[, sample_date := as.IDate(sample_date)]

precip <- 
  precip[sample_date > as.IDate("2011/09/30")]

setkey(precip, site, sample_date)

site_info <- 
  fread("data/site_descriptions.csv",
        select = c("site", "siteType", "treatment", "area_ha", "NAD8316N_easting_m", "NAD8316N_northing_m"))

site_info <- 
  cbind(site_info,
        site_info %>% 
          st_as_sf(coords = c("NAD8316N_easting_m", "NAD8316N_northing_m"), 
                   crs = 26916) %>% 
          st_transform(4326) %>% 
          st_coordinates())

site_info[, c("NAD8316N_easting_m", "NAD8316N_northing_m") := NULL]

setnames(site_info, c("siteType", "X", "Y"), c("study", "lon", "lat"))

study_types <- 
  c(Planting = "planting",
    `Ecological Processes` = "eco",
    `Paired Watershed` = "pws")
  
site_info[, `:=`(site = expand_site(site),
                 study = study_types[study],
                 treatment = factor(treatment, levels = c("Control", "Girdle", "Ash Cut")))]
setkey(site_info, site)

# No precip data for 111 & 139

temps <- 
  fread("data/temperature_daily.csv")

temps[, `:=`(site = expand_site(site),
             sample_date = as.IDate(sample_date), 
             gdd_10 = (max_temp_c + min_temp_c) / 2 - 10)]
setkey(temps, site, sample_date)

water_balance <- 
  temps[precip[daily_wl]]

# Define Variablies -------------------------------------------------------

water_balance[site_info, 
              `:=`(treatment = treatment,
                   treatment_period = set_treatment_period(sample_date, study = study),
                   lon = lon,
                   lat = lat)]

water_balance[, `:=`(field_season = as.character(year(sample_date)),
                     doy = yday(sample_date),
                     sample_month = month(sample_date),
                     sample_week = isoweek(sample_date + 1), # isoweek is Monday-Sunday, `+ 1` makes the week Sun-Sat
                     delta_wl = water_level_cm - shift(water_level_cm),
                     precip_2day_cm = precip_cm + shift(precip_cm, 1),
                     precip_3day_cm = precip_cm + shift(precip_cm, 2) + shift(precip_cm, 1))]

# Create Day and Week of Season
water_balance[, `:=`(dos = doy - yday(paste0(unique(year(sample_date)), "-05-15")),
                     wos = sample_week - isoweek(paste0(unique(year(sample_date)), "-05-16"))), # isoweek is Monday-Sunday, `-05-16` makes the week Sun-Sat
              by = .(field_season)]
water_balance[, scaled_dos := dos / 100]

# Add PET -----------------------------------------------------------------

# Aridity bias is something that should be considered if we see very dry future
# conditions. [@valiantzas-2013; @shahidian-2012]
# Using Hargreaves for now (with no coefficient adjustments because I'm not 
# worried about systematically under or over estimating). However [@hargreaves-2003]
# recommends using it for periods of 5 or more days.

# Calculate extrasolar radiation
extrasolar <- 
  unique(water_balance[, .(site, lat, doy)])

extrasolar <- 
  extrasolar[, .(doy, 
                 rad_MJ_m2 = extrat(doy, lat)$ExtraTerrestrialSolarRadiationDaily), 
             by = .(site, lat)][, .(site, doy, rad_MJ_m2)]

water_balance[extrasolar, 
              rad_MJ_m2 := rad_MJ_m2,
              on = c("site", "doy")]

lambda_MJ_kg <- 
  2.45

water_balance[, pet_hs_cm := 0.1 * 0.0023 * rad_MJ_m2 / lambda_MJ_kg * (mean_temp_c + 17.8) * sqrt(max_temp_c - min_temp_c)]
water_balance[, pet_hm_cm := 0.1 * 0.0023 * rad_MJ_m2 / lambda_MJ_kg * (mean_temp_c + 17.8) * (max_temp_c - min_temp_c - 0.00123 * precip_cm) ** 0.76]

# Remove periods for now -----------------------------------------

# Snowmelt
water_balance <- 
  water_balance[!(month(sample_date) < 5 | (month(sample_date) == 5 & mday(sample_date) < 15))]

# Bad couple of days
water_balance <- 
  water_balance[!(site == "135" & sample_date >= as.IDate("2015-10-28") & sample_date < as.IDate("2015-11-01"))]

# 2016 has too many missing models for a hierarchical GAM, 2019 is not ready yet
water_balance <- 
  water_balance[!(field_season %in% c("2016", "2019"))]

empty_seasons <- 
  water_balance[, .(N = sum(!(is.na(delta_wl) | is.na(precip_2day_cm) | is.na(max_temp_c)))), by = .(site, field_season)][N == 0, .(site, field_season)]

water_balance <- 
  water_balance[!empty_seasons, on = .(site, field_season)]

rm(empty_seasons)


# Define Precipitation Threshold ------------------------------------------

# Get precip threshold (right now it is the precipication value where the mean
# water level response is zero. This should probably be done with a GLM 

# Site-specific precip thresholds:
precip_thresholds <- 
  lm(delta_wl ~ precip_2day_cm:site, 
     data = water_balance[precip_2day_cm > 0]) %>% 
  coef() %>% 
  {-.[1] / .[-c(1)]}

precip_thresholds <- 
  data.table(site = sub(".*([0-9]{3})$", "\\1", names(precip_thresholds)),
             precip_threshold = precip_thresholds, 
             key = "site")

# Fail to reject the null that the thresholds are normally distributed. So
# if we use this function, we could simulate a threshold from a normal 
# distribution to generalize the model from any given site. Or we could run
# simulations using specified thresholds
shapiro.test(precip_thresholds$precip_threshold)

# ggplot(data = water_balance[precip_2day_cm > 0], 
#        aes(x = precip_2day_cm, y = delta_wl)) + 
#   geom_point(color = "gray50") + 
#   geom_hline(aes(yintercept = 0),
#              color = "gray10",
#              linetype = "dotted") +
#   geom_vline(data = precip_thresholds, 
#              aes(xintercept = precip_threshold), 
#              color = "red",
#              linetype = "dashed") +
#   geom_quantile(quantiles = 0.5,
#                 method = "rq",
#                 formula = "y~x") +
#   facet_wrap(~site) +
#   theme_minimal()

water_balance[precip_thresholds, 
              `:=`(precip_occurrence_threshold = as.numeric(precip_2day_cm >= precip_threshold),
                   precip_occurrence_raw = as.numeric(precip_cm > 0))]


# Identify Drivers --------------------------------------------------------

# Simple answer approach would be to use a 1 or two day lagged precip model
# More complex would be to use a moving average for the effect of precip
# In addition theres some autoregressive behavior of water levels and an AR 
# relationship between water level and response to precip

# Need to test for difference in month responses within a site

# The best predictors are likely going to be max temperature for days without
# precip and 2 day precip for days with precip

# Precip
{{ggplot(data = water_balance[precip_cm > 0, 
                            .(correlation = cor(delta_wl, 
                                                precip_cm, 
                                                use = "pairwise.complete.obs",
                                                method = "spearman")), 
                            by = .(site, sample_month)],
       aes(x = as.factor(sample_month), 
           y = correlation)) + 
  geom_boxplot() +
  ggtitle("Same Day Precip")} +
  {ggplot(data = water_balance[precip_2day_cm > 0, 
                               .(correlation = cor(delta_wl, 
                                                   precip_2day_cm, 
                                                   use = "pairwise.complete.obs",
                                                   method = "spearman")), 
                               by = .(site, sample_month)],
          aes(x = as.factor(sample_month), 
              y = correlation)) + 
      geom_boxplot() +
      ggtitle("2 Day Precip")} +
  {ggplot(data = water_balance[precip_3day_cm > 0, 
                              .(correlation = cor(delta_wl, 
                                                  precip_3day_cm, 
                                                  use = "pairwise.complete.obs",
                                                  method = "spearman")), 
                              by = .(site, sample_month)],
         aes(x = as.factor(sample_month), 
             y = correlation)) + 
  geom_boxplot() +
  ggtitle("3 Day Precip")}} *
  theme_minimal() *
  xlab("Month") *
  coord_cartesian(ylim = c(0, 1)) +
  plot_annotation(tag_levels = "A")


# Temperature
{{ggplot(data = water_balance[precip_cm == 0, 
                              .(correlation = cor(delta_wl, 
                                                  min_temp_c, 
                                                  use = "pairwise.complete.obs",
                                                  method = "spearman")), 
                              by = .(site, sample_month)],
         aes(x = as.factor(sample_month), 
             y = correlation)) + 
    geom_boxplot() +
    ggtitle("Minimum Daily Temp")} +
    {ggplot(data = water_balance[precip_cm == 0, 
                                 .(correlation = cor(delta_wl, 
                                                     mean_temp_c, 
                                                     use = "pairwise.complete.obs",
                                                     method = "spearman")), 
                                 by = .(site, sample_month)],
            aes(x = as.factor(sample_month), 
                y = correlation)) + 
        geom_boxplot() +
        ggtitle("Mean Daily Temp")} +
    {ggplot(data = water_balance[precip_cm == 0, 
                                 .(correlation = cor(delta_wl, 
                                                     max_temp_c, 
                                                     use = "pairwise.complete.obs",
                                                     method = "spearman")), 
                                 by = .(site, sample_month)],
            aes(x = as.factor(sample_month), 
                y = correlation)) + 
        geom_boxplot() +
        ggtitle("Maximum Daily Temp")}} *
  theme_minimal() *
  xlab("Month") *
  coord_cartesian(ylim = c(-1, 0.5)) +
  plot_annotation(tag_levels = "A")

ggplot(data = water_balance[, 
                            .(correlation = cor(delta_wl, 
                                                pet_cm, 
                                                use = "pairwise.complete.obs",
                                                method = "spearman")), 
                            by = .(site, sample_month)],
       aes(x = as.factor(sample_month), 
           y = correlation)) + 
  geom_boxplot() +
  ggtitle("HS PET") +
  theme_minimal() +
  xlab("Month") +
  coord_cartesian(ylim = c(-1, 0.5))

# GAM Model ---------------------------------------------------------------

# library(mgcv)
# test <- water_balance[site == "151" & field_season == 2012]
# gammod_test <- gam(delta_wl ~ s(max_temp_c, by = sample_month) + s(precip_2day_cm, by = sample_month), data = test)
# test[, predicted := predict(gammod_test, newdata = test)]
# ggplot(data = test, 
#        aes(x = delta_wl, 
#            y = predicted, 
#            color = as.factor(sample_month))) + 
#   geom_point() + 
#   geom_abline(intercept = 0, slope = 1) +
#   facet_wrap(~sample_month)



# library(tidyverse)

# Fit a single model with no random effects, but a different family. It looks
# like a good fit, but doesn't seem to work out at the site/field season level.
# I think this approach is worth exploring

gam(delta_wl ~ 
      s(dos, 
        bs = "cc") + 
      # s(precip_2day_cm, precip_cm) + 
      # s(max_temp_c) + 
      site_numeric*field_season, 
    data = .) %>% 
  broom::augment() %>% 
  group_by(site_numeric, 
           field_season) %>% 
  mutate(.fitted = cumsum(.fitted)) %>%
  ggplot(data = ., 
         aes(y = .fitted, 
             x = dos, 
             color = as.factor(site_numeric))) + 
  geom_line(stat = "summary", 
            fun.y = mean) + 
  facet_wrap(~field_season)

# Using bs = 'cc' is good for water level, but meaningless for delta_wl. If I 
# wanted to use 'cc' to get a good cycle (not necessarily the best option, eg
# low spring water levels could be exceeded after wet summer or fall)

# It looks like this is the model to start working from to improve. Start working
# with the training data!

# Need to look specifically at 2012 to try and capture extreme dry periods, and
# should identify a season for wet conditions. Perhaps training/testing will not
# be chronologically split but instead should be based on extreme wet years, or
# years that had large events

gamm4_mod <- 
  gamm4(delta_wl ~ s(dos) + s(pet_cm, by = precip_occurrence_threshold) + s(precip_2day_cm), 
        random = ~(dos | site/field_season) + (pet_cm | site/field_season) + (precip_2day_cm | site/field_season), 
        data = water_balance)

plot(gamm4_mod[[2]], scale = 0, pages = 1)

water_balance %>% 
  mutate(.fitted = predict(gamm4_mod[[2]], newdata = water_balance)) %>% 
  group_by(site, field_season) %>% 
  mutate(delta_wl = cumsum(replace_na(delta_wl, 0)),
         .fitted = cumsum(replace_na(.fitted, 0))) %>% 
  ggplot(data = .,
         aes(x = doy,
             y = .fitted)) +
  geom_line() +
  geom_line(aes(y = delta_wl),
            linetype = "dotted") +
  facet_grid(field_season ~ site, 
             scales = "free") +
  theme_minimal()

gamm4_mod_wl <- 
  gamm4(water_level_cm ~ s(dos, bs = "cc") + s(pet_cm, by = precip_occurrence_threshold) + s(precip_2day_cm), 
        random = ~(dos | site/field_season) + (pet_cm | site/field_season) + (precip_2day_cm | site/field_season), 
        data = water_balance)

water_balance %>% 
  mutate(.fitted = predict(gamm4_mod_wl[[2]], newdata = water_balance)) %>% 
  ggplot(data = .,
         aes(x = doy,
             y = .fitted)) +
  geom_line() +
  geom_line(aes(y = water_level_cm),
            linetype = "dotted") +
  facet_grid(field_season ~ site, 
             scales = "free") +
  theme_minimal()

# mod_gam <-
  gam(delta_wl ~
        s(doy, sp = 10) +
        s(max_temp_c, sp = 1) +
        s(precip_2day_cm, sp = 100),
      family = "scat",
      data = water_balance)
# plot(mod_gam)
# gam.check(mod_gam)

# Fit a single model using site and field season as random effects
# mod_gam <-  
  # gam(delta_wl ~
  #       s(doy) +
  #       s(max_temp_c,
  #         precip_occurrence,
  #         sp = 1000,
  #         by = sample_month) +
  #       s(precip_2day_cm,
  #         precip_occurrence,
  #         water_level_cm,
  #         sp = 1000) +
  #       ti(max_temp_c, precip_2day_cm) +
  #       s(site_numeric, field_season, bs = "re"),
  #     data = water_balance)
# 
# water_balance[, c("delta_wl_pred", "se_pred") := predict(mod_gam, newdata = water_balance, se.fit = TRUE)]

test_dat <- copy(water_balance)
test_dat[delta_wl >= 25, delta_wl:=NA_real_]
test_mod <- gam(delta_wl ~ max_temp_c + precip_2day_cm, family = "scat", data = test_dat) %T>% plot(pages = 1, scale = 0)
test_dat[, delta_wl_pred := predict(test_mod, newdata = test_dat)]
test_dat %>% group_by(site, field_season) %>% mutate(cum_wl = cumsum(delta_wl_pred)) %>% plot(cum_wl ~ doy, data = .)
test_dat %>% group_by(site, field_season) %>% mutate(cum_delta_wl_pred = cumsum(delta_wl_pred), cum_delta_wl = cumsum(replace_na(delta_wl, 0))) %>% ggplot(data = ., aes(x = doy, y = cum_delta_wl_pred)) + geom_line() + geom_line(aes(y = cum_delta_wl), col = "gray50") + facet_grid(site ~ field_season, scales = "free") + theme_bw()
# Looks like I have seasonal drawdown figured out, but not the rebound

# Fit separate models for each level of site:field_season 
gams <-
  water_balance %>%
  select(-sample_date) %>% 
  group_nest(site, field_season) %>%
  mutate(mod = map(data,
                   ~  gam(delta_wl ~
                            s(doy) +
                            s(max_temp_c, sp = 1000, by = sample_month) +
                            s(precip_2day_cm, sp = 1000, by = sample_month) +
                            ti(max_temp_c, precip_2day_cm),
                          # family = "scat",
                          data = .x))) %>%
  mutate(pred = map2(mod, data,
                     ~predict(.x,
                              newdata = .y,
                              se.fit = TRUE) %>%
                       as_tibble() %>%
                       set_names("delta_wl_pred", "pred_se")))

gam_models <- 
  gams %>%
  select(site, field_season, mod)

# This looks good. Just need to get starting levels better. Need to find 
# predictors for intial seasonal values. Also need to incorporate ARMA residual
# structure

gams <-
  gams %>% 
  select(-mod) %>% 
  unnest(cols = c(data, pred)) %>%
  mutate(sample_date = as.IDate(sample_date_numeric)) %>% 
  group_by(site, field_season) %>% 
  # mutate(predicted_water_level_cm = quantile(water_level_cm, probs = 0.95, na.rm = TRUE)) %>% 
  mutate(predicted_water_level_cm = ifelse(sample_date == min(sample_date[!is.na(water_level_cm)]),
                                           water_level_cm, NA_real_)) %>%
  fill(predicted_water_level_cm, .direction = "down") %>% # Uses observed starting value for each year
  # Uses median of starting values for each site
  # group_by(site) %>% 
  # mutate(predicted_water_level_cm = median(predicted_water_level_cm, na.rm = TRUE)) %>% 
  # group_by(site, field_season) %>% 
  mutate(delta_wl_pred = ifelse(sample_date == min(sample_date),
                                0, delta_wl_pred),
         predicted_water_level_cm = predicted_water_level_cm + cumsum(delta_wl_pred),
         sample_date = lubridate::date(as.IDate(sample_date)))


ggplot(data = gams,
       aes(x = delta_wl,
           y = delta_wl_pred,
           color = as.factor(precip_occurrence))) + 
  geom_point() +
  geom_abline(aes(intercept = 0,
                  slope = 1)) +
  facet_wrap(~site, 
             scales = "free")

ggplot(gams, 
       aes(x = delta_wl, 
           y = delta_wl_pred, 
           color = site)) + 
  geom_point() + 
  geom_abline(aes(slope = 1, 
                  intercept = 0)) + 
  facet_wrap(~sample_month,
             scale = "free")

ggplot(gams, 
       aes(x = delta_wl, 
           y = scale(delta_wl_pred - delta_wl), 
           color = site)) + 
  geom_point() +
  facet_wrap(~sample_month, scales = "free")


plot(density(gams$delta_wl_pred - gams$delta_wl, na.rm = TRUE))
qqnorm(gams$delta_wl_pred - gams$delta_wl); qqline(gams$delta_wl_pred - gams$delta_wl)



{ggplot(data = filter(gams, field_season == 2017),
       aes(x = sample_date)) +
  geom_line(aes(y = water_level_cm),
            linetype = "dotted") +
  geom_line(aes(y = predicted_water_level_cm)) +
  facet_wrap(~site,
             scales = "free") +
  theme_bw()}



# Potential GAM Models ----------------------------------------------------
# These could all be blocked by month
# s(doy) + s(max_temp_c) + s(precip_2day_cm)
# s(doy) + s(max_temp_c, precip_occurrence) + s(precip_2day_cm, water_level_cm)

# Potential Distrubitons: (found via lmomco)
# test_vals <- water_balance[!is.na(delta_wl), delta_wl]
# ggplot() + geom_density(aes(x = test_vals)) + geom_density(aes(x = rlmomco(1000, lmr2par(x = test_vals, type = "aep4"))), col = "red")
# aep4 <- works with raw data, negative values are okay
# exp
# emu, exp, gam, gep, gev, gld, glo, gno, gov, gpa, gum, kap, kmu, kur, lap, lmrq, ln3, nor, pe3, ray, revgum, rice, sla, st3, texp, wak, or wei.

# Water year might not be that useful for two reasons:
#   1. We don't have data from fall 2011 (I don't think, have to check), or 
#      the fall before paired watersheds
#   2. Not currently keeping data past the end of October because of 
#      inconsistencies with ice in last winter

water_balance[, water_year := ifelse(month(sample_date) >= 10, field_season +1, field_season)]
water_balance[, water_doy := julian(sample_date, origin = min(sample_date)), by = .(site, water_year)]

# water_balance[, precip_2day := frollsum(precip_cm, 2, align = "right"), by = .(site, field_season)]
# water_balance[, precip_lag := shift(precip_cm, 1, fill = NA_real_), by = .(site, field_season)]

# library(vlbuildr)
# vl_chart(data = water_balance) %>% 
#   vl_mark_line(tooltip = "data") %>% 
#   vl_encode_x("sample_date:T") %>% 
#   vl_encode_y("water_level_m:Q") %>% 
#   vl_encode_color("site:N") %>% 
#   vl_resolve_axis_x(how = "independent") %>% 
#   vl_facet_wrap(field = "water_year:Q")

# Create harmonic covariates <- include other covariates here
harmonics <- 
  water_balance[, .(seq_date = as.IDate(as.IDate("2012-01-01"):as.IDate("2019-12-31")))]

harmonics[, `:=`(no_harm = 1,
                 sin1 = 100 * sin(pi*as.integer(seq_date)/max(yday(seq_date))), 
                 sin2 = 100 * sin(2*pi*as.integer(seq_date)/max(yday(seq_date))), 
                 sin4 = 100 * sin(4*pi*as.integer(seq_date)/max(yday(seq_date))), 
                 cos1 = 100 * cos(pi*as.integer(seq_date)/max(yday(seq_date))), 
                 cos2 = 100 * cos(2*pi*as.integer(seq_date)/max(yday(seq_date))), 
                 cos4 = 100 * cos(4*pi*as.integer(seq_date)/max(yday(seq_date)))), 
          by = .(year(seq_date))]

setkey(harmonics, seq_date)

water_balance <- 
  water_balance[harmonics, 
                on = .(sample_date = seq_date),
                nomatch = 0]

control_train <- 
  water_balance[treatment_period == "Pre-treatment" | (treatment == "control" & field_season <= 2017), ]

control_test <- 
  water_balance[treatment == "control" & field_season > 2017, ]

# Need to test if Ash Cut & Girdle can be pooled like this
treatment_train <- 
  water_balance[treatment %in% c("girdle", "ash cut") & treatment_period == "Post-treatment" & field_season <= 2017]

treatment_test <- 
  water_balance[treatment %in% c("girdle", "ash cut") & field_season > 2017]

# Check that all datasets are unique, should be all zeroes
combn(x = c("control_train", "control_test", "treatment_train", "treatment_test"), 
      m = 2, 
      FUN = function(x){nrow(fintersect(get(x[1]), get(x[2]), all = TRUE))})

ggplot(data = control_train, 
       aes(x = yday(sample_date),
           y = water_level_cm, 
           color = as.factor(field_season))) + 
  geom_line() + 
  facet_wrap(~site)


ggplot(data = control_train, 
       aes(x = yday(sample_date),
           y = c(NA_real_, diff(water_level_cm)), 
           color = as.factor(field_season))) + 
  geom_line() + 
  facet_wrap(~site)

# Could try with PET in mm (which helps some, but does not completely solve the 
# problem). Could also try to use max temp C, which is of a larger scale. The
# appropriate way would be put them all on a standard scale and then weight the
# drivers within the model. PET may abe better than max temp, unless max temp is
# set to have an intercept of 0

control_train[, pet_mm := pet_cm / 10]
control_train[, site := as.factor(site)]
control_train[, field_season := as.factor(field_season)]

gamm4_mod <- 
  gamm4(delta_wl ~ s(scaled_dos) + t2(pet_cm, precip_2day_cm), 
        random = ~(scaled_dos | site/field_season) + (pet_cm | site/field_season) + (precip_2day_cm | site/field_season), 
        data = control_train)

plot(gamm4_mod[[2]], scale = 0, pages = 1)

library(tidyverse)

control_train %>% 
  mutate(.fitted = predict(gamm4_mod[[2]], newdata = control_train)) %>% 
  group_by(site, field_season) %>% 
  mutate(delta_wl = cumsum(replace_na(delta_wl, 0)),
         .fitted = cumsum(replace_na(.fitted, 0))) %>% 
  ggplot(aes(x = doy,
             y = .fitted)) +
  geom_line() +
  geom_line(aes(y = delta_wl),
            linetype = "dotted") +
  facet_grid(field_season ~ site, 
             scales = "free") +
  theme_minimal()

control_train[dos == 0]

# Model & extract residuals

control_mods <- 
  control_train[, .(mod = list(lm(water_level_cm ~ cbind(sin2, cos2, sin4, cos4)))),
                by = .(site, field_season)]

control_mods <- 
  control_mods[harmonics[, .(harmonics = list(data.table(seq_date, .SD))),
                         by = .(field_season = year(seq_date))], 
               on = .(field_season = field_season),
               nomatch = 0]
                                                
control_pred <- 
  control_mods[, list(sample_date = harmonics[[1]]$seq_date, 
                      predicted_wl = predict(mod[[1]], newdata = harmonics[[1]])), 
               by = .(site, field_season)]

control_pred <- 
  merge(control_pred,
        control_train[, .(site, sample_date, water_level_cm)],
        by = c("site", "sample_date"), all = TRUE)

ggplot(data = control_pred,
       aes(x = yday(sample_date), 
           y = predicted_wl, 
           color = as.factor(field_season))) +
  geom_line() +
  geom_line(aes(y = water_level_cm),
            linetype = "dashed") +
  # coord_cartesian(ylim = c(-1, 1)) +
  facet_wrap(~site, scales = "free")

# Look for autocorrelation in residuals (almost certain to be there)

control_train[, .(acf_1 = acf(water_level_m, plot = FALSE)$acf[4]), by = .(site, field_season)]
control_train[, .(list(lm(water_level_m ~ sample_date))), by = .(site, field_season)]
control_train[, diff(water_level_m), by = .(site, field_season)]

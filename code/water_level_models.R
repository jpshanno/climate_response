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

# Using final value from day, not mean
daily_wl <- 
  raw_wl[, 
         .(water_level_cm = 100*tail(water_level_m, 1)), 
         keyby = .(site, sample_date)]

precip <- 
  fread("data/precipitation_daily.csv",
        select = c(site = "character",
                   sample_date = "character",
                   precip_cm = "numeric",
                   filled_data = "logical",
                   fill_precip = "numeric",
                   fill_coefficient = "numeric"))

precip[, sample_date := as.IDate(sample_date)]
precip[, field_season := year(sample_date)]
precip[, cum_precip_cm := cumsum(precip_cm), 
       by = .(site, field_season)]

precip[, field_season := NULL]

precip <- 
  precip[sample_date > as.IDate("2011/09/30")]

setkey(precip, site, sample_date)

# Add snowmelt data to precip

snowmelt <- 
  fread("data/daily_snowmelt.csv", 
        select = c(site = "character",
                   sample_date = "character",
                   melt_cm = "numeric"))

snowmelt[, sample_date := as.IDate(sample_date)]

setkey(snowmelt, site, sample_date)

precip[snowmelt, 
       melt_cm := melt_cm]

precip[, melt_cm := nafill(melt_cm, "const", 0)]

precip[, total_I := precip_cm + melt_cm]

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

# write.csv(site_info, "tmp/site_info_with_coordinates.csv", row.names = FALSE)

# No precip data for 111 & 139

temps <- 
  fread("data/temperature_daily.csv")

temps[, `:=`(site = expand_site(site),
             sample_date = as.IDate(sample_date), 
             doy = yday(as.IDate(sample_date)),
             gdd_10 = (max_temp_c + min_temp_c) / 2 - 10,
             range_temp_c = ifelse(is.na(range_temp_c), max_temp_c - min_temp_c, range_temp_c))]
setkey(temps, site, sample_date)

temps[, gdd_10_base0 := pmax(gdd_10, 0)]

# Add PET -----------------------------------------------------------------

# Aridity bias is something that should be considered if we see very dry future
# conditions. [@valiantzas-2013; @shahidian-2012]
# Using Hargreaves for now (with no coefficient adjustments because I'm not 
# worried about systematically under or over estimating). However [@hargreaves-2003]
# recommends using it for periods of 5 or more days.

# Consider creating a calibration for PET using the White method from [@diamond-2018]

# Calculate extrasolar radiation
extrasolar <- 
  unique(site_info[, .(lat, doy = 1:365), by = .(site)])

extrasolar <- 
  extrasolar[, .(doy, 
                 rad_MJ_m2 = extrat(doy, lat)$ExtraTerrestrialSolarRadiationDaily), 
             by = .(site, lat)][, .(site, doy, rad_MJ_m2)]

temps[extrasolar, 
      rad_MJ_m2 := rad_MJ_m2,
      on = c("site", "doy")]

lambda_MJ_kg <- 
  2.45

# Original HS equation taken from [@hargreaves-2003]
temps[, pet_hs_cm := 0.1 * 0.0023 * rad_MJ_m2 / lambda_MJ_kg * (mean_temp_c + 17.8) * sqrt(max_temp_c - min_temp_c)]

# Need to fix this to add Precip to temps ---------------------------------

# Modified HS equation taken from [@droogers-2002]
temps[, pet_hm_cm := 0.1 * 0.0023 * rad_MJ_m2 / lambda_MJ_kg * (mean_temp_c + 17.8) * (max_temp_c - min_temp_c - 0.00123 * precip_cm) ** 0.76]


ggplot(data = temps, 
       aes(x = sample_date, 
           y = pet_hs_cm)) + 
  geom_line() + 
  facet_wrap(~year(sample_date), 
             scales = "free") + 
  geom_vline(aes(xintercept = as.IDate("2013-08-15")))


temps <- 
  temps[sample_date > as.IDate("2011/09/30")]

temps[, field_season := year(sample_date)]
temps[, `:=`(cum_gdd = cumsum(gdd_10),
             cum_gdd_base0 = cumsum(gdd_10_base0),
             cum_pet_cm = cumsum(pet_hs_cm)),
      by = .(site, field_season)]

temps[, `:=`(filled_data = NULL, 
             doy = NULL)]

met <- 
  precip[temps]

met[, `:=`(water_availability = cum_precip_cm - cum_pet_cm,
           p_pet = cum_precip_cm / cum_pet_cm)]

met[month(sample_date) %in% 5:9, 
    water_availability_mjjas := cumsum(pet_hs_cm) - cumsum(gdd_10), 
    by = .(site, field_season)]

met[, field_season := NULL]

water_balance <- 
  daily_wl[met]

# Water levels may start to rebound at the inflection point of cumulative GDD. 
# Look at the date of that inflection point vs the timing of rebound of water 
# levels.
ggplot(data = met, 
       aes(x = sample_date, 
           y = cum_gdd_base0,
           group = site)) + 
  geom_line() + 
  facet_wrap(~year(sample_date), 
             scales = "free") + 
  geom_vline(aes(xintercept = as.IDate("2013-08-15")))


# Define Variablies -------------------------------------------------------

water_balance[site_info, 
              `:=`(treatment = treatment,
                   treatment_period = set_treatment_period(sample_date, study = study))]

water_balance[, `:=`(field_season = year(sample_date),
                     field_season_fct = as.factor(year(sample_date)),
                     site_fct = as.factor(site),
                     doy = yday(sample_date),
                     sample_month = month(sample_date),
                     sample_week = isoweek(sample_date + 1))] # isoweek is Monday-Sunday, `+ 1` makes the week Sun-Sat

# Create Day and Week of Season
water_balance[, `:=`(dos = 1 + doy - yday(paste0(unique(year(sample_date)), "-05-15")),
                     wos = 1 + sample_week - isoweek(paste0(unique(year(sample_date)), "-05-16"))), # isoweek is Monday-Sunday, `-05-16` makes the week Sun-Sat
              by = .(field_season)]
water_balance[, `:=`(scaled_dos = dos / 100,
                     scaled_wos = wos / 100)]

# Add lagged water level and delta_wl

water_balance[, `:=`(lagged_wl = shift(water_level_cm, 1, type = "lag"),
                     delta_wl = water_level_cm - shift(water_level_cm),
                     lagged_precip_2day = shift(precip_cm + shift(precip_cm, 1)),
                     precip_2day_cm = precip_cm + shift(precip_cm, 1),
                     precip_3day_cm = precip_cm + shift(precip_cm, 1) + shift(precip_cm, -1)),
              by = .(site, field_season)]

# Remove periods for now -----------------------------------------

# 2016 has too many missing models for a hierarchical GAM, 2019 is not ready yet
water_balance <- 
  water_balance[!(field_season %in% c(2016, 2019))]

empty_seasons <- 
  water_balance[, .(N = sum(!(is.na(water_level_cm) | is.na(precip_2day_cm) | is.na(max_temp_c)))), 
                by = .(site, field_season)][N == 0, .(site, field_season)]

water_balance <- 
  water_balance[!empty_seasons, on = .(site, field_season)]

rm(empty_seasons)

water_balance[, site_year := interaction(site_fct, field_season_fct)]

# Define Precipitation Threshold ------------------------------------------

# Get precip threshold (right now it is the precipication value where the mean
# water level response is zero. This should probably be done with a GLM 

# Site-specific precip thresholds:
precip_thresholds <- 
  lm(delta_wl ~ precip_3day_cm:site, 
     data = water_balance[precip_3day_cm > 0]) %>% 
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

ggplot(data = water_balance[precip_3day_cm > 0],
       aes(x = precip_3day_cm, y = delta_wl)) +
  geom_point(color = "gray50") +
  geom_hline(aes(yintercept = 0),
             color = "gray10",
             linetype = "dotted") +
  geom_vline(data = precip_thresholds,
             aes(xintercept = precip_threshold),
             color = "red",
             linetype = "dashed") +
  geom_quantile(quantiles = 0.5,
                method = "rq",
                formula = "y~x") +
  facet_wrap(~site) +
  theme_minimal()

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


# Temperature & PET
{{{ggplot(data = water_balance[precip_cm == 0, 
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
        ggtitle("Maximum Daily Temp")}} / 
    {{ggplot(data = water_balance[, 
                                  .(correlation = cor(delta_wl, 
                                                      pet_hs_cm, 
                                                      use = "pairwise.complete.obs",
                                                      method = "spearman")), 
                                  by = .(site, sample_month)],
             aes(x = as.factor(sample_month), 
                 y = correlation)) + 
        geom_boxplot() +
        ggtitle("HS PET")} +
        {ggplot(data = water_balance[, 
                                     .(correlation = cor(delta_wl, 
                                                         pet_hm_cm, 
                                                         use = "pairwise.complete.obs",
                                                         method = "spearman")), 
                                     by = .(site, sample_month)],
                aes(x = as.factor(sample_month), 
                    y = correlation)) + 
            geom_boxplot() +
            ggtitle("HM PET")} + 
        plot_spacer()}} *
  theme_minimal() *
  xlab("Month") *
  coord_cartesian(ylim = c(-1, 0.5)) +
  plot_annotation(tag_levels = "A")


# Look at scaled values
scaled <- 
  water_balance[, lapply(.SD, function(x){if(!is.numeric(x)) return(x);as.numeric(scale(x))}), 
                .SDcols = patterns("temp_c|delta_wl|precip.*cm|pet|site$")]

scaled <- 
  melt(scaled, 
       id.vars = c("site", "delta_wl"), 
       variable.name = "driver")[!is.na(delta_wl) & !is.na(value)]

scaled_driver_summary <- 
  water_balance[, lapply(.SD, function(x){if(!is.numeric(x)) return(x);as.numeric(scale(x))}), 
                .SDcols = patterns("temp_c|delta_wl|precip.*cm|pet|site$")] %>% 
  melt(id.vars = c("site", "delta_wl"), 
       variable.name = "driver") %>% 
  .[!is.na(delta_wl) & !is.na(value)] %>% 
  .[, .(linear_slope = abs(coef(lm(delta_wl ~ value))[2]),
        abs_correlation = abs(cor(value, delta_wl, method = "spearman")),
        gam_r2 = summary(gam(delta_wl ~ s(value)))$r.sq), 
    by = .(driver, site)] %>% 
  .[, class := ifelse(grepl("precip", driver), "input", "loss")]

{{ggplot(aes(y = gam_r2, x = driver), 
         data = scaled_driver_summary) +
    geom_boxplot()} /
    {ggplot(aes(y = abs_correlation, x = driver), 
            data = scaled_driver_summary) +
        geom_boxplot()} /
    {ggplot(aes(y = linear_slope, x = driver), 
            data = scaled_driver_summary) +
        geom_boxplot()}} *
  aes(color = class) *
  plot_annotation(title = "Scaled Drivers")

unscaled_driver_summary <- 
  melt(water_balance, 
       id.vars = c("site", "delta_wl"), 
       measure.vars = patterns("temp_c|precip.*cm|pet"), 
       variable.name = "driver") %>% 
  .[!is.na(delta_wl) & !is.na(value)] %>% 
  .[,.(linear_slope = abs(coef(lm(delta_wl ~ value))[2]), 
       abs_correlation = abs(cor(value, delta_wl, method = "spearman")), 
       gam_r2 = summary(gam(delta_wl ~ s(value)))$r.sq), 
    by = .(driver, site)] %>% 
  .[, class := ifelse(grepl("precip", driver), "input", "loss")]

{ggplot(aes(y = gam_r2, x = driver), 
        data = unscaled_driver_summary) +
    geom_boxplot()} /
  {ggplot(aes(y = abs_correlation, x = driver), 
          data = unscaled_driver_summary) +
      geom_boxplot()} /
  {ggplot(aes(y = linear_slope, x = driver), 
          data = unscaled_driver_summary) +
      geom_boxplot()} *
  aes(color = class) *
  plot_annotation(title = "Unscaled Drivers")




# Potential GAM Models ----------------------------------------------------
# These could all be blocked by month
# s(doy) + s(max_temp_c) + s(precip_2day_cm)
# s(doy) + s(max_temp_c, precip_occurrence) + s(precip_2day_cm, water_level_cm)

# library(vlbuildr)
# vl_chart(data = water_balance) %>% 
#   vl_mark_line(tooltip = "data") %>% 
#   vl_encode_x("sample_date:T") %>% 
#   vl_encode_y("water_level_m:Q") %>% 
#   vl_encode_color("site:N") %>% 
#   vl_resolve_axis_x(how = "independent") %>% 
#   vl_facet_wrap(field = "water_year:Q")

# Split Data --------------------------------------------------------------

control_train <- 
  water_balance[treatment_period == "Pre-treatment" | 
                  (treatment == "Control" & field_season <= 2017), ]

control_test <- 
  water_balance[treatment == "Control" & 
                  field_season > 2017, ]

# Need to test if Ash Cut & Girdle can be pooled like this
treatment_train <- 
  water_balance[treatment %in% c("Girdle", "Ash Cut") & 
                  treatment_period == "Post-treatment" & 
                  field_season <= 2017]

treatment_test <- 
  water_balance[treatment %in% c("Girdle", "Ash Cut") & 
                  field_season > 2017]

# Check that all datasets are unique, should be all zeroes
combn(x = c("control_train", "control_test", "treatment_train", "treatment_test"), 
      m = 2, 
      FUN = function(x){nrow(fintersect(get(x[1]), get(x[2]), all = TRUE))})

ggplot(data = control_train, 
       aes(x = yday(sample_date),
           y = water_level_cm, 
           color = field_season_fct)) + 
  geom_line() + 
  facet_wrap(~site)


ggplot(data = control_train, 
       aes(x = dos,
           y = delta_wl, 
           color = field_season_fct)) + 
  geom_line() + 
  facet_wrap(~site)



# Model Selection ---------------------------------------------------------

# Should create smooth plots of different drivers (dos, wos, pet, precip, etc)
# for site/field season groups to evaluate the need for GS, GI, S, I models
# as in [@pedersen-2019, p. 31]

# Use t2(full = TRUE) to generate tensor products as they found oversmoothing of
# the global function when using te()

# So family = scat is a bad choice! The residual plots look good 

# Potential Model Improvements:
# - Treat P as a intervention in the time series sense with an autoregressive 
#   effect (at least in growing season)
# - Check for autocorrelation of residuals
# - Include random effects
# - Look at cumulative dry days as a driver
# - Improve wet year performance (77 2013)
# - Evaluate fit of final models under 
# - GDD looks to be the best driver for seasonal trend of water levels
# - Attempt to fit a seasonal trend using just GDD (or in combo with another 
#   driver) and use Arima() with external predictors to fit a model
# - Check out tensor product of PET and water stress
# - What about a two step model of Δ WL? First a logistic model to determine if
#   the water table moved up or down, and then a second model to determine the
#   magnitude? See "Variable Selection" in [@larsen-2015]
# - Try water year water_availability calculations to see if the extra snow helps
# - What about looking at a recession curve for the summer and the opposite for
#   the fall? It probably is not as a good a plan as the water balance, because it
#   will be making more assumptions about the hydrologic system staying stable

arma_orders <- 
  water_balance[, as.list(forecast::arimaorder(forecast::auto.arima(water_level_cm))), 
                by = .(site, field_season)]


week_knots <- 
  c(min(water_balance$sample_week) - 1, max(water_balance$sample_week) + 1)


# CHECKOUT THESE (include unused levels)

# Water stress seems to be have a hystersis loop for each precip event and
# for the entire season. If I can figure out how to model that loop then
# I could probably really improve the model. For the invididual events
# the loop could probably be improved by creating another precip variable
# to add into the water stress or the model. I have to think about where
# it would have the effect of removing the loops. It would be a decay 
# applied to each storm (likely an exponential decay, look at the 
# intervention analysis stuff from MA 5781). Perhaps there is a way to 
# include cumulative GDD to remove the seasonal hystersis, maybe by finding
# the inflection point of the cumulative GDD curve? It would mean fitting a
# second order polynomial to the daily GDD and then finding the root

# Attempting to find threshold for model to change
met[, .(gdd_10 = min(gdd_10)), 
    by = .(site, field_season, sample_week)][
      gdd_10 > 0, 
      .(first_week = min(sample_week),
        last_week = max(sample_week)), 
      by = .(site)]

# Can I use the hystersis loops (and changes in them) to say something 
# about changes in wetland processes after treatment?

# Need to find the trigger for water level reversal. It may be around when GDD
# peaks. If that's the case then 

# Generate climate time series
# Calculate GDD
# Fit GDD curve & find root
# Apply water level model that includes a threhold or switch for maximum GDD 
# timing. Would better PET estimates just make this a moot point as water stress
# would be more representative of the real water balance drivers?

ggplot(aes(y = water_level_cm,
           x = gdd_10, 
           color = as.factor(sample_month)), 
       data = water_balance[site == "140" & field_season == 2017, 
                            lapply(.SD, mean, na.rm = TRUE),
                            by = .(site, field_season, sample_month, sample_week),
                            .SDcols = c("water_level_cm", "gdd_10")]) + 
  geom_path() + 
  facet_grid(site ~ field_season, scales = "free")

  water_balance[, .(water_level_cm, doy, water_availability, gdd_10_7day = roll_mean(gdd_10, n = 7, fill = NA, na.rm = TRUE)), by = .(site, field_season)][, .(water_level_cm, doy,gdd_10_7day, water_availability, gdd_peak = ifelse(gdd_10_7day <= max(gdd_10_7day, na.rm = TRUE), "Pre", "Post")), by = .(site, field_season)] %>% ggplot(aes(x = doy, y = gdd_10_7day)) + geom_line() + facet_wrap(~site + field_season)

ggplot(aes(y = water_level_cm,
           x = water_availability, 
           color = ifelse(doy > 275, "Post", "Pre")), 
       data = water_balance[site %in% c("135", "152", "157", "111", "113")]) + 
  geom_path() + 
  facet_grid(site ~ field_season, scales = "free")

# Look at using MJJAS water stress; gets rid of some loops

ggplot(aes(y = water_level_cm,
           x = water_availability, 
           color = as.factor(sample_month)), 
       data = water_balance) + 
  geom_path() + 
  facet_grid(site ~ field_season, scales = "free")

ggplot(aes(y = water_level_cm,
           x = water_availability_mjjas, 
           color = field_season_fct), 
       data = water_balance) + 
  geom_path() + 
  facet_wrap(~site, scales = "free")

mod_water_availability_gs_siteyear <- 
  gam(water_level_cm ~ 
        s(water_availability) + 
        s(water_availability, 
          by = site_year, 
          bs = "fs"),
      control = list(nthreads = 4) , 
      drop.unused.levels = FALSE,
      data = control_train); beepr::beep(7)

mod_water_availability_gs_site <- 
  gam(water_level_cm ~ 
        s(water_availability_mjjas) + 
        s(water_availability_mjjas, 
          by = site_fct, 
          bs = "fs") +
        s(site_fct, 
          bs = "re"),
      control = list(nthreads = 4) , 
      drop.unused.levels = FALSE,
      data = control_train)


ggplot(control_train, 
       aes(x = doy, 
           y = water_level_cm)) + 
  geom_line() + 
  geom_line(linetype = "dotted", 
            aes(y = predict(test_mod, 
                            newdata = control_train))) + 
  facet_wrap(~field_season +site, 
             scales = "free")



ggplot(aes(y = water_level_cm,
           x = cum_gdd, 
           color = as.factor(sample_month)), 
       data = water_balance) + 
  geom_point() + 
  facet_grid(site ~ field_season, scales = "free")

ggplot(aes(y = water_level_cm,
           x = doy), 
       data = water_balance) + 
  geom_line() + 
  # geom_line(aes(y = water_availability2),
  #               linetype = "dashed") +
  geom_line(aes(y = water_availability),
            linetype = "dotted") +
  facet_grid(site ~field_season, scales = "free_y")

test_delta <- 
  gam(delta_wl ~ 
        s(gdd_10) +
        t2(gdd_10, 
           site_year, 
           bs = "fs", 
           full = TRUE), 
      data = control_train)

test_mod <- gam(water_level_cm ~ 
                  s(water_availability) + 
                  s(gdd_10) +
                  t2(water_availability, site_year, bs = "fs", full = TRUE) + 
                  s(site_year, bs = "re"), 
                  # s(precip_2day_cm), 
                data = control_train,
                drop.unused.levels = FALSE)
control_train[, pred_wl := predict(test_mod, newdata = .SD), 
              by = .(site, field_season)]
ggplot(data = control_train, 
       aes(x = doy)) + 
  geom_line(aes(y = water_level_cm),
            color = "gray60") + 
  geom_line(aes(y = pred_wl), linetype = "dotted") + 
  facet_wrap(~site + field_season, scales = "free_y")

control_test[, pred_wl := predict(test_mod, newdata = .SD), 
              by = .(site, field_season)]
ggplot(data = control_test, 
       aes(x = doy)) + 
  geom_line(aes(y = water_level_cm),
            color = "gray60") + 
  geom_line(aes(y = pred_wl - 460), linetype = "dotted") + 
  facet_wrap(~site + field_season, scales = "free_y")


delta_g <- 
  gam(delta_wl ~
        s(pet_hs_cm, by = sample_week) + 
        s(precip_2day_cm, by = sample_week),
      method = "REML", 
      data = control_train)

level_g_week <- 
  gam(water_level_cm ~ 
        te(wos, lagged_wl,
           bs = c("cc", "tp")) +
        s(wos,
          bs = c("cc")) +
        te(precip_2day_cm,
           pet_hs_cm),
      knots = list(wos = week_knots),
      method = "REML",
      data = control_train,
      control = list(nthreads = 4))

# level_g_te <- 
#   gam(water_level_cm ~ te(wos, pet_hs_cm, precip_2day_cm, lagged_wl, 
#                           bs = c("cc", "tp", "tp", "tp")), 
#       knots = list(wos = week_knots),
#       method = "REML",
#       data = control_train,
#       control = list(nthreads = 4))

level_g_gdd <- 
  gam(water_level_cm ~ 
        te(cum_gdd, lagged_wl) +
        s(cum_gdd) +
        te(precip_2day_cm,
           pet_hs_cm),
      method = "REML",
      data = control_train,
      control = list(nthreads = 4))

level_g_ws <- 
  gam(water_level_cm ~ 
        te(water_availability, lagged_wl) +
        s(water_availability) +
        te(precip_2day_cm,
           pet_hs_cm),
      method = "REML",
      data = control_train,
      control = list(nthreads = 4))

level_g_gamm <- 
  gamm(water_level_cm ~ 
        s(cum_gdd) +
        te(precip_2day_cm,
           pet_hs_cm),
      method = "REML",
      data = control_train,
      control = list(nthreads = 4),
      correlation = corARMA(form = ~ 1 | site_year, p = 1, q = 1))

# This should be allowed to peak above max for 1 day
predict_delta_wl <- 
  function(mod, df){
    # stopifnot(nrow(df) == length(wl))
    # stopifnot(any(is.na(d)))
    
    wl <- 
      df$water_level_cm
    
    max_wl <- 
      # max(wl, na.rm = TRUE)
      quantile(wl, probs = 0.99, na.rm = TRUE, names = FALSE)
    
    min_ind <- min(which(!is.na(wl)))
    
    wl_out <- 
      wl[min_ind]
    
    new_wl <- 
      numeric(nrow(df))
    
    new_wl[min_ind] <- wl_out
    
    for(i in (min_ind):nrow(df)){
      
      max_diff <- 
        max_wl - wl_out
      
      d <- 
        predict(mod, newdata = df[i, ])
      
      d_adj <- 
        min(c(d, max_diff))
      
      wl_out <- 
        wl_out + d_adj
      
      new_wl[i] <- wl_out
    }
    new_wl
  }


predict_lagged_wl <- 
  function(mod, df){
    stopifnot(all(c("water_level_cm", "precip_2day_cm", "lagged_wl", "sample_week", "pet_hs_cm") %in% names(control_train)))
    
    wl <- 
      df$water_level_cm
    
    max_wl <- 
      # max(wl, na.rm = TRUE)
      quantile(wl, probs = 0.99, na.rm = TRUE, names = FALSE)
    
    min_ind <- min(which(!is.na(wl)))
    
    wl_out <- 
      wl[min_ind]
    
    new_wl <- 
      numeric(nrow(df))
    
    new_wl[min_ind] <- wl_out
    
    for(i in (min_ind + 1):nrow(df)){
      
      new_data <- 
        df[i, ]
      
      new_data$lagged_wl <- 
        wl_out
      
      wl_out <- 
        predict(mod, newdata = new_data)
      
      new_wl[i] <- 
        min(c(wl_out, max_wl))
    }
    new_wl
  }


# Evaluate Models ---------------------------------------------------------

###### 
# Improve predict_delta_wl() to allow temporary exceedence of max_wl

# # Predict water level change
# control_train[, pred_delta_wl := predict(delta_g, newdata = control_train)]
# # Predict wetland water levels
# control_train[, .fitted_delta := predict_delta_wl(pred_delta_wl, water_level_cm), 
#               by = .(field_season, site)]

control_train[, .fitted_level_week := predict_lagged_wl(level_g_week, .SD), 
              by = .(site, field_season_fct)]

control_train[, .fitted_level_gdd := predict_lagged_wl(level_g_gdd, .SD), 
              by = .(site, field_season_fct)]

control_train[, .fitted_level_ws := predict_lagged_wl(level_g_ws, .SD), 
              by = .(site, field_season_fct)]

control_train[, .fitted_level_gamm := predict_lagged_wl(level_g_gamm$gam, .SD), 
              by = .(site, field_season_fct)]

ggplot(data = control_train,
       aes(x = dos)) +
  geom_line(aes(y = .fitted_level_week,
                color = "Level (Week)"),
            linetype = "dashed") +
  geom_line(aes(y = .fitted_level_gdd,
                color = "Level GDD"),
            linetype = "dashed") +
  geom_line(aes(y = .fitted_level_ws,
                color = "Level WS"),
            linetype = "dashed") +
  geom_line(aes(y = .fitted_level_gamm,
                color = "Level GAMM"),
            linetype = "dashed") +
  geom_line(aes(y = water_level_cm,
                color = "Observed"),
            linetype = "solid") +
  scale_color_manual(breaks = c("Level (Week)", "Level GDD", "Level WS", "Level GAMM", "Observed"),
                     values = c("red", "blue", "green", "orange", "gray60")) +
  facet_wrap(~field_season + site,
             scales = "free_y")

# Predict water level change
control_test[, .fitted_predict := predict(level_g_gdd, newdata = .SD)]

control_test[, .fitted_level_week := predict_lagged_wl(level_g_week, .SD), 
              by = .(site, field_season_fct)]

control_test[, .fitted_level_gdd := predict_lagged_wl(level_g_gdd, .SD), 
              by = .(site, field_season_fct)]

control_test[, .fitted_level_ws := predict_lagged_wl(level_g_ws, .SD), 
              by = .(site, field_season_fct)]

ggplot(data = control_test,
       aes(x = dos)) +
  # geom_line(aes(y = .fitted_predict,
  #               color = "Raw Prediction"),
  #           linetype = "dotted") +
  geom_line(aes(y = .fitted_level_week,
                color = "Level (Week)"),
            linetype = "dashed") +
  geom_line(aes(y = .fitted_level_gdd,
                color = "Level GDD"),
            linetype = "dashed") +
  geom_line(aes(y = .fitted_level_ws,
                color = "Level WS"),
            linetype = "dashed") +
  geom_line(aes(y = water_level_cm,
                color = "Observed"),
            linetype = "solid") +
  scale_color_manual(breaks = c("Level (Week)", "Level GDD", "Level WS", "Observed"),
                     values = c("red", "blue", "green", "gray60")) +
  facet_wrap(~field_season + site,
             scales = "free_y")


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

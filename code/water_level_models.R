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

# write.csv(site_info, "tmp/site_info_with_coordinates.csv", row.names = FALSE)

# No precip data for 111 & 139

temps <- 
  fread("data/temperature_daily.csv")

temps[, `:=`(site = expand_site(site),
             sample_date = as.IDate(sample_date), 
             gdd_10 = (max_temp_c + min_temp_c) / 2 - 10,
             range_temp_c = ifelse(is.na(range_temp_c), max_temp_c - min_temp_c, range_temp_c))]
setkey(temps, site, sample_date)

water_balance <- 
  temps[precip[daily_wl]]

# Define Variablies -------------------------------------------------------

water_balance[site_info, 
              `:=`(treatment = treatment,
                   treatment_period = set_treatment_period(sample_date, study = study),
                   lon = lon,
                   lat = lat)]

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
                     precip_3day_cm = precip_cm + shift(precip_cm, 2) + shift(precip_cm, 1)),
              by = .(site, field_season)]
              
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

# Original HS equation taken from [@hargreaves-2003]
water_balance[, pet_hs_cm := 0.1 * 0.0023 * rad_MJ_m2 / lambda_MJ_kg * (mean_temp_c + 17.8) * sqrt(max_temp_c - min_temp_c)]

# Modified HS equation taken from [@droogers-2002]
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
  water_balance[!(field_season %in% c(2016, 2019))]

empty_seasons <- 
  water_balance[, .(N = sum(!(is.na(delta_wl) | is.na(precip_2day_cm) | is.na(max_temp_c)))), 
                by = .(site, field_season)][N == 0, .(site, field_season)]

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

# Could try with PET in mm (which helps some, but does not completely solve the 
# problem). Could also try to use max temp C, which is of a larger scale. The
# appropriate way would be put them all on a standard scale and then weight the
# drivers within the model. PET may abe better than max temp, unless max temp is
# set to have an intercept of 0

# Use t2(full = TRUE) to generate tensor products as they found oversmoothing of
# the global function when using te()

# Looks like the water level model may overfit (even if the test data looks good)
# There's no observable treatment effect (which is interesting). The delta_wl 
# model may be better, but I kind of think I have to find a way to use precip as
# an intervention rather than as a predictor. I can capture the sharp spike, but
# not the rapid decline. Maybe I need a 'days since rain' predictor or something

week_knots <- 
  c(min(water_balance$scaled_wos) - 0.01, max(water_balance$scaled_wos) + 0.01)

mod_g <- 
  gam(water_level_cm ~ te(scaled_wos, pet_hs_cm, precip_2day_cm, lagged_wl, 
                          bs = c("cc", "tp", "tp", "tp")), 
      knots = list(scaled_wos = week_knots),
      method = "REML",
      data = control_train)

# Took about 20 minutes to fit, but better fit of extremes compared to Gaussian
# family. Should probably do scaled sample_week
mod_g_scat <- 
  gam(water_level_cm ~ te(scaled_wos, pet_hs_cm, precip_2day_cm, lagged_wl, 
                          bs = c("cc", "tp", "tp", "tp")), 
      knots = list(scaled_wos = week_knots),
      method = "REML",
      family = scat,
      data = control_train,
      control = list(nthreads = 4))

mod_g_dwl_full <- 
  gam(delta_wl ~ te(sample_week, pet_hs_cm, precip_2day_cm, water_level_cm), 
      method = "REML",
      data = control_train)

mod_g_dwl_only_pet <- 
  gam(delta_wl ~ te(sample_week, pet_hs_cm), 
      method = "REML",
      data = control_train)

gmm_full <- 
  gamm4(delta_wl ~ s(dos) + t2(pet_hs_cm, precip_2day_cm), 
        random = ~(dos | site/field_season) + (pet_hs_cm | site/field_season) + (precip_2day_cm | site/field_season), 
        data = control_train)

gmm_scaled_dos <- 
  gamm4(delta_wl ~ s(scaled_dos) + t2(pet_hs_cm, precip_2day_cm), 
        random = ~(scaled_dos | field_season_fct/site_fct) + 
          (pet_hs_cm | field_season_fct/site_fct) + 
          (precip_2day_cm | field_season_fct/site_fct), 
        data = control_train)

# This doesn't work quite right because it can't predict 2018 in the test data
gam_scaled_dos <- 
  gam(delta_wl ~ s(scaled_dos) + t2(pet_hs_cm, precip_2day_cm, by = field_season_fct) + s(site_fct, bs = "re"), 
      method = "REML", 
      drop.unused.levels = FALSE,
      data = control_train)

gmm_scaled_dos_arma <- 
  gamm(delta_wl ~ s(scaled_dos) + t2(pet_hs_cm, precip_2day_cm), 
       random = list(field_season_fct = ~scaled_dos + pet_hs_cm +precip_2day_cm,
                     site_fct = ~scaled_dos + pet_hs_cm +precip_2day_cm), 
       data = control_train,
       correlation = corARMA(form = ~scaled_dos|site/field_season, p = 1, q = 1))

gmm_full_cross <- 
  gamm4(delta_wl ~ t2(scaled_dos, pet_hs_cm, precip_2day_cm), 
        random = ~(scaled_dos | site/field_season) + (pet_hs_cm | site/field_season) + (precip_2day_cm | site/field_season), 
        data = control_train)

gmm_lagged_wl <- 
  gamm4(delta_wl ~ t2(lagged_wl, pet_hs_cm, precip_2day_cm), 
        random = ~(lagged_wl | site/field_season) + (pet_hs_cm | site/field_season) + (precip_2day_cm | site/field_season), 
        data = control_train)

gmm_drivers_only <- 
  gamm4(delta_wl ~ t2(pet_hs_cm, precip_2day_cm), 
        random = ~ (pet_hs_cm | site/field_season) + (precip_2day_cm | site/field_season), 
        data = control_train)

gmm_1day_precip <- 
  gamm4(delta_wl ~ t2(pet_hs_cm, precip_cm), 
        random = ~ (pet_hs_cm | site/field_season) + (precip_cm | site/field_season), 
        data = control_train)

par(mfcol = c(2, 3), ps = 24)
plot(gmm_full[[2]], cex = 3, scale = 0)
plot(gmm_scaled_dos[[2]], scale = 0)
plot(gmm_drivers_only[[2]], scale = 0)
par(mfcol = c(2, 3))
plot(gmm_full[[2]], scale = 0)
plot(gmm_scaled_dos[[2]], scale = 0)
plot.new()
plot(gmm_drivers_only[[2]], scale = 0)
par(mfcol = c(1, 1))


# Generating predicted water levels needs to account for maximum wetland water
# levels. Will have to figure out a way to represent maximum water level in the
# future wetlands

# This should be allowed to peak above max for 1 day
predict_wl <- 
  function(d, wl){
    stopifnot(length(d) == length(wl))
    # stopifnot(any(is.na(d)))
    
    max_wl <- 
      # max(wl, na.rm = TRUE)
      quantile(wl, probs = 0.99, na.rm = TRUE, names = FALSE)
    
    min_ind <- min(which(!is.na(wl)))
    
    wl_out <- 
      wl[min_ind]
    
    new_wl <- 
      numeric(length(d))
    
    new_wl[min_ind] <- wl_out
    
    for(i in (min_ind + 1):length(d)){
      
      max_diff <- 
        max_wl - wl_out
      
      d_adj <- 
        min(c(d[i], max_diff))
      
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
# Improve predict_wl() to allow temporary exceedence of max_wl

# Predict water level change
control_train[, pred_delta_wl := predict(mod_g_dwl_full, newdata = control_train)]
# Predict wetland water levels
control_train[, pred_water_level_cm := predict_wl(pred_delta_wl, water_level_cm), 
              by = .(field_season, site)]

ggplot(data = control_train,
       aes(x = dos)) +
  geom_line(aes(y = pred_water_level_cm)) +
  geom_line(aes(y = water_level_cm),
            linetype = "dotted") +
  facet_wrap(~field_season + site,
             scales = "free_y")

control_test[, pred_delta_wl := predict(mod_g_dwl, newdata = control_test)]
# Predict wetland water levels
control_test[, pred_water_level_cm := predict_wl(pred_delta_wl, water_level_cm), 
              by = .(field_season, site)]

ggplot(data = control_test,
       aes(x = dos)) +
  geom_line(aes(y = pred_water_level_cm)) +
  geom_line(aes(y = water_level_cm),
            linetype = "dotted") +
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

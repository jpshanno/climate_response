source("code/climate_packages.R")


# Functions ---------------------------------------------------------------

weighted_mean <- 
  function(x, w){
    if(sum(is.na(x)) == length(x)){
      return(NA_real_)
    }
    
    na_x <- 
      which(is.na(x))
    
    sum(x * w, na.rm = TRUE) / sum(ifelse(is.na(x), NA, w), na.rm = TRUE)
  }

sum_na <- 
  function(x, class.out = class(x)){
    if(sum(is.na(x)) == length(x)){
      res <- NA
      class(res) <- class.out
      return(res)
    }
    
    sum(x, na.rm = TRUE)
  }

mean_na <- 
  function(x, class.out = class(x)){
    if(sum(is.na(x)) == length(x)){
      res <- NA
      class(res) <- class.out
      return(res)
    }
    
    mean(x, na.rm = TRUE)
  }

median_na <- 
  function(x, class.out = class(x)){
    if(sum(is.na(x)) == length(x)){
      res <- NA
      class(res) <- class.out
      return(res)
    }
    
    median(x, na.rm = TRUE)
  }

max_na <- 
  function(x, class.out = class(x)){
    if(sum(is.na(x)) == length(x)){
      res <- NA
      class(res) <- class.out
      return(res)
    }
    
    max(x, na.rm = TRUE)
  }

min_na <- 
  function(x, class.out = class(x)){
    if(sum(is.na(x)) == length(x)){
      res <- NA
      class(res) <- class.out
      return(res)
    }
    
    min(x, na.rm = TRUE)
  }

# External Met Data -------------------------------------------------------

# Constant for PET calculation
lambda_MJ_kg <- 
  2.45

# This precip is NOT hourly precip (I think it's annual accumulated, have to
# check processing scripts)
ex_met <-
  cbind(fread("../Data/Raw/Downloaded/mesowest_met/WKFM4.2019-12-31.csv",
              select = c(1, 2, 4, 5, 6, 8, 9, 10, 11, 13),
              skip = 12,
              col.names = c("station_id", "sample_time", "air_temperature_c",
                            "relative_humidity", "wind_speed_m_s",
                            "wind_gust_m_s", "solar_rad_w_m2", "precip_mm",
                            "wind_peak_m_s", "dew_point_temperature_c")),
        setnames(data.table(t(fread("../Data/Raw/Downloaded/mesowest_met/WKFM4.2019-12-31.csv",
                                    skip = 6,
                                    nrows = 3,
                                    header = FALSE,
                                    sep = ":")[, 2])),
                 c("latitude", "longitude", "elevation_m")))

ex_met[, sample_time := sample_time - 60]

ex_met[, `:=`(sample_time = setattr(sample_time, "tzone", "EST"),
              elevation_m = elevation_m / 3.2808,
              solar_rad_MJ_m2_hr = solar_rad_w_m2 * 3600 * 1e-6,
              solar_rad_w_m2 = NULL)]

ex_met[, etr_cm_hr := 0.1 * water::hourlyET(data.frame(wind = wind_speed_m_s,
                                                       RH = relative_humidity,
                                                       temp = air_temperature_c,
                                                       radiation = solar_rad_MJ_m2_hr,
                                                       height = 6.1,
                                                       lat= latitude,
                                                       long = longitude,
                                                       elev = elevation_m),
                                            DOY = yday(sample_time),
                                            hours = hour(sample_time),
                                            ET = "ETr",
                                            long.z = longitude)]

ex_met[, etr_cm_hr := pmax(0, etr_cm_hr)]

# Get min and max temp by using the previous and next hour
# ex_met[, pet_hs_cm_hr := 0.1 * 0.0023 * solar_rad_MJ_m2_hr / lambda_MJ_kg * (air_temperature_c + 17.8) * sqrt(air_temperature_c + sd(c(shift(air_temperature_c, 1), air_temperature_c, shift(air_temperature_c, -1)), na.rm = TRUE) - (air_temperature_c - sd(c(shift(air_temperature_c, 1), air_temperature_c, shift(air_temperature_c, -1)), na.rm = TRUE))),
ex_met[, pet_hs_cm_hr := 0.1 * 0.0023 * solar_rad_MJ_m2_hr / lambda_MJ_kg * (air_temperature_c + 17.8) * sqrt(max_na(c(shift(air_temperature_c, 1), air_temperature_c, shift(air_temperature_c, -1))) - min_na(c(shift(air_temperature_c, 1), air_temperature_c, shift(air_temperature_c, -1)))),
       by = .(year(sample_time))]
ex_met[, pet_hs_cm_hr := pmax(0, pet_hs_cm_hr)]

# Add in Hargreaves PET

daily_ex_met <- 
  ex_met[, .(mean_temp_c = mean_na(air_temperature_c), 
             min_temp_c = min_na(air_temperature_c), 
             max_temp_c = max_na(air_temperature_c), 
             observed_rad_MJ_m2 = sum_na(solar_rad_MJ_m2_hr), 
             etr_cm_d = sum_na(etr_cm_hr),
             pet_hs_cm_summed = sum_na(pet_hs_cm_hr)), 
         by = .(sample_date = as.Date(sample_time, tz = "EST"))]

ha_coefs <- 
  hacal(lat = 46.440278,
        days = daily_ex_met$sample_date,
        rad_mea = daily_ex_met$observed_rad_MJ_m2,
        tmax = daily_ex_met$max_temp_c,
        tmin = daily_ex_met$min_temp_c)

daily_ex_met[, doy := yday(sample_date)]
daily_ex_met[, ext_rad_MJ_m2 := extrat(doy, radians(46.440278))$ExtraTerrestrialSolarRadiationDaily]
daily_ex_met[, ha_rad_MJ_m2 := ha(sample_date, 
                                  lat = 46.440278, 
                                  lon = -89.826667, 
                                  ext_rad_MJ_m2, 
                                  A = ha_coefs[["Ha"]], 
                                  B = ha_coefs[["Hb"]], 
                                  Tmax = max_temp_c, 
                                  Tmin = min_temp_c)]
daily_ex_met[, 
             `:=`(exrad_pet_hs_cm = 0.1 * 0.0023 * ext_rad_MJ_m2 / lambda_MJ_kg * (mean_temp_c + 17.8) * sqrt(max_temp_c - min_temp_c), 
                  harad_pet_hs_cm = 0.1 * 0.0023 * ha_rad_MJ_m2 / lambda_MJ_kg * (mean_temp_c + 17.8) * sqrt(max_temp_c - min_temp_c),
                  obsrad_pet_hs_cm = 0.1 * 0.0023 * observed_rad_MJ_m2 / lambda_MJ_kg * (mean_temp_c + 17.8) * sqrt(max_temp_c - min_temp_c))]



# Load Data ---------------------------------------------------------------

precip <- 
  fread("data/precipitation_daily.csv",
        select = c(site = "character",
                   sample_date = "Date",
                   filled_precip_cm = "numeric"),
        key = c("site", "sample_date"))

setnames(precip, "filled_precip_cm", "precip_cm")

# Fill in missing precip where possible (should be a weighted mean by each site)

precip_means <- 
  precip[, .(mean_precip_cm = mean_na(precip_cm)), by = .(sample_date)]

# Fill a few wayward NAs from 2013
precip_means[between(month(sample_date), 6, 10) & is.na(mean_precip_cm),
             mean_precip_cm := 0]

precip[precip_means,
       precip_cm := ifelse(!is.na(precip_cm), precip_cm, i.mean_precip_cm),
       on = "sample_date"]

hourly_precip <- 
  fread("data/precipitation_hourly.csv",
        select = c(site = "character",
                   sample_time = "character",
                   precip_cm = "numeric"))

hourly_precip[, sample_time := ymd_hms(sample_time, tz = "EST")]

setkey(hourly_precip, site, sample_time)

sub_hourly_precip <- 
  fread("data/precipitation_fifteen_minutes.csv",
        select = c(site = "character",
                   sample_time = "character",
                   precip_cm = "numeric"))

sub_hourly_precip[, sample_time := ymd_hms(sample_time, tz = "EST")]

setkey(sub_hourly_precip, site, sample_time)

melt <- 
  fread("data/daily_snowmelt.csv", 
        select = c(site = "character", 
                   sample_date = "Date", 
                   melt_cm = "numeric"),
        key = c("site", "sample_date"))

# Need melt for all sites. For not just using max predicted melt
melt <- 
  melt[, .(melt_cm = mean_na(melt_cm)), by = .(sample_date)]

# I haven't accounted for some melt periods having interception removed
precip[melt,
       `:=`(melt_cm = i.melt_cm,
            precip_cm = precip_cm + i.melt_cm),
       on = "sample_date"]

water_levels <- 
  fread("data/well_levels.csv",
        select = c(site = "character",
                   sample_time = "character",
                   water_level_cm = "numeric",
                   air_temperature_c = "numeric",
                   data_source = "character"))

water_levels[, sample_time := ymd_hms(sample_time, tz = "EST")]
water_levels[, sample_date := as.Date(sample_time, tz = "EST")]
water_levels[, sample_year := year(sample_date)]


setkey(water_levels, "site", "sample_date")

water_balance <- 
  water_levels[precip,
               precip_cm := i.precip_cm]

water_balance[hourly_precip,
              h_precip_cm := i.precip_cm,
              on = c("site", "sample_time")]

water_balance[sub_hourly_precip,
              sh_precip_cm := i.precip_cm,
              on = c("site", "sample_time")]

water_balance[, 
             Ds_cm := shift(water_level_cm, -1) - water_level_cm,
             by = .(site, sample_year)]



site_info <- 
  fread("data/site_info.csv", 
        select = c("site", "study", "treatment"), 
        colClasses = "character",
        key = "site")

water_balance[site_info,
             `:=`(treatment = i.treatment,
                  study = i.study)]

water_balance[, treatment_period := set_treatment_period(sample_date, study = study)]

# Shorthand for control sites and pre-treatment treated sites
water_balance[, site_status := ifelse(treatment_period == "Pre-treatment" | site == "Control",
                                      "Control",
                                      treatment)]

# Set growing season
water_balance[, season := fifelse(month(sample_date) %in% 6:7 |
                                    (month(sample_date) == 5 & mday(sample_date) > 15) |
                                    (month(sample_date) == 9 & mday(sample_date) <= 15),
                                  "growing",
                                  "dormant")]

daily_water_balance <- 
  water_balance[,.(site_status = first(site_status),
                   season = first(season),
                   treatment = first(treatment),
                   i_water_level_cm = first(water_level_cm),
                   f_water_level_cm = last(water_level_cm),
                   min_water_level_cm = min_na(water_level_cm),
                   median_water_level_cm = median_na(water_level_cm),
                   max_water_level_cm = max_na(water_level_cm),
                   Dl_cm = max_na(water_level_cm) - min_na(water_level_cm),
                   Dl_signed_cm = sign(which.max(water_level_cm) - which.min(water_level_cm)) * (max_na(water_level_cm) - min_na(water_level_cm)),
                   Dl_time_range_s = (which.max(water_level_cm) - which.min(water_level_cm)) * 900,
                   gross_precip_cm = first(precip_cm),
                   observed_precip_cm = sum_na(h_precip_cm),
                   Ds_cm = last(water_level_cm) - first(water_level_cm)),
                by = .(site, sample_date)]

daily_water_balance[, best_precip_cm := fcoalesce(observed_precip_cm, gross_precip_cm)]

interception <- 
  daily_water_balance[best_precip_cm > 0,
                      .(treatment = first(treatment),
                        mod = list(lmrob(Ds_cm ~ best_precip_cm,
                                         data = .SD,
                                         setting = "KS2014"))),
                      by = .(site, site_status, season)]

interception[, c("intercept", "slope") := map_dfr(mod, coef)]
interception[, i_cm := -intercept / slope]

interception[, i_cm := ifelse(is.na(i_cm), mean_na(i_cm), i_cm),
          by = .(site_status)]

interception[i_cm < 0, i_cm := 0]


daily_water_balance[interception, 
                    `:=`(i_cm = i.i_cm,
                         net_precip_cm = pmax(0, best_precip_cm - i.i_cm)),
                    on = c("site", "site_status", "season")]

daily_water_balance[, precip_intensity_cm_hr := net_precip_cm / (Dl_time_range_s / 3600)]

# ESy ------------------------------------------------------------------
# To consider:
# Combining multi-day precip events
# -Removing days with preceeding precipitation- Checked this, no discernible patterns
# Adding an upper limit to ESy derived as the ratio of watershed:wetland area
# Could include dormant season data as well

sy_dat <- 
  daily_water_balance[net_precip_cm > 0, 
                      .(site,
                        season,
                        sample_date,
                        water_level_cm = min_water_level_cm,
                        i_cm = first(i_cm),
                        Ds_cm,
                        Dl_signed_cm,
                        Dl_time_range_s,
                        net_precip_cm,
                        precip_intensity_cm_hr)]

sy_dat[, sy := net_precip_cm / Dl_signed_cm]

# Less than 1.1 to remove problem points from bad precip estimates
# I was clipping out sy > 1.1 because there were some bad precip estimates. I
# have found and corrected most of those and only have ~2 questionable esimates
# of sy now. Some storms still show no actual response (or a muted responsed) 
# and so Dwt goes down over the course of the day

# Low intensity storms have extremely low Sy, even at high water levels
sy_dat <-
  sy_dat[24*precip_intensity_cm_hr > i_cm]

# Remove April (dominated by snowmelt not rain)
sy_dat <- 
  sy_dat[month(sample_date) > 4]

# Artificially high sy values (may be from well bottoming out)] 
sy_dat <- 
  sy_dat[sy < 2]

# Artificially high sy values (may be from well bottoming out)] 
sy_dat <- 
  sy_dat[!(water_level_cm < -60 & sy > 0.3)]

# High water levels at two sites with no outflow
sy_dat <- 
  sy_dat[!(site %in% c("119", "156") & water_level_cm > 10)]

ggplot(sy_dat,
       aes(x = water_level_cm,
           y = sy)) +
  geom_point() + 
  geom_smooth(method = nlrob,
              method.args = list(start = list(b = 0.2, m = 1.25, c = -50),
                                 algorithm = "port",
                                 maxit = 100,
                                 control = nls.control(maxiter = 100),
                                 lower = list(b = 0, m = 0, c =-100000)),
              se = FALSE,
              formula = y ~ b + m ** (x - c),
              color = "blue") +
  geom_smooth(method = robustbase::glmrob,
              method.args = list(family = Gamma(inverse)),
              se = FALSE,
              formula = y ~ x + I(x^2),
            linetype = "dashed",
            color = "red") +
  facet_wrap(~site, scales = "free")


# Large storms result in Hortonian overland flow, which likely pulses streamflow,
# resulting in lower delta WL. This means Sy is lower than would be predicted for
# the same amount of precip at a lower intensity. This could potentially be 
# remedied by looking at storm events

sy_mods <- 
  sy_dat[, 
         .(mod = list(glmrob(sy ~ water_level_cm,
                              family = Gamma(inverse))),
           mod_quad = list(glmrob(sy ~ water_level_cm + I(water_level_cm^2),
                                  family = Gamma(inverse))),
           mod_nl = list(nlrob(sy ~ b + m ** (water_level_cm - c),
                               data = .SD,
                               start = list(b = 0.2, m = 1.25, c = -50),
                               algorithm = "port",
                               maxit = 50,
                               lower = list(b = min_na(.SD$sy), m = 0, c =-100000),
                               control = nls.control(warnOnly = TRUE, maxiter = 100)))),
         by = .(site)]


sy_dat[site %in% sy_mods$site, 
       `:=`(pred_sy = predict(sy_mods[.BY[[1]], mod[[1]], on = "site"], 
                          newdata = .SD,
                          type = "response"),
            quad_sy = predict(sy_mods[.BY[[1]], mod_quad[[1]], on = "site"], 
                              newdata = .SD,
                              type = "response")
            , nl_sy = predict(sy_mods[.BY[[1]], mod_nl[[1]], on = "site"],
                          newdata = .SD,
                          type = "response")), 
       by = .(site)]

ggplot(sy_dat,
       aes(x = water_level_cm,
           y = sy)) +
  geom_point() + 
  geom_line(aes(y = pred_sy),
            color = "blue") +
  geom_line(aes(y = quad_sy),
            linetype = "dashed",
            color = "blue") +
  geom_line(aes(y = nl_sy),
            linetype = "dotted",
            color = "blue") +
  facet_wrap(~site, scales = "free")

ggplot(CJ(site = unique(sy_mods$site),
          water_level_cm = seq(-150, 50, by = 1)
          )[, 
            `:=`(uni_sy = ifelse(pmin(1, predict(sy_mods[.BY[[1]], 
                                                         mod[[1]], 
                                                         on = "site"], 
                                                 newdata = .SD, 
                                                 type = "response")) > 0,
                                 pmin(1, predict(sy_mods[.BY[[1]], 
                                                         mod[[1]], 
                                                         on = "site"], 
                                                 newdata = .SD, 
                                                 type = "response")),
                                 1),
                 quad_sy = ifelse(pmin(1, predict(sy_mods[.BY[[1]], 
                                                          mod_quad[[1]], 
                                                          on = "site"], 
                                                  newdata = .SD, 
                                                  type = "response")) > 0,
                                  pmin(1, predict(sy_mods[.BY[[1]], 
                                                          mod_quad[[1]], 
                                                          on = "site"], 
                                                  newdata = .SD, 
                                                  type = "response")),
                                  1),
                 nl_sy = pmin(1, predict(sy_mods[.BY[[1]],
                                                 mod_nl[[1]],
                                                 on = "site"],
                                         newdata = .SD))), 
            by = .(site)],
       aes(x = water_level_cm)) +
  geom_line(aes(y = nl_sy),
            linetype = "dotted",
            color = "blue") +
  geom_line(aes(y = quad_sy),
            linetype = "dashed",
            color = "blue") +
  geom_line(aes(y = uni_sy),
            linetype = "solid",
            color = "blue") +
  facet_wrap(~site, scales = "free")

# It may be worth using the non-linear models which have an asymptote at low
# water levels. The quadratic form could result in extremely low estimates of
# sy when predicting outside the water level range of sy dat

# Daily water balance from  -------------------------------
# Approach from McLaughlin, 2019 doesn't quite work. Nor did McLaughlin, 2014

daily_water_balance[daily_ex_met,
                    `:=`(exrad_pet_hs_cm = i.exrad_pet_hs_cm,
                         harad_pet_hs_cm = i.harad_pet_hs_cm,
                         obsrad_pet_hs_cm = i.obsrad_pet_hs_cm,
                         etr_cm_d = i.etr_cm_d,
                         observed_rad_MJ_m2 = i.observed_rad_MJ_m2),
                    on = c("sample_date")]

daily_water_balance[site %in% sy_mods$site, 
                    sy := predict(sy_mods[.BY[[1]], mod[[1]], on = "site"], 
                                  newdata = .SD[, .(water_level_cm = i_water_level_cm)],
                                  type = "response"), 
                    by = .(site)]

# Deal with exponential increase at asymptotes by setting everything after the
# asymptote to 1
daily_water_balance[, sy := pmin(1, sy)]
daily_water_balance[sy < 0, sy := 1]


daily_water_balance[, Ds_actual := sy * Ds_cm]

daily_water_balance[, 
                    net_flow := Ds_actual - observed_precip_cm + obsrad_pet_hs_cm - net_precip_cm,
                    by = .(site, sample_date)]

ggplot(daily_water_balance,
       aes(x = i_water_level_cm,
           y = net_flow)) +
  geom_point() +
  facet_wrap(~site, 
             scales = "free")

daily_water_balance[, sample_year := year(sample_date)]

# Did it rain?
daily_water_balance[, dry_status := as.numeric(!(net_precip_cm > 0))]

# If status is na for some reason (no precip records) then consider it a dry day
# This may be a mistake where rainfall happens in what is expected to be a 
# recession
daily_water_balance[is.na(dry_status), dry_status := 1]

# Create a unique ID for each wet/dry period
# Keep the numeric rleid to be able to combined periods more easily
daily_water_balance[, 
                    precip_period_count := rleid(dry_status),
                    by = .(site, sample_year)]

daily_water_balance[, 
            response_id := paste(site, sample_year, precip_period_count, sep = "-")]

# Combine dry period and preceeding wet period
daily_water_balance[, storm_id := response_id]

daily_water_balance[dry_status == 0,
                    storm_id := paste(site, sample_year, precip_period_count + 1, sep = "-"),
                    by = .(site, sample_year)]

# Count the number of dry days in each dry period
daily_water_balance[, 
            dry_days := cumsum(dry_status),
            by = .(response_id)]

# Add negative dry days for multiday storms where 0 is the last day of the storm
daily_water_balance[dry_status == 0, 
                    dry_days := (-(.N - 1)):0,
                    by = .(storm_id)]

# Get peak water level from the storm period
# There could be outliers produced if there is a large storm one day, followed
# by a small storm the next day. Most of the recession could occur during a
# 'non-dry' day. A multiday storm could also have a large second day event that
# is not quite a new max, so I can't just say that the day with the max should
# be considered day 0. Going to compromise with max water level from the last
# precip day in a storm period. That's what .SD[dry_status == 0, last(max_water_level_cm)]
# does

daily_water_balance[, storm_peak := .SD[dry_days == 0, max_water_level_cm],
            by = .(storm_id)]

daily_water_balance[, storm_rise := .SD[dry_days == 0, sy * Dl_cm],
                    by = .(storm_id)]

# Change in daily water level from the storm peak. 
daily_water_balance[,
                    delta_Ds_cm := f_water_level_cm - storm_peak]

daily_water_balance[, delta_Ds_actual := sy * delta_Ds_cm]

# Normalized Ds within a storm period
daily_water_balance[,
                    `:=`(recession_cm_day = delta_Ds_cm / dry_days,
                         recession_actual_day = delta_Ds_actual / dry_days)]

# Add up precip amount
daily_water_balance[, 
            storm_precip_cm := sum_na(net_precip_cm),
            by = .(storm_id)]

# Predict day 1 draw down rate from precip size, then fit a log decay to the
# median daily drawdown rate for a season. Could consider using a true
# intervention approach from time-series analysis. I can predict delta_Ds_cm
# when day == 1, using delta_Ds_cm * sy ~ log(storm_precip_cm), or delta_Ds_cm ~
# sy {family = gaussian(inverse)} Dormant season recession_cm_day never quite
# approaches the median or mean daily Ds. It probably should go to recession
# constant (delta_Ds_cm ~ dry_days), which can be approximated as the mean
# recession_cm_day on day 7

ggplot(daily_water_balance[between(dry_days, 1, 8) & site_status == "Control" & storm_precip_cm > i_cm],
       aes(x = dry_days,
           y = recession_actual_day,
           color = season)) +
  geom_point() + 
  geom_line(aes(group = storm_id),
            alpha = 0.3) +
  facet_wrap(~site, scales = "free")

# It's a geometric reduction
ggplot(daily_water_balance[between(dry_days, 1, 8) & site_status == "Control" & storm_precip_cm > i_cm],
       aes(x = dry_days,
           y = recession_actual_day)) +
  geom_point(aes(color = season)) + 
  geom_line(aes(group = storm_id,
                color = season),
            alpha = 0.3) +
  geom_function(fun = ~-2*(1.5**(-(.x-1)))) +
  facet_wrap(~site, scales = "free")

# Looks to perform better:
ggplot(daily_water_balance[dry_days == 1 & site_status == "Control" & storm_precip_cm > i_cm],
       aes(x = storm_precip_cm,
           y = delta_Ds_actual)) +
  geom_point() +
  geom_smooth(method = lmrob, 
              formula = y ~ log(x)) +
  facet_wrap(~site, scales = "free")


# Water Budget ------------------------------------------------------------

# Inflows | Outflows
# --------|---------
# Precip  | Interception
# GW      | Precip Recession
# Melt    | ET
#         | Q (SW & GW)

# Interception Model

# Melt Model

# Working in Water Level Space
# = fifelse(month(sample_date) %in% 6:8 |
#             (month(sample_date) == 5 & mday(sample_date) > 15) |
#             (month(sample_date) == 9 & mday(sample_date) <= 15),
#           "growing",
#           "dormant")
water_budget <- 
  daily_water_balance[between(month(sample_date), 6, 10), 
                      .(sample_date, 
                        season,
                        site_status,
                        dry_days,
                        sy, 
                        gross_precip_cm = gross_precip_cm / sy,
                        net_precip_cm = net_precip_cm / sy,
                        best_precip_cm = best_precip_cm / sy,
                        water_level_cm = i_water_level_cm, 
                        Dl_signed_cm,
                        Ds_cm = Ds_cm,
                        Ds_less_pet = Ds_cm + harad_pet_hs_cm / sy,
                        site_flow = Ds_cm + harad_pet_hs_cm / sy - net_precip_cm / sy,
                        precip_loss = Ds_cm - fifelse(net_precip_cm > 0, Dl_signed_cm, 0),
                        storm_recession_cm = shift(delta_Ds_cm, -1),
                        obsrad_pet_hs_cm = obsrad_pet_hs_cm / sy,
                        harad_pet_hs_cm = harad_pet_hs_cm / sy),
                      by = .(site, sample_year)]

water_budget[net_precip_cm == 0, storm_recession_cm := NA_real_]


# Net-flow Model ----------------------------------------------------------
# Need to perform this model and apply it before doing the PET model, then I
# should just be able to do one PET model that more accurately captures only PET
# ggplot(data = water_budget[!is.na(site_flow),
#                            .(site_flow = median_na(site_flow), nobs = .N),
#                            by = .(site, water_level_cm = round(water_level_cm, 0))],
#        aes(x = water_level_cm,
#            y = site_flow)) +
#   geom_point(alpha = 0.4) +
#   geom_smooth(method = lmrob, method.args = list(setting = "KS2014")) +
#   facet_wrap(~ site,
#              scales = "free")
# 
# # There are arguments for using dormant, growing season, or separate models. I 
# # think I am going to use dormant. You could argue that upland contributions 
# # could be reduced for a given water level during the growing season relative to
# # the dormant season. But I think that the hydraulics are likely the control in
# # this case and higher water levels means increased flow into the site, even if
# # upland transpiration is still in effect.
# 
# # flow_mods <- 
# #   water_budget[,
# #                .(mod= list(lmrob(site_flow ~ water_level_cm,
# #                                  data = .SD,
# #                                  # weights = abs(1 / (mean_na(water_level_cm) - water_level_cm)),
# #                                  setting = "KS2014"))),
# #                keyby = .(site)]
# # 
# # flow_mods[, c("intercept", "slope") := map_dfr(mod, coef)]
# # flow_mods[, slope_p := map(mod, ~broom::tidy(.x)$p.value[2])]
# # flow_mods[, intercept_p := map(mod, ~broom::tidy(.x)$p.value[1])]
# # flow_mods[, x_intercept := -intercept / slope]
# # 
# # water_budget[flow_mods,
# #              net_flow_cm := i.intercept + i.slope * water_level_cm, 
# #              on = "site"]
# 
# flow_mods <-
  # water_budget[!is.na(site_flow),
  #              .(site_flow = median_na(site_flow), nobs = .N),
  #              by = .(site, water_level_cm = round(water_level_cm, 0))] %>%
#   .[,
#     .(mod = list(nls(site_flow ~ b - (p**(water_level_cm)),
#                      data = .SD,
#                      start = list(b = 1, p = 1.06),
#                      algorithm = "port",
#                      # maxit = 100,
#                      lower = list(b = -20, p = 1),
#                      upper = list(b = 20, p = 2),
#                      weights = 1/abs(mean_na(pmax(-max_na(water_level_cm), water_level_cm)) - pmax(-max_na(water_level_cm), water_level_cm))^2))),
#     by = .(site)]
# 
# flow_mods[, c("beta0", "pow_base") := map_dfr(mod, coef)]
# flow_mods[, slope_p := map(mod, ~broom::tidy(.x)$p.value[2])]
# flow_mods[, intercept_p := map(mod, ~broom::tidy(.x)$p.value[1])]
# flow_mods[, x_intercept := logb(beta0, pow_base)]
# 
# water_budget[flow_mods,
#              net_flow_cm := i.beta0 - (i.pow_base^water_level_cm),
#              on = "site"]
# 
# water_budget[, Ds_cm := raw_Ds_cm - net_flow_cm]
# water_budget[, Dl_signed_cm := Dl_signed_cm - net_flow_cm]
# water_budget[flow_mods,
#              threshold := i.x_intercept,
#              on = c("site")]
# 
# ggplot(data = water_budget[dry_days > 1 & net_precip_cm == 0],
#        aes(x = water_level_cm,
#            y = site_flow)) +
#   geom_point() +
#   geom_line(data = water_budget[!is.na(net_flow_cm),
#                                 .(net_flow_cm = median_na(net_flow_cm), nobs = .N),
#                                 by = .(site, water_level_cm = round(water_level_cm, 0))],
#             alpha = 0.4,
#             aes(y = net_flow_cm)) +
#   facet_wrap(~ site,
#              scales = "free")
# 
# # Check Residuals
# ggplot(data = water_budget[dry_days > 1 & net_precip_cm == 0],
#        aes(x = water_level_cm,
#            y = site_flow - net_flow_cm)) + 
#   geom_vline(aes(xintercept = threshold)) +
#   geom_point(alpha = 0.4) +
#   facet_wrap(~ site, 
#              scales = "free")



# Ds Models

water_budget[, doy_decimal := yday(sample_date) / 366]
water_budget[daily_ex_met, 
             gdd := (i.max_temp_c + i.min_temp_c) / 2 - 10,
             on = "sample_date"]
water_budget[, cum_gdd := cumsum(gdd),
             by = .(site, sample_year)]
water_budget[, `:=`(cum_pet_cm = cumsum(ifelse(!is.na(harad_pet_hs_cm), harad_pet_hs_cm, 0)),
                    cum_precip_cm = cumsum(ifelse(!is.na(net_precip_cm), net_precip_cm, 0))),
             by = .(site, sample_year)]
water_budget[, `:=`(water_availability = (cum_precip_cm - cum_pet_cm) / 100,
                    p_pet = cum_precip_cm / cum_pet_cm)]


wb_mods <- 
  water_budget[, .(dat = list(.SD)), keyby = .(site_status)]


# Precip rise mod ---------------------------------------------------------
# For future model this could be improved by using quantile regression to get at
# the change variance of response to storm size
ggplot(water_budget[net_precip_cm > 0 & Dl_signed_cm > 0],
       aes(x = net_precip_cm,
           y = Dl_signed_cm,
           color = site)) +
  geom_point() +
  geom_smooth(method = lmrob, 
              formula = y ~ x,
              method.args = list(setting = "KS2014"),
              color = "blue") +
  facet_wrap(~site_status,
             scales = "free")

# wb_mods[, 
#         precip_rise := map(dat,
#                            ~lmrob(Dl_signed_cm ~ net_precip_cm,
#                                   data = .x[net_precip_cm > 0 & Dl_signed_cm > 0],
#                                   setting = "KS2014"))]

wb_mods[,
        precip_rise := map(dat,
                           ~lmer(Dl_signed_cm ~ net_precip_cm + (net_precip_cm | site),
                                  data = .x[net_precip_cm > 0 & Dl_signed_cm > 0]))]

# Residual Check
ggplot(data = water_budget[net_precip_cm > 0 & Dl_signed_cm > 0, 
                           .(net_precip_cm, 
                             std_resid = as.numeric(scale(predict(wb_mods[.BY[[1]], precip_rise[[1]]], newdata = .SD) - Ds_cm))), 
                           by = .(site_status)],
       aes(x = net_precip_cm,
           y = std_resid)) + 
  geom_point(alpha = 0.4) +
  facet_wrap(~site_status, 
             scales = "free")


# PET Model ---------------------------------------------------------------
# There are periods where precip is NA, but there is precip. So by specifying
# net_precip_cm == 0 it removes those points
# Could improve precip model by looking at Ds_cm/Dl_cm and identifying days with
# precip

# Single PET Model after Net Flow Adjustment
# More difference bewteen seasons than between high and low conditions
ggplot(data = water_budget[dry_days > 3 & net_precip_cm == 0 & site_status == "Control"],
       aes(x = water_availability,
           y = Ds_cm)) +
  geom_point(alpha = 0.4,
             aes(color = harad_pet_hs_cm)) +
  geom_smooth(method = lmrob,
              formula = y ~ x + I(x^2),
              se = FALSE,
              method.args = list(setting = "KS2014")) +
  facet_wrap(~ site,
             scales = "free")

ggplot(data = water_budget[dry_days > 3 & net_precip_cm == 0 & site_status == "Control"],
       aes(x = harad_pet_hs_cm,
           y = Ds_cm)) +
  geom_point(alpha = 0.4,
             aes(color = water_availability)) +
  geom_smooth(method = lmrob,
              formula = y ~ I(x^0.5),
              se = FALSE,
              method.args = list(setting = "KS2014")) +
  facet_wrap(~ site,
             scales = "free")

# Need to fine-tune this model
wb_mods[, 
        outflow_pet_single := map(dat,
                                  ~lmer(Ds_cm ~ I(harad_pet_hs_cm^0.5) + I(harad_pet_hs_cm^0.5):cos(pi * water_availability) + (I(harad_pet_hs_cm^0.5) + I(harad_pet_hs_cm^0.5):cos(pi * water_availability) || site),
                                         data = .x[dry_days > 3 & net_precip_cm == 0]))]

# Check Fit
ggplot(data = water_budget[dry_days > 3 & net_precip_cm == 0,
                           .(harad_pet_hs_cm,
                             water_availability,
                             sample_year,
                             site,
                             Ds_cm,
                             pred_Ds = predict(wb_mods[.BY[[1]], outflow_pet_single[[1]]], newdata = .SD)),
                           by = .(site_status)],
       aes(x = water_availability,
           y = Ds_cm)) +
  geom_point(alpha = 0.4) +
  geom_smooth(aes(y = pred_Ds)) +
  facet_wrap(~site,
             scales = "free")

ggplot(data = water_budget[dry_days > 3 & net_precip_cm == 0,
                           .(harad_pet_hs_cm,
                             water_availability,
                             sample_year,
                             site,
                             Ds_cm,
                             pred_Ds = predict(wb_mods[.BY[[1]], outflow_pet_single[[1]]], newdata = .SD)),
                           by = .(site_status)],
       aes(x = harad_pet_hs_cm,
           y = Ds_cm)) +
  geom_point(alpha = 0.4) +
  geom_smooth(aes(y = pred_Ds)) +
  facet_wrap(~site,
             scales = "free")

# Residual Check
ggplot(data = water_budget[dry_days > 3 & net_precip_cm == 0,
                           .(obsrad_pet_hs_cm,
                             site,
                             std_resid = as.numeric(Ds_cm - predict(wb_mods[.BY[[1]], outflow_pet_single[[1]]], newdata = .SD))),
                           by = .(site_status)],
       aes(x = obsrad_pet_hs_cm,
           y = std_resid)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm") +
  facet_wrap(~site,
             scales = "free")

ggplot(data = water_budget[dry_days > 3 & net_precip_cm == 0,
                           .(obsrad_pet_hs_cm,
                             site,
                             sample_date,
                             std_resid = as.numeric(scale(Ds_cm - predict(wb_mods[.BY[[1]], outflow_pet_single[[1]]], newdata = .SD)))),
                           by = .(site_status)],
       aes(x = sample_date - years(year(sample_date)),
           y = std_resid,
           color = as.factor(year(sample_date)))) +
  geom_point(alpha = 0.4) +
  facet_wrap(~site,
             scales = "free")



# Apply Models ------------------------------------------------------------

water_budget[, 
             `:=`(pred_pet_single = predict(wb_mods[CJ(.BY[[1]]), outflow_pet_single[[1]]],
                                            re.form = ~0,
                                             newdata = .SD),
                  pred_precip_rise = predict(wb_mods[CJ(.BY[[1]]), precip_rise[[1]]],
                                             re.form = ~0,
                                             newdata = .SD)),
             by = .(site_status)]

# Could add Markov component to Ds, it should rise if it was rising previously

water_budget[, c("pred_in", "pred_out", "pred_Ds_cm", "precip_recession_cm", "streamflow_cm") := list(0, 0, 0, 0, 0)]

# Apply Constraints to Predictions
water_budget[!(net_precip_cm > 0),
             pred_precip_rise := 0]

# # Add precip recession (estimated & probably varies by season)
# water_budget[,
#              precip_recession_cm := -reduce(shift(pred_precip_rise, 1:2, fill = 0), `+`) / 3,
#              by = .(site, sample_year)]
# water_budget[water_level_cm < 0,
#              precip_recession_cm := 0]

# water_budget[,
#              precip_recession_cm := -shift(pred_precip_rise, 1, fill = 0) / 3]
# 
# water_budget[water_level_cm < 0,
#              precip_recession_cm := 0]

# # Add streamflow
# water_budget[water_level_cm >= threshold,
#              streamflow_cm := predict(lmrob(site_flow ~ water_level_cm,
#                                          data = .SD,
#                                          setting = "KS2014"),
#                                    newdata = .SD),
#              by = .(site)]

# Sum up in and out components
water_budget[,
             pred_in := pred_precip_rise]

water_budget[, 
             pred_out := pred_pet_single]

# Get Predicted Ds
water_budget[, pred_Ds_cm := pred_in + pred_out]
water_budget[, residual_flow := pred_Ds_cm - Ds_cm]

# Adjust for net flow

ggplot(data = water_budget[,
                           .(residual_flow = median_na(residual_flow), nobs = .N),
                           by = .(site, water_level_cm = round(water_level_cm, 0))],
       aes(x = water_level_cm,
           y = residual_flow)) + 
  geom_point(alpha = 0.4) +
  geom_smooth(method = nls,
              method.args = list(start = list(m = 1.25, c = -20),
                                 control = nls.control(maxiter = 100,
                                                       warnOnly = TRUE)),
              se = FALSE,
              formula = y ~ m ** (x - c),
              color = "blue") +
  facet_wrap(~ site, 
             scales = "free")


# Use this if not using flow_mods above
net_flow_mods <-
  water_budget[,
               .(residual_flow = median_na(residual_flow), nobs = .N),
               by = .(site, water_level_cm = round(water_level_cm, 0))]

net_flow_mods <-
  water_budget[, 
                .(mod = list(nls(residual_flow ~ m ** (water_level_cm + c),
                             start = list(m = 1.25, c = 20),
                             control = nls.control(maxiter = 100,
                                                   warnOnly = TRUE)))),
                keyby = .(site)]

net_flow_mods[, c("m", "c") := map_dfr(mod, coef)]

net_flow_mods[, treatment := fcase(site %in% c("009", "053", "077", "139", "156"), "Ash Cut",
                                   site %in% c("006", "119", "140", "151"), "Girdle",
                                   default = "Control")]

net_flow_mods[, morphology := fcase(site %in% c("006", "113", "135", "140", "152"), "flowthrough",
                                    site %in% c("009", "053", "077", "111", "139", "151", "157"), "diffuse",
                                    site %in% c("119", "156"), "closed")]

water_budget[, pred_flow := predict(net_flow_mods[.BY[[1]], mod[[1]]],
                                            newdata = .SD),
             by = .(site)]

ggplot(data = water_budget,
       aes(x = water_level_cm,
           y = residual_flow)) + 
  geom_point(alpha = 0.4) +
  geom_line(aes(y = pred_flow),
            color = "blue") +
  facet_wrap(~ site, 
             scales = "free")

water_budget[, Ds_resid := (pred_in + pred_out - pred_flow) - Ds_cm]


# Need to improve dormant drawdown (probably has to do with streamflow)
split(water_budget[!is.na(net_precip_cm)],
      by = "site_status") %>% 
  map(~wrap_elements(ggplot(.x,
                            aes(x = Ds_cm,
                                y = pred_Ds_cm)) + 
                       geom_point(aes(color = season)) +
                       geom_abline() +
                       geom_smooth(method = lmrob,
                                   formula = y ~ x,
                                   method.args = list(setting = "KS2014")) +
                       facet_zoom(xy = between(Ds_cm, -2, 2),
                                  horizontal = FALSE))) %>% 
  reduce(`+`)

ggplot(water_budget[!is.na(net_precip_cm)],
       aes(x = Ds_cm,
           y = pred_Ds_cm,
           color = site_status)) + 
  geom_point() +
  geom_abline() +
  geom_smooth(method = lmrob,
              formula = y ~ x,
              method.args = list(setting = "KS2014")) +
  facet_wrap(~site, 
             scales = "free")

split(water_budget[!is.na(net_precip_cm)],
      by = "site_status") %>% 
  map(~wrap_elements(ggplot(.x,
              aes(x = raw_Ds_cm,
                  y = Ds_resid)) + 
        geom_point() +
        facet_zoom(xy = between(raw_Ds_cm, -2, 2),
                   horizontal = FALSE))) %>% 
  reduce(`+`)

# Metrics

water_budget[, 
             .(rmse = caret::RMSE((pred_in + pred_out - pred_flow), Ds_cm, na.rm = TRUE),
               r2 = caret::R2((pred_in + pred_out - pred_flow), Ds_cm, na.rm = TRUE),
               mad = median_na((pred_in + pred_out - pred_flow) / Ds_cm)),
             by = .(site_status, season)]

water_budget[between(Ds_cm, -2, 2), 
             .(rmse = caret::RMSE((pred_in + pred_out - pred_flow), Ds_cm, na.rm = TRUE),
               r2 = caret::R2((pred_in + pred_out - pred_flow), Ds_cm, na.rm = TRUE),
               mad = median_na((pred_in + pred_out - pred_flow) / Ds_cm)),
             by = .(site_status, season)]
# Error Sources

ggplot(water_budget[!is.na(net_precip_cm)],
       aes(x = obsrad_pet_hs_cm,
           y = Ds_resid,
           color = season)) + 
  geom_point() +
  geom_smooth(method = lm) +
  facet_wrap(~site,
             scales = "free")

ggplot(water_budget[!is.na(net_precip_cm)],
       aes(x = net_flow_cm,
           y = Ds_resid)) + 
  geom_point() +
  geom_smooth(method = lm) +
  facet_wrap(~site,
             scales = "free")

ggplot(water_budget[!is.na(net_precip_cm)],
       aes(x = water_level_cm,
           y = Ds_resid)) + 
  geom_point() +
  geom_smooth(method = lm) +
  facet_wrap(~site,
             scales = "free")

ggplot(water_budget[!is.na(net_precip_cm) & precip_recession_cm > 0],
       aes(x = precip_recession_cm,
           y = Ds_resid)) + 
  geom_point() +
  geom_smooth(method = lm) +
  facet_wrap(~site,
             scales = "free")

ggplot(water_budget[net_precip_cm > 0],
       aes(x = net_precip_cm,
           y = Ds_resid)) + 
  geom_point() +
  geom_smooth(method = lm) +
  facet_wrap(~site,
             scales = "free")


# Apply Predictions

predictions <- 
  water_budget[water_budget[!is.na(water_level_cm), 
                            .(start_date = pmax(min_na(sample_date), 
                                                ymd(paste0(.BY[[2]], "0615"))),
                              end_date = pmin(max_na(sample_date),
                                              ymd(paste0(.BY[[2]], "0930")))), 
                            by = .(site, sample_year)], on = c("site", "sample_year")][between(sample_date, start_date, end_date)]

predictions[net_flow_mods,
            `:=`(m = i.m,
                 c = i.c),
            on = "site"]

predict_wl <-
  function(init, x, m, c){
    m <- unique(m)
    c <- unique(c)
     
    stopifnot(length(m) == 1)
    stopifnot(length(c) == 1)
    
    wl <- 
      numeric(length(x))
    
    wl[1] <- 
      init
    
    for(i in 1:(length(x)-1)){
      flow <- m ** (wl[i] + c)
      
      flow <- 
        ifelse(wl[i] > (-c),
               pmin(flow, wl[i] + c),
               flow)
      
      wl[i+1] <- wl[i] + x[i] - flow
    }
    
    wl
  }

predictions[, 
            pred_wl := predict_wl(first(water_level_cm), pred_Ds_cm, m, c),
            by = .(site, sample_year)]

ggplot(predictions[sample_year == 2014],
       aes(x = sample_date,
           color = site_status)) +
  geom_line(aes(y = water_level_cm)) +
  geom_line(aes(y = pred_wl),
            linetype = "dotted") +
  facet_wrap(~site, 
             scales = "free") + 
  

ggplot(predictions,
       aes(x = yday(sample_date),
           color = factor(sample_year))) +
  geom_line(aes(y = pred_wl - water_level_cm)) +
  facet_wrap(~site, 
             scales = "free")

ggplot(predictions,
       aes(x = yday(sample_date),
           color = site)) +
  geom_line(aes(y = pred_wl - water_level_cm)) +
  facet_wrap(~factor(sample_year), 
             scales = "free")

# Coefficients

precip_coefs <- 
  rbind(wb_mods[, map_dfr(precip_rise, 
                          ~data.table(site = rownames(ranef(.x)$site), 
                                      fixef(.x) + ranef(.x)$site)), 
                by = .(site_status)],
        wb_mods[, map_dfr(precip_rise, 
                          ~data.table(site = "base", 
                                      fixef(.x)[1],
                                      fixef(.x)[2])), 
                by = .(site_status)],
        use.names = FALSE)

setnames(precip_coefs, c("site_status", "site", "intercept", "slope"))

pet_coefs <- 
  rbind(wb_mods[, map_dfr(outflow_pet_single, 
                          ~data.table(site = rownames(ranef(.x)$site), 
                                      fixef(.x) + ranef(.x)$site)), 
                by = .(site_status)],
        wb_mods[, map_dfr(outflow_pet_single, 
                          ~as.data.table(c(list(site = "base"), 
                                           as.list(fixef(.x))))), 
                by = .(site_status)],
        use.names = FALSE)

pet_coefs <- 
  flow_mods[, .(site, beta0, pow_base)]

setnames(pet_coefs, c("site_status", "site", "intercept", "slope_pet", "slope_water", "slope_inter"))

# Just pool them for now. May be less standard deviation with ash cut, but have
# to test that to make sure it's not just 156 skewing the results
starting_values <- 
  water_budget[format(sample_date, "%m%d") == "0601" & !is.na(water_level_cm),
               .(site, site_status, sample_year, water_level_cm)]

# Using predict() rather than the coefficients means I can't as easily represent
# non-observed combinations




# Save Everything ---------------------------------------------------------

saveRDS(net_flow_mods, "tmp/models/net_flow_mods.rds")
saveRDS(interception, "tmp/models/interception.rds")
saveRDS(sy_mods, "tmp/models/sy_mods.rds")
saveRDS(starting_values, "tmp/models/starting_values.rds")
saveRDS(wb_mods, "tmp/models/wb_mods.rds")
saveRDS(water_budget[, unique(.SD[, .(site, site_status)])], "tmp/models/site_status_combinations.rds")


full_predict <-
  function(site_status, 
           hydrology_model = NULL, 
           initial.wl = NA,
           weather, 
           components = FALSE,
           random.effects = NULL){
    
    # weather: data.frame with doy, gross_precip_cm, min_temp_c, max_temp_c, 
    # mean_temp_c, water_availability
    
    stopifnot(is.logical(components))
    
    site_status <- 
      unique(site_status)
    
    stopifnot(length(site_status) == 1)
    
    if(is.null(hydrology_model)){
      
      hydrology_model <- 
        sample(net_flow_mods[treatment == site_status, site], 1)
      
    } else {
      
      hydrology_model <- 
        unique(hydrology_model)
      
      stopifnot(length(hydrology_model) == 1)
      
    }
    
    if(is.na(initial.wl)){
      initial.wl <- 
        sample(starting_values[site == hydrology_model, water_level_cm], 1)
    }
    
    precip_rise_mod <- 
      wb_mods[site_status, 
              precip_rise[[1]]]
      
    pet_mod <- 
      wb_mods[site_status, 
              outflow_pet_single[[1]]]
    
    # Can't use predict for flow if using the power function
    flow_mod <- 
      net_flow_mods[hydrology_model]
    
    slope <- flow_mod$m
    offset <- flow_mod$c
    
    # Get ESy Mod
    esy_mod <- 
      sy_mods[hydrology_model, mod[[1]]]
    
    dat <- 
      copy(weather)
    
    if(any(is.na(dat$gross_precip_cm))){
      message("Some missing precip values were filled in as 0.")
      
      dat[, gross_precip_cm := nafill(gross_precip_cm,
                                      type = "const",
                                      fill = 0)]
    }
    
    dat[, `:=`(site = hydrology_model,
               site_status = site_status)]
    
    dat[, season := fifelse(between(doy, 135, 258), 
                            "growing",
                            "dormant")]
    
    dat[interception, 
        interception_cm := i.i_cm,
        on = c("site", "site_status", "season")]
    
    # Generated for Wakefield (could be done spatially with GLM)
    ha_coefs <- 
      c(Ha = 0.1587725,
        Hb = -1.7921571,
        Hr2 = 0.7908223)

    # Predict Radation & PET
    dat[, ET_rad_MJ_m2 := extrat(doy, radians(46.440278))$ExtraTerrestrialSolarRadiationDaily]
    dat[, ha_rad_MJ_m2 := ha(as.Date(doy, origin = "2018-12-31"),  # non-leap year date from doy (leap-year vs not does not change radiation estimate)
                             lat = 46.440278,
                             lon = -89.826667,
                             ET_rad_MJ_m2,
                             A = ha_coefs[["Ha"]],
                             B = ha_coefs[["Hb"]],
                             Tmax = max_temp_c,
                             Tmin = min_temp_c)]

    dat[, harad_pet_hs_cm := 0.1 * 0.0023 * ha_rad_MJ_m2 / lambda_MJ_kg * (mean_temp_c + 17.8) * sqrt(max_temp_c - min_temp_c)]

    # Could add Markov component to Ds, it should rise if it was rising previously
    
    
    
    wl <- pet <- esy <- precip_rise <- flow <- 
      numeric(nrow(dat))
    
    wl[1] <- 
      initial.wl
    
    i_water_availability <- 
      0
    
    for(i in 1:(nrow(dat)-1)){

      esy[i] <- 
        predict(esy_mod,
                type = "response",
                newdata = data.frame(water_level_cm = wl[i]))
      
      esy[i] <- 
        pmin(1, esy[i])
      
      # There probably are not any asymptotes in the univariate gamma model  &
      # definitely none in the NL model. Leaving this in case I try the quadratic
      # gamma model again
      
      esy[i] <- 
        ifelse(esy[i] < 0, 1, esy[i])
      
      i_pet <- 
        dat[i, harad_pet_hs_cm / esy[i]]
        
      i_precip <- 
        dat[i, pmax(0, gross_precip_cm - interception_cm) / esy[i]]
      
      i_water_availability <- 
        i_water_availability + ((i_precip - i_pet) / 100)
      
      pet[i] <- 
        predict(pet_mod,
                re.form = random.effects,
                newdata = data.frame(site = hydrology_model,
                                     harad_pet_hs_cm = i_pet,
                                     water_availability = i_water_availability))
      
      precip_rise[i] <- 
        predict(precip_rise_mod, 
                re.form = random.effects,
                newdata = data.frame(site = hydrology_model,
                                     net_precip_cm = i_precip))
      
      # Set rise to zero for days with no precip (model has an intercept)
      precip_rise[i] <- 
        fifelse(i_precip > 0, precip_rise[i], 0)
      
      flow[i] <- 
        slope ** (wl[i] + offset)
      
      flow[i] <- 
        ifelse(wl[i] > (-offset),
               pmin(flow[i], wl[i] + offset),
               flow[i])
      
      wl[i+1] <- wl[i] + precip_rise[i] + pet[i] - flow[i]
    }
    
    if(components){
      return(data.table(hydrology_model, site_status, doy = dat$doy, wl, interception = dat$interception, precip_rise, pet, flow, esy))
    }
    wl
  }


# Test Full Prediction ----------------------------------------------------

water_budget[daily_ex_met,
            `:=`(min_temp_c = i.min_temp_c,
                 max_temp_c = i.max_temp_c,
                 mean_temp_c = i.mean_temp_c),
            on = "sample_date"]
  
water_budget[, raw_precip := gross_precip_cm * sy]

water_budget[, 
             full_pred := full_predict(site_status = site_status, 
                                       hydrology_model = site,
                                       initial.wl = first(water_level_cm),
                                       weather = .SD[, .(doy = yday(sample_date),
                                                         gross_precip_cm = raw_precip,
                                                         min_temp_c,
                                                         max_temp_c,
                                                         mean_temp_c)]),
             by = .(site, sample_year)]

YEAR <- 2014
ggplot(water_budget[sample_year == YEAR],
       aes(x = sample_date,
           color = site_status)) +
  geom_line(aes(y = water_level_cm)) +
  geom_line(data = predictions[sample_year == YEAR], 
            aes(y = pred_wl),
            linetype = "dotted",
            color = "gray50") +
  geom_line(aes(y = full_pred),
            linetype = "dotted") +
  facet_wrap(~site, 
             scales = "free")

ggplot(water_budget,
       aes(x = yday(sample_date),
           color = factor(sample_year))) +
  geom_line(aes(y = full_pred - water_level_cm)) +
  facet_wrap(~site, 
             scales = "free")

ggplot(water_budget,
       aes(x = yday(sample_date),
           color = site)) +
  geom_line(aes(y = full_pred - water_level_cm)) +
  facet_wrap(~factor(sample_year), 
             scales = "free")


# Test Treatment Simulation

test_dat <- 
  water_budget[site == "140" & sample_year == 2014, 
              .(site,
                sample_date,
                sample_year,
                water_level_cm,
                control_pred = full_predict(site_status = "Control", 
                                        hydrology_model = site,
                                        initial.wl = first(water_level_cm),
                                        weather = .SD[, .(doy = yday(sample_date),
                                                          gross_precip_cm = raw_precip,
                                                          min_temp_c,
                                                          max_temp_c,
                                                          mean_temp_c)]),
              girdle_pred = full_predict(site_status = "Girdle", 
                                           hydrology_model = site,
                                           initial.wl = first(water_level_cm),
                                           weather = .SD[, .(doy = yday(sample_date),
                                                             gross_precip_cm = raw_precip,
                                                             min_temp_c,
                                                             max_temp_c,
                                                             mean_temp_c)]),
              cut_pred = full_predict(site_status = "Ash Cut", 
                                         hydrology_model = NULL,
                                         initial.wl = first(water_level_cm),
                                         weather = .SD[, .(doy = yday(sample_date),
                                                           gross_precip_cm = raw_precip,
                                                           min_temp_c,
                                                           max_temp_c,
                                                           mean_temp_c)]))]

ggplot(test_dat[sample_year == 2014],
       aes(x = sample_date)) +
  geom_line(aes(y = control_pred),
            linetype = "solid") +
  geom_line(aes(y = girdle_pred),
            linetype = "dashed") +
  geom_line(aes(y = cut_pred),
            linetype = "dotted") +
  geom_line(aes(y = water_level_cm),
            color = "gray40") +
  ggthemes::theme_few()


treatment_simulations <- 
  water_budget[site != "006", 
               .(sample_date,
                 water_level_cm,
                 control_pred = full_predict(site_status = "Control", 
                                             hydrology_model = site,
                                             initial.wl = first(water_level_cm),
                                             weather = .SD[, .(doy = yday(sample_date),
                                                               gross_precip_cm = raw_precip,
                                                               min_temp_c,
                                                               max_temp_c,
                                                               mean_temp_c)]),
                 treatment_pred = full_predict(site_status = site_status, 
                                            hydrology_model = site,
                                            initial.wl = first(water_level_cm),
                                            weather = .SD[, .(doy = yday(sample_date),
                                                              gross_precip_cm = raw_precip,
                                                              min_temp_c,
                                                              max_temp_c,
                                                              mean_temp_c)])),
               by = .(site, sample_year)]


ggplot(treatment_simulations[sample_year == 2015],
       aes(x = sample_date)) +
  geom_line(aes(y = control_pred),
            linetype = "dotted",
            color = "blue") +
  geom_line(aes(y = treatment_pred),
            linetype = "solid",
            color = "blue") + 
  geom_line(aes(y = water_level_cm),
            color = "gray40") +
  ggthemes::theme_few() +
  facet_wrap(~site,
             scales = "free")

# Simulations -------------------------------------------------------------

weather_sequences <- 
  predictions[predictions[!is.na(gross_precip_cm) & !is.na(obsrad_pet_hs_cm),
                          .(.N),
                          by = .(site, sample_year)][N == 108][, sequence_id := 1:.N],
              .(sequence_id, 
                doy = yday(sample_date), 
                net_precip_cm,
                obsrad_pet_hs_cm),
              nomatch = NULL,
              on = c("site", "sample_year")]




flow_mods[flow_mods[water_budget,
                    on = "site"][water_level_cm >= threshold, 
                                 .(streamflow_mod = list(lmrob(site_flow ~ water_level_cm,
                                                               data = .SD,
                                                               setting = "KS2014"))),
                                 by = .(site)],
          streamflow_mod := i.streamflow_mod,
          on = "site"]

hydrology_models <- 
  flow_mods[, .(hydrology_id = site,
                netflow_mod = mod,
                streamflow_mod,
                x_intercept)]



wb_mods

# Select Weather Sequence
# Select Starting Value
# Select Hydrology Model (threshold, net flow & streamflow)
# Start with initial water level and work through weather sequence
# More effective if done bayesian like

# Simulate water level

# Simulate Future ---------------------------------------------------------

future_sequences <- 
  fread("output/tabular/weather_sequences_gfdl-cm3_2070-2099.csv")

future_sequences[, mean_temp_c := (tmin + tmax) / 2]

future_sequences <- 
  future_sequences[between(month(sample_date), 6, 10)]

future_sequences <- 
  water_budget[, unique(.SD[, .(site, site_status)])][, 
               .(dat = list(future_sequences)), 
               by = .(site, site_status)
               ][, .(dat = list(future_sequences)), 
                 by = .(site, site_status)
                 ][, dat[[1]], 
                   by = .(site, site_status)]

set.seed(1234)
future_simulations <- 
  split(future_sequences, 
        by = c("site", "site_status", "seq_ID")) %>% 
  mclapply(function(x){
    x[,
      c("hydrology_model", "site_status_mod", "doy_mod", "wl", 
        "interception", "precip_rise", "pet", "flow", "esy") :=
        full_predict(hydrology_model = x$site,
                 site_status = x$site_status,
                 weather = .SD[, 
                             .(doy, 
                               gross_precip_cm = precip / 10,
                               min_temp_c = tmin,
                               max_temp_c = tmax,
                               mean_temp_c)],
                 components = TRUE,
                 random.effects = NULL)]
    },
    mc.cores = 4) %>% 
  rbindlist()



# # First Test of Sequences
# future_simulations <- 
#   future_sequences[,
#                    rbindlist(map(c("Control", "Ash Cut", "Girdle"),
#                                  ~full_predict(site_status = .x,
#                                                weather = .SD[, 
#                                                              .(doy, 
#                                                                gross_precip_cm = precip / 10,
#                                                                min_temp_c = tmin,
#                                                                max_temp_c = tmax,
#                                                                mean_temp_c)],
#                                                components = TRUE,
#                                                random.effects = NA))),
#                    by = .(seq_ID)]

fwrite(future_simulations, "tmp/test_future_simulations.csv")

RPushbullet::pbPost("note",
                    "Simulations Complete")

treatment_simulations <- 
  water_budget[site != "006", 
               .(sample_date,
                 water_level_cm,
                 control_pred = full_predict(site_status = "Control", 
                                             hydrology_model = site,
                                             initial.wl = first(water_level_cm),
                                             weather = .SD[, .(doy = yday(sample_date),
                                                               gross_precip_cm = raw_precip,
                                                               min_temp_c,
                                                               max_temp_c,
                                                               mean_temp_c)]),
                 treatment_pred = full_predict(site_status = site_status, 
                                               hydrology_model = site,
                                               initial.wl = first(water_level_cm),
                                               weather = .SD[, .(doy = yday(sample_date),
                                                                 gross_precip_cm = raw_precip,
                                                                 min_temp_c,
                                                                 max_temp_c,
                                                                 mean_temp_c)])),
               by = .(site, sample_year)]
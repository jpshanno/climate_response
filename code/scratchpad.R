source("code/climate_packages.R")


# Functions ---------------------------------------------------------------

safe_lm <- 
  possibly(lm, NA)

safe_coef <- 
  possibly(coef, c(`(Intercept)` = NA_real_, gross_precip_cm = NA_real_))

extract_changepoint <- 
  function(model, n = 1){
    cp <- 
      paste0("cp_", n)
    
    res <- mcp::fixef(model)
    res <- subset(res, name == cp, drop = TRUE)
    res[["mean"]]
  }

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

calculate_detrended_g <- 
  function(data){
    data <- copy(data)
    
    recharge_period <- 
      0:7
    
    # nested <-
    #   map_dfr(unique(data$sample_date)[-c(length(unique(data$sample_date)))],
    #            # c(list(c(0, 1, 2)), replicate(length(unique(data$sample_date))-2, -1:1, simplify = FALSE)),
    #            ~data.table(sample_date = .x,
    #                        dat = list(data[hour(sample_time) %in% recharge_period & sample_date %in% c(.x + -1:1)])))
    # 
    # # nested[sample_date %in% as.Date(paste0("2018-08-", 10:23)), dat] %>%
    # #   map(~ggplot(data = .x,
    # #               aes(x = sample_time, y = corrected_compensated_level_cm)) +
    # #         ggtitle(unique(.x$sample_date)) +
    # #         geom_point() +
    # #         geom_smooth(method = "lm", se = FALSE)) %>%
    # #   reduce(`+`)
    # 
    # nested[, `:=`(raw_detrend_mod = map(dat, function(x){lm(raw_compensated_level_cm ~ sample_time, data = x)}),
    #               corrected_detrend_mod = map(dat, function(x){lm(corrected_compensated_level_cm ~ sample_time, data = x)}))]
    # 
    # nested[, `:=`(raw_detrend_m_cm_s = map_dbl(raw_detrend_mod,
    #                                    ~coef(.x)[["sample_time"]]),
    #               raw_detrend_b_cm_s = map_dbl(raw_detrend_mod,
    #                                    ~coef(.x)[["(Intercept)"]]),
    #               corrected_detrend_m_cm_s = map_dbl(corrected_detrend_mod,
    #                                          ~coef(.x)[["sample_time"]]),
    #               corrected_detrend_b_cm_s = map_dbl(corrected_detrend_mod,
    #                                          ~coef(.x)[["(Intercept)"]]))]
    # 
    # data[nested,
    #      `:=`(raw_detrend_m_cm_s = i.raw_detrend_m_cm_s,
    #           raw_detrend_b_cm_s = i.raw_detrend_b_cm_s,
    #           corrected_detrend_m_cm_s = i.corrected_detrend_m_cm_s,
    #           corrected_detrend_b_cm_s = i.corrected_detrend_b_cm_s),
    #      on = "sample_date"]
    # 
    # data[,`:=`(raw_detrended_cm = raw_compensated_level_cm - (raw_detrend_b_cm_s + raw_detrend_m_cm_s * as.numeric(sample_time)),
    #            corrected_detrended_cm = corrected_compensated_level_cm - (corrected_detrend_b_cm_s + corrected_detrend_m_cm_s * as.numeric(sample_time)))]
    
    # daily detrend:
    nested <-
      data[hour(sample_time) %in% recharge_period, 
           .(corrected_i_cm = first(water_level_cm)),
           by = .(sample_date)]
    
    nested[, `:=`(corrected_f_cm = shift(corrected_i_cm, -1))]
    
    nested[, `:=`(corrected_detrend_m_cm_s = (corrected_f_cm - corrected_i_cm) / 86400)]
    
    nested[, `:=`(corrected_detrend_b_cm_s =corrected_i_cm)]
    
    data[nested,
         `:=`(corrected_detrend_m_cm_s = i.corrected_detrend_m_cm_s,
              corrected_detrend_b_cm_s = i.corrected_detrend_b_cm_s),
         on = "sample_date"]
    
    data[,`:=`(corrected_detrended_cm = water_level_cm - (corrected_detrend_b_cm_s + corrected_detrend_m_cm_s * as.numeric(sample_time - min(sample_time)))),
         by = .(sample_date)]
    
    data[!(hour(sample_time) %in% recharge_period),
         `:=`(corrected_detrend_m_cm_s = NA_real_,
              corrected_detrend_b_cm_s = NA_real_)]
    
    data[, `:=`(corrected_detrend_m_cm_s = zoo::na.spline(corrected_detrend_m_cm_s), 
                corrected_detrend_b_cm_s = zoo::na.spline(corrected_detrend_b_cm_s))]
    
    data[, `:=`(corrected_dwtDT_dt = three_point_slope(sample_time, corrected_detrended_cm))]
    
    # Gamma_dat <- 
    #   data[, .(sample_time,
    #            sample_date,
    #            raw_compensated_level_cm,
    #            corrected_compensated_level_cm,
    #            raw_detrended_cm,
    #            corrected_detrended_cm,
    #            raw_dwtDT_dt,
    #            corrected_dwtDT_dt)]
    # 
    # Gamma_dat <- 
    #   Gamma_dat[hour(sample_time) %in% recharge_period]
    
    # Need to create Gamma(wt) for each day separately. Then need to find a 
    # way to smoothly move between the two days to create a function. If it were
    # just slope it would be easy to use linear interpretation between the two 
    # recharge period models. But it isn't, need to think about how to transition
    # slope and intercept
    # The assumption of an approximately linear relationship for Gamma(wt) 
    # between recharge periods does not hold
    
    # Gamma_dat <- 
    #   map_dfr(unique(Gamma_dat$sample_date)[-c(length(unique(Gamma_dat$sample_date)))], 
    #           ~data.table(sample_date = .x, 
    #                       data = list(Gamma_dat[sample_date %in% c(.x)])))
    # 
    # Gamma_dat[sample_date %in% as.Date(paste0("2018-08-", 10:23)), data] %>%
    #   map(~ggplot(data = .x,
    #               aes(x = corrected_detrended_cm, y = corrected_dwtDT_dt)) +
    #         ggtitle(unique(.x$sample_date)) +
    #         geom_point(aes(color = as.factor(sample_date)),
    #                    show.legend = FALSE) +
    #         geom_smooth(method = "lm", se = FALSE, formula = "y ~ x")) %>%
    #   reduce(`+`)
    # 
    # Gamma_dat[sample_date %in% as.Date(paste0("2018-08-", 10:23)), data] %>%
    #   map(~ggplot(data = .x,
    #               aes(x = sample_time, y = three_point_slope(sample_time, corrected_dwtDT_dt))) +
    #         ggtitle(unique(.x$sample_date)) +
    #         geom_point(aes(color = as.factor(sample_date)),
    #                    show.legend = FALSE) +
    #         geom_smooth(method = "lm", se = FALSE, formula = "y ~ x")) %>%
    #   reduce(`+`)
    
    # Gamma_dat[sample_date %in% as.Date(paste0("2018-08-", 10:23)), data][1:13] %>% 
    #   rbindlist() %>% 
    #   ggplot(aes(x = corrected_detrended_cm, 
    #              y = corrected_dwtDT_dt, 
    #              color = as.factor(sample_date))) +
    #   geom_point() + scale_color_viridis_d() +
    #   geom_smooth(method = "lm", se = FALSE)
    
    # For daily data need to remove days with no valid observations
    # null_dates <- 
    #   Gamma_dat[, .(sample_date, 
    #                 nobs = map_dbl(data, nrow), 
    #                 na_dwtDT_dt = map_dbl(data, ~sum(is.na(.x$corrected_dwtDT_dt))), 
    #                 na_WTDT = map_dbl(data, ~sum(is.na(.x$corrected_detrended_cm))))]
    # 
    # null_dates <- 
    #   null_dates[nobs == na_dwtDT_dt | nobs == na_WTDT, sample_date]
    #     
    # Gamma_wt <- 
    #   Gamma_dat[!(sample_date %in% null_dates), 
    #             c(as.list(coef(lm(corrected_dwtDT_dt ~ corrected_detrended_cm, data = data[[1]]))),
    #               as.list(coef(lm(raw_dwtDT_dt ~ raw_detrended_cm, data = data[[1]])))), 
    #             by = .(sample_date)]
    # 
    # setnames(Gamma_wt, c("sample_date", "corrected_b", "corrected_m", "raw_b", "raw_m"))
    # 
    # data[Gamma_wt,
    #      `:=`(raw_gamma_m = i.raw_m,
    #           raw_gamma_b = i.raw_b,
    #           corrected_gamma_m = i.corrected_m,
    #           corrected_gamma_b = i.corrected_b),
    #      on = "sample_date"]
    
    # data[!(hour(sample_time) %in% recharge_period),
    #      `:=`(raw_gamma_m = NA_real_,
    #           raw_gamma_b = NA_real_,
    #           corrected_gamma_m = NA_real_,
    #           corrected_gamma_b = NA_real_)]
    # 
    # data[, `:=`(raw_gamma_m = zoo::na.spline(raw_gamma_m),
    #             raw_gamma_b = zoo::na.spline(raw_gamma_b), 
    #             corrected_gamma_m = zoo::na.spline(corrected_gamma_m), 
    #             corrected_gamma_b = zoo::na.spline(corrected_gamma_b))]
    # 
    # data[,`:=`(raw_Gamma_cm_s = raw_gamma_b + raw_gamma_m * raw_detrended_cm,
    #            corrected_Gamma_cm_s = corrected_gamma_b + corrected_gamma_m * corrected_detrended_cm)]
    
    data[, `:=`(corrected_gamma_m = corrected_dwtDT_dt + three_point_slope(sample_time, corrected_dwtDT_dt))]
    
    data[!(hour(sample_time) %in% recharge_period),
         `:=`(corrected_Gamma_cm_s = NA_real_)]
    
    data[!is.na(corrected_gamma_m) & hour(sample_time) %in% recharge_period,
         corrected_Gamma_cm_s := fitted(lm(corrected_gamma_m ~ sample_time)),
         by = .(sample_date)]
    
    # data[, `:=`(raw_gamma_m = zoo::na.approx(raw_gamma_m, rule = 2),
    #             corrected_gamma_m = zoo::na.approx(corrected_gamma_m, rule = 2))]
    
    # data[,`:=`(raw_Gamma_cm_s = raw_gamma_m * as.numeric(sample_time),
    #            corrected_Gamma_cm_s = corrected_gamma_m * as.numeric(sample_time))]
    
    data[, `:=`(corrected_Gamma_cm_s = zoo::na.spline(corrected_Gamma_cm_s))]
    
    data[, `:=`(corrected_net_in_cm_s = sy * (corrected_Gamma_cm_s + corrected_detrend_m_cm_s))]
    
    data[]
  }

calculate_delta_s <- 
  function(data){
    
    data <- copy(data)
    
    data[, `:=`(corrected_dwt_cm_s = water_level_cm - shift(water_level_cm))]
    
    # data[, `:=`(corrected_dwt_cm_s = rolling_slope(sample_time, water_level_cm, 3))]
    
    # data <- 
    #   smooth_data(data, "raw_dwt_cm_s", "corrected_dwt_cm_s", n = 13)
    
    data[] 
    # data[, `:=`(raw_dwt_cm_s = (shift(raw_compensated_level_cm, -1) - raw_compensated_level_cm) / 900,
    #             corrected_dwt_cm_s = (shift(corrected_compensated_level_cm, -1) - corrected_compensated_level_cm) / 900)]
  }

# The decline of water levels following rain is steeper than anticpated by the 
# basic detrending procedure. This means that for 1-2 days after rain there is
# a faster decline in water levels than would be anticipated. A simple/post-hoc
# solution may be to use a left-aligned window to detrend (which will cause 
# problems for days with preceding rainfall). A more mechanistic solution would
# be to develop a recession curve (which would likely vary with water level) and
# apply that in addition to the detrending slope. Or perhaps there is an event-
# based recession curve, but also a recession curve that can be applied to non-
# precipitation events, which is just how quickly the wetland drains, which 
# would certainly vary with water level.

# @watras-2017 was able to apply the @loheideii-2008 method because of the 
# small magnitude of seasonal water level fluctations (an the corresponding 
# small responses to preciptation)

# We have an outflow stream! Its really calculating net loss, not just ET. May
# need a PQ ratio for these wetlands to subtract out surface losses

calculate_et <- 
  function(data){
    data <- copy(data)
    data[, `:=`(corrected_et_cm_s = corrected_net_in_cm_s - sy * corrected_dwt_cm_s)]
    
    # data <- 
    #   smooth_data(data, "raw_et_cm_s", "corrected_et_cm_s", n = 13)
    
    data[]
    
  }

adjust_water_balance <- 
  function(data){
    data <- 
      copy(data)
    
    data[corrected_et_cm_s < 0,
         `:=`(corrected_net_in_cm_s = corrected_net_in_cm_s - corrected_et_cm_s,
              corrected_et_cm_s = corrected_et_cm_s - corrected_et_cm_s)]
    
    data[corrected_net_in_cm_s < 0,
         `:=`(corrected_et_cm_s = corrected_et_cm_s - corrected_net_in_cm_s,
              corrected_net_in_cm_s = corrected_net_in_cm_s - corrected_net_in_cm_s)]
    
    data[]
  }

three_point_slope <- 
  function(x, y){
    if(inherits(x, "POSIXt")){
      x <- as.numeric(x)
    }
    
    m1 <- 
      (shift(y, -1) - y) / (shift(x, -1) - x)
    
    m2 <- 
      (y - shift(y, 1)) / (x - shift(x, 1))
    
    (m1 + m2) / 2
  }

rolling_slope <- 
  function(x, y, n){
    if(inherits(x, "POSIXt")){
      x <- as.numeric(x)
    }
    
    if(n %% 2 == 0){
      stop("n must be odd")
    }
    
    n0 <- (n-1)/2
    
    m1 <- 
      (shift(y, -n0) - y) / (shift(x, -n0) - x)
    
    m2 <- 
      (y - shift(y, n0)) / (x - shift(x, n0))
    
    (m1 + m2) / 2
  }





# External Met Data -------------------------------------------------------

# Constant for PET calculation
lambda_MJ_kg <- 
  2.45

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
water_balance[, season := fifelse(month(sample_date) %in% 6:8 |
                                    (month(sample_date) == 5 & mday(sample_date) > 15) |
                                    (month(sample_date) == 9 & mday(sample_date) <= 15),
                                  "growing",
                                  "dormant")]

daily_water_balance <- 
  water_balance[,.(site_status = first(site_status),
                   season = first(season),
                   treatment = first(treatment),
                   i_water_level_cm = first(water_level_cm),
                   min_water_level_cm = min_na(water_level_cm),
                   median_water_level_cm = median_na(water_level_cm),
                   max_water_level_cm = max_na(water_level_cm),
                   Dl_cm = max_na(water_level_cm) - min_na(water_level_cm),
                   Dl_signed_cm = sign(which.max(water_level_cm) - which.min(water_level_cm)) * (max_na(water_level_cm) - min_na(water_level_cm)),
                   Dl_time_range_s = (which.max(water_level_cm) - which.min(water_level_cm)) * 900,
                   gross_precip_cm = first(precip_cm),
                   Ds_cm = last(water_level_cm) - first(water_level_cm)),
                by = .(site, sample_date)]

interception <- 
  daily_water_balance[,.(treatment = first(treatment),
                         mod = list(safe_lm(Ds_cm ~ gross_precip_cm))),
                      by = .(site, site_status, season)]

interception[, c("intercept", "slope") := map_dfr(mod, safe_coef)]
interception[, i_cm := -intercept / slope]

interception[, i_cm := ifelse(is.na(i_cm), mean_na(i_cm), i_cm),
          by = .(site_status)]

interception[i_cm < 0, i_cm := 0]


daily_water_balance[interception, 
                    `:=`(i_cm = i.i_cm,
                         net_precip_cm = pmax(0, gross_precip_cm - i.i_cm)),
                    on = c("site", "site_status", "season")]

daily_water_balance[, precip_intensity_cm_hr := net_precip_cm / (Dl_time_range_s / 3600)]

# McLaughlin ESy ----------------------------------------------------------

rr_dat <- 
  water_balance[, 
                .(site,
                  sample_time,
                  sample_year,
                  site_status,
                  season,
                  water_level_cm,
                  Ds_cm, 
                  gross_precip_cm = sh_precip_cm)]
                
# Set precip status as wet/dry
rr_dat[, status := fifelse(gross_precip_cm > 0, "wet", "dry")]

# If precip is NA set storm to 'dry' to calculate run lengths correctly
rr_dat[is.na(gross_precip_cm), status := "dry"]

# Get run lengths of wet/dry
rr_dat[, 
       precip_run_length := rep(rle(status)$lengths, times = rle(status)$lengths),
       by = .(site, sample_year)]

# Define storms by precip separated by 4 hours of dry
rr_dat[, storm := fifelse(status == "dry" & precip_run_length > 32,
                          "dry", 
                          "storm")]

# Give storms an id
rr_dat[, storm_id := paste(.BY[[1]], .BY[[2]], rleid(storm), sep = "-"), 
       by = .(site, sample_year)]

# Calculate running cumulative gross precip for each storm
rr_dat[, cumulative_gross_precip_cm := cumsum(gross_precip_cm),
       by = .(storm_id)]

# Add interception
rr_dat[interception,
       i_cm := i.i_cm,
       on = c("site", "site_status", "season")]

# Calculate cumulative net precip by removing interception losses
rr_dat[, cumulative_net_precip_cm := cumulative_gross_precip_cm - i_cm]

# Remove negative net precip values
rr_dat[, cumulative_net_precip_cm := pmax(0, cumulative_net_precip_cm)]

# Change periods of cumulative net precip to 'dry'
rr_dat[cumulative_net_precip_cm == 0,
       storm := "dry"]

# Recalculate storm ID based on net_precip storms
rr_dat[, storm_id := paste(.BY[[1]], .BY[[2]], rleid(storm), sep = "-"), 
       by = .(site, sample_year)]

# Calculate cumulative change in water level
rr_dat[, cumulative_Ds_cm := cumsum(Ds_cm),
       by = .(storm_id)]

# Calculate elapsed time and total storm length
rr_dat[, elapsed_time_hr := as.numeric(difftime(sample_time, min(sample_time), units = "hours")),
       by = .(storm_id)]
rr_dat[, storm_length_hr := max(elapsed_time_hr),
       by = .(storm_id)]

# Add preceding dry period length
rr_dat[, preceding_dry_hr := shift(storm_length_hr, 1)]
rr_dat[, preceding_dry_hr := first(preceding_dry_hr), 
       by = .(storm_id)]

# Calculate Esy
rr_dat[, esy := cumulative_net_precip_cm / cumulative_Ds_cm]

# Calculate precip intensity 
rr_dat[, precip_intensity := cumulative_net_precip_cm / elapsed_time_hr]

# There is a declining impact of precip for long precipitation events. The
# responses of ESy to water level looks like a recession curve within a single 
# storm
ggplot(rr_dat[storm == "storm", first(.SD, n = 4 * 3), by = .(storm_id)],
       aes(x = water_level_cm,
           y = esy)) +
  geom_point() + 
  facet_wrap(~site,
             scales = "free") +
  scale_y_log10()


ggplot(rr_dat[storm == "storm", last(.SD), by = .(storm_id)],
       aes(x = water_level_cm,
           y = esy)) +
  geom_point() + 
  facet_wrap(~site,
             scales = "free")

# Non-linear Sy -----------------------------------------------------------

# Bad precip estimates
# sy_dat <-
#   sy_dat[!(site == "006" & sy > 4 |
#            site == "111" & sample_date %in% ymd(c("2015-08-09", "2017-07-06", "2018-07-12")) |
#            site == "009" & sy > 1.5)]
# 
# 
# nl_sy <- 
#   sy_dat[24 * precip_intensity_cm_hr > i_cm]
# 
# nl_sy[,
#        init_b := .SD[water_level_cm <= quantile(water_level_cm, 0.1, na.rm = TRUE),
#                      mean_na(sy)],
#        by = .(site)]
# 
# nl_sy[,
#        init_c := .SD[sy == max_na(sy),
#                        mean_na(water_level_cm)],
#        by = .(site)]
# 
# nl_sy[,
#        init_m := 1 / init_c]
# 
# nl_mods <-
#   nl_sy[,
#          .(mod = list(nls(sy ~ b + m ** (water_level_cm - c),
#                             start = list(b = unique(init_b), m = 1.25,  c = unique(init_c)),
#                             # maxit = 100,
#                             control = nls.control(warnOnly = TRUE),
#                             na.action = na.exclude,
#                             data = .SD))),
#          by = .(site)]
# #
# nl_mods[, c("b", "m", "c") := map_dfr(mod, coef)]
# 
# nl_sy[nl_mods, pred_sy := i.b + i.m ** (water_level_cm - i.c), on = "site"]
# 
# ggplot(nl_sy,
#        aes(x = water_level_cm, y = sy)) +
#   geom_point() +
#   geom_line(aes(y = pred_sy), color = "red") +
#   geom_smooth(method = glmrob,
#               formula = y ~ x + I(x^2),
#               se = FALSE,
#               method.args = list(family = Gamma(inverse))) +
#   facet_wrap(~site, scales = "free")

# GLM Sy ------------------------------------------------------------------
# To consider:
# Combining multi-day precip events
# -Removing days with preceeding precipitation- Checked this, no discernible patterns
# Adding an upper limit to ESy derived as the ratio of watershed:wetland area

sy_dat <- 
  daily_water_balance[net_precip_cm > 0 & season == "growing", 
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

# This stopped working (way down at the flow adjustment models) when I removed
# low intensity storms before I removed extreme high water levels. Probably
# indicates some sort of issue
sy_dat <-
  sy_dat[24*precip_intensity_cm_hr > i_cm]


sy_dat <- 
  sy_dat[!(site %in% c("009", "077", "157") & water_level_cm < -50) & # Artificially high sy values (may be from well bottoming out)
           !(site %in% c("006", "009") & sy > 3)] # Strange high values

# Clip some extremely high water levels
# sy_dat[, upper_x := median_na(water_level_cm) + 1.28 * mad(water_level_cm, na.rm = TRUE),
#        by = .(site)]
# sy_dat <-
#   sy_dat[water_level_cm <= upper_x]

# Sy needs to be adjusted for precip intensity, which is derived as net precip
# divided by the time difference between min and max water levels. This makes 
# sense where small storms may not show any response in water table that 
# significantly deviates from the normal sine pattern, resulting in negative
# sy values

ggplot(sy_dat,
       aes(x = precip_intensity_cm_hr,
           y = sy)) +
  geom_point() +
  facet_wrap(~site, scales = "free")


# Rather than adjust for intensity, we can remove low-intensity/small storms and
# it strengthens the exponential/gamma fit 
ggplot(sy_dat[sy > 0],
       aes(x = water_level_cm,
           y = sy)) +
  geom_point() + 
  geom_point(data = sy_dat[sy > 0 & 24*precip_intensity_cm_hr <= i_cm],
             color = "red") + 
  geom_smooth(method = glm,
              method.args = list(family = Gamma(inverse)),
              se = FALSE, 
              formula = y ~ x + I(x^2)) + 
  geom_smooth(data = sy_dat[24*precip_intensity_cm_hr > i_cm], 
              method = glm, 
              method.args = list(family = Gamma(inverse)), 
              se = FALSE, 
              linetype = "dotted",
              formula = y ~ x + I(x^2)) + 
  facet_wrap(~site, scales = "free")

ggplot(sy_dat,
       aes(x = water_level_cm,
           y = sy)) +
  geom_point() + 
  geom_smooth(method = robustbase::glmrob,
              method.args = list(family = Gamma(inverse)),
              se = FALSE,
              formula = y ~ x,
            color = "blue") +
  geom_smooth(method = robustbase::glmrob,
              method.args = list(family = Gamma(inverse)),
              se = FALSE,
              formula = y ~ x + I(x^2),
            linetype = "dashed",
            color = "blue") +
  facet_wrap(~site, scales = "free")


# Large storms result in Hortonian overland flow, which likely pulses streamflow,
# resulting in lower delta WL. This means Sy is lower than would be predicted for
# the same amount of precip at a lower intensity. This could potentially be 
# remedied by looking at storm events

# GLM performs better than NL approach as measured by RMSE of daily ET vs PET
# Additionally quadratic GLM works better than linear

# RMSE of ET, PET for GLMs:
  # With outliers of et greater than 0.5 removed:
  # site   uni_rmse  quad_rmse
  # 006 0.10449259 0.09904485
  # 009 0.09186648 0.09197557
  # 053 0.09173676 0.09185025
  # 077 0.09910674 0.09971839
  # 111 0.09228345 0.09273795
  # 113 0.11032220 0.11062702
  # 119 0.10024449 0.08804028
  # 135 0.10593098 0.10539962
  # 139 0.09841121 0.09712307
  # 140 0.11733474 0.11720316
  # 151 0.10111083 0.09981528
  # 152 0.11195858 0.11091667
  # 156 0.09877201 0.08868542
  # 157 0.09595529 0.09685958

  # Without outliers removed:
  # site   uni_rmse  quad_rmse
  # 006 1.04136006 0.12778260
  # 009 0.12049691 0.12106574
  # 053 0.12071625 0.10568511
  # 077 0.17703904 0.18321488
  # 111 0.10690248 0.09273795
  # 113 0.28465484 0.11225210
  # 119 0.80264890 0.08804028
  # 135 0.34114639 0.14128561
  # 139 0.25580621 0.12745262
  # 140 0.39746485 0.15926680
  # 151 0.37738213 0.10801099
  # 152 0.11967333 0.11240844
  # 156 0.86780250 0.09194840
  # 157 0.09960029 0.10404289

# R^2 of square-root transformed PET and calculated ET
# flow_eval[, .(uni_rmse = cor(et_uni_cm_d, sqrt(pet_cm_d), 
#                              use = "pairwise.complete",
#                              method = "spearman")^2, 
#               quad_rmse = cor(et_quad_cm_d, pet_cm_d, 
#                               use = "pairwise.complete",
#                               method = "spearman")^2),
#           by = .(site)]

    # site   uni_rmse  quad_rmse
    # 006 0.03701899 0.04797261
    # 009 0.23565416 0.23575935
    # 053 0.15123796 0.16354319
    # 077 0.07540678 0.07455430
    # 111 0.19127962 0.19395045
    # 113 0.09179230 0.10021615
    # 119 0.15100786 0.19678138
    # 135 0.13064624 0.12077897
    # 139 0.08497514 0.07668538
    # 140 0.09795278 0.10035826
    # 151 0.11715590 0.11919963
    # 152 0.17469229 0.18107918
    # 156 0.21186710 0.22485594
    # 157 0.17138341 0.15422581

# Fit Models
# Visually the high model and the quad model show the best forms for et_index ~ water_level_cm
# Need to evaulate RMSE of ET:PET of low/no flow water levels to tell which 
# accutaly performs better. Also comparing the ET:PET > 1 with streamflow at the
# paired watersheds would be good
# ggplot(daily_flows, 
#        aes(x = water_level_cm, y = et_index)) + 
#   geom_point() + 
#   facet_wrap(~site, scales = "free") + 
#   scale_y_log10() + 
#   labs(title = "uni sy", 
#        caption = "Ratio of ET:PET as a function of water level. When ET:PET > 1,
#        or greater than it's usual median value (or maybe there's a linear trend
#        there), we can assume streamflow is the excess." )
#   geom_hline(aes(yintercept = 1))


safe_nls <- 
  function(data, site){
    if(site %in% c("119", "156")){
      return(lmrob(sy ~ log(200 + water_level_cm),
                   data = data))
    }
    nlrob(sy ~ b + m ** (water_level_cm - c),
          data = data,
          start = list(b = 0, m = 1.25, c = -50),
          control = nls.control(warnOnly = TRUE))
  }

sy_mods <- 
  sy_dat[, # between(net_precip_cm, 0.5, 1), 
         .(# mod = list(glmrob(sy ~ water_level_cm,
           #                    family = Gamma(inverse))),
           mod_quad = list(glmrob(sy ~ water_level_cm + I(water_level_cm^2),
                                  family = Gamma(inverse)))
# CAN"T RUN NL WITH 119 & 156
           , mod_nl = list(safe_nls(.SD, .BY[[1]]))
           # ,mod_high = list(glmrob(sy ~ water_level_cm,
           #                         family = Gamma(inverse),
           #                        data = .SD[24*precip_intensity_cm_hr >= i_cm])),
           # mod_low = list(glmrob(sy ~ water_level_cm,
           #                       family = Gamma(inverse),
           #                       data = .SD[24*precip_intensity_cm_hr < i_cm]))
           # mod_ratio = list(glmrob(sy_ratio ~ water_level_cm,
           #                         family = Gamma(inverse)))
           ),
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
                           type = "response")
            # , high_sy = predict(sy_mods[.BY[[1]], mod_high[[1]], on = "site"],
            #                   newdata = .SD,
            #                   type = "response"),
            # low_sy = predict(sy_mods[.BY[[1]], mod_low[[1]], on = "site"],
            #                   newdata = .SD,
            #                   type = "response")
            ), 
       by = .(site)]

# sy_dat[, `:=`(pred_sy = pmin(1, pred_sy),
#               quad_sy = pmin(1, quad_sy))]

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
          )[, pred_sy := pmin(1, abs(predict(sy_mods[.BY[[1]], 
                                         mod_quad[[1]], 
                                         on = "site"], 
                                 newdata = .SD, 
                                 type = "response"))), 
            by = .(site)],
       aes(x = water_level_cm)) +
  geom_line(aes(y = pred_sy),
            linetype = "dashed",
            color = "blue") +
  facet_wrap(~site, scales = "free")

# It may be worth using the non-linear models which have an asymptote at low
# water levels. The quadratic form could result in extremely low estimates of
# sy when predicting outside the water level range of sy dat

# Daily water balance from McLaughlin, 2019 -------------------------------
# water_balance[interception,
#               `:=`(i_cm = i.i_cm,
#                    net_precip_cm = pmax(0, precip_cm - i.i_cm)),
#               on = c("site", "site_status", "season")]

daily_water_balance[daily_ex_met,
                    `:=`(exrad_pet_hs_cm = i.exrad_pet_hs_cm,
                         harad_pet_hs_cm = i.harad_pet_hs_cm,
                         obsrad_pet_hs_cm = i.obsrad_pet_hs_cm,
                         etr_cm_d = i.etr_cm_d,
                         observed_rad_MJ_m2 = i.observed_rad_MJ_m2),
                    on = c("sample_date")]

daily_water_balance[site %in% sy_mods$site, 
                    sy := predict(sy_mods[.BY[[1]], mod_quad[[1]], on = "site"], 
                                  newdata = .SD[, .(water_level_cm = i_water_level_cm)],
                                  type = "response"), 
                    by = .(site)]

# Deal with exponential increase at asymptotes
daily_water_balance[, sy := pmin(1, abs(sy))]


daily_water_balance[, 
                    net_flow := Ds_cm * sy - net_precip_cm + obsrad_pet_hs_cm,
                    by = .(site, sample_date)]

ggplot(daily_water_balance,
       aes(x = i_water_level_cm,
           y = net_flow)) +
  geom_point() +
  facet_wrap(~site, 
             scales = "free")


daily_water_balance[, Ds_actual := sy * Ds_cm]

daily_water_mods <- 
  daily_water_balance[, .(mod = list(lmrob(actual_Ds_cm ~ observed_r + net_precip_cm))), by = .(site)]
daily_water_balance[, pred_Ds := predict(daily_water_mods[site == .BY[[1]], mod[[1]]], newdata = .SD) / sy, by = .(site)]


# Separate Flows ----------------------------------------------------------

water_balance[site %in% sy_mods$site, 
              sy := predict(sy_mods[.BY[[1]], mod_nl[[1]], on = "site"], 
                            newdata = .SD,
                            type = "response"), 
              by = .(site)]

# Deal with exponential increase at asymptotes
water_balance[, sy := pmin(1, abs(sy))]

# # Add asymptote for sy at low water levels:
# water_balance[sy_dat[,
#                      .(lower_thresh = quantile(water_level_cm, 0.1, na.rm = TRUE),
#                        b = .SD[water_level_cm <= quantile(water_level_cm, 0.1, na.rm = TRUE),
#                                mean_na(sy)]),
#                      by = .(site)],
#               sy := fifelse(water_level_cm < i.lower_thresh,
#                             i.b,
#                             sy),
#               on = "site"]

dat <- 
  water_balance[, .(data = list(.SD)),
                by = .(site, sample_year)]


dat[, data := map(data, calculate_detrended_g)]
dat[, data := map(data, calculate_delta_s)]
dat[, data := map(data, calculate_et)]
dat[, data := map(data, adjust_water_balance)]


flows <- 
  dat[, data[[1]], by = .(site, sample_year)]

hourly_flows <- 
  flows[, 
        .(sy = mean_na(sy),
          water_level_cm = mean_na(water_level_cm),
          precip_cm = first(precip_cm),
          Ds_cm = sum(Ds_cm),
          i_water_level_cm = first(water_level_cm),
          f_water_level_cm = last(water_level_cm),
          corrected_gamma_m = mean_na(corrected_gamma_m),
          corrected_Gamma_cm_s = mean_na(corrected_Gamma_cm_s),
          g_cm_hr = sum_na(corrected_net_in_cm_s),
          et_cm_hr = sum_na(corrected_et_cm_s)),
        by = .(site, sample_time = floor_date(sample_time, "hour"))]

# Check that water balance and water level change match
ggplot(hourly_flows[, .(site, 
                        sample_time, 
                        sy, 
                        Ds_cm = Ds_cm * sy,
                        WB = (g_cm_hr - et_cm_hr))], 
       aes(x = Ds_cm, y = WB)) + 
  geom_abline() + 
  geom_point() + 
  facet_wrap(~site, scales = "free")


hourly_flows[ex_met, 
             `:=`(etr_cm_hr = i.etr_cm_hr,
                  pet_hs_cm_hr = i.pet_hs_cm_hr),
             on = "sample_time"]

hourly_flows[, et_resid := et_cm_hr - etr_cm_hr]
hourly_flows[, et_index := et_cm_hr/etr_cm_hr]


# There is probably a good way to integrate the Sy over the day

daily_flows <- 
  hourly_flows[, 
               .(sy = median_na(sy),
                 water_level_cm = min_na(water_level_cm), 
                 Ds_cm = sum_na(Ds_cm),
                 precip_cm = first(precip_cm),
                 etr_cm_d = sum_na(etr_cm_hr), 
                 g_cm_d = sum_na(g_cm_hr),
                 et_cm_d = sum_na(et_cm_hr)), 
               by = .(site, sample_date = as.Date(sample_time, tz = "EST"))]

# Check that water balance and water level change match
ggplot(daily_flows[, .(site, 
                       sample_date, 
                       sy, 
                       Ds_cm = Ds_cm * sy,
                       WB = g_cm_d - et_cm_d)], 
       aes(x = Ds_cm, y = WB)) + 
  geom_abline() + 
  geom_point() + 
  facet_wrap(~site, scales = "free")

daily_flows[water_balance, 
            `:=`(treatment = i.treatment, 
                 treatment_period = i.treatment_period,
                 site_status = i.site_status,
                 season = i.season), 
            on = c("site", "sample_date")]

daily_flows[daily_ex_met,
  `:=`(obsrad_pet_hs_cm = i.obsrad_pet_hs_cm,
       pet_hs_cm_summed = i.pet_hs_cm_summed,
       harad_pet_hs_cm = i.harad_pet_hs_cm,
       exrad_pet_hs_cm = i.exrad_pet_hs_cm,
       mean_temp_c = i.mean_temp_c,
       min_temp_c = i.min_temp_c,
       max_temp_c = i.max_temp_c,
       temp_range_c = i.max_temp_c - i.min_temp_c,
       ha_rad_MJ_m2 = i.ha_rad_MJ_m2,
       observed_rad_MJ_m2 = i.observed_rad_MJ_m2),
  on = "sample_date"]

daily_flows[daily_water_balance,
            `:=`(i_cm = i.i_cm,
                 net_precip_cm = i.net_precip_cm),
            on = c("site", "sample_date")]

daily_flows[, `:=`(et_resid = et_cm_d - obsrad_pet_hs_cm,
                   et_index = et_cm_d / obsrad_pet_hs_cm)]

ggplot(daily_flows[site_status == "Control"],
       aes(x = obsrad_pet_hs_cm,
           y = et_cm_d)) + 
  geom_point() + 
  facet_wrap(~site, scales = "free")

ggplot(daily_flows[site_status == "Control"], 
       aes(x = water_level_cm, 
           y = et_resid)) + 
  geom_point() + 
  facet_wrap(~site, scales = "free")

ggplot(daily_flows[site_status == "Control"],
       aes(x = water_level_cm,
           y = et_cm_d)) + 
  geom_point() + 
  facet_wrap(~site, scales = "free")

# Low Water Level PET -----------------------------------------------------

no_flows <- 
  daily_flows[site_status == "Control"]

no_flows[, wl_threshold := quantile(water_level_cm, 0.1, na.rm = TRUE), 
         by = .(site)]

no_flows <- 
  no_flows[water_level_cm <= wl_threshold]

ggplot(no_flows, 
       aes(x = water_level_cm, 
           y = et_resid,
           color = season)) + 
  geom_point() + 
  facet_wrap(~site, scales = "free")

no_flows[, scaled_et := as.numeric(scale(et_cm_d)),
         by = .(site, season)]

while(nrow(no_flows[scaled_et > qnorm(0.995)]) > 0){
  no_flows <- 
    no_flows[abs(scaled_et) <= qnorm(0.995)]
  
  no_flows[, scaled_et := as.numeric(scale(et_cm_d)),
           by = .(site, season)]
}

ggplot(no_flows, 
       aes(x = obsrad_pet_hs_cm, 
           y = et_cm_d,
           color = season)) + 
  geom_point() + 
  geom_smooth(method = lmrob) +
  facet_wrap(~site, scales = "free")

et_mods <- 
  no_flows[, 
           .(mod = list(lm(et_cm_d ~ obsrad_pet_hs_cm))), 
           keyby = .(site, season)]


et_mods[, `:=`(slope = map_dbl(mod, ~coef(.x)[[1]]), 
                                r_squared = map_dbl(mod, ~summary(.x)$adj.r.squared),
                                p_value = map_dbl(mod, ~(summary(.x)$coefficients)[4]),
                                rmedianse_cm = map_dbl(mod, ~sqrt(median(resid(.x)^2))),
                                rmedianse_pct = map_dbl(mod, ~100 * sqrt(median((resid(.x)/fitted(.x))^2))))][]

et_mods <- 
  et_mods[slope != 0]


daily_flows[et_mods[, .(site, season)], 
            pred_et := predict(et_mods[CJ(.BY[[1]], .BY[[2]]), mod[[1]]],
                               newdata = .SD),
            by = .(site, season),
            on = c("site", "season")]

daily_flows[, pred_resid := et_cm_d - pred_et]

ggplot(daily_flows,
       aes(x = obsrad_pet_hs_cm,
           y = pred_et,
           color = season)) +
  geom_point() +
  geom_abline() +
  facet_wrap(~site,
             scales = "free")

ggplot(daily_flows,
       aes(x = et_cm_d,
           y = pred_et,
           color = season)) +
  geom_point() +
  facet_wrap(~site,
             scales = "free")

ggplot(daily_flows,
       aes(x = water_level_cm,
           y = pred_resid,
           color = season)) +
  geom_point() +
  facet_wrap(~site,
             scales = "free")

# Adjust Flows ------------------------------------------------------------
# The nls() models fit below are not final! They need to be improved.

ggplot(daily_flows[site_status == "Control" &
                     !(site == "053" & water_level_cm > 0) &
                     !(site == "151" & water_level_cm > 20) & 
                     !(site == "135" & water_level_cm > 40) &
                     !(site == "152" & water_level_cm > 26) &
                     !(site == "140" & water_level_cm > 26)], 
       aes(x = water_level_cm, y = et_resid)) + 
  geom_point() + 
  geom_smooth(method = nlrob, formula = y ~ b + m * exp(c*x), 
              se = FALSE, 
              method.args = list(start = list(b = -0.1, m = 0.2, c = 0.1), 
                                 maxit = 200,
                                 control = nls.control(maxiter = 1500, minFactor = 2^(-32)))) +
  facet_wrap(~site, scales = "free")

# Remove some outliers. Some of the removals are based on the Gamma glms for sy.
# They have an asymptotic relationship at the maximum water level used in the 
# sy function calculation. Perhaps the exponential function would be better 
# suited. 

flow_adjust_mods <-
  daily_flows[!(site == "111" & water_level_cm > 18) &
                !(site == "135" & water_level_cm > 25) &
                !(site == "139" & water_level_cm > 15) &
                !(site == "140" & water_level_cm > 29) &
                !(site == "151" & water_level_cm > 15) &
                !(site == "152" & water_level_cm > 20)]

# flow_adjust_mods <-
#   copy(daily_flows)

flow_adjust_mods[, scaled_resid := as.numeric(scale(pred_resid)),
                 by = .(site)]

flow_adjust_mods <- 
  flow_adjust_mods[abs(scaled_resid) <= qnorm(0.995)]

ggplot(flow_adjust_mods, 
       aes(x = water_level_cm, y = pred_resid)) + 
  geom_point() + 
  geom_smooth(method = nlrob, 
              formula = y ~ b + m * exp(c*x),
              se = FALSE, 
              method.args = list(start = list(b = 0, m = 0.1, c = 0.2))) +
  facet_wrap(~site, scales = "free")

# Need to get nlrob working here

# Figure out initial estimates


# inits <- 
#   fread(text = 
#         "site, b, m, c
#         006, -0.1, 0.1, 0.2
#         009, -0.1, 0.1, 0.2
#         053, -0.1, 1, 0.1
#         077, -0.1, 0.1, 0.2
#         111, -0.1, 0.1, 0.2
#         113, -0.1, 0.1, 0.15
#         119, -0.1, 0.1, 0.2
#         135, -0.1, 0.1, 0.2
#         139, -0.1, 0.1, 0.2
#         140, -0.1, 0.04, 0.2
#         151, -0.1, 0.1, 0.2
#         152, -0.1, 0.1, 0.2
#         156, -0.1, 0.1, 0.2
#         157, -0.1, 0.1, 0.2",
#         colClasses = c(site = "character"),
#         key = "site")

inits <- 
  fread(text = 
          "site, b, m, c
        006, -0.1, 0.3, 0.2
        009, -0.1, 0.1, 0.2
        053, -0.1, 0.1, 0.2
        077, -0.1, 0.1, 0.2
        111, -0.1, 0.1, 0.2
        113, -0.1, 0.1, 0.2
        119, -0.1, 0.2, 0.03
        135, -0.1, 0.1, 0.2
        139, -0.1, 0.1, 0.2
        140, -0.1, 0.1, 0.2
        151, -0.1, 0.1, 0.2
        152, -0.1, 0.1, 0.2
        156, -0.1, 0.1, 0.2
        157, -0.1, 0.1, 0.2",
        colClasses = c(site = "character"),
        key = "site")

# 119 look like offset y ~ log(x + 150) fits okay

flow_adjust_mods <- 
  flow_adjust_mods[, .(dat = list(.SD)), keyby = .(site)]

# flow_adjust_mods[,
#                  mod := map2(site,
#                              dat,
#                              ~nlrob(et_resid ~ b + m * exp(c*water_level_cm),
#                                     data = .y,
#                                     start = inits[.x, .(b, m, c)],
#                                     na.action = na.exclude,
#                                     maxit = 500,
#                                     control = nls.control(warnOnly = TRUE, maxiter = 1500, minFactor = 2^(-32))))]

flow_adjust_mods[, init := list(list(b = 0, m = 0.1, c= 0.2))]
flow_adjust_mods[site == "156",
                 init := list(list(b = -0.012, m = 0.05, c = 0.04))]

flow_adjust_mods[,
                 mod := map2(dat, init,
                            ~nlrob(pred_resid ~ b + m * exp(c*water_level_cm),
                                   data = .x,
                                   start = .y,
                                   na.action = na.exclude))]

# flow_adjust_mods[,
#                  mod := map(dat,
#                             ~lmrob(et_resid ~ water_level_cm,
#                                    data = .x[water_level_cm < 0]))]

# flow_adjust_mods[site %in% c("119", "156"),
#                  mod := map(dat,
#                             ~lmrob(et_resid ~ log(150 + water_level_cm),
#                                    data = .x))]

daily_flows[site %in% unique(flow_adjust_mods$site),
            flow_cm := predict(flow_adjust_mods[.BY[[1]], mod[[1]]],
                               newdata = .SD),
            by = .(site)]

ggplot(daily_flows[site_status == "Control"],
       aes(x = water_level_cm, 
           y = flow_cm)) + 
  geom_point() + 
  facet_wrap(~site, scales = "free")

daily_flows[, aet_cm_d := et_cm_d - flow_cm]

daily_flows[, sample_month := month(sample_date)]

ggplot(daily_flows[net_precip_cm == 0 & site_status == "Control" & abs(aet_cm_d) < 1.5],
       aes(x = obsrad_pet_hs_cm,
           y = aet_cm_d)) + 
  geom_point() + 
  geom_abline() +
  geom_smooth(method = lmrob,
              formula = y ~ 0 + x,
              se = FALSE) +
  facet_wrap(~site, scales = "free")

# R-squared of ET estimation
daily_flows[site_status == "Control", 
            .(mod = list(lmrob(aet_cm_d ~ 0 + obsrad_pet_hs_cm))), 
            by = .(site)][, `:=`(slope = map_dbl(mod, ~coef(.x)[[1]]), 
                                 r_squared = map_dbl(mod, ~summary(.x)$adj.r.squared),
                                 p_value = map_dbl(mod, ~(summary(.x)$coefficients)[4]),
                                 rmedianse_cm = map_dbl(mod, ~sqrt(median(resid(.x)^2))),
                                 rmedianse_pct = map_dbl(mod, ~100 * sqrt(median((resid(.x)/fitted(.x))^2))))][]

# Removing Sy from all low intensity storms (harad for resid calcuation. Rerun with obsrad)
# site         mod     slope r_squared       p_value   rmse_cm  rmse_pct
# 009 <lmrob[23]> 1.0468900 0.8075816  1.816937e-60 0.1517488  2.733062
# 053 <lmrob[23]> 0.7229603 0.7679042  3.822383e-33 0.1486921 18.708524
# 077 <lmrob[23]> 0.8961368 0.6121542  9.104044e-27 0.3261004  5.657638
# 111 <lmrob[23]> 0.8305789 0.7343130  1.637171e-86 0.2093083  2.357422
# 113 <lmrob[23]> 0.9326998 0.8234367 4.068339e-139 0.2098048  3.859579
# 119 <lmrob[23]> 1.1374668 0.7416650  1.103770e-44 0.2151114  1.802269
# 135 <lmrob[23]> 0.8733030 0.7535042 4.053265e-121 1.0305379  6.601660
# 139 <lmrob[23]> 1.2030189 0.7758504  1.013224e-24 0.2289345  3.522902
# 140 <lmrob[23]> 1.0052801 0.8731178  4.770734e-80 0.1931104  4.473873
# 151 <lmrob[23]> 1.0016261 0.6794711  1.336575e-35 6.5952237 43.204789
# 152 <lmrob[23]> 0.9053358 0.8835780 7.342773e-179 0.1295101  2.166442
# 156 <lmrob[23]> 1.0722171 0.8106684  6.811327e-60 0.1843767  2.289339
# 157 <lmrob[23]> 0.7809795 0.6439957  8.888347e-83 0.2089187  2.952013

# Removing only negative Sy
# site         mod     slope r_squared       p_value   rmse_cm  rmse_pct
# 009 <lmrob[23]> 1.0762209 0.8139123  4.611856e-66 0.1587956  3.527756
# 053 <lmrob[23]> 0.8123119 0.7845350  3.553563e-40 0.1456170 18.855250
# 077 <lmrob[23]> 1.0496294 0.7384602  1.277528e-40 0.2407469  4.326870
# 111 <lmrob[23]> 0.8717127 0.7658312  4.112514e-92 0.2053895  2.240297
# 113 <lmrob[23]> 0.7447574 0.8255592 3.271574e-136 0.1312721  3.421390
# 119 <lmrob[23]> 1.1289198 0.6791274  1.116394e-38 0.2304309  1.856575
# 135 <lmrob[23]> 0.9778028 0.8108758 2.796589e-149 0.3391158  2.419453
# 139 <lmrob[23]> 1.1896962 0.7935207  1.672269e-27 0.2082007  3.322835
# 140 <lmrob[23]> 0.7330046 0.8626001  4.489006e-80 0.1787243  4.544009
# 151 <lmrob[23]> 0.9955179 0.7437711  4.851038e-35 0.2272365  5.416622
# 152 <lmrob[23]> 0.9341568 0.8948362 9.721541e-190 0.1228548  2.269515
# 156 <lmrob[23]> 1.1766212 0.8359266  7.862597e-66 0.1912364  2.016008
# 157 <lmrob[23]> 0.8869471 0.6754919  8.550353e-90 0.2120567  2.226756



# Multi-Driver Outlfow Models ---------------------------------------------

daily_flows[, water_level_m := 0.01 * water_level_cm]

multi_mods <- 
  daily_flows[site_status == "Control" & !(site %in% c("119", "156")), 
              .(mod = list(nls(et_cm_d ~ b + mpet * obsrad_pet_hs_cm + mwl * exp(c*water_level_m), 
                               start = list(b = 0, mpet = 1, mwl = 0.2, c = 0.1),
                               data = .SD))),
              keyby = .(site)]


multi_mods2 <-
  daily_flows[precip_cm == 0 & !(site %in% c("119", "156")),
              .(mod = list(nls(et_cm_d ~ b + mpet * obsrad_pet_hs_cm + mwl ** (water_level_cm - c),
                               start = list(b = 0, mpet = 1, mwl = 1.1, c = 5),
                               na.action = na.exclude,
                               data = .SD))),
              keyby = .(site)]

daily_flows[site %in% unique(multi_mods$site), 
            multi_et := predict(multi_mods[.BY[[1]], mod[[1]]], newdata = .SD),
            by = .(site)]

ggplot(daily_flows[net_precip_cm == 0 & site_status == "Control"],
       aes(x = et_cm_d,
           y = multi_et)) + 
  geom_point(aes(color = season)) +
  geom_abline() +
  facet_wrap(~site,
             scales = "free")

daily_flows[precip_cm == 0, .(rmse = caret::RMSE(multi_et))]

# Check Flow Against Flume Stage ------------------------------------------

flume <- 
  fread("../Data/Raw/Manually_Entered/well_and_flume_measurements.csv",
        select = c(site = "numeric", location = "character", 
                   sample_date = "character", manual_level = "numeric", 
                   value_unit = "character"))

flume[, sample_date := mdy(sample_date)]
flume[, site := sprintf("%03d", site)]
flume <- flume[location == "fl"]
flume[value_unit == "ft", manual_level := manual_level * 25.4]

daily_flows[flume, 
            flume_stage_cm := i.manual_level,
            on = c("site", "sample_date")]

## LOOK AT 2013 FLUME STAGE
daily_flows[, flume_stage_cm := ifelse(year(sample_date) == 2013, 
                                       0.1 * flume_stage_cm, 
                                       flume_stage_cm)]


daily_flows[, offset_flow := flow_cm - min_na(flow_cm), by = .(site)]

ggplot(daily_flows[!is.na(flume_stage_cm) & flow_cm < 1],
       aes(x = flow_cm,
           y = flume_stage_cm,
           color = site_status)) +
  geom_point() +
  facet_wrap(~ site,
             scales = "free")

daily_flows[!is.na(flume_stage_cm), 
            .(mod = list(lmrob(flume_stage_cm ~ 0 + offset_flow))), 
            by = .(site)][, `:=`(slope = map_dbl(mod, ~coef(.x)[[1]]), 
                                 r_squared = map_dbl(mod, ~summary(.x)$adj.r.squared),
                                 p_value = map_dbl(mod, ~(summary(.x)$coefficients)[4]),
                                 rmedianse_cm = map_dbl(mod, ~sqrt(median(resid(.x)^2))),
                                 rmedianse_pct = map_dbl(mod, ~100 * sqrt(median((resid(.x)/fitted(.x))^2))))][]


# Net Inflow Drivers ------------------------------------------------------
# It looks like the simple driver of precip will work. need to also look at
# a recession curve for precip to add longer-term precip impact. Then also need
# to consider spring inflow (probably need better representation of melt)


daily_flows[, cumulative_precip_cm := cumsum(precip_cm),
            by = .(site, year(sample_date))]

daily_flows[, cumulative_pet_cm := cumsum(obsrad_pet_hs_cm),
            by = .(site, year(sample_date))]

daily_flows[, water_availability := cumulative_precip_cm - cumulative_pet_cm]
daily_flows[, water_index := cumulative_pet_cm / cumulative_precip_cm]


ggplot(daily_flows[daily_flows[, .(upper_g = mean_na(g_cm_d) + 1.96 * sd(g_cm_d, na.rm = TRUE)), by = .(site)], on = "site"][g_cm_d < upper_g],
       aes(x = precip_cm,
           y = g_cm_d)) +
  geom_point() + 
  facet_wrap(~site,
             scales = "free")


# Recession Curves --------------------------------------------------------


daily_flows[, precip_status := as.numeric(!(precip_cm > 0))]
daily_flows[is.na(precip_status), precip_status := 1]
daily_flows[, 
            precip_id := paste(.BY[[1]], .BY[[2]], rleid(precip_status), sep = "-"),
            by = .(site, sample_year = year(sample_date))]
daily_flows[, 
            dry_days := cumsum(precip_status),
            by = .(precip_id)]
daily_flows[,
            delta_g := g_cm_d - first(g_cm_d),
            by = .(precip_id)]

daily_flows[daily_water_balance,
            max_water_level_cm := i.max_water_level_cm,
            on = c("site", "sample_date")]

daily_flows[, storm_peak := shift(max_water_level_cm),
            by = .(site, year(sample_date))]

daily_flows[,
            delta_Ds := (water_level_cm - first(storm_peak)),
            by = .(precip_id)]

daily_flows[, 
            storm_precip_cm := sum_na(net_precip_cm),
            by = .(precip_id)]

daily_flows[, 
            storm_precip_cm := ifelse(precip_status == 0, NA, shift(storm_precip_cm, 1))]

daily_flows[, storm_precip_cm := first(storm_precip_cm), by = .(precip_id)]

ggplot(daily_flows[between(dry_days, 1, 10)],
       aes(x = dry_days,
           y = delta_Ds)) +
  geom_line(aes(group = precip_id,
                color = season)) +
  geom_quantile(quantiles = 0.5, 
                data = daily_flows[between(dry_days, 1, 7)],
                size = rel(2),
                color = "black") +
  facet_wrap(~site, scales = "free")

ggplot(daily_flows[(daily_flows[, between(max_na(dry_days), 3, 5) & season == "dormant", by = .(precip_id)]$V1)],
       aes(x = dry_days,
           y = delta_Ds)) +
  geom_line(aes(group = precip_id)) +
  facet_wrap(~site, scales = "free")

recessions <- 
  daily_flows[daily_flows[, between(max_na(dry_days), 3, 7), by = .(precip_id)][(V1)],
              on = "precip_id",
              nomatch = NULL][!is.na(delta_Ds)]

recession_ks <- 
  recessions[, 
             .(k = coef(lm(delta_Ds ~ 0 + dry_days))[[1]],
               season = first(season),
               site_status = first(site_status),
               sy = first(sy),
               water_level_cm = first(water_level_cm)),
             by = .(site, precip_id)]

ggplot(recession_ks[season == "growing"], 
       aes(x = water_level_cm, 
           y = k, 
           color = site_status)) + 
  geom_point() + 
  ggtitle("Growing Season Recession By Water Level") +
  facet_wrap(~site, 
             scales = "free")

ggplot(recession_ks[season == "dormant"], 
       aes(x = water_level_cm, 
           y = k, 
           color = site_status)) + 
  geom_point() + 
  ggtitle("Dormant Season Recession By Water Level") +
  facet_wrap(~site, 
             scales = "free")

ggplot(recession_ks[site_status == "Control"], 
       aes(x = water_level_cm, 
           y = k, 
           color = season)) + 
  geom_point() + 
  ggtitle("Control Period Recession By Water Level & Season") +
  facet_wrap(~site, 
             scales = "free")

ggplot(recession_ks[site_status == "Control"], 
       aes(x = site, 
           y = k, 
           fill = season)) + 
  geom_boxplot() 

seasonal_k <- 
  recession_ks[, 
               .(mod = list(lm(k ~ water_level_cm))),
               by = .(site, season)]

seasonal_k <- 
  seasonal_k[, map_dfr(mod, broom::tidy),
           by = .(site, season)]

seasonal_k[p.value > 0.05 & term == "water_level_cm",
           estimate := 0]

seasonal_k <- 
  dcast(seasonal_k,
      site + season ~ term,
      value.var = "estimate")

setnames(seasonal_k, c("site", "season", "intercept", "slope"))



# Soylu -------------------------------------------------------------------



# Basic Water Balance -----------------------------------------------------

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

# Get peak water level from the storm period
# There could be outliers produced if there is a large storm one day, followed
# by a small storm the next day. Most of the recession could occur during a
# 'non-dry' day. A multiday storm could also have a large second day event that
# is not quite a new max, so I can't just say that the day with the max should
# be considered day 0. Going to compromise with max water level from the last
# precip day in a storm period. That's what .SD[dry_status == 0, last(max_water_level_cm)]
# does

daily_water_balance[, storm_peak := .SD[dry_status == 0, last(max_water_level_cm)],
            by = .(storm_id)]

# Change in daily water level from the storm peak. 
daily_water_balance[,
                    delta_Ds_cm := median_water_level_cm - storm_peak]

# Normalized Ds within a storm period
daily_water_balance[,
                    recession_cm_day := delta_Ds_cm / dry_days]

# Add up precip amount
daily_water_balance[, 
            storm_precip_cm := sum_na(net_precip_cm),
            by = .(storm_id)]

# # Predict day 1 draw down rate from precip size, then fit a log decay to the median
# daily drawdown rate for a season. Could consider using a true intervention 
# approach from time-series analysis. I can predict delta_Ds_cm when day == 1, 
# using delta_Ds_cm * sy ~ log(storm_precip_cm), or delta_Ds_cm ~ sy {family = gaussian(inverse)}

# Dormant season recession_cm_day never quite approaches the median or mean 
# daily Ds. It probably should go to recession constant (delta_Ds_cm ~ dry_days),
# which can be approximated as the mean recession_cm_day on day 7
ggplot(daily_water_balance[between(dry_days, 1, 8) & site_status == "Control" & storm_precip_cm > i_cm],
       aes(x = dry_days,
           y = recession_cm_day,
           color = season)) +
  geom_point() + 
  geom_line(aes(group = response_id),
            alpha = 0.3) +
  facet_wrap(~site, scales = "free")

ggplot(daily_water_balance[dry_days == 1 & site_status == "Control" & storm_precip_cm > i_cm],
       aes(x = storm_precip_cm,
           y = delta_Ds_cm * sy)) +
  geom_point() +
  geom_smooth(method = lmrob, 
              formula = y ~ log(x)) +
  facet_wrap(~site, scales = "free")

ggplot(daily_water_balance[dry_days == 1 & site_status == "Control" & storm_precip_cm > i_cm],
       aes(x = sy,
           y = delta_Ds_cm)) +
  geom_point() +
  geom_smooth(method = glmrob, 
              formula = y ~ x,
              method.args = list(family = gaussian(inverse))) +
  facet_wrap(~site, scales = "free")

recession_mods <- 
  daily_water_balance[, 
              .(dat = list(.SD),
              by = .(site)]

# WORK IN REAL SPACE. USE DS_ACTUAL RATHER THAN DIVIDING PRECIP AND PET BY SY. 
# I WOULD RATHER HAVE LOW SY VALUES LEAD TO A ZERO, RATHER THAN INFINITY


data %>% 
  nest_by(site) %>% 
  model_initial_Ds_response(glmrob(Ds_cm ~ sy, family = Gamma(inverse)),
                            dry_days ==1,
                            output = initial_ds_hat) %>% 
  model_recession(recession_cm_day ~ b + m * log(dry_days),
                  data = 
                    dry_days = c(1, 7),
                    recession_cm_day = c(initial_ds_hat, mean(recession_cm_day | dry_days == 7)),
                  output = recession_function)


precip / sy = Ds
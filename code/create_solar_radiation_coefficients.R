# Generate coefficients for Hargreaves Solar Radiation Calculations
# Data downloaded from https://download.synopticdata.com on 2020-12-03

source('code/climate_packages.R')

meso_files <- 
  list.files("data/mesowest_met/",
             full.names = TRUE)

ex_met <- 
  data.table(file = meso_files)

ex_met[, station := substr(basename(file), 1, 5)]

ex_met[, dat := lapply(file, 
                       function(x){
                         cbind(fread(x,
                                     select = 2:4,
                                     skip = 12,
                                     col.names = c("sample_time", "air_temp_c", "solar_rad_w_m2")),
                               setnames(data.table(t(fread(x,
                                                           skip = 6,
                                                           nrows = 3,
                                                           header = FALSE,
                                                           sep = ":")[, 2])),
                                        c("latitude", "longitude", "elevation_ft")))
                         })]

# It looks like the downloaded files do not apply DST at these stations. Assuming
# that all times are standard, not daylight then.
tzones <- 
  fread(text = 
          "station,tz
          WKFM4,Etc/GMT+6
          PIEM4,Etc/GMT+5
          KTNM4,Etc/GMT+5
          WMTM4,Etc/GMT+6
          GDNW3,Etc/GMT+6
          F0708,Etc/GMT+6
          BPLM4,Etc/GMT+5
          WINM4,Etc/GMT+5")

# Add time zones to each station
ex_met[tzones, 
       tz := i.tz,
       on = "station"]

# Convert time to local and convert elevation & solar radiation units
ex_met[, dat := map2(dat, tz, 
                     ~.x[, `:=`(sample_time = setattr(sample_time, "tzone", .y),
                                elevation_m = elevation_ft / 3.2808,
                                solar_rad_MJ_m2_hr = solar_rad_w_m2 * 3600 * 1e-6,
                                solar_rad_w_m2 = NULL,
                                elevation_ft = NULL)])]

# Summarize to daily data
# We only want to use complete days, so na.rm = FALSE
ex_met[, daily_dat := map2(dat, tz,
                           ~.x[, .(latitude = first(latitude),
                                   longitude = first(longitude),
                                   mean_temp_c = mean(air_temp_c), 
                                   min_temp_c = min(air_temp_c), 
                                   max_temp_c = max(air_temp_c), 
                                   observed_rad_MJ_m2 = sum(solar_rad_MJ_m2_hr)), 
                               by = .(sample_date = as.Date(sample_time, tz = .y))])]

# Remove days where max temp is not greater than min temp
# This likely indicates an error
ex_met[, daily_dat := map(daily_dat, 
                          ~.x[max_temp_c > min_temp_c])]

# Clip off partial days from the beginning & end of the data
ex_met[, daily_dat := map(daily_dat, 
                          ~.x[between(sample_date, 
                                      ymd("2011-10-01"),
                                      ymd("2020-09-30"))])]
    

# Compare Solar Radiation Models ------------------------------------------

# Calculate clear sky transmissivity for BC data
ex_met[, clear_sky_t := map_dbl(daily_dat,
                                ~cst(RefRad = .x$observed_rad_MJ_m2,
                                     days = .x$sample_date,
                                     lat = first(.x$latitude)))]

# Calculate Hargreaves Coefs
ex_met[, c("Ha", "Hb", "Hr2") := map_dfr(daily_dat,
                                         ~hacal(lat = first(.x$latitude), 
                                                days = .x$sample_date,
                                                rad_mea = .x$observed_rad_MJ_m2,
                                                tmax = .x$max_temp_c,
                                                tmin = .x$min_temp_c))]

# Calculate Bristow-Campbell Coefs
ex_met[, BCb := map2_dbl(daily_dat, clear_sky_t,
                         ~bccal(lat = first(.x$latitude), 
                                days = .x$sample_date,
                                rad_mea = .x$observed_rad_MJ_m2,
                                Tmax = .x$max_temp_c,
                                Tmin = .x$min_temp_c,
                                tal = .y))]

# Calculate solar radiation using Mahmood-Hubbard model
ex_met[, daily_dat := map(daily_dat,
                          ~.x[, mh_rad := mh(days = sample_date, 
                                             lat = first(latitude),
                                             Tmax = max_temp_c,
                                             Tmin = min_temp_c)])]

# Calculate solar radiation using Mahmood-Hubbard model
ex_met[, daily_dat := pmap(list(daily_dat, Ha, Hb),
                          function(X, Ha, Hb){
                            set(x = X, 
                                j = "ha_rad",
                                value = ha(days = X$sample_date,
                                           lat = first(X$latitude),
                                           lon = first(X$longitude),
                                           A = Ha,
                                           B = Hb,
                                           Tmax = X$max_temp_c,
                                           Tmin = X$min_temp_c))
                            })]

# Calculate solar radiation using Bristow-Campbell
ex_met[, daily_dat := pmap(list(daily_dat, BCb, clear_sky_t),
                           function(X, B, Ta){
                             set(x = X, 
                                 j = "bc_rad",
                                 value = bc(days = X$sample_date,
                                            lat = first(X$latitude),
                                            BCb = B,
                                            Tmax = X$max_temp_c,
                                            Tmin = X$min_temp_c,
                                            tal = Ta))
                           })]

# Unnest daily data
daily_ex_met <- 
  ex_met[, daily_dat[[1]], by = .(station)]

# Evaluate Solar Rad Models
daily_ex_met[, 
             lapply(.SD, modeval, measured = observed_rad_MJ_m2, stat = "RRMSE"),
             by = .(station),
             .SDcols = c("mh_rad", "ha_rad", "bc_rad")]

# Evaluate Solar Rad Models Seasonally
daily_ex_met[, season := fcase(month(sample_date) %in% 3:5, "mam",
                               month(sample_date) %in% 6:8, "jja",
                               month(sample_date) %in% 9:11, "son", 
                               default = "djf")]

daily_ex_met[, 
             lapply(.SD, modeval, measured = observed_rad_MJ_m2, stat = "RRMSE"),
             by = .(station, season),
             .SDcols = c("mh_rad", "ha_rad", "bc_rad")]

# All of the models perform fairly well for many of the measure, but when
# looking at signed-error HA really outperforms the other models. Is this just
# because it is a linear regression fitting A & B, so OLS dictates that is has a
# very small signed error? Tried about applying seasonal coefficients, but that
# has worse performance for BC and seemingly no impact for HA. 

# Plotting the data reveals more about the seasonal fit of the models. Based on 
# the approximately equal performance measures. Bristow-Campbell shows a better
# linear fit to the 1:1 line, but does does not perform as well at low solar rad
# (djf, son). The non-linear relationship between obs ~ measured for HC, and to 
# a less extent MH, leads to me use BC as my final method for calculating daily
# solar radiation. This could potentially be imporved if I downloaded an extra-
# terrestrial raditation dataset rather than calculating it. I did attempt to 
# apply an additional correction using the slope and intercept from modeval for 
# HA and MC. It can be done: 
# corrected_rad = calc_rad + slope * (calc_rad - intercept / slope)
# which rotates the values to get closer to 1:1. It still has a lot of error and
# is probably not worth the additional complexity

ggplot(daily_ex_met, 
       aes(x = observed_rad_MJ_m2, 
           y = bc_rad)) + 
  geom_point(alpha = 0.3) + 
  geom_abline() + 
  facet_grid(station~season)

solrad_coefs <- 
  ex_met[, .(station_name = station, 
             lat = map_dbl(dat, ~first(.x$latitude)), 
             lon = map_dbl(dat, ~first(.x$longitude)), 
             B = BCb, 
             C = 2)]

fwrite(solrad_coefs,
       "output/tabular/bristow-campbell_solar_radiation_coefficients.csv")

fwrite(daily_ex_met[, .(station_name = station,
                        sample_date,
                        tmean_c = mean_temp_c,
                        tmin_c = min_temp_c,
                        tmax_c = max_temp_c,
                        solrad_MJ_m2 = observed_rad_MJ_m2)],
       "output/tabular/mesowest_daily_temperature_and_radiation.csv")

# b_grid <- 
#   CJ(lat = seq(min(solrad_coefs$lat), max(solrad_coefs$lat), by = 0.01),
#      lon = seq(min(solrad_coefs$lon), max(solrad_coefs$lon), by = 0.01))[
#        , B := predict(glm(B ~ lat + lon,
#                           family = Gamma(inverse),
#                           data = ex_met[, .(station, lat = map_dbl(dat, ~first(.x$latitude)), lon = map_dbl(dat, ~first(.x$longitude)), B = BCb, C = 2)]),
#                       type = "response", 
#                       newdata = .SD)]
# 
# leaflet::leaflet(data = solrad_coefs) %>% 
#   leaflet::addTiles()  %>% 
#   # addCircles(data = b_grid, 
#   #            lng = ~lon, 
#   #            lat = ~lat, 
#   #            color = ~colorNumeric(palette = "viridis", domain = b_grid[, B])(B)) %>%
#   leaflet::addCircleMarkers(lng = ~lon,
#                             lat = ~lat, 
#                             color = ~colorNumeric(palette = "viridis", 
#                                                   domain = b_grid[,B])(B))


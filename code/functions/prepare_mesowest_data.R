##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author Joe Shannon
##' @export
prepare_mesowest_data <- 
  function(dir, start.date, end.date) {
    
    meso_files <- 
      list.files(dir,
                 full.names = TRUE)
    
    ex_met <- 
      data.table(file = meso_files)
    
    ex_met[, station := substr(basename(file), 1, 5)]
    
    ex_met[, dat := lapply(file, 
                           function(x){
                             cbind(fread(x,
                                         select = 2:4,
                                         skip = 12,
                                         col.names = c("sample_time", "air_temp_c", "solar_rad_w_m2")))
                           })]
    
    ex_met[, c("lat", "lon", "elevation_ft") := rbindlist(lapply(file,
                                                                   function(.x){data.table(t(fread(.x,
                                                                                       skip = 6,
                                                                                       nrows = 3,
                                                                                       header = FALSE,
                                                                                       sep = ":")[, 2]))}))]
    
    # Convert Elevation to meters
    ex_met[, `:=`(elevation_m = elevation_ft / 3.2808,
                  elevation_ft = NULL)]
    
    # Add time zones to each station
    # It looks like the downloaded files do not apply DST at these stations. Assuming
    # that all times are standard, not daylight then.
    ex_met[, tz := get_gmt_tzone(lat = lat,
                                 lon = lon)]
    
    # Convert time to local and convert solar radiation units
    ex_met[, dat := map2(dat, tz, 
                         ~.x[, `:=`(sample_time = setattr(sample_time, "tzone", .y),
                                    solar_rad_MJ_m2_hr = solar_rad_w_m2 * 3600 * 1e-6,
                                    solar_rad_w_m2 = NULL)])]
    
    # Summarize to daily data
    # We only want to use complete days, so na.rm = FALSE
    ex_met[, daily_dat := map2(dat, tz,
                               ~.x[, .(mean_temp_c = mean(air_temp_c), 
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
                                          start.date,
                                          end.date)])]
    
    # Unnest daily data
    daily_ex_met <- 
      ex_met[, daily_dat[[1]], 
             by = .(station, lat, lon)]
   
    daily_ex_met[, doy := yday(sample_date)]
    
    solrad <- 
      daily_ex_met[, 
                   .(doy = 1:366,
                     et_solrad_MJ_m2 = extrat(1:366, radians(first(lat)))$ExtraTerrestrialSolarRadiationDaily), 
                   by = .(station)]
    
    daily_ex_met[solrad, 
                 et_solrad_MJ_m2 := i.et_solrad_MJ_m2, 
                 on = c("doy", "station")]
    
    # Match to standard column names for project
    daily_ex_met[, .(station_name = station,
                     sample_date,
                     lat,
                     lon,
                     tmin_c = min_temp_c,
                     tmax_c = max_temp_c,
                     solrad_MJ_m2 = observed_rad_MJ_m2,
                     et_solrad_MJ_m2)]
}

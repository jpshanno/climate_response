##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param ncei.path
##' @param ghcnd.units
##' @return
##' @author Joe Shannon
##' @export
prepare_swg_data <- 
  function(ncei.paths, 
           ghcnd.units,
           additional.stations = NULL,
           force.download = FALSE){

    ghcnd_stations <- 
      setDT(read.fortran("data/ncei_data/ghcnd/station_inventory.txt",
                   format = c( "A11", "F9", "F10", "F7", "X1","A2",
                               "X1","A30", "X1", "A3", "X1", "A3", "X1", "A5"), 
                   col.names = c("ID", "LAT", "LON", "ELEV", "ST", "NAME","GSN", "HCN", "WMOID"),
                   comment.char=""))
    
    ghcnd_stations[, station_name := gsub("_$", "", gsub("[^A-z09]{1,}", "_", tolower(NAME)))]
    
    if(!is.null(additional.stations)){
      
      fs::dir_create("data/ncei_extra")
      
      additional.names <- 
        ghcnd_stations[ID %in% additional.stations, 
                       station_name]
      
      
      urls <- 
        paste0("https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/all/",
               additional.stations,
               ".dly")
      
      filenames <- 
        fs::path("data/ncei_extra", additional.names, ext = "dly")
        
      missing_files <- 
        !file.exists(filenames)
        
      if(force.download){
        missing_files <- 
          rep(TRUE, missing_files)
      }
      
      if(sum(missing_files) != 0){
        download.file(url = urls[missing_files],
                      destfile = filenames[missing_files],
                      method = "libcurl")
      }
      
      ncei.paths <- 
        c(ncei.paths, filenames)
      
      }
    
    ghcnd_data <- 
      lapply(ncei.paths,
             prepare_ghcnd_data,
             ghcnd.units = ghcnd.units) %>% 
      rbindlist(fill = TRUE)
    
    # Keep 1979-12-31 to avoid NAs when calculating lags with format_swg
    ghcnd_data <- 
      ghcnd_data[between(sample_date, as.Date("1979-12-31"), as.Date("2009-12-31")),
                 .(station_name, station_id, sample_date, tmin_c, tmax_c, precip_cm = prcp_cm)]
    
    ghcnd_data[ghcnd_stations,
               `:=`(lon = i.LON, 
                    lat = i.LAT, 
                    elevation_m = ELEV), 
               on = c("station_id" = "ID")]
    
    expanded <- 
      CJ(station_name = unique(ghcnd_data$station_name),
         sample_date = unique(ghcnd_data$sample_date))
    
    expanded <- 
      expanded[ghcnd_data[, .SD[1, .(station_id, lon, lat, elevation_m)], by = .(station_name)],
                 on = "station_name"]
    
    ghcnd_data <- 
      ghcnd_data[expanded,
                 on = intersect(names(ghcnd_data), names(expanded))]
    
    missing_days <- 
      ghcnd_data[, c(days_in_month = .N, map(.SD, ~sum(is.na(.x)))),
                 by = .(station_name, sample_year = year(sample_date), sample_month = month(sample_date)),
                 .SDcols = c("precip_cm", "tmin_c", "tmax_c")]
    
    stations_to_use <- 
      missing_days[, .(n_months = sum(precip_cm <= 0.75 * days_in_month & 
                                        tmin_c <= 0.75 * days_in_month & 
                                        tmax_c <= 0.75 * days_in_month)),
                   by = .(station_name)][n_months > 336, .(station_name)]
    
    ghcnd_data <- 
      ghcnd_data[stations_to_use,
                 on = "station_name"]
  
    format_swg_data(ghcnd_data)
      
  }

format_swg_data <- 
  function(data){
    
    dat_names <- 
      names(data)
    
    stopifnot("station_name" %in% dat_names,
              "sample_date" %in% dat_names,
              "precip_cm" %in% dat_names,
              "tmax_c" %in% dat_names,
              "tmin_c" %in% dat_names)
    
    data <- 
      copy(data)
  
    data[, `:=`(sample_month = month(sample_date),
                sample_year = year(sample_date),
                sample_season = as.climate_season(sample_date, TRUE))]
    
    # Use NA_real_ not 0 because the precip amount model only deals with the
    # amount, so the relationship between the covariates and precip amount would
    # get muddled if we included 0s
    data[, `:=`(precip_occur = as.numeric(precip_cm > 0.01),
                precip_amount = ifelse(precip_cm < 0.01, NA_real_, precip_cm))]
    
    data[, `:=`(tmin_l1 = shift(tmin_c, 1),
                tmax_l1 = shift(tmax_c, 1),
                precip_occur_l1 = shift(precip_occur, 1)),
         by = .(station_name)]
    
    # After lagging variables remove 1979-12-31 to avoid messing up trend &
    # harmonics
    data <- 
      data[sample_date != as.Date("1979-12-31")]
    
    data[, temperature_trend := seq(-1, 1, length.out = .N),
         by = .(station_name)]
    
    data[, `:=`(harmonic_cos = cos((2*pi*(1:.N / .N))),
                harmonic_sin = sin((2*pi*(1:.N / .N)))),
         by = .(station_name, sample_year)]
    
    
    # Site & Regional Monthly & Seasonal Temperatures
    data[, `:=`(site_monthly_tmin_c = mean_na(tmin_c), 
                site_monthly_tmax_c = mean_na(tmax_c)),
         by = .(station_name, sample_year, sample_month)]

    data[, `:=`(regional_monthly_tmin_c = mean_na(site_monthly_tmin_c),
                regional_monthly_tmax_c = mean_na(site_monthly_tmax_c)),
         by = .(sample_year, sample_month)]
    
    data[, `:=`(regional_seasonal_tmin_c = mean_na(regional_monthly_tmin_c),
                regional_seasonal_tmax_c = mean_na(regional_monthly_tmax_c)),
         by = .(sample_year, sample_season)]
    
    # Precip
    # sum precip by site, year, month
    data[, site_total_monthly_precip_cm := sum_na(precip_cm),
         by = .(station_name, sample_year, sample_month)]

    # mean of summed precip by year, month
    data[, regional_mean_monthly_precip_cm := mean_na(site_total_monthly_precip_cm),
         by = .(sample_year, sample_month)]
    
    # sum of mean_monthly_precip by year, season
    # Needs to be divided by .N because regional_mean_monthly_precip_cm is repeated across sites
    data[, regional_total_seasonal_precip_cm := sum_na(regional_mean_monthly_precip_cm) / .N,
         by = .(sample_year, sample_season)]
    
   
    dcast(data[, 
               .(station_name,
                 lat,
                 lon,
                 elevation_m,
                 sample_date, 
                 sample_month,
                 sample_season, 
                 tmin_c,
                 tmax_c,
                 tmin_l1,
                 tmax_l1,
                 harmonic_cos, 
                 harmonic_sin, 
                 precip_occur,
                 precip_occur_l1,
                 precip_amount,
                 temperature_trend,
                 regional_precip_cm = regional_total_seasonal_precip_cm,
                 regional_tmin_c = regional_seasonal_tmin_c,
                 regional_tmax_c = regional_seasonal_tmax_c)], 
          ... ~ sample_season,
          value.var = c("regional_tmin_c", "regional_tmax_c", "regional_precip_cm"), 
          fill = 0)

  }
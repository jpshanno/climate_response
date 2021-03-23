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
prepare_ghcnd_swg_data <- 
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
    
    # Use IDW to fill-in missing data (do it prior to removing unused stations
    # because they still have useful data for the IDW backfilling)
    
    station_weights <- 
      calculate_weights(ghcnd_data)
    
    idw_data <- 
      map_dfr(colnames(station_weights), 
              ~ghcnd_data[station_name != .x, 
                          c(station_name = .x,
                            map(.SD, 
                                function(X){
                                  weighted.mean(x = X, 
                                                w = station_weights[.x, station_name],
                                                na.rm = TRUE)})), 
                          by = .(sample_date), 
                          .SDcols = c("tmin_c", "tmax_c", "precip_cm")])
    
    ghcnd_data <- 
      ghcnd_data[stations_to_use,
                 on = "station_name"]
  
    ghcnd_data <- 
      ghcnd_data[idw_data,
                 on = c("station_name", "sample_date"),
                 `:=`(tmin_c = fcoalesce(tmin_c, i.tmin_c),
                      tmax_c = fcoalesce(tmax_c, i.tmax_c),
                      precip_cm = fcoalesce(precip_cm, i.precip_cm))]
    
    swg_format_data(ghcnd_data)
      
  }


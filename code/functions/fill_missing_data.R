##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author Joe Shannon
##' @export
fill_missing_data <- 
  function(data,
           additional.data = NULL,
           source) {

    stopifnot(length(source) == 1)
    
    out <- 
      switch(source,
             "ncei" = fill_ncei_idw(data,
                                    additional.data = additional.data))
    
    out
    
}

fill_ncei_idw <- 
  function(data,
           additional.data){
    
    data <- copy(data)
    
    # Caclulate Distances and Weights
    
    station_weights <- 
      calculate_weights(rbind(data, additional.data, fill = TRUE))
    
    # Fill Missing Precip Data ------------------------------------------------
    # Tested with GLM -> predicted precip -> weighted mean of predicted precip and
    # saw slight decline in model metrics
    # Create nested datatable to performed IDW for each station using all other 
    # stations. This requires nesting all data except for a given station
    idw <- 
      data[, 
           .(daily_dat = list(set(x = copy(data[station_name != .BY[[1]]]),
                                  j = "site",
                                  value = .BY[[1]]))), 
           by = .(station_name)]
    
    idw[, daily_dat := 
          lapply(daily_dat,
                 function(x){
                   x$weights <- 
                     station_weights[, first(x$site)][x$station_name]
                   x
                 })]
    
    idw[, daily_dat := 
          lapply(daily_dat, 
                 function(x){
                   x[, .(idw_precip_cm = weighted.mean(precip_cm, weights, na.rm = TRUE)), 
                     by = .(sample_date)]})]
    
    data[idw[, daily_dat[[1]], by = .(station_name)],
         idw_precip_cm := i.idw_precip_cm,
         on = c("station_name", "sample_date")]
   
    # Fill Missing Temperatures -----------------------------------------------
    
    temperature_dat <- 
      rbind(data[, .(station_name, sample_date = as.Date(sample_date), lat, lon, tmin_c, tmax_c)],
            additional.data[, .(station_name, sample_date, lat, lon, tmin_c, tmax_c)])
    
    idw_temp <- 
      data[, 
           .(daily_dat = list(set(x = copy(temperature_dat[station_name != .BY[[1]]]),
                                  j = "site",
                                  value = .BY[[1]]))), 
           by = .(station_name)]

    idw_temp[, daily_dat := lapply(daily_dat, 
                                   function(x){
                                     x$weights <- 
                                       station_weights[, first(x$site)][x$station_name]
                                     x})]
    
    idw_temp[, daily_dat := lapply(daily_dat, 
                                   function(x){x[, lapply(.SD, 
                                                          weighted.mean, 
                                                          w = weights,
                                                          na.rm = TRUE), 
                                                 by = .(sample_date),
                                                 .SDcols = patterns("t(min|max)_c")]})]
    
    setnames(data,
             c("precip_cm", "tmin_c", "tmax_c"),
             c("obs_precip_cm", "obs_tmin_c", "obs_tmax_c"))
    
    data[idw_temp[, daily_dat[[1]], by = .(station_name)],
         `:=`(idw_tmin_c = i.tmin_c,
              idw_tmax_c = i.tmax_c),
         on = c("station_name", "sample_date")]
    
    data[, `:=`(precip_cm = fcoalesce(obs_precip_cm, idw_precip_cm),
                tmin_c = fcoalesce(obs_tmin_c, idw_tmin_c),
                tmax_c = fcoalesce(obs_tmax_c, idw_tmax_c))]
    
    # Set logicals if final data is filled or original
    data[, .(station_name,
             station_id,
             lat,
             lon,
             sample_date,
             sample_year,
             water_year,
             dowy,
             precip_cm = fcoalesce(obs_precip_cm, idw_precip_cm),
             tmin_c = fcoalesce(obs_tmin_c, idw_tmin_c),
             tmax_c = fcoalesce(obs_tmax_c, idw_tmax_c),
             filled_precip = is.na(obs_precip_cm),
             filled_tmin = is.na(obs_tmin_c),
             filled_tmax = is.na(obs_tmax_c))]
    
  }

calculate_weights <- 
  function(data){
    
    data <- 
      unique(data[, .(station_name, lon, lat)])
    
    station_distances <- 
      geodist(data,
              measure = "geodesic") / 1000
    
    # Set row and column names to station names for easy lookup
    dimnames(station_distances) <- 
      list(data$station_name,
           data$station_name)
    
    # Cacluate IDW weights
    ifelse(station_distances == 0,
           0,
           1 / station_distances^3)
  }
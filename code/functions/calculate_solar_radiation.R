##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author Joe Shannon
##' @export
calculate_solar_radiation <- 
  function(data,
           coefs) {

  # lat in bc() can only be a single value, so this has to be done by each 
  # station name
    
  data[, 
       solrad_MJ_m2 := bc(days = sample_date, 
                          lat = first(lat),
                          BCb = coefs[, mean(BC_B)],
                          Tmax = tmax_c,
                          Tmin = tmin_c,
                          tal = coefs[, mean(clear_sky_t)]),
       by = .(station_name)][]
    
  }

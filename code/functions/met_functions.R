##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param lambda_MJ_kg
##' @return
##' @author Joe Shannon
##' @export
calculate_hargreaves_pet <- 
  function(data, lambda.MJ.kg) {
    
    data[, pet_cm := 0.1 * 0.0023 * solrad_MJ_m2 / lambda.MJ.kg * (tmean_c + 17.8) * sqrt(tmax_c - tmin_c)][]
    
  }

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

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author Joe Shannon
##' @export
calculate_mean_temp <-
  function(data) {
    
    data[, tmean_c := (tmin_c - tmax_c)/2][]
    
  }
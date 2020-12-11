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

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author Joe Shannon
##' @export
calculate_water_availability <- 
  function(data) {
  
  data[, water_availability_cm := cumsum(precip_cm - pet_cm),
       by = .(station_name,
              water_year)]
  
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author Joe Shannon
##' @export
calculate_solrad_coefs <- 
  function(data) {
    
    stopifnot(all(c("sample_date", "lat", "solrad_MJ_m2", "tmax_c", "tmin_c") %in% names(data)))
    
    data <- 
      copy(data)
    
    data <- 
      data[,
           .(dat = list(.SD)),
           by = .(station_name)]
    
    data[, 
         clear_sky_t := map_dbl(dat,
                                ~cst(RefRad = .x$solrad_MJ_m2,
                                     days = .x$sample_date,
                                     lat = radians(first(.x$lat)),
                                     extraT = .x$et_solrad_MJ_m2))]
    
    data[,
         mod := 
           map2(dat, clear_sky_t,
                ~nls(solrad_MJ_m2 ~ .y * (1-exp(-B*(tmax_c - tmin_c)**C)) * et_solrad_MJ_m2,
                     start = list(B = 0.06, C = 2),
                     data = .x))]
    
    data[, c("BC_B", "BC_C") := map_dfr(mod, coef)]
    
    data[, .(station_name, 
             BC_A = clear_sky_t, 
             BC_B,
             BC_C)]

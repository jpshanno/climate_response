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
    # 0.014375 = 0.0023/0.16, taken from Hargreaves & Allen 2003 (eq. 8 & Ks =
    # 0.16 from precedeing paragraph)
    # Eq 3 from Hargreaves & Allen (2003), using Rs as calculated from 
    # Bristow-Campbell.
    data[, pet_cm := 0.1 * 0.0135 * solrad_MJ_m2 / lambda.MJ.kg * (tmean_c + 17.8)][]
    
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
    
    # Calculate daily external solar radiation by DOY & lat to avoid duplicate 
    # calculations in met data
    solrad <- 
      data[, 
           .(doy = 1:366,
             et_solrad_MJ_m2 = extrat(1:366, radians(first(lat)))$ExtraTerrestrialSolarRadiationDaily), 
           by = .(station_name)]
    
    has_doy <- 
      "doy" %in% names(data)
    
    if(!has_doy){
      data[, doy := yday(sample_date)]  
    }
    
    
    data[solrad,
         et_solrad_MJ_m2 := i.et_solrad_MJ_m2,
         on = c("station_name", "doy")]
    
    # lat in bc() can only be a single value, so this has to be done by each 
    # station name
    data[, 
         solrad_MJ_m2 := sirad::bc(days = sample_date, 
                            lat = first(lat),
                            BCb = coefs[, mean(BC_B)],
                            BCc = coefs[, mean(BC_C)],
                            Tmax = tmax_c,
                            Tmin = tmin_c,
                            extraT = et_solrad_MJ_m2,
                            tal = coefs[, mean(BC_A)]),
         by = .(station_name)]
    
    # Check that only missing values are from days where bc_tmin is greater than
    # tmax, this seems to happen for winter days that are significantly warmer
    # than the following day, likely a result of a relaxation of cold air inputs
    # from the north for a day. Without removing these bc() returns NaN
    if(any(is.na(data$solrad_MJ_m2))){
      # Set tmin_c to average of tmin_c(t, t+1) according to Bristow-Campbell
      data[, bc_tmin := (tmin_c + shift(tmin_c, -1)) / 2,
           by = .(station_name)]
      
      # Set last value in each series to last value (mean of last two observations)
      data[, bc_tmin := nafill(bc_tmin, "locf")]
      
      # If all errors are not explained by bc_tmin > tmax_c then throw an error
      if(!data[is.na(solrad_MJ_m2), all(bc_tmin > tmax_c)]){
        stop("There are unexplained NA values from solar radiation calculation 
             using bc().")
      }
      
      # Calculate solar radiation for errors by setting tdiff = 0, which sets
      # solrad to 0. Doing it this way rather than directly setting solard := 0
      # in case I want to change my approach in the future. Tried using tdiff = 
      # tmax - tmin, but that ended up with sigificantly higher solrad and PET 
      # than any of the other stations for the error days (2015-02-07)
      
      data[,
           solrad_MJ_m2 := ifelse(!is.na(solrad_MJ_m2),
                                  solrad_MJ_m2,
                                  coefs[, mean(BC_A)] * (1-exp((-coefs[, mean(BC_B)]*(0)**coefs[, mean(BC_C)]) / mean((tmax_c - tmin_c)))) * et_solrad_MJ_m2),
           by = .(station_name, month(sample_date))]

      data[, bc_tmin := NULL]
    }

    if(!has_doy){
      data[, doy := NULL]
    }
    
    data[, et_solrad_MJ_m2 := NULL]
    
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
    
    data[, tmean_c := (tmin_c + tmax_c)/2][]
    
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
    
    # Set tmin_c to average of tmin_c(t, t+1) according to Britow-Campbell
    data[, bc_tmin := (tmin_c + shift(tmin_c, -1)) / 2,
         by = .(station_name)]
    
    # Set last value in each series to last value (mean of last two observations)
    data[, bc_tmin := nafill(bc_tmin, "locf")]
    
    # Drop days where the new bc_tmin is greater than tmax, this seems to happen
    # for winter days that are significantly warmer than the following day,
    # likely a result of a relaxation of cold air inputs from the north for a
    # day
    
    data <- 
      data[bc_tmin < tmax_c]
    
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
                ~nls(solrad_MJ_m2 ~ .y * (1-exp((-B*tdiff**C) / mnth_tdiff_c)) * et_solrad_MJ_m2,
                     start = list(B = 0.06, C = 2),
                     data = .x[, .(solrad_MJ_m2, tmax_c, bc_tmin, et_solrad_MJ_m2, 
                                   tdiff = tmax_c - bc_tmin, mnth_tdiff_c = mean(tmax_c - bc_tmin)),
                               by = .(month(sample_date))]))]
    
    data[, c("BC_B", "BC_C") := map_dfr(mod, coef)]
    
    data[, .(station_name, 
             BC_A = clear_sky_t, 
             BC_B,
             BC_C)]
  }

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @param data
##' @param precip.col
##' @param pet.col
##' @return
##' @author Joe Shannon
##' @export
calculate_water_availability <- 
  function(data){
    
    data[, 
         water_availability_cm := cumsum(total_input_cm - pet_cm),
         by = .(station_name)]
    
    data[, ytd_water_availability_cm := (water_availability_cm - first(water_availability_cm)) - max(water_availability_cm - first(water_availability_cm)),
         by = .(station_name, sample_year)]

    data    
  }
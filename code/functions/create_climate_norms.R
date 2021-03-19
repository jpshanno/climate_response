##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param norms.file,
##' @param solrad.coefs
##' @return
##' @author Joe Shannon
##' @export
##' 
prepare_climate_norms <- 
  function(norms.file,
           ncei.path,
           ghcnd.units,
           solrad.coefs) {

    
    norms <- 
      fread(norms.file, 
            na.strings = c("", "-7777"))
    
    norms <- 
      norms[,
            .(station_name = "bergland_dam",
              sample_date =as.Date(paste0("2020-", DATE)),
              doy = yday(as.Date(paste0("2020-", DATE))),
              lat = LATITUDE,
              lon = LONGITUDE,
              tmax_c = (0.1 * `DLY-TMAX-NORMAL` - 32) / 1.8,
              tmean_c = (0.1 * `DLY-TAVG-NORMAL` - 32) / 1.8,
              tmin_c = (0.1 * `DLY-TMIN-NORMAL` - 32) / 1.8,
              precip_cm = `DLY-PRCP-50PCTL` / 39.37)]
    
    ex_solrad <- 
      unique(norms[, .(station_name, lat, doy)])
    
    ex_solrad[, et_solrad := extrat(doy, first(radians(lat)))$ExtraTerrestrialSolarRadiationDaily,
              by = .(station_name)]
    
    norms[ex_solrad, 
          et_solrad := i.et_solrad,
          on = .(station_name)]
    
    norms[, solrad_MJ_m2 := sirad::bc(days = sample_date, 
                                      lat = first(lat),
                                      BCb = solrad.coefs$BC_B,
                                      extraT = et_solrad,
                                      Tmax = tmax_c,
                                      Tmin = tmin_c,
                                      BCc = solrad.coefs$BC_C,
                                      tal = solrad.coefs$BC_A),
          by = .(station_name)]
    
    calculate_hargreaves_pet(norms, 2.45)
    
    norms[, dowy := as.dowy(sample_date, 11)]
    
    norms[, sample_date := NULL]
    
    norms
  }
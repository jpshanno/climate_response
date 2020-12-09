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
    
    data <- 
      copy(data)
    
    data <- 
      data[,
           .(dat = list(.SD)),
           by = .(station_name)]
    
    data[, clear_sky_t := map_dbl(dat,
                                  ~cst(RefRad = .x$solrad_MJ_m2,
                                       days = .x$sample_date,
                                       lat = first(.x$lat)))]
    
    data[, BC_B := map2_dbl(dat, clear_sky_t,
                           ~bccal(lat = first(.x$lat), 
                                  days = .x$sample_date,
                                  rad_mea = .x$solrad_MJ_m2,
                                  Tmax = .x$tmax_c,
                                  Tmin = .x$tmin_c,
                                  tal = .y))]
    
    data[, .(station_name, clear_sky_t, BC_B)]
}

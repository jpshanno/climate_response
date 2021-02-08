##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param data
##' @param external.met
##' @return
##' @author Joe Shannon
##' @export
merge_external_met <- function(data, external.met) {

  water_budget <- 
    copy(data)
  
  wghts <- 
    calculate_weights(rbind(unique(data[, .(station_name = site, lat, lon)]), 
                            unique(external.met[, .(station_name, lat, lon)])))
  
  setkey(external.met, "sample_date")
  
  water_budget[, c("pet_cm", 
                   "tmin_c",
                   "tmax_c",
                   "water_availability_cm",
                   "ytd_water_availability_cm") := 
                 map2_dfr(site, sample_date, 
                          function(site, ex_date){
                            
                            ex_dat <- 
                              external.met[J(ex_date)]
                            
                            ex_dat[, weight := wghts[site, station_name]]
                            
                            ex_dat[, .(pet_cm = weighted.mean(pet_cm, weight),
                                       tmin_c = weighted.mean(tmin_c, weight),
                                       tmax_c = weighted.mean(tmax_c, weight),
                                       water_availability_cm = weighted.mean(water_availability_cm, weight),
                                       ytd_water_availability_cm = weighted.mean(ytd_water_availability_cm, weight))]
                          })
  ]
  
  water_budget
  
}

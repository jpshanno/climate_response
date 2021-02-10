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
  
  external.met <- 
    copy(external.met)
  
  wghts <- 
    calculate_weights(rbind(unique(data[, .(station_name = site, lat, lon)]), 
                            unique(external.met[, .(station_name, lat, lon)])))
  
  setkey(external.met, "sample_date")
  
  water_budget[, c("pet_cm", 
                   "tmin_c",
                   "tmax_c",
                   "melt_cm",
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
                                       melt_cm = weighted.mean(melt_cm, weight),
                                       water_availability_cm = weighted.mean(water_availability_cm, weight),
                                       ytd_water_availability_cm = weighted.mean(ytd_water_availability_cm, weight))]
                          })
  ]
  
  external.met[, `:=`(l1_pet_cm = shift(pet_cm, 1),
                      l1_melt_cm = shift(melt_cm, 1),
                      l1_water_availability_cm = shift(water_availability_cm, 1),
                      l1_ytd_water_availability_cm = shift(ytd_water_availability_cm, 1)),
               by = .(station_name)]
  
  water_budget[, c("l1_pet_cm", 
                   "l1_melt_cm",
                   "l1_water_availability_cm",
                   "l1_ytd_water_availability_cm") := 
                 map2_dfr(site, sample_date, 
                          function(site, ex_date){
                            
                            ex_dat <- 
                              external.met[J(ex_date)]
                            
                            ex_dat[, weight := wghts[site, station_name]]
                            
                            ex_dat[, .(l1_pet_cm = weighted.mean(l1_pet_cm, weight, na.rm = TRUE),
                                       l1_melt_cm = weighted.mean(l1_melt_cm, weight, na.rm = TRUE),
                                       l1_water_availability_cm = weighted.mean(l1_water_availability_cm, weight, na.rm = TRUE),
                                       l1_ytd_water_availability_cm = weighted.mean(l1_ytd_water_availability_cm, weight, na.rm = TRUE))]
                          })
  ]
  
  water_budget
  
}

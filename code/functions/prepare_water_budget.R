##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param data
##' @param external.met
##' @param sy.mods
##' @return
##' @author Joe Shannon
##' @export
prepare_water_budget <- 
  function(data, 
           external.met, 
           interception.mods,
           sy.mods) {

    data <- 
      copy(data)
    
    data[, sy := sy.mods[.BY[[1]], f_predict[[1]]](wl_min_cm),
         by = .(site)]
    
    data[, net_precip_cm := best_precip_cm - interception.mods[.BY[[1]], f_predict[[1]]](best_precip_cm),
         by = .(site_status)]
    
    water_budget <- 
      data[, 
           .(site, 
             lat,
             lon,
             doy,
             dowy,
             sample_date, 
             sample_year,
             water_year,
             season,
             site_status,
             sy, 
             wl_initial_cm,
             wl_median_cm,
             obs_precip_cm = obs_precip_cm,
             net_precip_cm = net_precip_cm,
             best_precip_cm = best_precip_cm,
             melt_cm = melt_cm,
             Dl_signed_cm = Ds_cm * sy,
             Ds_cm = Ds_cm * sy)]
    
    wghts <- 
      calculate_weights(rbind(unique(data[, .(station_name = site, lat, lon)]), 
                              unique(external.met[, .(station_name, lat, lon)])))
    
    setkey(external.met, "sample_date")
    
    water_budget[, c("pet_cm", 
                     "tmin_c",
                     "tmax_c",
                     "water_availability_cm") := 
                   map2_dfr(site, sample_date, 
                            function(site, ex_date){
                              
                              ex_dat <- 
                                external.met[J(ex_date)]
                              
                              ex_dat[, weight := wghts[site, station_name]]
                              
                              ex_dat[, .(pet_cm = weighted.mean(pet_cm, weight),
                                         tmin_c = weighted.mean(tmin_c, weight),
                                         tmax_c = weighted.mean(tmax_c, weight),
                                         water_availability_cm = weighted.mean(water_availability_cm, weight))]
                            })
    ]

    water_budget[, `:=`(pet_cm = pet_cm)]
    
    water_budget[, net_flow_cm := Ds_cm + pet_cm - net_precip_cm - melt_cm]
    
    water_budget
    
  }

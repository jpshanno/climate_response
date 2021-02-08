##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author Joe Shannon
##' @export
prepare_water_levels <- 
  function(path) {

  water_levels <- 
    fread(path,
          select = c(site = "character",
                     sample_time = "character",
                     water_level_cm = "numeric"))
  
  water_levels[, `:=`(sample_time = ymd_hms(sample_time, tz = "EST"),
                      sample_date = as.Date(sample_time, tz = "EST"),
                      sample_year = year(sample_time))]
  
  setkey(water_levels, "site", "sample_date")
  
  daily_water_levels <- 
    water_levels[,.(doy = yday(.BY[[2]]),
                    doy_decimal = yday(.BY[[2]]) / 366,
                    sample_year = first(sample_year),
                    wl_initial_cm = first(water_level_cm),
                    wl_final_cm = last(water_level_cm),
                    wl_min_cm = min_na(water_level_cm),
                    wl_median_cm = median_na(water_level_cm),
                    wl_max_cm = max_na(water_level_cm),
                    Dl_cm = max_na(water_level_cm) - min_na(water_level_cm),
                    Dl_signed_cm = sign(which.max(water_level_cm) - which.min(water_level_cm)) * (max_na(water_level_cm) - min_na(water_level_cm)),
                    Dl_time_range_s = abs((which.max(water_level_cm) - which.min(water_level_cm)) * 900),
                    Ds_cm = last(water_level_cm) - first(water_level_cm)),
                 keyby = .(site, sample_date)]
  
  site_info <- 
    fread(text = 
            "site,lon,lat,study,treatment,treatment_date,morphology
             006,-89.19036,46.32392,planting,Girdle,2014-03-31,flowthrough
             009,-89.71207,46.41624,eco,Ash Cut,2014-03-31,threshold
             053,-89.53043,46.61028,pws,Ash Cut,2015-03-31,flowthrough
             077,-88.86684,46.71778,eco,Ash Cut,2014-03-31,threshold
             111,-89.71369,46.57636,planting,Control,2014-03-31,threshold
             113,-89.80990,46.35090,pws,Control,2015-03-31,flowthrough
             119,-89.59433,46.28794,eco,Girdle,2014-03-31,closed
             135,-88.87838,46.75382,eco,Control,2014-03-31,flowthrough
             139,-88.95577,46.59619,planting,Ash Cut,2014-03-31,threshold
             140,-89.61962,46.43279,eco,Girdle,2014-03-31,flowthrough
             151,-88.89502,46.67307,eco,Girdle,2014-03-31,threshold
             152,-89.71468,46.41277,eco,Control,2014-03-31,flowthrough
             156,-89.58431,46.28327,eco,Ash Cut,2014-03-31,closed
             157,-89.58479,46.28173,eco,Control,2014-03-31,threshold",
          select = c(site = "character",
                     lon = "numeric",
                     lat = "numeric",
                     study = "character",
                     treatment = "character",
                     treatment_date = "Date",
                     morphology = "character"),
          key = "site")
  
  
  daily_water_levels[site_info, 
                     `:=`(lon = i.lon,
                          lat = i.lat,
                          treatment = i.treatment,
                          site_status = ifelse(sample_date > i.treatment_date & i.treatment != "Control",
                                               "Treated",
                                               "Control"),
                          morphology = i.morphology,
                          season = fifelse(month(sample_date) %in% 6:7 |
                                             (month(sample_date) == 5 & mday(sample_date) > 15) |
                                             (month(sample_date) == 9 & mday(sample_date) <= 15),
                                           "growing",
                                           "dormant"))]
  
  daily_water_levels[, `:=`(water_year = as.water_year(sample_date, 11),
                            dowy = as.dowy(sample_date, 11))]
  
  daily_water_levels
  
  }

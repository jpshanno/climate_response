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
    water_levels[,.(sample_year = first(sample_year),
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
            "site,study,treatment,treatment_date,morphology
             006,planting,Girdle,2014-03-31,flowthrough
             009,eco,Ash Cut,2014-03-31,threshold
             053,pws,Ash Cut,2015-03-31,flowthrough
             077,eco,Ash Cut,2014-03-31,threshold
             111,planting,Control,2014-03-31,threshold
             113,pws,Control,2015-03-31,flowthrough
             119,eco,Girdle,2014-03-31,closed
             135,eco,Control,2014-03-31,flowthrough
             139,planting,Ash Cut,2014-03-31,threshold
             140,eco,Girdle,2014-03-31,flowthrough
             151,eco,Girdle,2014-03-31,threshold
             152,eco,Control,2014-03-31,flowthrough
             156,eco,Ash Cut,2014-03-31,closed
             157,eco,Control,2014-03-31,threshold",
          select = c(site = "character",
                     study = "character",
                     treatment = "character",
                     treatment_date = "Date",
                     morphology = "character"),
          key = "site")
  
  
  daily_water_levels[site_info, 
                     `:=`(treatment = i.treatment,
                          site_status = ifelse(sample_date > i.treatment_date,
                                               "Treated",
                                               "Control"),
                          morphology = i.morphology,
                          season = fifelse(month(sample_date) %in% 6:7 |
                                             (month(sample_date) == 5 & mday(sample_date) > 15) |
                                             (month(sample_date) == 9 & mday(sample_date) <= 15),
                                           "growing",
                                           "dormant"))]
  
  daily_water_levels
  
  }

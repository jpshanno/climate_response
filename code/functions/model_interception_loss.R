##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param data
##' @return
##' @author Joe Shannon
##' @export
model_interception_loss <- 
  function(data,
           precip.column) {

    interception <- 
      data[best_precip_cm > 0,
                          .(treatment = first(treatment),
                            mod = list(lmrob(Ds_cm ~ best_precip_cm,
                                             data = .SD,
                                             setting = "KS2014"))),
                          by = .(site, site_status, season)]
    
    interception[, c("intercept", "slope") := map_dfr(mod, coef)]
    interception[, i_cm := -intercept / slope]
    
    interception[, i_cm := ifelse(is.na(i_cm), mean_na(i_cm), i_cm),
                 by = .(site_status)]
    
    interception[i_cm < 0, i_cm := 0]
    
    interception
    
    # daily_water_balance[interception, 
    #                     `:=`(i_cm = i.i_cm,
    #                          net_precip_cm = pmax(0, best_precip_cm - i.i_cm)),
    #                     on = c("site", "site_status", "season")]
    # 
    # daily_water_balance[, precip_intensity_cm_hr := net_precip_cm / (Dl_time_range_s / 3600)]

}

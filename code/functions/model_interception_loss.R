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
           y,
           x,
           g) {

    form <- 
      make_formula(y, x)
    
    interception <- 
      data[data[[x[1]]] > 0,
           .(treatment = first(treatment),
             mod = list(lmrob(form,
                              data = .SD,
                              setting = "KS2014"))),
           keyby = g]
    
    interception[, c("intercept", "slope", "interaction_slope") := map_dfr(mod, coef)]

    # Doing this using split inside := because doing it by group with env = .SD
    # is only using the last row of values

    prediction_functions <- 
      split(interception, 
            by = "site_status") %>% 
      map(~as.function(list(doy = NULL, 
                            substitute(-intercept / (slope + interaction_slope * cos(2*pi*((doy / 366)))), 
                                       env = .x))))
        
    interception[, f_predict := prediction_functions]
    
    # Add Observed Precip Range
    interception[data[, .(precip_range = list(range(best_precip_cm, na.rm = TRUE))),
                      by = .(site_status)],
                 precip_range_cm := i.precip_range,
                 on = c("site_status")]
    
    interception
    
    # daily_water_balance[interception, 
    #                     `:=`(i_cm = i.i_cm,
    #                          net_precip_cm = pmax(0, best_precip_cm - i.i_cm)),
    #                     on = c("site", "site_status", "season")]
    # 
    # daily_water_balance[, precip_intensity_cm_hr := net_precip_cm / (Dl_time_range_s / 3600)]

}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @return
##' @author Joe Shannon
##' @export
calculate_net_precip <- 
  function(data, intercept.models, precip.col) {

  stopifnot("doy_decimal" %in% names(data))
    
  data[, i_cm := intercept.models[.BY[[1]], f_predict[[1]]](doy_decimal),
       by = .(site_status)]

  set(data, 
      j = "net_precip_cm",
      value = pmax(0, data[[precip.col]] - data$i_cm))
  
  data[]
}

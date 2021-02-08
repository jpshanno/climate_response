##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author Joe Shannon
##' @export
model_precip_rise <- 
  function(data) {

  # Precip models use net precip, so treatment period is not a necessary factor
    
  mods <- 
    data[net_precip_cm > 0 & Dl_signed_cm > 0, 
         .(mod = list(lmrob(Dl_signed_cm ~ 0 + net_precip_cm,
                            setting = "KS2014"))),
         keyby = .(site)]
  
  mods[, 
       c("intercept", "slope") := map_dfr(mod, coef), 
       by = .(site)]
  
  prediction_functions <- 
    split(mods, 
          by = "site") %>% 
    map(~as.function(list(net.precip = NULL, 
                          substitute({ifelse(net.precip == 0,
                                             0,
                                             slope * net.precip)}, 
                                     env = .x))))
  
  mods[, 
       f_predict := prediction_functions]
  
}

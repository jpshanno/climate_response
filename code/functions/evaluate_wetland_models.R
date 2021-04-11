##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param data
##' @param model.params
##' @return
##' @author Joe Shannon
##' @export
evaluate_wetland_models <- 
  function(data) {

  data[!is.na(wl_hat + wl_initial_cm), 
       cbind(site_status = site_status[1],
             as.data.table(t(hydroGOF::gof(wl_hat, wl_initial_cm)))),
       by = .(site, water_year)][]

}

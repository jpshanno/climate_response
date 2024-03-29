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
predict_test_period <- function(data, model.params) {

  data[, 
       cbind(sample_date,
             wl_initial_cm,
             site_status = site_status[1],
             wetland_model(data = .SD,
                           params = model.params[CJ(.BY[[1]], .SD$site_status[1]), params[[1]]])),
       by = .(site, water_year)][]

}


predict_test_period_future_forest <- function(data, model.params) {
  
  data[, 
       cbind(sample_date,
             wl_initial_cm,
             site_status = "Future Forested",
             wetland_model(data = .SD,
                           params = model.params[CJ(.BY[[1]], "Future Forested"), params[[1]]])),
       by = .(site, water_year)][]
  
}

predict_train_period <- function(data, model.params) {
  
  lapply(data,
    function(x) {
      x[, 
           cbind(sample_date,
                 wl_initial_cm,
                 site_status = site_status[1],
                 wetland_model(data = .SD,
                               params = model.params[CJ(.BY[[1]], .SD$site_status[1]), params[[1]]])),
           by = .(site, water_year)][]
    }
  )
  
  
}

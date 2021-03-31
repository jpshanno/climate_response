##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param control_optimization
##' @param treated_optimization
##' @return
##' @author Joe Shannon
##' @export
concatenate_model_params <- 
  function(...) {
    
    params <- 
      list(...)
    
    stopifnot(all(names(params) %in% c("Control", "Treated")))
    
    params[["Future Forested"]] <- 
      params[["Control"]][fread("data/overstory_ash_percentage.csv", 
                                colClasses = c("character", "numeric", "character")), 
                          .(site, params, future.forest.change = i.percent.basalarea.as.ash), 
                          on = "site"]
    params <- 
      rbindlist(params,
                idcol = "site_status",
                fill = TRUE)
    
    params <- 
      params[, .(site, site_status, params, future.forest.change)]
    
    setkey(params, site, site_status)
    
    params
    
  }

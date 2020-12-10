##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author Joe Shannon
##' @export
append_study_precip <- 
  function(data,
           precip.path,
           snow.path) {

    precip <- 
      fread(precip.path,
            select = c(site = "character",
                       sample_date = "Date",
                       obs_precip_cm = "numeric",
                       idw_precip_cm = "numeric"),
            key = c("site", "sample_date"))
    
    melt <- 
      fread(snow.path, 
            select = c(site = "character", 
                       sample_date = "Date", 
                       melt_cm = "numeric"),
            key = c("site", "sample_date"))
    
    # Need melt for all sites. For not just using max predicted melt
    melt <- 
      melt[, .(melt_cm = mean_na(melt_cm)), by = .(sample_date)]
    
    data[precip,
         `:=`(obs_precip_cm = i.obs_precip_cm,
              idw_precip_cm = i.idw_precip_cm)]
    
    data[melt,
        melt := i.melt_cm,
        on = "sample_date"]
    
    data[, best_precip_cm := fcoalesce(obs_precip_cm, idw_precip_cm)]
    
    data
}

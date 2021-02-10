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
           precip.path) {

    precip <- 
      fread(precip.path,
            select = c(site = "character",
                       sample_date = "Date",
                       obs_precip_cm = "numeric",
                       idw_precip_cm = "numeric"),
            key = c("site", "sample_date"))
    
    # Add Lagged Data  
    # Adding lagged data here because daily water levels only covers periods 
    # with water level data, which means the start of each season of data would
    # have an NA if the lag was performed in that dataset
    precip[, 
           `:=`(l1_obs_precip_cm = shift(obs_precip_cm, 1),
                l1_idw_precip_cm = shift(idw_precip_cm, 1)),
           by = .(site)]
    
    data[precip,
         `:=`(obs_precip_cm = i.obs_precip_cm,
              idw_precip_cm = i.idw_precip_cm,
              l1_obs_precip_cm = i.l1_obs_precip_cm,
              l1_idw_precip_cm = i.l1_idw_precip_cm)]
    
    data[, best_precip_cm := fcoalesce(obs_precip_cm, idw_precip_cm)]
    data[, l1_best_precip_cm := fcoalesce(l1_obs_precip_cm, l1_idw_precip_cm)]
    
    data
}

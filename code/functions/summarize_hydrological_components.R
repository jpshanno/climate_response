summarize_hydrological_components <- function(files, variable, dummy) {
  
  dat <- vector("list", length = length(files))
  
  for(.x in files) {
      dat[[.x]] <- vroom::vroom(.x, n_max)
      
      data.table::setDT(dat[[.x]])

      dat[[.x]] <- dat[[.x]][
        j = data.table::as.data.table(ggdist::median_hdci(.SD[[variable]], .width = 0.67, na.rm = TRUE)),
        by = .(site_status, gcm, scenario, simulation_date)
      ]
      
      gc()
      
      cat(.x, "complete \n")
      
  }
  
  rbindlist(dat)
  
}


summarize_pet <- function(files, dummy) {
  
  dat <- vector("list", length = length(files))
  
  for(.x in seq_along(files)) {
    dat[[.x]] <- vroom::vroom(
      file = files[.x],
      col_select = c(site, site_status, gcm, scenario, simulation_date, pet_hat, gradient))
    
    data.table::setDT(dat[[.x]])
    
    dat[[.x]] <- dat[[.x]][
      j = .(
        aet_hat = median(filled_rolling_mean(fifelse(gradient == 0, 0, abs(pet_hat) / gradient)), na.rm= TRUE),
        aet_effect = median(filled_rolling_mean(fifelse(gradient == 0, 0, abs(pet_hat))), na.rm = TRUE)
      ),
      by = .(site, site_status, gcm, scenario, simulation_date)
    ]
    
    dat[[.x]] <- dat[[.x]][
      j = ggdist::median_hdci(.SD, aet_hat, aet_effect, .width = 0.67, na.rm = TRUE),
      by = .(site_status, gcm, scenario, simulation_date)
    ]
    
    gc()
    
    cat(files[.x], "complete \n")
    
  }
  
  rbindlist(dat)
  
}

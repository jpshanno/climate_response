summarize_hydrological_components <- function(files, variable, dummy) {
  
  dat <- vector("list", length = length(files))
  
  for(.x in files) {
      dat[[.x]] <- vroom::vroom(.x, n_max = 365*10000*2)
      
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

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param data
##' @return
##' @author Joe Shannon
##' @export
select_training_data <- 
  function(data) {
    
    # Prefer 2012 becuase most sites have data, it's the control period, and it
    # had a strong drawdown, which means good ESy training data. But less good
    # melt data
    training_years <- 
      data[!is.na(wl_initial_cm) & site_status == "Control", 
           .(training_year = min(water_year)),
           by = .(site)]
    
    # When 2012 data doesn't exist we should choose the year with the greatest 
    # range of water level variation (excluding melt)
    training_years[, 
                  training_year := 
                    fifelse(training_year == 2012, 
                            2012,
                            data[!is.na(wl_initial_cm) & site == .BY[[1]] & site_status == "Control" & month(sample_date) %in% 7:11, 
                                         .(wl_range = diff(range(wl_initial_cm))), 
                                         by = .(water_year)][wl_range == max(wl_range), water_year]), 
      by = .(site)]
    
    data[training_years,
         on = .(site, water_year = training_year)]
    
  }

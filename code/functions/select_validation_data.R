##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param data
##' @return
##' @author Joe Shannon
##' @export
select_validation_years <- 
  function(data, 
           date.col,
           n.in.group,
           groups) {
    
    set.seed(1234)
    
    data[, .sample_year := year(get(date.col))]
    
    # If number of length(unique(years)) = 1 sample does not work as expected!
    # Simple work around is rep each value twice (since I'm already doing unique)
    # this does not affect the probability of any given year being selected, it 
    # goest from 1 of each year value to two of each years values
    
    samples <- 
      data[, 
           .(.sample_year = as.numeric(sample(x = rep(unique(.sample_year), 2),
                                             size = n.in.group))),
           by = groups]
    
    data <- 
      data[samples,
           on = c(groups, ".sample_year")]
    
    data[, .sample_year := NULL]
    
    data
    
  }

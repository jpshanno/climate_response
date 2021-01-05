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
    # goes from 1 of each year value to two of each years values. This was a 
    # problem previously, but is not longer a problem now that I remove site/status
    # where there is only one year of data. If there is only one year I would 
    # rather have it as a training year set than a validation set to ensure that
    # all necessary models exist. If validation is an issue those years could be
    # split, but that would likely limit data too much for model fitting and 
    # require updating the validation approach. There is only one site/status 
    # combination (139 Control) with this issue.
    
    samples <- 
      data[, 
           .(.sample_year = as.numeric(sample(x = rep(unique(.sample_year), 2),
                                             size = n.in.group))),
           by = groups]
    
    # Identify site/status combinations that only have a single year of data
    single_year_combinations <- 
      data[, .(.n_years = uniqueN(.sample_year)), 
           by = groups
      ][.n_years == 1, 
        ..groups]
    
    # Remove single year combinations from consideration as validation data (see
    # note above). I was removing this from the data before sampling, which 
    # likely changed the sampled years. This led to an unexpected error when 
    # setting method = DASvar in model_meteo_flow.R. I'm not sure what was 
    # causing that warning (and it went away if you tried to run the lines again
    # this the console). But I decided to take this approach and be able to 
    # leave the method specification.
    samples <- 
      samples[!single_year_combinations,
              on = groups]
    
    data <- 
      data[samples,
           on = c(groups, ".sample_year")]
    
    data[, .sample_year := NULL]
    
    data
    
  }

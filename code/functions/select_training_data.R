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
    
    set.seed(1234)
    
    
    # new Plan:
    # Control Sites: 
    #   2012 and two other random years
    # Treatment Sites:
    #   Control Status: Even years even months, odd years odd months
    #   Treatment Status: Random three years from post-treatment
    #   The only difficulty here is that Sy calculations may be more difficult
    #   Because the water levels never get so low. To solve this create the 
    #   Control models and then use the posterior distributions of the Esy coefs
    #   as priors in the treatment models.
    
    data[, .sample_year := year(sample_date)]

    # Identify site/status combinations that only have a single year of data
    single_year_combinations <- 
      data[, .(.n_years = uniqueN(.sample_year)), 
           by = .(site, site_status)
      ][.n_years == 1, 
        .(site, site_status)]    
    
    
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
      data[treatment == "Control" & .sample_year != 2012, 
           .(.sample_year = as.numeric(sample(x = unique(.sample_year),
                                              replace = FALSE,
                                              size = 2))),
           by = .(site)]
    
    # Identify site/status combinations that only have a single year of data
    single_year_combinations <- 
      data[, .(.n_years = uniqueN(.sample_year)), 
           by = .(site, site_status)
      ][.n_years == 1, 
        .(site, site_status)]
    
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

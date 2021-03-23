##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author Joe Shannon
##' @export
swg_single_site <-   
  function(con.data, swg.models = NULL, simulation.dates, n.workers = 1, n.simulations, seed = 1234, solar.coefs = NULL){
  
  # For some reason this does not work when used with .SD inside a nested 
  # data.table. It hangs on exit. 
  
  set.seed(seed)
  
  if(nrow(con.data) < 365){
    stop("con.data must be more than 1 year of data (30 is nice)")
  }
  
  if(is.null(swg.models)){
    swg.models <- 
      swg_build_models(con.data)
  }
    
    station_name <- 
      con.data$station_name[1]
    
    lat <- 
      con.data$lat[1]
    
    # Should calculated PET etc when donig format_swg
    con.data <- 
      copy(con.data) %>% 
      calculate_mean_temp() %>% 
      calculate_solar_radiation(coefs = solar.coefs,
                                return.vector = FALSE) %>% 
      calculate_hargreaves_pet(lambda.MJ.kg = 2.45)
    
    
    mean_snowfall <- 
      calculate_snowmelt(con.data, 
                         cn.params = c(0.4, 2.9), 
                         return.mean.snowfall = TRUE)$mean_snowfall_mm
    
    swg.coefs <-
      flatten(swg.models[1, .(occur_coefs,
                           amount_coefs,
                           tmin_coefs,
                           tmax_coefs,
                           regional_monthly_sd_tmin_c,
                           regional_monthly_sd_tmax_c,
                           amount_shape)])
    
  sim_data <- 
    data.table(simulation_date = simulation.dates,
               simulation_month = month(simulation.dates),
               simulation_season = as.climate_season(simulation.dates, TRUE),
               simulation_year = year(simulation.dates))
  
  # Shift harmonic pair to align them with curves starting on 1/1
  harmonics <- 
    data.table(temp_date = seq(as.Date(paste0(year(simulation.dates[1]), "-01-01")), 
                               as.Date(paste0(year(tail(simulation.dates, 1)), "-12-31")),
                               by = "days"))
  
  harmonics[, `:=`(harmonic_cos = cos((2*pi*((1:.N) / .N))),
                   harmonic_sin = sin((2*pi*((1:.N) / .N)))),
            by = .(year(temp_date))]
  
  sim_data[harmonics, 
           `:=`(harmonic_cos = i.harmonic_cos,
                harmonic_sin = i.harmonic_sin),
           on = c("simulation_date" = "temp_date")]
  
  # Not worrying about linear trends within the model construction period
  # because I am not simulating a given year. (Though for testing the model it
  # could help reduce bias between the study period and the SWG)
  sim_data[, temperature_trend := 0]
  
  days_per_sim <- 
    length(simulation.dates)
  
  total_n <- 
    days_per_sim * n.simulations
  
  on.exit(plan(sequential))

  plan(multisession,
       workers = n.workers)
  # 
  # simulation_chunks <- 
  #   split(seq_len(n.simulations), paste0("chunk_", seq_len(n.simulations) %% 4))
  
  chunk_size <- 
    n.simulations / n.workers
  
  simulations <- 
    # future_map_dfr(simulation_chunks,
                   # .options = furrr_options(seed = TRUE),
                   # .id = "simulation_id",
    future_map(seq_len(n.workers), 
               .options = furrr_options(seed = TRUE),
               ~{
                 replicate(n = chunk_size,
                           simplify = FALSE,
                           {
                             simulation_season <- 
                               sim_data$simulation_season
                             
                             simulation_month <- 
                               sim_data$simulation_month
                             
                             sim_occur <- 
                               sim_amount <- 
                               sim_tmax <- 
                               sim_tmin <- 
                               numeric(days_per_sim)
                             
                             min_t_difference <- 0.5
                             
                             init_conditions <- 
                               con.data[format(simulation.dates[1] - 1, "%m%d") == format(sample_date, "%m%d"),
                                        .(mean_precip_occur = mean(precip_occur, na.rm = TRUE), 
                                          mean_precip_amount = mean(precip_amount, na.rm = TRUE),
                                          sd_precip_amount = sd(precip_amount, na.rm = TRUE), 
                                          mean_tmax_c = mean(tmax_c, na.rm = TRUE),
                                          sd_tmax_c = sd(tmax_c, na.rm = TRUE), 
                                          mean_tmin_c = mean(tmin_c, na.rm = TRUE),
                                          sd_tmin_c = sd(tmin_c, na.rm = TRUE))]
                             
                             # initialize with a sample from binomial with probability equal to 
                             # the observed probability of occurrence for 31 Dec (all years)
                             sim_occur[1] <- sim_occur_l1 <- 
                               rbinom(n=1, size=1, prob = init_conditions$mean_precip_occur)
                             
                             # initialize temperatures with the climatological mean of 31 Dec
                             initial_tdiff <- -10
                             
                             while(initial_tdiff < min_t_difference){
                               
                               sim_tmax[1] <- sim_tmax_l1 <- 
                                 rnorm(1, mean = init_conditions$mean_tmax_c, sd = init_conditions$sd_tmax_c)
                               
                               sim_tmin[1] <- sim_tmin_l1 <- 
                                 rnorm(1, mean = init_conditions$mean_tmin_c, sd = init_conditions$sd_tmin_c)
                               
                               initial_tdiff <- sim_tmax[1] - sim_tmin[1]
                             }
                             
                             # Sample a conditioning year at random each simulation
                             
                             conditioning <- 
                               con.data[year(sample_date) == sample(unique(year(sample_date)), 1)]
                             
                             regional_precip_cm <- 
                               as.list(conditioning[, map_dbl(.SD, ~.x[.x!=0][1]), 
                                                    .SDcols = patterns("regional_precip_cm")]) %>% 
                               set_names(~stringr::str_extract(.x, "[a-z]{3}$")) %>% 
                               imap(~ifelse(simulation_season == .y, .x, 0))
                             
                             regional_tmin_c <- 
                               as.list(conditioning[, map_dbl(.SD, ~.x[.x!=0][1]), 
                                                    .SDcols = patterns("regional_tmin_c")]) %>% 
                               set_names(~stringr::str_extract(.x, "[a-z]{3}$")) %>% 
                               imap(~ifelse(simulation_season == .y, .x, 0))
                             
                             regional_tmax_c <- 
                               as.list(conditioning[, map_dbl(.SD, ~.x[.x!=0][1]), 
                                                    .SDcols = patterns("regional_tmax_c")]) %>% 
                               set_names(~stringr::str_extract(.x, "[a-z]{3}$")) %>% 
                               imap(~ifelse(simulation_season == .y, .x, 0))
                             
                             # Scale depends on simulation dates & conditioning covariates
                             amount_scale <- 
                               gamma.scale(coefs = swg.coefs$amount_coefs, 
                                           X = rbind(1, 
                                                     sim_data$harmonic_cos, 
                                                     sim_data$harmonic_sin, 
                                                     regional_precip_cm$djf,
                                                     regional_precip_cm$mam,
                                                     regional_precip_cm$jja,
                                                     regional_precip_cm$son),
                                           shape = swg.coefs$amount_shape)
                             
                             monthly_tmin_sd <- 
                               swg.coefs$regional_monthly_sd_tmin_c
                             
                             monthly_tmax_sd <- 
                               swg.coefs$regional_monthly_sd_tmax_c
                             
                             # Reduce subsets within loop to speed things up
                             regional_precip_cm_djf <- regional_precip_cm$djf
                             regional_precip_cm_mam <- regional_precip_cm$mam
                             regional_precip_cm_jja <- regional_precip_cm$jja
                             regional_precip_cm_son <- regional_precip_cm$son
                             
                             regional_tmin_c_djf <- regional_tmin_c$djf
                             regional_tmin_c_mam <- regional_tmin_c$mam
                             regional_tmin_c_jja <- regional_tmin_c$jja
                             regional_tmin_c_son <- regional_tmin_c$son
                             
                             regional_tmax_c_djf <- regional_tmax_c$djf
                             regional_tmax_c_mam <- regional_tmax_c$mam
                             regional_tmax_c_jja <- regional_tmax_c$jja
                             regional_tmax_c_son <- regional_tmax_c$son
                             
                             for(i in 2:days_per_sim){
                               
                               # A note from Verdin:
                               # Set very small amounts equal to zero and then set those days to have zero 
                               # occurence. Also set days with zero occurence to zero amount. Why are amounts
                               # always simulated? Faster to simulate and remove than to have an IF statement
                               # dependent on occrence? If amount is less than 0.1 should it actually be 
                               # re-simulated within the SWG so that your occurance probabilities are not 
                               # different than what you intended? Does this have a measurable impact? It 
                               # probably could in areas with many small storms. -- It did have an impact
                               # for my SWG, so I put in a while() loop to make sure it all amounts are
                               # greater than 0.01. The while loop makes the number of days with precip 
                               # slightly better aligned, but it overestimates precipitation totals at 
                               # the monthly and annual scales
                               
                               
                               while(sim_amount[i] < 0.01){
                                 # probit has variance 1 by definition, so we use a copula
                                 # to transform a mean zero, variance unity random sample
                                 sim_amount[i] <- qgamma(pnorm(rnorm(n = 1, mean = 0, sd = 1)), 
                                                         shape= swg.coefs$amount_shape, 
                                                         scale= amount_scale[i])
                                 
                                 # occurrence is mean function + rnorm(mean=0, sd=1)
                                 # covariates for mean function are POCC, ct, st
                                 # initial sim_occur is overwritten when i == 1, 
                                 # which is okay because it is using the initial values from outside the
                                 # while loop to provide starting conditions for the true i == 1 inside the
                                 # for loop
                                 sim_occur[i] <- 
                                   as.numeric((sum(swg.coefs$occur_coefs * 
                                                     c(1, 
                                                       sim_occur_l1, 
                                                       sim_data$harmonic_cos[i],
                                                       sim_data$harmonic_sin[i],
                                                       regional_precip_cm_djf[i],
                                                       regional_precip_cm_mam[i],
                                                       regional_precip_cm_jja[i],
                                                       regional_precip_cm_son[i])) + rnorm(1, 0, 1)) > 0)
                               }
                               
                               sim_t_difference <- -10
                               
                               while(sim_t_difference < min_t_difference){
                                 
                                 # max is mean function + rnorm TMAX.sd according to month
                                 # covariates for mean function are PMN, PMX, ct, st, OCC, Rt
                                 sim_tmax[i] <- 
                                   sum(swg.coefs$tmax_coefs *
                                         c(1,
                                           sim_tmin_l1,
                                           sim_tmax_l1,
                                           sim_data$harmonic_cos[i],
                                           sim_data$harmonic_sin[i],
                                           sim_occur[i],
                                           sim_data$temperature_trend[i],
                                           regional_tmin_c_djf[i],
                                           regional_tmin_c_mam[i],
                                           regional_tmin_c_jja[i],
                                           regional_tmin_c_son[i],
                                           regional_tmax_c_djf[i],
                                           regional_tmax_c_mam[i],
                                           regional_tmax_c_jja[i],
                                           regional_tmax_c_son[i])) + 
                                   rnorm(n=1, mean=0, sd = monthly_tmin_sd[simulation_month[i]])
                                 
                                 
                                 # min is mean function + rnorm TMIN.sd according to month
                                 # covariates for mean function are PMN, PMX, ct, st, OCC, Rt
                                 sim_tmin[i] <- 
                                   sum(swg.coefs$tmin_coefs *
                                         c(1,
                                           sim_tmin_l1,
                                           sim_tmax_l1,
                                           sim_data$harmonic_cos[i],
                                           sim_data$harmonic_sin[i],
                                           sim_occur[i],
                                           sim_data$temperature_trend[i],
                                           regional_tmin_c_djf[i],
                                           regional_tmin_c_mam[i],
                                           regional_tmin_c_jja[i],
                                           regional_tmin_c_son[i],
                                           regional_tmax_c_djf[i],
                                           regional_tmax_c_mam[i],
                                           regional_tmax_c_jja[i],
                                           regional_tmax_c_son[i])) + 
                                   rnorm(n = 1, mean = 0, sd = monthly_tmax_sd[simulation_month[i]])
                                 
                                 
                                 sim_t_difference <- 
                                   sim_tmax[i] - sim_tmin[i]
                               }
                               
                               sim_occur_l1 <- sim_occur[i]
                               sim_tmin_l1 <- sim_tmin[i]
                               sim_tmax_l1 <- sim_tmax[i]
                             }
                             
                             dat <- 
                               data.table(station_name = station_name,
                                          lat = lat, 
                                          sample_date = simulation.dates,
                                          simulation_month = simulation_month,
                                          simulation_season = simulation_season,
                                          precip_cm = sim_amount * sim_occur,
                                          tmax_c = sim_tmax,
                                          tmin_c = sim_tmin,
                                          tmean_c = (sim_tmax + sim_tmin) / 2)
                             
                             if(!is.null(solar.coefs)){
                               dat <- 
                                 dat %>% 
                                 calculate_solar_radiation(coefs = solar.coefs,
                                                           stop.on.error = FALSE,
                                                           return.vector = FALSE) %>% 
                                 calculate_hargreaves_pet(lambda.MJ.kg = 2.45) %>% 
                                 calculate_snowmelt(mean.snowfall = mean_snowfall)
                             }
                             
                             setnames(dat, "sample_date", "simulation_date")
                             dat[, c("station_name", "lat") := NULL][]
                           })
               })
  
  RPushbullet::pbPost("note", "Simulations Complete")
  # simulations <- copy(simulations)
  
  simulations <- rbindlist(flatten(simulations))
  simulations[, simulation_id := rep(seq_len(n.simulations), each = days_per_sim)]
  setcolorder(simulations, "simulation_id")
  
  # simulations[, precip_cm := precip_cm * precip_occur]
  # simulations[, tmean_c := (tmin_c + tmax_c)/2]
  
  # if(!is.null(solar.coefs)){
  #   
  #   simulations[, `:=`(station_name = con.data$station_name[1],
  #                      lat = con.data$lat[1],
  #                      sample_date = simulation_date)]
  #   
  #   simulations[, solrad_MJ_m2 := calculate_solar_radiation(.SD,
  #                                                           coefs = solar.coefs,
  #                                                           stop.on.error = FALSE,
  #                                                           return.vector = TRUE),
  #               by = .(simulation_id)]
  #   
  #   simulations[, pet_cm := calculate_hargreaves_pet(simulations,
  #                                                    lambda.MJ.kg = 2.45,
  #                                                    return.vector = TRUE)]
  #   
  #   # The CN models (at least the hydrologic component) needs actual dates to 
  #   # run and they cannot be duplicated. It also should have a warm up period 
  #   # (this may not be strictly necessary with the CN snowmelt modeule but I 
  #   # haven't had a chance to test it). So I am giving fake years to the sample
  #   # dates.
  #   simulations[, sample_date := sample_date + years(simulation_id)]
  #   
  #   RPushbullet::pbPost("note", "Simulating Snow")
  #   
  #   # Having trouble getting parallel to speed this up at all. It keeps being
  #   # slower future_map_dfr. Even when manually chunking the data as for the SWG
  #   simulations <- 
  #     calculate_snowmelt(simulations)
  #   
  #   simulations[, c("station_name", "lat", "sample_date") := NULL]
  # }
  # 
  # RPushbullet::pbPost("note", "Snow Simulation Complete")
  
  simulations[]

}



##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author Joe Shannon
##' @export
gamma.scale <- 
  function(fit = NULL, coefs = NULL, X, shape= NULL){
    
    if(!is.null(fit)){
      if(!is.null(coefs) | !is.null(shape)){
        message("When fit is supplied coefs and shape are superceded by values extracted from the fit.")
      }
      coefs <- coef(fit)[1:3]
      shape <- MASS::gamma.shape(fit)$alpha
    }
    
    if(nrow(X) != length(coefs)){
      # stop("X must be supplied as a P x N matrix with P rows corresponding to intercept = 1, a harmonic pair. N should be the length of the desired predictions.")
      stop("nrow(X) must equal length(coefs)")
    }    
    
    exp(colSums(coefs*X, na.rm=T))/shape
    
  }


##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author Joe Shannon
##' @export
swg_build_models <-   
  function(data){
    
    # Keeping model building nested so it's easier if I want to eventually build
    # multisite SWG
    
    dat_names <- 
      names(data)

    swg_models <- 
      data[, 
           .(dat = list(.SD)),
           keyby = .(station_name)]
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    # Using na.action = na.exclude keeps the NAs in place when you call residuals()
    swg_models[, occur_mod := map(dat,
                                  ~glm(precip_occur ~ precip_occur_l1 + harmonic_cos + harmonic_sin + regional_precip_cm_djf + regional_precip_cm_mam + regional_precip_cm_jja + regional_precip_cm_son,
                                       na.action = na.exclude,
                                       data = .x,
                                       family = binomial(probit)))]
    
    swg_models[, occur_resids := map(occur_mod, residuals, type = "working")]
    swg_models[, occur_coefs := map(occur_mod, coef)]
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    swg_models[, amount_mod := map(dat,
                                   ~glm(precip_amount ~ harmonic_cos + harmonic_sin + regional_precip_cm_djf + regional_precip_cm_mam + regional_precip_cm_jja + regional_precip_cm_son,
                                        na.action = na.exclude,
                                        data = .x,
                                        family = Gamma(link = log)))]
    
    swg_models[, amount_coefs := map(amount_mod, coef)]
    swg_models[, amount_shape := map_dbl(amount_mod, ~MASS::gamma.shape(.x)$alpha)]
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    swg_models[, tmin_mod := map(dat,
                                 ~lm(tmin_c ~ tmin_l1 + tmax_l1 + harmonic_cos + harmonic_sin + precip_occur + temperature_trend + regional_tmin_c_djf + regional_tmin_c_mam + regional_tmin_c_jja + regional_tmin_c_son + regional_tmax_c_djf + regional_tmax_c_mam + regional_tmax_c_jja + regional_tmax_c_son, 
                                     data = .x,
                                     na.action = na.exclude))]
    
    swg_models[, tmin_resids := map(tmin_mod, residuals, type = "working")]
    swg_models[, tmin_coefs := map(tmin_mod, coef)]
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    swg_models[, tmax_mod := map(dat,
                                 ~lm(tmax_c ~ tmin_l1 + tmax_l1 + harmonic_cos + harmonic_sin + precip_occur + temperature_trend + regional_tmin_c_djf + regional_tmin_c_mam + regional_tmin_c_jja + regional_tmin_c_son + regional_tmax_c_djf + regional_tmax_c_mam + regional_tmax_c_jja + regional_tmax_c_son, 
                                     data = .x,
                                     na.action = na.exclude))]
    
    swg_models[, tmax_resids := map(tmax_mod, residuals, type = "working")]
    swg_models[, tmax_coefs := map(tmax_mod, coef)]
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    # Temperature standard deviations
    swg_models[,
               regional_monthly_sd_tmin_c := map2(dat, tmin_resids,
                                                  ~tapply(.y, .x$sample_month, sd, na.rm = TRUE))]
    
    swg_models[,
               regional_monthly_sd_tmax_c := map2(dat, tmax_resids,
                                                  ~tapply(.y, .x$sample_month, sd, na.rm = TRUE))]
    
    swg_models[]
  }

swg_format_data <- 
  function(data){
    
    dat_names <- 
      names(data)
    
    stopifnot("station_name" %in% dat_names,
              "sample_date" %in% dat_names,
              "precip_cm" %in% dat_names,
              "tmax_c" %in% dat_names,
              "tmin_c" %in% dat_names)
    
    data <- 
      copy(data)
    
    data[, `:=`(sample_month = month(sample_date),
                sample_year = year(sample_date),
                sample_season = as.climate_season(sample_date, TRUE))]
    
    # Use NA_real_ not 0 because the precip amount model only deals with the
    # amount, so the relationship between the covariates and precip amount would
    # get muddled if we included 0s
    data[, `:=`(precip_occur = as.numeric(precip_cm > 0.01),
                precip_amount = ifelse(precip_cm < 0.01, NA_real_, precip_cm))]
    
    data[, `:=`(tmin_l1 = shift(tmin_c, 1),
                tmax_l1 = shift(tmax_c, 1),
                precip_occur_l1 = shift(precip_occur, 1)),
         by = .(station_name)]
    
    # After lagging variables remove 1979-12-31 to avoid messing up trend &
    # harmonics (loca data doesn't have 1979-12-31, but the single NA is fine)
    data <- 
      data[sample_date != as.Date("1979-12-31")]
    
    data[, temperature_trend := seq(-1, 1, length.out = .N),
         by = .(station_name)]
    
    data[, `:=`(harmonic_cos = cos((2*pi*(1:.N / .N))),
                harmonic_sin = sin((2*pi*(1:.N / .N)))),
         by = .(station_name, sample_year)]
    
    # Site & Regional Monthly & Seasonal Temperatures
    data[, `:=`(site_monthly_tmin_c = mean_na(tmin_c), 
                site_monthly_tmax_c = mean_na(tmax_c)),
         by = .(station_name, sample_year, sample_month)]
    
    data[, `:=`(regional_monthly_tmin_c = mean_na(site_monthly_tmin_c),
                regional_monthly_tmax_c = mean_na(site_monthly_tmax_c)),
         by = .(sample_year, sample_month)]
    
    data[, `:=`(regional_seasonal_tmin_c = mean_na(regional_monthly_tmin_c),
                regional_seasonal_tmax_c = mean_na(regional_monthly_tmax_c)),
         by = .(sample_year, sample_season)]
    
    # ggplot(data, aes(x = yday(sample_date), 
    #                  group = sample_year)) + 
    #   geom_line(aes(y = tmin_c,
    #                 group = interaction(sample_year, station_name, drop = TRUE)),
    #             color = 'lightblue', 
    #             alpha = 0.4) +
    #   geom_line(aes(y = tmax_c,
    #                 group = interaction(sample_year, station_name, drop = TRUE)), 
    #             color = 'red4', alpha = 0.4) + 
    #   geom_line(aes(y = regional_monthly_tmin_c), 
    #             color = 'blue') + 
    #   geom_line(aes(y = regional_monthly_tmax_c), 
    #             color = 'red') + 
    #   geom_line(aes(y = regional_seasonal_tmin_c)) +
    #   geom_line(aes(y = regional_seasonal_tmax_c))
    
    # Precip
    # sum precip by site, year, month
    data[, site_total_monthly_precip_cm := sum_na(precip_cm),
         by = .(station_name, sample_year, sample_month)]
    
    # mean of summed precip by year, month
    data[, regional_mean_monthly_precip_cm := mean_na(site_total_monthly_precip_cm),
         by = .(sample_year, sample_month)]
    
    # sum of mean_monthly_precip by year, season
    # The .SD argument inside sum_na takes just the first row from each month
    # in a year because the monthy mean data is repeated daily in data. This 
    # could be done more straightforward by creating separate monthly mean and
    # season total objects then joining them back to data. I chose to do it this
    # because I wanted to have the option to use monthly data rather than seasonal
    # data. With this approach just a tweak to the dcast() call would make 
    # monthly columns.
    data[, regional_total_seasonal_precip_cm := sum_na(.SD[, .(precip = regional_mean_monthly_precip_cm[1]), by = .(sample_month)]$precip),
         by = .(sample_year, sample_season)]
    
    # Drop unnecesary columns & rename before dcast
    
    data[, c("site_monthly_tmax_c",
             "regional_monthly_tmax_c",
             "regional_mean_monthly_precip_cm",
             "site_total_monthly_precip_cm"):= NULL]
    
    setnames(data,
             c("regional_total_seasonal_precip_cm", "regional_seasonal_tmin_c", "regional_seasonal_tmax_c"),
             c("regional_precip_cm", "regional_tmin_c", "regional_tmax_c"))
    
    dcast(data, 
          ... ~ sample_season,
          value.var = c("regional_tmin_c", "regional_tmax_c", "regional_precip_cm"), 
          fill = 0)
    
  }
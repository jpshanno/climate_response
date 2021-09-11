##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param n
##' @param ...
##' @return
##' @author Joe Shannon
##' @export
run_simulations <- function(n, ...) {

  # Need to include a check if climate variables (specifically water_availabilty)
  # is outside of the observed range. There is a good chance that it could drop
  # below the observed values, which is a problem for the meteo_flow models. 
  # In order to include a harmonic term for water availability in the model I 
  # had to convert water_availability_cm to water_availability_decimal, which 
  # divides observed water level by the diff(range(x)). I think if this occurs,
  # I will just use the minimum observed water availability in place of the 
  # simulated value and then make sure to provide some sort of QA column or notes
  # on affected simulations
  
  NULL

}

##' .. content for \description{} (no empty lines) ..
##'
##' Does not allow for unspecified site/status. Could expand to include fixed
##' effects general models
##'
##' @title
##' 
##' @param site.id
##' @param lat
##' @param site.status
##' @param intercept.models
##' @param morpho.flow.models
##' @param meteo.flow.models
##' @param precip.models
##' @param pet.models
##' @param esy.models
##' @param solrad.coefs
##' @param initial.wl
##' @param initial.water.availability
##' @param weather
##' @param water.levels.only
##' @return
##' @author Joe Shannon
##' @export
simulate_water_levels <- 
  function(site.id,
           lat,
           site.status, 
           intercept.models,
           morpho.flow.models,
           meteo.flow.models,
           precip.models,
           pet.models,
           esy.models,
           solrad.coefs,
           initial.wl, 
           initial.water.availability,
           weather, 
           water.levels.only = FALSE) {
    
    # Need to figure out how to generalized solar radiation coefficients
    # Perhaps by looking at how well the climate simulations match each met 
    # station and then use the coefficients from the station that has the best
    # correlation?
    
    # weather: data.frame with doy, precip_cm, min_temp_c, max_temp_c, water_availability_cm
    # mean_temp_c, water_availability
    
    weather <- copy(weather)
    
    # Checks
    
    stopifnot(all(c("doy", "precip_cm", "tmin_c", "tmax_c") %in% names(weather)))
    
    stopifnot(is.logical(water.levels.only))
    
    stopifnot(!is.na(initial.wl))
    
    site.status <- 
      unique(site.status)
    
    stopifnot(length(site.status) == 1)
    
    site.id <- 
      unique(site.id)
    
    stopifnot(length(site.id) == 1)
    
    lat <- 
      unique(lat)
    
    stopifnot(length(lat) == 1)
    
    # Get Functions
    f_interception <- 
      intercept.models[site.status, f_predict[[1]]]
    
    f_morpho <- 
      morpho.flow.models[site.id, f_predict[[1]]]
      
    f_meteo <- 
      meteo.flow.models[CJ(site.id, site.status), f_predict[[1]]]
    
    f_esy <- 
      esy.models[site.id, f_predict[[1]]]
    
    f_pet <- 
      pet.models[CJ(site.id, site.status), f_predict[[1]]]
    
    f_rise <- 
      precip.models[site.id, f_predict[[1]]]
    
    # Get Observed Ranges
    # pet_range <- 
    #   meteo.flow.models[CJ(site.id, site.status),
    #                     pet_range_cm[[1]]]
    # 
    # water_availability_range <- 
    #   meteo.flow.models[CJ(site.id, site.status), 
    #                     water_availability_range_cm[[1]]]
    # 
    # precip_range <- 
    #   intercept.models[site.status, 
    #                    precip_range_cm[[1]]]
    #   
    # water_level_range <- 
    #   morpho.flow.models[site.id, 
    #                      wl_range_cm[[1]]]

    # Need to make model carry over previous year's final wl + winter precip
    # if(is.na(initial.wl)){
    #   initial.wl <- 
    #     sample(starting_values[site == hydrology_model, water_level_cm], 1)
    # }
    
    if(any(is.na(weather$precip_cm))){
      message("Some missing precip values were filled in as 0.")
      
      weather[, precip_cm := nafill(precip_cm,
                                    type = "const",
                                    fill = 0)]
    }


# Calculate Solar Radiation & PET -----------------------------------------

    # bc() needs sample date but just converts it back to doy, so year doesn't
    # matter
    # calculate_solar_radiation uses station_name (just for grouping) and lat
    weather[, `:=`(lat = lat,
                   station_name = site.id,
                   sample_date = as.Date(doy, origin = "2020-01-01"))]
    
    # calculate_solar_radiation() calculate_mean_temp() and
    # calculate_hargreaves_pet() were designed to work in a pipe so they return
    # modified data as a side effect. This may change in the future
    calculate_solar_radiation(weather, solrad.coefs)
    calculate_mean_temp(weather)
    calculate_hargreaves_pet(weather, lambda.MJ.kg = 2.45)
    
# Calculate Interception ---------------------------------------------------

    weather[, interception_cm := f_interception(doy)]
    weather[, net_precip_cm := pmax(0, precip_cm - interception_cm)]
    
    # Could add Markov component to Ds, it should rise if it was rising previously
    
    wl <- Ds <- esy <- rain_rise <- morpho_flow <- meteo_flow <- water_availability <- drawdown <- slow_flow <- 
      numeric(nrow(weather))
    
    wl[1] <- initial.wl
    water_availability[1] <- initial.water.availability
    dry_days <- 0
    last_precip <- 0
    
    for(i in 1:(nrow(weather)-1)){
      
      esy[i] <- 
        f_esy(wl[i])
      
      morpho_flow[i] <- 
        f_morpho(wl[i])
      
      meteo_flow[i] <- 
        f_meteo(pet = weather[i, pet_cm], water_availability[i]) 
      
      rain_rise[i] <-
        f_rise(weather[i, net_precip_cm])
      
      drawdown[i] <-
        f_pet(weather[i, pet_cm])

      if(weather[i, net_precip_cm] > 0){
        dry_days <- 0
        last_precip <- weather[i, net_precip_cm]
        slow_flow[i] <- 0
      } else {
        dry_days <- dry_days + 1
        slow_flow[i] <- last_precip / 3^dry_days
      }
      
      # drawdown[i] <- 
      #   weather[i, pet_cm]
      
      # if(rain_rise[i] > 0){
      #   drawdown[i] <-
      #     0.5 * drawdown[i]*0.5
      # }
      
      # HAVE TO EXCLUDE MELT for now
      Ds[i] <- 
        (morpho_flow[i] - meteo_flow[i]) - drawdown[i] + rain_rise[i] + 0
      
      water_availability[i+1] <- 
        water_availability[i] + (weather[i, precip_cm] - weather[i, pet_cm])
      
      wl[i+1] <- wl[i] + Ds[i] / esy[i]
    }
    
    data.table(doy = weather$doy,
               wl_hat = wl,
               Ds_hat = Ds,
               esy_hat = esy,
               P_hat = rain_rise,
               ET_hat = drawdown,
               G_morph_hat = morpho_flow,
               G_meteo_hat = meteo_flow,
               wa_hat = water_availability)
    
    # pred_budget <-
    #   data.table(doy = weather$doy,
    #              g_precip = weather$precip_cm,
    #              pred_i = weather$interception_cm,
    #              pred_wl = wl,
    #              pred_Ds = Ds,
    #              pred_esy = esy,
    #              rain_rise,
    #              drawdown,
    #              morpho_flow,
    #              meteo_flow,
    #              pred_water_availability = water_availability)
    # 
    # ggplot(test_dat,
    #        aes(x = doy)) +
    #   geom_line(aes(y = wl_initial_cm)) +
    #   geom_line(data = pred_budget,
    #             aes(y = pred_wl),
    #             color = "red",
    #             linetype = "dashed") +
    #   geom_col(data = pred_budget,
    #            aes(y = morpho_flow),
    #            fill = "green",
    #            alpha = 0.25) +
    #   theme_bw()
    # 
    # ggplot(test_dat,
    #        aes(x = doy)) +
    #   geom_line(aes(y = Ds_cm)) +
    #   geom_line(data = pred_budget,
    #             aes(y = pred_Ds),
    #             color = "red",
    #             linetype = "dashed") +
    #   theme_bw()
    # 
    # ggplot(test_dat,
    #        aes(x = doy)) +
    #   geom_line(aes(y = net_flow_cm)) +
    #   geom_line(data = pred_budget,
    #             aes(y = (morpho_flow - meteo_flow)),
    #             color = "red",
    #             linetype = "dashed") +
    #   theme_bw()
    # 
    # ggplot(test_dat,
    #        aes(x = doy)) +
    #   geom_line(aes(y = pet_cm)) +
    #   geom_line(data = pred_budget[drawdown != 0],
    #             aes(y = drawdown),
    #             color = "red",
    #             linetype = "dashed") +
    #   theme_bw()
    # 
    # if(components){
    #   return(data.table(hydrology_model, site_status, doy = dat$doy, wl, interception = dat$interception, precip_rise, pet, flow, esy))
    # }
    # 
    # # Check for out of observed range values
    # 
    # 
    # wl
  }

set.seed(1234)
source("climate_packages.R")
library(parallel)

# Constant for PET calculation
lambda_MJ_kg <- 
  2.45

net_flow_mods <- 
  readRDS("net_flow_mods.rds")
interception <- 
  readRDS("interception.rds")
sy_mods <- 
  readRDS("sy_mods.rds")
starting_values <- 
  readRDS("starting_values.rds")
wb_mods <- 
  readRDS("wb_mods.rds")
water_budget <- 
  readRDS("site_status_combinations.rds")

full_predict <-
  function(site_status, 
           hydrology_model = NULL, 
           initial.wl = NA,
           weather, 
           components = FALSE,
           random.effects = NULL){
    
    # weather: data.frame with doy, gross_precip_cm, min_temp_c, max_temp_c, 
    # mean_temp_c, water_availability
    
    stopifnot(is.logical(components))
    
    site_status <- 
      unique(site_status)
    
    stopifnot(length(site_status) == 1)
    
    if(is.null(hydrology_model)){
      
      hydrology_model <- 
        sample(net_flow_mods[treatment == site_status, site], 1)
      
    } else {
      
      hydrology_model <- 
        unique(hydrology_model)
      
      stopifnot(length(hydrology_model) == 1)
      
    }
    
    if(is.na(initial.wl)){
      initial.wl <- 
        sample(starting_values[site == hydrology_model, water_level_cm], 1)
    }
    
    precip_rise_mod <- 
      wb_mods[site_status, 
              precip_rise[[1]]]
    
    pet_mod <- 
      wb_mods[site_status, 
              outflow_pet_single[[1]]]
    
    # Can't use predict for flow if using the power function
    flow_mod <- 
      net_flow_mods[hydrology_model]
    
    slope <- flow_mod$m
    offset <- flow_mod$c
    
    # Get ESy Mod
    esy_mod <- 
      sy_mods[hydrology_model, mod[[1]]]
    
    dat <- 
      copy(weather)
    
    if(any(is.na(dat$gross_precip_cm))){
      message("Some missing precip values were filled in as 0.")
      
      dat[, gross_precip_cm := nafill(gross_precip_cm,
                                      type = "const",
                                      fill = 0)]
    }
    
    dat[, `:=`(site = hydrology_model,
               site_status = site_status)]
    
    dat[, season := fifelse(between(doy, 135, 258), 
                            "growing",
                            "dormant")]
    
    dat[interception, 
        interception_cm := i.i_cm,
        on = c("site", "site_status", "season")]
    
    # Generated for Wakefield (could be done spatially with GLM)
    ha_coefs <- 
      c(Ha = 0.1587725,
        Hb = -1.7921571,
        Hr2 = 0.7908223)
    
    # Predict Radation & PET
    dat[, ET_rad_MJ_m2 := extrat(doy, radians(46.440278))$ExtraTerrestrialSolarRadiationDaily]
    dat[, ha_rad_MJ_m2 := ha(as.Date(doy, origin = "2018-12-31"),  # non-leap year date from doy (leap-year vs not does not change radiation estimate)
                             lat = 46.440278,
                             lon = -89.826667,
                             ET_rad_MJ_m2,
                             A = ha_coefs[["Ha"]],
                             B = ha_coefs[["Hb"]],
                             Tmax = max_temp_c,
                             Tmin = min_temp_c)]
    
    dat[, harad_pet_hs_cm := 0.1 * 0.0023 * ha_rad_MJ_m2 / lambda_MJ_kg * (mean_temp_c + 17.8) * sqrt(max_temp_c - min_temp_c)]
    
    # Could add Markov component to Ds, it should rise if it was rising previously
    
    
    
    wl <- pet <- esy <- precip_rise <- flow <- 
      numeric(nrow(dat))
    
    wl[1] <- 
      initial.wl
    
    i_water_availability <- 
      0
    
    for(i in 1:(nrow(dat)-1)){
      
      esy[i] <- 
        predict(esy_mod,
                type = "response",
                newdata = data.frame(water_level_cm = wl[i]))
      
      esy[i] <- 
        pmin(1, esy[i])
      
      # There probably are not any asymptotes in the univariate gamma model  &
      # definitely none in the NL model. Leaving this in case I try the quadratic
      # gamma model again
      
      esy[i] <- 
        ifelse(esy[i] < 0, 1, esy[i])
      
      i_pet <- 
        dat[i, harad_pet_hs_cm / esy[i]]
      
      i_precip <- 
        dat[i, pmax(0, gross_precip_cm - interception_cm) / esy[i]]
      
      i_water_availability <- 
        i_water_availability + ((i_precip - i_pet) / 100)
      
      pet[i] <- 
        predict(pet_mod,
                re.form = random.effects,
                newdata = data.frame(site = hydrology_model,
                                     harad_pet_hs_cm = i_pet,
                                     water_availability = i_water_availability))
      
      precip_rise[i] <- 
        predict(precip_rise_mod, 
                re.form = random.effects,
                newdata = data.frame(site = hydrology_model,
                                     net_precip_cm = i_precip))
      
      # Set rise to zero for days with no precip (model has an intercept)
      precip_rise[i] <- 
        fifelse(i_precip > 0, precip_rise[i], 0)
      
      flow[i] <- 
        slope ** (wl[i] + offset)
      
      flow[i] <- 
        ifelse(wl[i] > (-offset),
               pmin(flow[i], wl[i] + offset),
               flow[i])
      
      wl[i+1] <- wl[i] + precip_rise[i] + pet[i] - flow[i]
    }
    
    if(components){
      return(data.table(hydrology_model, site_status, doy = dat$doy, wl, interception = dat$interception, precip_rise, pet, flow, esy))
    }
    wl
  }

future_sequences <- 
  fread("weather_sequences_gfdl-cm3_2070-2099.csv")

future_sequences[, mean_temp_c := (tmin + tmax) / 2]

future_sequences <- 
  future_sequences[between(month(sample_date), 6, 10)]

future_sequences <-
  CJ(site_status = c("Control", "Ash Cut", "Girdle"),
     site = c("077", "151", "156"))[, .(dat = list(future_sequences)), by = .(site, site_status)]

future_sequences <- 
  future_sequences[, dat[[1]], by = .(site, site_status)]

future_sequences[starting_values[, .(water_level_cm = median(water_level_cm)), by = .(site)],
                 start_level_cm := i.water_level_cm,
                 on = "site"]

cat("\n  Running Simulations \n")

set.seed(1234)
simulations <- 
  split(future_sequences, 
        by = c("site", "site_status", "seq_ID")) %>%
  mclapply(function(x){
    x[,
      c("hydrology_model", "site_status_mod", "doy_mod", "wl", 
        "interception", "precip_rise", "pet", "flow", "esy") :=
        full_predict(hydrology_model = x$site,
                     site_status = x$site_status,
                     initial.wl = x$start_level_cm,
                     weather = .SD[, 
                                   .(doy, 
                                     gross_precip_cm = precip / 10,
                                     min_temp_c = tmin,
                                     max_temp_c = tmax,
                                     mean_temp_c)],
                     components = TRUE,
                     random.effects = NA)]
  },
  mc.cores = 16) %>% 
  rbindlist()

saveRDS(simulations, "simulations_fixed_effects_2070-2099.rds")
fwrite(simulations, "simulations_fixed_effects_2070-2099.csv")
RPushbullet::pbPost("note", "Future Simulations Complete")


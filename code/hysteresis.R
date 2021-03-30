# Hysteresis --------------------------------------------------------------

# This looks like the best approach. There seems to be little chance of runaway
# water level declines because the function is built on WA not an ESY derived
# from water levels. 

# Get just the declining portion of water availability (should probably consider
# doing this as high sprint water level to lowest water level, possible with the
# detrended annual levels). For now the max and min available water serves as
# a good proxy. This could probably be slightly improved by getting just the 
# first local minimum, or try to remove any rewetting periods (but that may 
# reduce too much data). 
# The only serious improvement that seems to be necessary is to have slightly
# better interannual alignment of the water availability. It looks like using
# the previous year's water deficit/surplus may be enough to align the curves. 
# I tried aligning by adjusting to previous year's water defecit/surplus, but 
# that doesn't quite do it. It aligned the curve better if you set 
# limits on the adjustment, but I would have to figure out what controlled those
# limits. I am not messing with it becuase it may be that the adjustment has to do with how
# much more water is necessary to bring the water level up (which is dependent
# on the slope of the recession curve at the time), or with the melt equation 
# calibration.

source("code/load_project.R")
tar_load(external_met)
tar_load(water_budget)

SITE <- 
  "157"

kenton <- 
  external_met[station_name == "kenton"]

dat <- 
  kenton[water_year == 2012 & format(sample_date, "%m%d") != "0229",
           .(sample_date, 
             water_year,
             pet_cm = -pet_cm,
             # ytd_pet_cm = cumsum(-pet_cm),
             rain_cm,
             # ytd_p_cm = cumsum(rain_cm),
             melt_cm,
             # ytd_m_cm = cumsum(melt_cm),
             total_input_cm,
             # ytd_input_cm = cumsum(total_input_cm),
             tmax_c,
             tmin_c,
             tmean_c)
         ][water_budget[site == SITE],
           `:=`(site = i.site,
                wl_initial_cm = i.wl_initial_cm),
           on = "sample_date"]

dat[, site := unique(na.omit(site))]

# Drop anomalous melt
dat[between(melt_cm, 0, 0.1) & month(sample_date) %in% 6:9, 
    `:=`(total_input_cm = total_input_cm - melt_cm,
         melt_cm = 0)]

# Drop anomalous pet
dat[tmax_c <= 0,
    pet_cm := 0]

# Calculate YTD after dropping above points
dat[, 
    `:=`(ytd_pet_cm = cumsum(pet_cm),
         ytd_p_cm = cumsum(rain_cm),
         ytd_m_cm = cumsum(melt_cm),
         ytd_input_cm = cumsum(rain_cm + melt_cm)),
    by = .(water_year)]

dat[, ytd_water_balance := ytd_pet_cm + ytd_p_cm + ytd_m_cm]

# Quadratic ---------------------------------------------------------------

quad <- 
  function(wa, b0, b1, b2){
    b0 + b1*wa + b2*wa^2
  }

quad_prime <- 
  function(wa, b1, b2){
    b1 + b2 * 2 * wa
  }

# wd is water deficit from previous year

# Initial Parameters:
# LOOK CLOSE wl_initial_cm is -8 not -max.wl (8 is max.wl for 152), probably
# should run it with max.wl
max.wl <- 
  density(water_budget[site == SITE, na.omit(wl_initial_cm)])$x[which.max(density(water_budget[site == SITE, na.omit(wl_initial_cm)])$y)]

(init_mod <- 
  robustbase::nlrob(wl_initial_cm ~ quad(ytd_water_balance, 
                                       b0 = b0, b1 = b1, b2 = b2),
                  data = dat[1:which.min(wl_initial_cm), .(wl_initial_cm = wl_initial_cm - 8, ytd_water_balance = ytd_water_balance - max(ytd_water_balance))],
                  na.action = na.exclude,
                  maxit = 50,
                  start = list(b0 = 8, b1 = 1, b2 = -1)))

ggplot(data = dat[1:which.min(wl_initial_cm), 
                  .(wl_initial_cm = wl_initial_cm - 8, 
                    ytd_water_balance = ytd_water_balance - max(ytd_water_balance))],
       aes(x = ytd_water_balance,
           y = wl_initial_cm)) +
  geom_path() +
  geom_function(fun = ~predict(init_mod, newdata = data.frame(ytd_water_balance = .x)),
                color = palette()[4]) 

b0 <- coef(init_mod)[["b0"]]
b1 <- coef(init_mod)[["b1"]]
b2 <- coef(init_mod)[["b2"]]

PET <- dat$pet_cm
P <- dat$rain_cm
M <- dat$melt_cm
D <- dat$sample_date
WL <- dat$wl_initial_cm
WA <- dat$ytd_water_balance - max(dat$ytd_water_balance)

quad_curve_wl <- 
  function(PET, P, M, max.wl, b0, b1, b2){
    
    # Length of weather vector
    n <- 
      length(PET)
    
    # Create empty vectors
    wa_hat <- 
      wl_hat <- 
      gradient <- 
      q_hat <- 
      g_hat <- 
      max_wa <- 
      numeric(n)
    
    # Initialize model at full water level
    wl_hat[1] <- 
      max.wl
    
    # Initialize water availability so that max annual water availability = 0
    # (max_wa could be a single value rather than a vector, but was turned into
    # a vector to make debugging easier)
    max_wa[1:365] <- 
      -max(cumsum(PET + P + M)[1:365])
    
    wa_hat[1] <- 
      -max(cumsum(PET + P + M)[1:365])
    
    # Loop through weather data
    for(t in 2:n){
      
      # Reset WA after each year so that each annual max water availability = 0
      if((t-1) %% 365 == 0){
        
        # Calculate maximum water availability for the next year of data
        max_wa[t:(t+364)] <-
          max(cumsum(PET[t:(t+364)] + P[t:(t+364)] + M[t:(t+364)]))

        # Reset wa_hat
        wa_hat[t] <-
          -max_wa[t]
        
        # Reset wl_hat
        wl_hat[t] <- 
          max.wl
        
        next
        
      } else {
        
        # If not start of year adjust water availability according to drivers
        wa_hat[t] <- 
          wa_hat[t-1] + PET[t] + M[t] + P[t]        
      }
      
      # Calculate gradient of drawdown 
      gradient[t] <- 
        quad_prime(wa_hat[t], b1 = b1, b2 = b2)
      
      # If gradient is negative (water level is above b0 from model) assume
      # that there is flow into the system
      if(gradient[t] <= 0){
        
        # Calculate g_hat as the water level difference between where water
        # availability leads to gradient == 0 & the modeled water level based
        # on water availability. Divided by 4 right now because it was too big
        # an input otherwise, this should probably be an optimized parameter
        # To improve the model I need to put in a system to 'bank' excess water
        # availablility and allow some of it to be drawn down more slowly. It 
        # may be something as simple as increasing the maximum water level w/o
        # increasing b0, so that there can be greater WL rise from g_hat before
        # Q draws WL down. But better would be some system that keeps WA elevated
        # for some amount of time after snowmelt
        g_hat[t] <- 
          (quad((-b1/(2*b2)), b0, b1, b2) - quad(wa_hat[t], b0, b1, b2))/4
        
        # Add G to water previous water level
        wl_hat[t] <- 
          wl_hat[t-1] + g_hat[t]
        
        # If gradient is positive:
      } else {
        
        # Use net input to determine if water level increases or decreases
        if((P[t] + PET[t]) <= 0){
          
          # Water level drawdown = PET, if P <= PET then it can be assumed to be
          # less than interception (not necessarily true, but works as a 
          # simplifying assumption)
          wl_hat[t] <-
            wl_hat[t-1] + (PET[t]) * gradient[t]
        
          } else {

            # Water rise is P - PET, which fits better, but need to rethink my
            # justification for this
            wl_hat[t] <-
              wl_hat[t-1] + (P[t] - PET[t]) * gradient[t]
            
          }
        
      }
      
      # Directly add melt to water level. This should probably have some sort of
      # multiplier
      wl_hat[t] <-
        wl_hat[t] + M[t]
      
      # If WL is above spill point threshold then lose some to streamflow. 
      # This could probably be improved using the morphology models to determine
      # streamflow
      
      if(wl_hat[t-1] > max.wl){
        
        q_hat[t] <- 
          (wl_hat[t-1] - max.wl)
        
        wl_hat[t] <- 
          wl_hat[t] - q_hat[t]
      }
      
    }
    
    data.table(wl_hat, wa_hat, gradient, q_hat, g_hat, max_wa)
    
  }

# This looks good
lplot(wl_hat ~ D); lines(WL ~ D, col = 'gray40')
# lplot(wa_hat ~ D)
# lplot(q_hat ~ D)
lplot(g_hat ~ D)
lplot(s_hat ~ D)
# lplot(gradient ~ D)
# hydroGOF::gof(wl_hat[!is.na(WL)], WL[!is.na(WL)])
# lplot(c + cumsum(gradient * P + gradient * PET))
# c + last(cumsum(gradient * P + gradient * PET))

# Optimize ----------------------------------------------------------------

quad_opt <- 
  function(wobs, met, ...){
    
    mod_dat <- 
      data.frame(wl = wobs)
    
    mod_dat$ytd_water_balance <- 
      met$rain_cm + met$pet_cm + met$melt_cm
    
    mod_dat <- 
      mod_dat[1:which.min(wobs), ]
    
    mod_dat$ytd_water_balance <- 
      mod_dat$ytd_water_balance - max(mod_dat$ytd_water_balance)
    
    ss_mod <- 
      robustbase::nlrob(wl ~ quad(ytd_water_balance, 
                                  b0 = b0, b1 = b1, b2 = b2),
                        data = mod_dat,
                        na.action = na.exclude,
                        maxit = 100,
                        control = nls.control(maxiter = 100),
                        start = list(b0 = 8, b1 = 1, b2 = -1))
    
    opt_start <- 
      list(b0 = coef(ss_mod)[["b0"]],
           b1 = coef(ss_mod)[["b1"]],
           b2 = coef(ss_mod)[["b2"]])
    
    nse <- 
      function(params, wobs, met, max.wl, wa.init){
        
        wl_hat <- 
          quad_curve_wl(PET = met$pet_cm,
                        P = met$rain_cm, 
                        M = met$melt_cm, 
                        max.wl = max.wl,
                        b0 = params[["b0"]],
                        b1 = params[["b1"]],
                        b2 = params[["b2"]])$wl_hat
        
        hydroGOF::NSE(wl_hat[!is.na(wobs)], wobs[!is.na(wobs)])
        # dtw::dtw(wl_hat[!is.na(wobs)], wobs[!is.na(wobs)])$distance
      }
    
    optim(par = opt_start,
          fn = nse,
          gr = NULL,
          met = met,
          wobs = wobs,
          ...)
  }

# Initial Parameters from nlrob for 152, 2012
# Could make it self-start
test_opt <- 
  quad_opt(wobs = dat$wl_initial_cm, 
           met = dat[, .(pet_cm, rain_cm, melt_cm)],
           max.wl = max.wl,
           control = list(fnscale = -1))

lplot(wl_initial_cm ~ sample_date,
      data = dat, 
      col = "gray40")
lines(y = quad_curve_wl(dat$pet_cm, dat$rain_cm, dat$melt_cm, 
                    max.wl = max.wl,
                    b0 = test_opt$par[["b0"]],
                    b1 = test_opt$par[["b1"]],
                    b2 = test_opt$par[["b2"]])$wl_hat, x= dat$sample_date)
# 2014 & 2018 need work
mod_dat <- 
  kenton[between(water_year, 2012, 2019),
         .(sample_date, 
           water_year,
           pet_cm = -pet_cm,
           rain_cm,
           melt_cm,
           total_input_cm,
           tmax_c,
           tmin_c,
           tmean_c)
  ][water_budget[site == "077",
                 .(sample_date, wl_initial_cm, idw_precip_cm)],
    `:=`(wl_initial_cm = i.wl_initial_cm,
         idw_precip_cm = i.idw_precip_cm),
    on = "sample_date"]

# Drop Leap year
mod_dat <- 
  mod_dat[format(sample_date, "%m%d") != "0229"]

# Drop anomalous melt
mod_dat[between(melt_cm, 0, 0.1) & month(sample_date) %in% 6:9, 
         `:=`(total_input_cm = total_input_cm - melt_cm,
              melt_cm = 0)]

# Drop anomalous pet
mod_dat[tmax_c <= 0,
         pet_cm := 0]

# Calculate YTD after dropping above points
mod_dat[, 
       `:=`(ytd_pet_cm = cumsum(pet_cm),
            ytd_p_cm = cumsum(rain_cm),
            ytd_m_cm = cumsum(melt_cm),
            ytd_input_cm = cumsum(rain_cm + melt_cm)),
       by = .(water_year)]

mod_dat[, 
        `:=`(ytd_water_balance_cm = ytd_pet_cm + ytd_p_cm + ytd_m_cm,
             water_balance_cm = cumsum(pet_cm + rain_cm + melt_cm))]

# Split Data
tr_dat <- 
  mod_dat[water_year %in% c(2012)]

tr_params <- 
  quad_opt(wobs = tr_dat$wl_initial_cm, 
           met = tr_dat[, .(pet_cm, rain_cm, melt_cm)],
           max.wl = density(mod_dat$wl_initial_cm, na.rm = TRUE)$x[which.max(density(mod_dat$wl_initial_cm, na.rm = TRUE)$y)],
           control = list(fnscale = -1))

mod_dat[, c("wl_hat_cm", "wa_hat_cm", "gradient", "q_hat", "g_hat", "max_wa") := 
          quad_curve_wl(PET = pet_cm,
                        P = fcoalesce(idw_precip_cm, rain_cm), 
                        M = melt_cm, 
                        max.wl = density(mod_dat$wl_initial_cm, na.rm = TRUE)$x[which.max(density(mod_dat$wl_initial_cm, na.rm = TRUE)$y)],
                        b0 = tr_params$par[["b0"]],
                        b1 = tr_params$par[["b1"]],
                        b2 = tr_params$par[["b2"]])]

mod_dat[, wl_hat_annual := quad_curve_wl(PET = pet_cm,
                                         P = fcoalesce(idw_precip_cm, rain_cm),
                                         M = melt_cm,
                                         max.wl = density(mod_dat$wl_initial_cm, na.rm = TRUE)$x[which.max(density(mod_dat$wl_initial_cm, na.rm = TRUE)$y)],
                                         b0 = tr_params$par[["b0"]],
                                         b1 = tr_params$par[["b1"]],
                                         b2 = tr_params$par[["b2"]]),
        by = .(water_year)]

ggplot(mod_dat,
       aes(x = sample_date)) +
  geom_line(aes(y = wl_initial_cm),
            col = "gray40") +
  geom_line(aes(y = wl_hat_cm)) +
  geom_line(aes(y = wl_hat_annual),
            linetype = 'dotted',
            color = 'blue') +
  facet_wrap(~water_year,
             scales = "free")

ggplot(mod_dat,
       aes(x = sample_date)) +
  geom_line(aes(y = ytd_water_balance_cm),
            col = "gray40") +
  geom_line(aes(y = wa_hat_cm)) +
  facet_wrap(~water_year,
             scales = "free")

ggplot(mod_dat,
       aes(x = sample_date)) +
  geom_line(aes(y = gradient)) +
  facet_wrap(~water_year,
             scales = "free")

ggplot(mod_dat,
       aes(x = sample_date)) +
  geom_line(aes(y = g_hat-q_hat)) +
  facet_wrap(~water_year,
             scales = "free")

tr_params
hydroGOF::gof(sim = mod_dat[!is.na(wl_initial_cm), wl_hat_cm], 
              obs = mod_dat[!is.na(wl_initial_cm), wl_initial_cm])


# Fit All Models ----------------------------------------------------------

site_info <- 
  fread(text = 
          "site,lon,lat,study,treatment,treatment_date,morphology
             006,-89.19036,46.32392,planting,Girdle,2013-10-31,flowthrough
             009,-89.71207,46.41624,eco,Ash Cut,2013-10-31,threshold
             053,-89.53043,46.61028,pws,Ash Cut,2014-10-31,flowthrough
             077,-88.86684,46.71778,eco,Ash Cut,2013-10-31,threshold
             111,-89.71369,46.57636,planting,Control,2013-10-31,threshold
             113,-89.80990,46.35090,pws,Control,2014-10-31,flowthrough
             119,-89.59433,46.28794,eco,Girdle,2013-10-31,closed
             135,-88.87838,46.75382,eco,Control,2013-10-31,flowthrough
             139,-88.95577,46.59619,planting,Ash Cut,2013-10-31,threshold
             140,-89.61962,46.43279,eco,Girdle,2013-10-31,flowthrough
             151,-88.89502,46.67307,eco,Girdle,2013-10-31,threshold
             152,-89.71468,46.41277,eco,Control,2013-10-31,flowthrough
             156,-89.58431,46.28327,eco,Ash Cut,2013-10-31,closed
             157,-89.58479,46.28173,eco,Control,2013-10-31,threshold",
        select = c(site = "character",
                   lon = "numeric",
                   lat = "numeric",
                   study = "character",
                   treatment = "character",
                   treatment_date = "Date",
                   morphology = "character"),
        key = "site")


full_dat <- 
  water_budget[CJ(site = unique(water_budget[, site]), 
                  sample_date = kenton[between(water_year, 2012, 2019), sample_date])
  ][kenton[between(water_year, 2012, 2019)],
    .(sample_date, 
      site, 
      wl_initial_cm, 
      idw_precip_cm,
      water_year = i.water_year,
      pet_cm = -i.pet_cm,
      rain_cm = i.rain_cm,
      melt_cm = i.melt_cm,
      total_input_cm = i.total_input_cm,
      tmax_c = i.tmax_c,
      tmin_c = i.tmin_c,
      tmean_c = i.tmean_c),
    on = "sample_date"]

setkey(full_dat, site, sample_date)

# Set Site status
full_dat[site_info,
         site_status := fifelse(sample_date > i.treatment_date & i.treatment != "Control",
                                "Treated",
                                "Control"),
         on = "site"]

# Drop Leap year
full_dat <- 
  full_dat[format(sample_date, "%m%d") != "0229"]

# Drop anomalous melt
full_dat[between(melt_cm, 0, 0.1) & month(sample_date) %in% 6:9, 
        `:=`(total_input_cm = total_input_cm - melt_cm,
             melt_cm = 0)]

# Drop anomalous pet
full_dat[tmax_c <= 0,
        pet_cm := 0]

# Calculate YTD after dropping above points
full_dat[, 
        `:=`(ytd_pet_cm = cumsum(pet_cm),
             ytd_p_cm = cumsum(rain_cm),
             ytd_m_cm = cumsum(melt_cm),
             ytd_input_cm = cumsum(rain_cm + melt_cm)),
        by = .(water_year)]

full_dat[, 
        `:=`(ytd_water_balance_cm = ytd_pet_cm + ytd_p_cm + ytd_m_cm,
             water_balance_cm = cumsum(pet_cm + rain_cm + melt_cm))]

full_dat[,
         max_wl_cm := possibly(function(x){density(x, na.rm = TRUE)$x[which.max(density(x, na.rm = TRUE)$y)]}, 
                               NA_real_)(wl_initial_cm),
         by = .(site, site_status)]

mod_dat <- 
  full_dat[, .(max_wl_cm = first(max_wl_cm),
               training_dat = list(.SD[water_year == ifelse(.BY[[2]] == "Control", 2012, 2015)]),
               valid_dat = list(.SD[water_year != ifelse(.BY[[2]] == "Control", 2012, 2015)])),
           by = .(site, site_status)
           ][map_dbl(training_dat, ~sum(!is.na(.x$wl_initial_cm))) > 0]

mod_dat[, params := map2(training_dat, max_wl_cm,
                        ~possibly(quad_opt, otherwise = NULL)(wobs = .x$wl_initial_cm, 
                                  met = .x[, .(pet_cm, rain_cm, melt_cm)],
                                  max.wl = .y,
                                  control = list(fnscale = -1)))]

mod_dat[, training_nse := map_dbl(params, pluck, "value")]

mod_dat[, valid_dat := pmap(list(valid_dat, params, max_wl_cm),
                            function(valid_dat, params, max_wl_cm){
                              valid_dat[, c("wl_hat_cm", "wa_hat_cm", "gradient", "q_hat", "g_hat", "max_wa") := 
                                        quad_curve_wl(PET = pet_cm,
                                                      P = fcoalesce(idw_precip_cm, rain_cm), 
                                                      M = melt_cm, 
                                                      max.wl = first(max_wl_cm),
                                                      b0 = params$par[["b0"]],
                                                      b1 = params$par[["b1"]],
                                                      b2 = params$par[["b2"]])]
                            })]

mod_dat[, validation_nse := map_dbl(valid_dat,
                                    ~hydroGOF::NSE(sim = .x[!is.na(wl_initial_cm), wl_hat_cm], 
                                                  obs = .x[!is.na(wl_initial_cm), wl_initial_cm]))]

valid_dat <- 
  mod_dat[, valid_dat[[1]],
          by = .(site, site_status)]

ggplot(valid_dat[water_year == 2017],
       aes(x = sample_date)) + 
  geom_line(aes(y = wl_hat_cm)) + 
  geom_line(aes(y = wl_initial_cm),
            color = "gray40") +
  facet_wrap(~site + site_status,
             scales = "free")

# Unused Curves -----------------------------------------------------------

pow <- 
  function(wa, b0, b1, b2){
    b0 + b1*(wa) - b2**(wa)
  }

pow_prime <- 
  function(wa, b1, b2){
    b1 - b2^(wa) * log(b2)
  }

alt_gompertz <- 
  function(wa, b, c, d, e){
    c + (d-c)*(1-exp(-exp(b*(wa - e))))
  }

gompertz <- 
  function(wa, b, c, d, e){
    c + (d-c)*(exp(-exp(-b*(wa - e))))
  }

sig_prime <- 
  function(x, b, c, d, e){
    (d - c) * ((exp(-(exp(b * (x - e))))) * ((exp(b * (x - e))) * b))
  }
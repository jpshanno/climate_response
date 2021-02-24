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
  "152"

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
    b1 + b2 * 2 * pmin(0, wa)
  }

# wd is water deficit from previous year

# Initial Parameters:
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
max.wl <- 
  density(water_budget[site == SITE, na.omit(wl_initial_cm)])$x[which.max(density(water_budget[site == SITE, na.omit(wl_initial_cm)])$y)]
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

            # Water rise is P - PET, which fit better, but need to rethink my
            # justification for this, or try without it
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
    
    # This looks good
    # lplot(wl_hat ~ D); lines(WL ~ D, col = 'gray40')
    # lplot(wa_hat ~ D)
    # lplot(q_hat ~ D)
    # lplot(g_hat ~ D)
    # lplot(gradient ~ D)
    # hydroGOF::gof(wl_hat[!is.na(WL)], WL[!is.na(WL)])
    # lplot(c + cumsum(gradient * P + gradient * PET))
    # c + last(cumsum(gradient * P + gradient * PET))
    
    data.table(wl_hat, wa_hat, gradient, q_hat, g_hat, max_wa)
    
  }

# Optimize ----------------------------------------------------------------

quad_opt <- 
  function(start, ...){
    
    
    nse <- 
      function(params, wobs, met, max.wl){
        
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
    
    optim(par = start,
          fn = nse,
          gr = NULL,
          ...)
  }

# Initial Parameters from nlrob for 152, 2012
# Could make it self-start
test_opt <- 
  quad_opt(start = list(b0 = coef(init_mod)[["b0"]],
                        b1 = coef(init_mod)[["b1"]], 
                        b2 = coef(init_mod)[["b2"]]),
           wobs = dat$wl_initial_cm, 
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
  ][water_budget[site == "135",
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
  mod_dat[water_year %in% c(2012, 2013)]

tr_params <- 
  quad_opt(start = list(b0 = coef(init_mod)[["b0"]],
                        b1 = coef(init_mod)[["b1"]],
                        b2 = coef(init_mod)[["b2"]]),
           wobs = tr_dat$wl_initial_cm, 
           met = tr_dat[, .(pet_cm, rain_cm, melt_cm)],
           max.wl = density(mod_dat$wl_initial_cm, na.rm = TRUE)$x[which.max(density(mod_dat$wl_initial_cm, na.rm = TRUE)$y)],
           control = list(fnscale = -1))
  # pow_opt(start = list(b0 = b0, b1 = b1, b2 = b2, I = 0.5),
  #         wobs = tr_dat$wl_initial_cm, 
  #         met = tr_dat[, .(pet_cm, rain_cm, melt_cm)],
  #         control = list(fnscale = -1),
  #         method = "L-BFGS-B",
  #         lower = c(0, -Inf, 0, 0),
  #         upper = c(Inf, Inf, 2, Inf))

mod_dat[, c("wl_hat_cm", "wa_hat_cm", "gradient", "q_hat", "g_hat", "max_wa") := 
          quad_curve_wl(PET = pet_cm,
                        P = fcoalesce(idw_precip_cm, rain_cm), 
                        M = melt_cm, 
                        max.wl = density(mod_dat$wl_initial_cm, na.rm = TRUE)$x[which.max(density(mod_dat$wl_initial_cm, na.rm = TRUE)$y)],
                        b0 = tr_params$par[["b0"]],
                        b1 = tr_params$par[["b1"]],
                        b2 = tr_params$par[["b2"]])]

ggplot(mod_dat,
       aes(x = sample_date)) +
  geom_line(aes(y = wl_initial_cm),
            col = "gray40") +
  geom_line(aes(y = wl_hat_cm)) +
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







for(i in 1:n){
  if(i == i){
    lmax <- x[1]
    lmaxi <- 1
  }
}




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
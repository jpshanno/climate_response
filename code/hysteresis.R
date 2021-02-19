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

# Drawdown curves are best calculated from peak melt or peak water avilability
# depending on the site. * sites have water level issues to resolve
# 006 - WA
# 009 - WA
# 053 - WA *
# 077 - WA
# 111 - WA
# 113 - WA *
# 119 - WA *
# 135 - WA
# 139 - WA
# 140 - WA
# 151 - WA
# 152 - WA
# 156 - WA *
# 157 - melt

kenton <- 
  external_met[station_name == "kenton"]

dat <- 
  water_budget[site == "152" & sample_year == 2012,
               .(sample_date, wl_initial_cm, best_precip_cm)
  ][kenton[sample_year == 2012,
           .(sample_date, 
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
             tmean_c)],
    on = "sample_date"]

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
    by = .(year(sample_date))]

# Fit this loop on a single site/year, then possibly attempt to optimize for
# m & b once you know the loop works
# f(PET) <- b - m**(-ytdPET)
# for t in 1:N
# ytdPET[t] <- 
#  ytdPET[t-1] + PET[t]
# If [P[t] == 0] 
# Then
#   WL[t] <- WL[t-1] - PET[t] * f`(PET[t])
# Else
#   WL[t] <- WL[t-1] + P[t] * f`(PET[t])
#   ytdPET[t] <- ytdPET[t-1] - (?P[t] * f`(PET[t]) or just P[t] possible)

plot(wl_initial_cm ~ ytd_pet_cm,
     data = dat, type = 'l')

# This will give a starting point, but it will likely be a slightly high estimate
# of b
nls(wl_initial_cm ~ c - b**(ytd_pet_cm), 
    data = dat[1:which.min(ytd_input_cm + ytd_pet_cm)], 
    na.action = na.exclude,
    start = list(c = 20, b = 0.8))



f_prime <- 
  function(x, b = b){
    -b**(x)*log(b)
  }

plot(wl_initial_cm ~ ytd_pet_cm,
     data = dat, type = 'l')
lines(na.omit(dat$wl_initial_cm)[1] + 
        cumsum(f_prime(dat[, ytd_pet_cm], 0.94)*
                 dat[order(ytd_pet_cm), pet_cm]) ~ dat[, ytd_pet_cm], 
      lty = 2)


# Power Curve -------------------------------------------------------------

pow <- 
  function(wa, b0, b1, b2){
    b0 + b1*(wa) - b2**(wa)
  }

pow_prime <- 
  function(wa, b1, b2){
    b1 - b2^(wa) * log(b2)
  }

# wd is water deficit from previous year

# Initial Parameters:
robustbase::nlrob(wl_initial_cm ~ pow(water_balance, 
                                      b0 = b0, b1 = b1, b2 = b2),
                  dat[1:which.min(wl_initial_cm), .(wl_initial_cm, water_balance = wd + ytd_p_cm + ytd_pet_cm + ytd_m_cm)],
                  na.action = na.exclude,
                  maxit = 50,
                  start = list(b0 = 20, b1 = 1, b2 = 0.9))

# b0 <- 24.9487008
# b1 <- 0.7650523
# b2 <- 0.8731924

power_curve_wl <- 
  function(PET, P, M, I, b0, b1, b2, wd){
    
    n <- 
      length(PET)
    
    wa_hat <- wl_hat <- gradient <- 
      numeric(n)
    
    wl_hat[1] <- b0
    
    wa_hat[1] <- wd
    
    # ind_min <-
    #   which.min(cumsum(P + PET))
    
    for(t in 2:n){
      wa_hat[t] <- 
        wa_hat[t-1] + PET[t] + M[t] + P[t]
      
      gradient[t] <- 
        pow_prime(wa_hat[t], b1 = b1, b2 = b2)
      
      if(P[t] <= abs(PET[t])){
        wl_hat[t] <-
          wl_hat[t-1] + (PET[t]) * gradient[t]
      } else {
        wl_hat[t] <-
          wl_hat[t-1] + (P[t] - I*PET[t]) * gradient[t]

        # if(t >= ind_min){
        #   wa_hat[t] <-
        #     pmin(0, YTD[t] + P[t])
        # }
      }
      
      wl_hat[t] <-
        wl_hat[t] + M[t]
      
      if(wl_hat[t-1] > b0){
        wl_hat[t] <- 
          wl_hat[t] - (wl_hat[t-1] - b0)/2
      }
      
      # Reset YTD
      # if(n %% 365 == 0){
      #   YTD[t] <- 
      #     c + (cumsum(gradient[(t-365):t] * (P[(t-365):t] + PET[(t-365):t])))
      # }
    }
    
    # This looks good
    # lplot(wl_hat ~ D); lines(WL ~ D, col = 'gray40')
    # hydroGOF::gof(wl_hat[!is.na(WL)], WL[!is.na(WL)])
    # lplot(c + cumsum(gradient * P + gradient * PET))
    # c + last(cumsum(gradient * P + gradient * PET))
    
    data.table(wl_hat, wa_hat, gradient)
    
  }

# tr_pet <- 

lines(power_curve_wl(PET, P, M, 0.25, 20, 0.932, -10.36) ~ D, col = 'green', lty = 3)



PET2013 <- kenton[sample_year == 2014, pet_cm]
P2013 <- kenton[sample_year == 2014, rain_cm]
M2013 <- kenton[sample_year == 2014, melt_cm]
D2013 <- kenton[sample_year == 2014, sample_date]
WD2012 <- -kenton[sample_year == 2013, pmin(0, last(water_availability_cm) - first(water_availability_cm))]

lplot(wl_initial_cm ~ sample_date,
      data = water_budget[site == "135" & sample_year == 2014], 
      col = 'gray40')
lines(wl_level(PET = PET2013, 
               P = P2013, 
               M = M2013, 
               I = 0.28, 
               int = 20, 
               base = 0.93, 
               i.YTD = WD2012) ~ D2013)

ggplot(water_budget[site == "135"]) + 
  aes(x = ytd_pet_cm, y = wl_initial_cm, color = factor(sample_year)) + geom_point()



mle2(wl_initial_cm ~ dnorm(wl_level(PET = pet_cm, P = rain_cm, M = melt_cm,
                                 I = Bi, int = 20, base = B1, i.YTD = 17.6)),
     start = list(Bi = 0.25, B1 = 0.935),
     data = dat[!is.na(wl_initial_cm)]) %>% summary()


LL <- 
  function(y, pet, p, m, wd, i, b0, b1, wd){
    -sum(dt(y - wl_level(PET = pet, P = p, M = m,
                         I = i, int = b0, base = b1, 
                         i.YTD = wd),
            df = 3,
            log = TRUE), na.rm = TRUE)
  }

with(dat, LL(wl_initial_cm, pet_cm, rain_cm, melt_cm, 17.6, 0.25, 20, 0.935, 17.6))

-sum(dt(dat$wl_initial_cm - wl_level(PET = dat$pet_cm, P = dat$rain_cm, M = dat$melt_cm,
                                     I = 0.25, int = 20, base = 0.935, i.YTD = 17.6),
        df = 3,
        log = TRUE), na.rm = TRUE)

mle2(LL(y = dat$wl_initial_cm, 
        pet = dat$pet_cm, 
        p = dat$rain_cm, 
        m = dat$melt_cm, 
        i = interception, 
        b0 = 0, 
        b1 = base,
        wd = 17.6),
     start = c(interception = 0.25, B0 = 20, base = 0.935))


# Sigmoidal Curve ---------------------------------------------------------

# robustbase::nlrob(wl_initial_cm ~ c + (d-c)*(1-exp(-exp(b*((ytd_p_cm + ytd_pet_cm + ytd_m_cm) - e)))),
#                   dat[1:278],
#                   na.action = na.exclude,
#                   maxit = 50,
#                   start = list(c = 20, d = -70, b = -0.25, e = -25))
# ggplot(dat, aes(x = ytd_p_cm + ytd_pet_cm, y = wl_initial_cm)) + geom_path() +
#   geom_function(fun = ~.x*sig_prime(.x, -0.1599038, 17.8725671, -127.2214549, -21.2614489), col = 'blue')

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

# From NLS on 2012 at site == 135 without melt (as opposed to the commented section above)
uasym <- 17.637899; lasym <- -126.562612; inflection <- -33.180073; rate <- -0.161178

plot(wl_initial_cm ~ ytd_pet_cm,
     data = dat, type = 'l')
lines(na.omit(dat$wl_initial_cm)[1] + 
        cumsum(f_prime(dat[, ytd_pet_cm], 0.94)*
                 dat[order(ytd_pet_cm), pet_cm]) ~ dat[, ytd_pet_cm], 
      lty = 2)

# Power Curve
sigmoid_wl <- 
  function(PET, P, M, uasym, lasym, rate, inflection, wd){
    
    n <- 
      length(PET)
    
    wl_hat <- gradient <- wa_hat <- 
      numeric(n)
    
    wl_hat[1] <- uasym
    
    wa_hat[1] <- wd
    
    ind_min <-
      which.min(cumsum(P + PET))
    
    for(t in 2:n){
      wa_hat[t] <- 
        wa_hat[t-1] + PET[t] + M[t] + P[t]
      
      gradient[t] <- 
        sig_prime(wa_hat[t], b = rate, c = uasym, d = lasym, e = inflection)
      
      if(P[t] == 0){
        wl_hat[t] <-
          wl_hat[t-1] + (PET[t]) * gradient[t]
      } else {
        wl_hat[t] <-
          wl_hat[t-1] + (P[t] + PET[t]) * gradient[t]
        
        # wa_hat[t] <- 
        #   wa_hat[t] + P[t]
        # 
        if(t >= ind_min){
          wl_hat[t] <-
            wl_hat[t] - (PET[t]) * gradient[t]
        }
      }
      
      wl_hat[t] <-
        wl_hat[t] + M[t]
      
      if(wl_hat[t-1] > uasym){
        wl_hat[t] <- 
          wl_hat[t] - (wl_hat[t-1] - uasym)/2
        
        # wa_hat[t] <- 
        #   wa_hat[t] - (wl_hat[t-1] - uasym)/2
      }
      
    }
    
    # This looks good
    # lplot(wl_hat ~ D); lines(WL ~ D, col = 'gray40')
    # hydroGOF::KGE(wl_hat[!is.na(WL)], WL[!is.na(WL)])
    # lplot(c + cumsum(gradient * P + gradient * PET))
    # c + last(cumsum(gradient * P + gradient * PET))
    
    wl_hat
    
  }

sig_wl <- 
  sigmoid_wl(PET = PET,
               P = P, 
               M = M, 
               uasym = uasym, 
               lasym = lasym, 
               rate = rate,
               inflection = inflection,
               wd = -10)

lplot(sig_wl ~ D); lines(WL ~ D, col = 'gray40'); lines(wl_hat ~ D, col = 'blue')
# lplot(sig_wl$gradient ~ D)

hydroGOF::gof(sig_wl[!is.na(dat$wl_initial_cm)], dat[!is.na(wl_initial_cm), wl_initial_cm])
hydroGOF::NSE(sig_wl[!is.na(dat$wl_initial_cm)], dat[!is.na(wl_initial_cm), wl_initial_cm])

# Optimize ----------------------------------------------------------------


sig_max_nse <- 
  function(x, wobs, met){
    
    wl_hat <- 
      sigmoid_wl(PET = met$pet_cm,
                 P = met$rain_cm, 
                 M = met$melt_cm, 
                 uasym = x[1], 
                 lasym = x[2], 
                 rate = x[3],
                 inflection = x[4],
                 wd = -10)
      
    hydroGOF::NSE(wl_hat[!is.na(wobs)], wobs[!is.na(wobs)])
    
  }

optim(par = c(uasym, lasym, rate, inflection),
      fn = sig_max_nse,
      gr = NULL,
      wobs = dat$wl_initial_cm, 
      met = dat[, .(pet_cm, rain_cm, melt_cm)],
      control = list(fnscale = -1))

lplot(sigmoid_wl(PET = PET,
                   P = P, 
                   M = M, 
                   uasym = 16.4034215, 
                   lasym = -121.8338205, 
                   rate = -0.1449116,
                   inflection = -31.7902495,
                   wd = -10) ~ D); lines(WL ~ D, col = 'gray40')

pow_max_nse <- 
  function(params, wobs, met){
    
    wl_hat <- 
      power_curve_wl(PET = met$pet_cm,
                     P = met$rain_cm, 
                     M = met$melt_cm, 
                     b0 = params[["b0"]],
                     b1 = params[["b1"]],
                     b2 = params[["b2"]],
                     I = params[["I"]],
                     wd = -10)$wl_hat
    
    hydroGOF::NSE(wl_hat[!is.na(wobs)], wobs[!is.na(wobs)])
    
  }

optim(par = list(b0 = b0, b1 = b1, b2 = b2, I = 0.5),
      fn = pow_max_nse,
      gr = NULL,
      wobs = dat$wl_initial_cm, 
      met = dat[, .(pet_cm, rain_cm, melt_cm)],
      control = list(fnscale = -1),
      method = "L-BFGS-B",
      lower = c(-Inf, -Inf, 0.7, 0))

lplot(power_curve_wl(PET = PET,
                     P = P, 
                     M = M, 
                     I = 0.93788,
                     b0 = 24.7347820,
                     b1 = 2.4423119,
                     b2 = 0.8784166,
                     wd = -10) ~ D); lines(WL ~ D, col = 'gray40')

# Optimized coefficients from testing with 135 & 2012
cbind(
  hydroGOF::gof(sigmoid_wl(PET = PET,
                           P = P, 
                           M = M, 
                           uasym = 16.4034215, 
                           lasym = -121.8338205, 
                           rate = -0.1449116,
                           inflection = -31.7902495,
                           wd = -10)[!is.na(WL)],
                WL[!is.na(WL)]),
  hydroGOF::gof(power_curve_wl(PET = PET,
                               P = P, 
                               M = M, 
                               I = 0.93788,
                               b0 = 24.7347820,
                               b1 = 2.4423119,
                               b2 = 0.8784166,
                               wd = -10)[!is.na(WL)],
                WL[!is.na(WL)])
)

pow_opt <- 
  function(start, wobs, met, ...){
    
    
    nse <- 
      function(params, wobs, met){
        
        wl_hat <- 
          power_curve_wl(PET = met$pet_cm,
                         P = met$rain_cm, 
                         M = met$melt_cm, 
                         b0 = params[["b0"]],
                         b1 = params[["b1"]],
                         b2 = params[["b2"]],
                         I = params[["I"]],
                         wd = -10)$wl_hat
        
        hydroGOF::NSE(wl_hat[!is.na(wobs)], wobs[!is.na(wobs)])
        
      }
    
    optim(par = start,
          fn = nse,
          gr = NULL,
          wobs = wobs, 
          met = met,
          ...)
  }

pow_opt(start = list(b0 = b0, b1 = b1, b2 = b2, I = 0.5),
        wobs = dat$wl_initial_cm, 
        met = dat[, .(pet_cm, rain_cm, melt_cm)],
        control = list(fnscale = -1),
        method = "L-BFGS-B",
        lower = c(-Inf, -Inf, 0.7, 0))

mod_dat <- 
  kenton[sample_year %in% unique(water_budget[site == "135", sample_year]),
         .(sample_date, 
           sample_year,
           pet_cm = -pet_cm,
           rain_cm,
           melt_cm,
           total_input_cm,
           tmax_c,
           tmin_c,
           tmean_c)
  ][water_budget[site == "135",
                 .(sample_date, wl_initial_cm, best_precip_cm)],
    `:=`(wl_initial_cm = i.wl_initial_cm,
         best_precip_cm = i.best_precip_cm),
    on = "sample_date"]

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
       by = .(year(sample_date))]

mod_dat[, 
        `:=`(ytd_water_balance_cm = ytd_pet_cm + ytd_p_cm + ytd_m_cm,
             water_balance_cm = cumsum(pet_cm + rain_cm + melt_cm))]

# Split Data
tr_dat <- 
  mod_dat[sample_year %in% c(2012, 2015)]

tr_params <- 
  pow_opt(start = list(b0 = b0, b1 = b1, b2 = b2, I = 0.5),
          wobs = tr_dat$wl_initial_cm, 
          met = tr_dat[, .(pet_cm, rain_cm, melt_cm)],
          control = list(fnscale = -1),
          method = "L-BFGS-B",
          lower = c(0, -Inf, 0, 0),
          upper = c(Inf, Inf, 2, Inf))

mod_dat[, c("wl_hat_cm", "wa_hat_cm", "gradient") := power_curve_wl(PET = pet_cm,
                                      P = rain_cm, 
                                      M = melt_cm, 
                                      I = tr_params$par[["I"]],
                                      b0 = tr_params$par[["b0"]],
                                      b1 = tr_params$par[["b1"]],
                                      b2 = tr_params$par[["b2"]],
                                      wd = -10)]

ggplot(mod_dat,
       aes(x = sample_date)) +
  geom_line(aes(y = wl_initial_cm),
            col = "gray40") +
  geom_line(aes(y = wl_hat_cm)) +
  facet_wrap(~sample_year,
             scales = "free")

hydroGOF::gof(sim = mod_dat[!is.na(wl_initial_cm), wl_hat_cm], 
              obs = mod_dat[!is.na(wl_initial_cm), wl_initial_cm])

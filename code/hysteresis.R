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

SITE <- 
  "135"
align_var <- 
  # "melt_cm"
  "water_availability_cm"

drying_dates <- 
  kenton[, .SD[
    # Index of maximum water availability (melt period)
    which.max(get(align_var)):
      # Index of max WA + Index of min WA that occurs after max WA
      (which.max(get(align_var)) + which.min(.SD[which.max(get(align_var)):.N, water_availability_cm]))], 
    by = .(sample_year)][, sample_date]

drying_curves <- 
  water_budget[site == SITE & sample_date %in% drying_dates]
# , 
#                .SD[
#                  # Index of maximum water availability (melt period)
#                  which.max(water_availability_cm):
#                    # Index of max WA + Index of min WA that occurs after max WA
#                    (which.max(water_availability_cm) + which.min(.SD[which.max(water_availability_cm):.N, water_availability_cm]))], 
#                by = .(sample_year)]

# The offset to align the curves lay in the selection of the start period for
# the drawdown. When I offset each curve to the water availability at the start
# of the drawdown, the curves sync much better. Still some work to identify the
# best start ands top points for the drawdown curve to get the best curve fit.
# It also seems like just using water avilability will work too
drying_curves <- 
  drying_curves[drying_curves[, .SD[1, .(sample_date, water_availability_cm)], 
                              by = .(sample_year)], 
                on = "sample_year", 
                .(wl_initial_cm, sample_year, sample_date, site_status, pet_cm, Ds_cm,
                  water_availability_cm = water_availability_cm - i.water_availability_cm)]

# Power Curve
# drying_mod <-
#   robustbase::nlrob(wl_initial_cm ~ b + m1*water_availability_cm - m2**water_availability_cm,
#                     start = list(b = 20,
#                                  m1 = 1,
#                                  m2 = 0.86),
#                     na.action = na.exclude,
#                     maxit = 50,
#                     data = drying_curves[sample_year != 2015])
# 
# hyst_function <-
#   deriv(bquote(.(b) + .(m1)*water.availability - .(m2)**water.availability,
#                where = list(b = coef(drying_mod)[["b"]],
#                             m1 = coef(drying_mod)[["m1"]],
#                             m2 = coef(drying_mod)[["m2"]])),
#         namevec = "water.availability",
#         func = TRUE)

# Sigmoid Curve

# 2015 does not have a good identification of drawdown curve
drying_mod <-
  robustbase::nlrob(wl_initial_cm ~ c + (d-c)*(1-exp(-exp(b*(water_availability_cm - e)))),
                    drying_curves[sample_year != 2015],
                    na.action = na.exclude,
                    start = list(c = 20, d = -70, b = -0.25, e = -25))

hyst_function <-
  deriv(bquote(.(c) + (.(d) - .(c))*(1-exp(-exp(.(b)*(water.availability - .(e))))),
               where = as.list(coef(drying_mod))),
        namevec = "water.availability",
        func = TRUE)


ggplot(drying_curves, 
       aes(x = water_availability_cm, 
           y = wl_initial_cm)) + 
  geom_point(aes(color = factor(sample_year))) + 
  geom_function(fun = hyst_function)

D <- kenton[between(sample_year, 2007, 2019), sample_date]
PET <- kenton[between(sample_year, 2007, 2019), pet_cm]
P <- kenton[between(sample_year, 2007, 2019), total_input_cm]
R <- kenton[between(sample_year, 2007, 2019), rain_cm]
M <- kenton[between(sample_year, 2007, 2019), melt_cm]

# Drawdown should be controlled by this function:
plot(x = drying_curves[sample_year != 2015, wl_initial_cm],
     y = attr(hyst_function(drying_curves[sample_year!= 2015, water_availability_cm]), "gradient")[,1],
     main = "Drying Curve Gradient",
     xlab = "Water Level (cm)",
     ylab = expression(delta/delta~X))

gradient_mod <- 
  mcp::mcp(list(gradient ~ wl_initial_cm, ~ 0 + wl_initial_cm),
      data = drying_curves[sample_year != 2015 & !is.na(wl_initial_cm), 
                           .(wl_initial_cm,
                             gradient = attr(hyst_function(water_availability_cm), "gradient")[,1])],
      cores = 1)

plot(gradient_mod, 
     q_predict = TRUE)

grad_coefs <- 
  mcp::fixef(gradient_mod)

# gradient_function <- 
#   as.function(list(water.level = NULL,
#                    bquote(.(int_1) + .(cp_1) * .(wl_initial_cm_1) + (water.level - .(cp_1)) * .(wl_initial_cm_2),
#                           where = set_names(as.list(grad_coefs$mean), grad_coefs$name))))

gradient_function <-
  as.function(list(water.level = NULL,
                   bquote((.(cp_1) >= water.level) * (.(int_1) + water.level * .(wl_initial_cm_1)) +
                            (.(cp_1) < water.level) * (.(int_1) + .(cp_1) * .(wl_initial_cm_1) + (water.level - .(cp_1)) * .(wl_initial_cm_2)),
                          where = set_names(as.list(grad_coefs$mean), grad_coefs$name))))

for(t in 1:(length(P) - 1)){
  
  # Set up on first iteration
  if(t == 1){
    
    # X is water calculated water availability
    # Y is predicted water level
    # y_hat is water level output from model
    # gradient is the derivative of the model slope 
    # dWA is the change in X
    X <- Y <- y_hat <- gradient <- dX <- 
      numeric(length(P))
    
    # Start the predictions at water yield == 0
    Y[t] <- 
      hyst_function(X[t])
    
    # gradient[t] <- 
    #   gradient_function(Y[t])
    # 
    # Y[t+1] <-
    #   X[t] + (gradient[t]*P[t] - PET[t])
    
    # Move onto t = 2
    # next
  }
  
  # Save slope of gradient at t
  gradient[t] <- 
    gradient_function(Y[t])
  
  # Calculate dX at step t
  dX[t] <- 
    P[t] - PET[t]
  
  X[t + 1] <- 
    X[t] + dX[t]
  
  if(format(D[t], "%m%d") == "1231"){
    X[t + 1] <- 
      (P[t+1] - PET[t+1]) - max(cumsum(P[(t+1):(t+365)] - PET[(t+1):(t+365)]))
    # e + log(-log(-((Y[t]-c)/(d-c) - 1)))/b
  }
  
  
  if(dX[t] < 0){
    Y[t+1] <-
      # hyst_function(X[t])
      # Y[t] - pmax(PET[t], PET[t] * gradient[t]) #+ pmin(dX[t], dX[t] * gradient[t])
      Y[t] + pmin(0, (hyst_function(X[t]) - ifelse(t > 1, hyst_function(X[t-1]), 0)))
  } else {
    Y[t+1] <-
      Y[t] + gradient[t] * P[t]
  }
 
}


ggplot(water_budget[site == SITE], 
       aes(y = wl_initial_cm, x = sample_date)) +
  geom_line(color = 'gray40') +
  geom_line(data = data.table(sample_date = D, 
                              wl_initial_cm = Y,
                              sample_year = year(D))[sample_year %between% c(2012, 2019)], 
            color = "blue",
            linetype = "dotted") + 
  facet_wrap(~sample_year, scales = "free")

ggplot(kenton[sample_year %between% c(2012, 2019)], 
       aes(y = ytd_water_availability_cm, x = sample_date)) +
  geom_line(color = 'gray40') +
  geom_line(data = data.table(sample_date = D, 
                              ytd_water_availability_cm = X,
                              sample_year = year(D))[sample_year %between% c(2012, 2019)], 
            color = "blue",
            linetype = "dotted") + 
  facet_wrap(~sample_year, scales = "free_x")











# Hyst Optimization -------------------------------------------------------

# Drawdown is power function for ytd_pet_cm
# Recovery is the slope of the power function when rain falls

# Set melt < 0.01 == 0
# Need to investigate CN parameters a bit more


dat <- 
  water_budget[site == "135" & sample_year == 2012,
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

# Power Curve
power_curve_wl <- 
  function(PET, P, M, I, c, b, wd){
    
    n <- 
      length(PET)
    
    YTD <- wl_hat <- gradient <- 
      numeric(n)
    
    wl_hat[1] <- c
    
    YTD[1] <- wd
    
    ind_min <-
      which.min(cumsum(P + PET))
    
    for(t in 2:n){
      YTD[t] <- 
          YTD[t-1] + PET[t] + M[t]
      
      gradient[t] <- 
        f_prime(YTD[t], b)
      
      if(P[t] <= I){
        wl_hat[t] <-
          wl_hat[t-1] + (PET[t]) * gradient[t]
      } else {
        wl_hat[t] <-
          wl_hat[t-1] + (P[t] - PET[t]) * gradient[t]

        if(t >= ind_min){
          YTD[t] <-
            pmin(0, YTD[t] + P[t])
        }
      }
      
      wl_hat[t] <-
        wl_hat[t] + M[t]
      
      if(wl_hat[t-1] > int){
        wl_hat[t] <- 
          wl_hat[t] - (wl_hat[t-1] - int)/2
      }
      
      # Reset YTD
      if(n %% 365 == 0){
        YTD[t] <- 
          c + (cumsum(gradient[(t-365):t] * (P[(t-365):t] + PET[(t-365):t])))
      }
    }
    
    # This looks good
    # lplot(wl_hat ~ D); lines(WL ~ D, col = 'gray40')
    # hydroGOF::KGE(wl_hat[!is.na(WL)], WL[!is.na(WL)])
    # lplot(c + cumsum(gradient * P + gradient * PET))
    # c + last(cumsum(gradient * P + gradient * PET))
    
    wl_hat
    
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


max_nse <- 
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
      fn = max_nse,
      gr = NULL,
      wobs = dat$wl_initial_cm, 
      met = dat[, .(pet_cm, rain_cm, melt_cm)],
      control = list(fnscale = -1))



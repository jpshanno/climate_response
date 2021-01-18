# Data from climate change impacts

# WL is ARIMA(1,1,0), but Ds is Arima(0,0,0), so it's easier to just do the 
# difference and then model with lm or some other regular regression or brms
library(targets)
source("code/load_project.R")
library(TSA)


tar_load(water_budget)

dat <- 
  water_budget[site == "152" & sample_year == 2017]

# Prewhitening example to look at ccf of pet and wl_initial:
# https://online.stat.psu.edu/stat510/lesson/9/9.1
# pet has an ARIMA structure:
forecast::auto.arima(dat[!is.na(wl_initial_cm), pet_cm])

prewhite_wl <- 
  filter(x = dat[!is.na(wl_initial_cm), wl_initial_cm],
         filter = c(coef(forecast::auto.arima(dat[!is.na(wl_initial_cm), pet_cm]))),
         method = "con",
         sides = 1)

ccf(x = resid(forecast::auto.arima(dat[!is.na(wl_initial_cm), pet_cm])), 
    y = prewhite_wl,
    na.action = na.pass)

# precip has an ARIMA structure
forecast::auto.arima(dat[!is.na(wl_initial_cm), best_precip_cm], allowmean = FALSE)

# But so we can prewhiten the prewhiteded wl_initial using the ARIMA structure for precip
pre_prewhite_wl <- 
  filter(prewhite_wl,
         filter = coef(forecast::auto.arima(dat[!is.na(wl_initial_cm), best_precip_cm], allowmean = FALSE)),
         method = "con",
         sides = 1)

ccf(x = dat[!is.na(wl_initial_cm), best_precip_cm], 
    y = filter(prewhite_wl,
               filter = coef(forecast::auto.arima(dat[!is.na(wl_initial_cm), best_precip_cm], allowmean = FALSE)),
               method = "con",
               sides = 1),
    na.action = na.pass)

# My final models predict Ds which shows now ARIMA structure, because it is the
# differenced water level and water level has an ARIMA(0,1,0) structure after the
# effects of PET and precip are accounted for
auto.arima()

ccf(x = resid(forecast::auto.arima(dat[!is.na(wl_initial_cm), pet_cm])), 
    y = filter(x = dat[!is.na(wl_initial_cm), wl_initial_cm],
               filter = c(coef(forecast::auto.arima(dat[!is.na(wl_initial_cm), pet_cm]))),
               method = "con",
               sides = 1),
    na.action = na.pass)

  
dat[, c("prev_precip", "prev2_precip", "prev3_precip", "prev4_precip") := map(1:4, ~shift(best_precip_cm, .x))]
dat[, c("prev_pet", "prev2_pet", "prev3_pet", "prev4_pet") := map(1:4, ~shift(pet_cm, .x))]
dat[, c("diff_precip", "diff2_precip", "diff3_precip", "diff4_precip",
        "diff_pet", "diff2_pet", "diff3_pet", "diff4_pet") := map(.SD, ~c(0, diff(.x))),
    .SDcols = c("prev_precip", "prev2_precip", "prev3_precip", "prev4_precip",
                "prev_pet", "prev2_pet", "prev3_pet", "prev4_pet")]
dat[, prev_wl := shift(wl_initial_cm, 1)]
dat[, prev2_wl := shift(wl_initial_cm, 2)]
dat[, prev_ds := shift(Ds_cm, 1)]
dat[, prev2_ds := shift(Ds_cm, 2)]

big_lm <- 
  lm(wl_initial_cm ~ 0 + offset(prev_wl) + best_precip_cm + prev_precip + prev2_precip + prev3_precip + pet_cm + prev_pet, data = dat)

dat[, wl_lm := predict(big_lm, newdata = .SD)]


plot(wl_lm ~ doy, data = dat[!is.na(wl_initial_cm)], type = 'l')
lines(wl_initial_cm ~ doy, data = dat[!is.na(wl_initial_cm)], col = "blue")

ccf(dat[!is.na(wl_initial_cm), wl_initial_cm],
    dat[!is.na(wl_initial_cm), best_precip_cm])

ccf(dat[!is.na(wl_initial_cm), wl_initial_cm],
    dat[!is.na(wl_initial_cm), pet_cm])

mod_ts <- 
  arimax(dat[!is.na(wl_initial_cm), wl_initial_cm], 
         order = c(1,1,0)
        , xtransf = dat[!is.na(wl_initial_cm),
                        .(best_precip_cm)],
        transfer=list(c(0,0)))

dat[!is.na(wl_initial_cm), wl_ts := wl_initial_cm - as.numeric(resid(mod_ts))]

plot(wl_initial_cm ~ doy, data = dat[!is.na(wl_initial_cm)], col = 'gray20', type = 'l')
# lines(wl_lm ~ doy, data = dat[!is.na(wl_initial_cm)], col = "blue")
lines(wl_ts ~ doy, data = dat, col = "red", lty = "dashed")
lines(wl_ar ~ doy, data = dat, col = "green", lty = "dotted")

dat[!is.na(wl_initial_cm), 
    wl_ar := 
      # Last Value
      shift(wl_initial_cm, 1) - 
      # AR1 impact of differenced time series
      as.numeric(filter({(shift(wl_initial_cm, 2, 0) - shift(wl_initial_cm, 1, 0))},
                        coef(mod_ts)[["ar1"]],
                        "con",
                        sides = 1))]

# Intervention calculations adapted from https://rstudio-pubs-static.s3.amazonaws.com/274357_701459df7e59463cb3a67e7f7ce3c64c.html
dat[!is.na(wl_initial_cm), p_occ := best_precip_cm > 0 | shift(best_precip_cm, 1, 0) > 0]
dat[!is.na(wl_initial_cm), pp_occ := shift(best_precip_cm, 1, 0) > 0 | shift(best_precip_cm, 2, 0) > 0]
dat[!is.na(wl_initial_cm), p_sign := ifelse(sum(p_occ, pp_occ) == 2, -1, 1)]

# Still some errors on the intervention calcualtions. I think it has to do with
# periods where there are multip precip days in a row. The filters probably have
# to be nested or something

dat[!is.na(wl_initial_cm), 
    wl_pred := 
      # Last Value
      shift(wl_initial_cm, 1) - 
      # AR1 impact of differenced time series
      filter({(shift(wl_initial_cm, 2, 0) - shift(wl_initial_cm, 1, 0))},
             coef(mod_ts)[["ar1"]],
             "con",
             sides = 1) -
      # Sign flips for the points where p_occ and pp_occ are both true. Has to 
      # do with 2 negative signs
      coef(mod_ts)[["best_precip_cm-MA0"]] * (filter(best_precip_cm, c(-1, 0.1645), "con", sides = 1) - filter(shift(best_precip_cm, 1, 0), c(-1, 0.1645), "con", sides = 1))]
# coef(mod_ts)[["best_precip_cm-MA0"]] * ((shift(best_precip_cm, 1) - best_precip_cm) - (shift(best_precip_cm, 2) - shift(best_precip_cm, 1)) * coef(mod_ts)[["ar1"]])]
dat[, wl_pred := as.numeric(wl_pred)]

plot(wl_ts ~ wl_pred, data = dat); abline(0, 1)
plot(wl_ts ~ doy, data = dat, type = "l"); lines(wl_pred ~ doy, data = dat, col = "blue", lty = "dashed")
lm(wl_resid ~ best_precip_cm, 
   data = dat[, .(wl_resid = wl_ts - wl_pred,
                  best_precip_cm = ifelse(shift(best_precip_cm, 1) == 0 & 
                                            shift(best_precip_cm, 2) == 0 & 
                                            shift(best_precip_cm, 3) == 0 & 
                                            best_precip_cm > 0,
                                          shift(best_precip_cm, 0), 
                                          NA_real_))])

t0  0.4709
t1 -0.07747
t2 -0.3934
t3 -0.311
t4  0

mean P = 0.2664371
 
mod_ds <- 
  arimax(dat[!is.na(Ds_cm), Ds_cm], 
         order = c(1,0,0), 
         include.mean = FALSE,
         xtransf = dat[!is.na(Ds_cm),
                       best_precip_cm],
         transfer = list(c(1,0)))

dat[!is.na(Ds_cm), 
    ds_ts := as.numeric(fitted(mod_ds))]
dat[!is.na(Ds_cm), ds_ar := shift(Ds_cm, 1) * coef(mod_ds)[["ar1"]]]
dat[, ds_pred := 
      shift(Ds_cm, 1) * coef(mod_ds)[["ar1"]] +
      best_precip_cm * coef(mod_ds)[["T1-MA0"]]
      (best_precip_cm - shift(best_precip_cm, 1)) * coef(mod_ds)[["ar1"]] + filter(best_precip_cm,filter= coef(mod_ds)[["T1-AR1"]],method='recursive',side=1)*(coef(mod_ds)[["T1-MA0"]])]

plot(ds_ts ~ ds_pred, data = dat); abline(0, 1)
plot(ds_ts ~ Ds_cm, data = dat); abline(0, 1)
plot(ds_ts ~ doy, data = dat, type = "l"); lines(ds_pred ~ doy, data = dat, col = "blue", lty = "dashed")




library(MASS)
set.seed(1)
N = ts(mvrnorm(50, mu=c(0,0), Sigma=matrix(c(1,0.56,0.56,1), ncol=2), 
               empirical=TRUE), frequency=12)
f = Arima(N[,1], order=c(2,0,0), xreg=N[,2], include.constant=FALSE)
n=f$x-f$coef[3]*f$xreg #finds the AR error part of the model
n_l1=lag(n,-1)
n_l2=lag(n,-2)
n_fit=n_l1*f$coef[1]+f$coef[2]*n_l2 #finds the fitted AR errors
pred=n_fit+f$coef[3]*f$xreg #Gives correct fitted values  


mod_ts <- 
  auto.arima(dat$wl_initial_cm, 
             max.q = 0,
             seasonal = FALSE, 
             xreg = as.matrix(dat[, .(prev_precip, prev2_precip, prev3_precip, prev_pet, prev2_pet)]))

eta <- uschange[,"Consumption"] - regress_fit
by_hand_fits3 <- regress_fit + lm_ar_fit$coef["ar1"]*lag(as.vector(eta))
plot(by_hand_fits3, lm_ar_fit$fitted)
summary(by_hand_fits3 - lm_ar_fit$fitted)

coef_ts <- 
  coef(mod_ts)

regress_fit <- 
  dat[!is.na(wl_initial_cm), prev_precip] * coef_ts[["prev_precip"]] +
  dat[!is.na(wl_initial_cm), prev2_precip] * coef_ts[["prev2_precip"]] +
  dat[!is.na(wl_initial_cm), prev3_precip] * coef_ts[["prev3_precip"]] +
  dat[!is.na(wl_initial_cm), prev_pet] * coef_ts[["prev_pet"]] +
  dat[!is.na(wl_initial_cm), prev2_pet] * coef_ts[["prev2_pet"]]



by_hand_fits <- regress_fit + 

plot(by_hand_fits ~ fitted(mod_ts)); abline(0, 1)

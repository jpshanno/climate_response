source("code/load_project.R")
tar_load(external_met)
tar_load(daily_water_levels)

# Plot all water levels
ggplot(daily_water_levels, 
       aes(x = sample_date,
           y = wl_initial_cm)) + 
  geom_line() + 
  facet_grid(site ~ sample_year,
             scales = "free") + 
  scale_x_date(date_breaks = "1 month",
               date_labels = "%m",
               expand = expansion()) + 
  theme_bw()

# Plot all water availabilities

ggplot(external_met[external_met[format(sample_date, "%m%d") == "0415", 
                                 .(station_name, 
                                   sample_year,
                                   init_wa = water_availability_cm)], 
                    on = c("station_name", "sample_year")
                    ][sample_year %in% 2012:2020, 
                      .(sample_date, 
                        water_availability_cm = water_availability_cm - init_wa), 
                      by = .(station_name, sample_year)], 
       aes(x = sample_date, y = water_availability_cm)) + 
  geom_line() + 
  geom_hline(aes(yintercept = 0)) +
  facet_grid(station_name ~ sample_year, scales = "free") + 
  scale_x_date(date_breaks = "1 month", date_labels = "%m", expand = expansion()) + 
  theme_bw()

ggplot(external_met, 
       aes(x = sample_date, y = water_availability_cm)) + 
  geom_line() + 
  geom_point(data = external_met[format(sample_date, "%m%d") == "0415"], 
             color = "red") + 
  facet_wrap(~station_name, ncol = 1) + 
  scale_x_date(breaks = as.Date(paste0(2000:2020, "-04-15")))

# WA Trends

peak_dates <- 
  external_met[, .(sos = .SD[which.max(resid(lm(water_availability_cm ~ sample_date, 
                                                        data = .SD,
                                                        weights = c(10000, rep(1, .N-2), 10000)))), sample_date],
                   eos = as.Date(paste0(.BY[[2]], "-12-31")),
                   max_wa_cm = .SD[which.max(resid(lm(water_availability_cm ~ sample_date, 
                                                      data = .SD,
                                                      weights = c(10000, rep(1, .N-2), 10000)))), water_availability_cm]),
               by = .(station_name, sample_year)]

wa_trends <- 
  external_met[peak_dates, 
               on = c("station_name", "sample_year")
  ][sample_date == sos | sample_date == eos,
    .(trend_x = sos,
      trend_y = max_wa_cm,
      trend_adj = diff(water_availability_cm) / as.numeric(diff(sample_date, units = "days"))),
    by = .(station_name, sample_year)]

external_met[wa_trends,
             `:=`(trend_y = i.trend_y,
                  trend_adj = as.numeric(difftime(i.trend_x, sample_date, units = "days")) * i.trend_adj),
             on = c("station_name", "sample_year")]

external_met[, detrend_wa_cm := (water_availability_cm + trend_adj) - trend_y]

# Rather than use the maximum annual water water availability, calculate the max
# springtime water avaialbility as above and use that for normalizing. Need to
# simulate how that would look for long-term drying

external_met[peak_dates,
             seasonal_max_wa_cm := i.max_wa_cm,
             on = c("station_name", "sample_year")]

# Using the seasonal max approach from above (ie spring max water level) show 
# divergent curves in Ds/normalized wa because it allows for some positive 
# values in normalized_wa_cm

external_met[,
             normalized_wa_cm := (water_availability_cm - first(water_availability_cm)) - max(water_availability_cm - first(water_availability_cm)),
             by = .(station_name, sample_year)]

daily_water_levels[, Ds_continuous := diff_na(wl_initial_cm),
                   by = .(site)]

daily_water_levels[, Ds_cumulative := ytd_sum(Ds_continuous),
                   by = .(site)]

external_met[,
             ytd_pet_cm := cumsum(pet_cm),
             by = .(station_name, sample_year)]

daily_water_levels[external_met[station_name == "kenton"],
                   `:=`(normalized_wa_cm = i.normalized_wa_cm,
                        detrend_wa_cm = i.detrend_wa_cm,
                        water_availability_cm = i.water_availability_cm,
                        pet_cm = i.pet_cm,
                        ytd_pet_cm = i.ytd_pet_cm),
                   on = c("sample_date")]

daily_water_levels[daily_water_levels[!is.na(Ds_cumulative), 
                                      .(sos = min(sample_date)), 
                                      by = .(site, sample_year)], 
                   sos := i.sos,
                   on = c("site", "sample_year")]

daily_water_levels[sample_date >= sos, 
                   wa_cumulative := water_availability_cm - first(water_availability_cm),
                   by = .(site, sample_year)]

daily_water_levels[, esy_wa_emp := Ds_cumulative / normalized_wa_cm]

daily_water_levels[, threshold := coef(lm(wl_l1_cm ~ Ds_cumulative,
                                          data = .SD))[[1]],
                   by = .(site)]


# Big Bayesian Model ------------------------------------------------------

# library(brms)
options(mc.cores = 4)

# Look for relationship comparing PET:Drawdown on days with precip vs days without

# daily_water_levels[, `:=`(wl_l1_cm = shift(wl_initial_cm, 1),
#                           pet_l1_cm = shift(pet_cm, 1),
#                           precip_l1_cm = shift(best_precip_cm, 1),
#                           precip_l2_cm = shift(best_precip_cm, 2),
#                           wa_l1_cm = shift(normalized_wa_cm)),
#                    by = .(site, sample_year)]

brm_test <- 
  water_budget[site == "135" & sample_year %in% c(2012, 2013),
               .(wlObs = wl_initial_cm,
                 wlL1 = shift(wl_initial_cm, 1),
                 precip = l1_total_input_cm,
                 pet = l1_pet_cm,
                 melt = l1_melt_cm,
                 rain = l1_rain_cm,
                 wa = l1_ytd_water_availability_cm),
               by = .(sample_year)][!is.na(wlObs)]
  # daily_water_levels[site == "152" & sample_year %in% c(2012, 2013) & !is.na(wl_initial_cm), 
  #                    .(wlObs = wl_initial_cm,
  #                      wlL1 = wl_l1_cm,
  #                      pet = pet_l1_cm,
  #                      wa = wa_l1_cm,
  #                      precip = precip_l1_cm)]

# Prediction loops are made to do i-1 for drivers, so I'm not using the 
# pre-lagged predictors in val_dat. This needs to be switched but lagged met 
# predictors need to be calculated in external_met and then joined to water_balance
# so that there is not an NA at the beginning of each year of data.

val_dat <- 
  water_budget[site == "135" & sample_year %in% 2012:2014,
               .(sample_date, 
                 wlObs = wl_initial_cm)][min(which(!is.na(wlObs))):.N][
                   external_met[station_name == "kenton" & 
                                  sample_year %in% 2012:2014 & 
                                  sample_date >= as.Date("2012-04-01")], 
                   .(sample_date, 
                     sample_year = year(sample_date), 
                     wlObs, 
                     precip = total_input_cm, 
                     pet = pet_cm, 
                     melt = melt_cm,
                     rain = rain_cm,
                     wa = ytd_water_availability_cm,
                     raw_wa = water_availability_cm), 
                   on = c("sample_date")]
  # daily_water_levels[site == "152" & sample_year %in% c(2012:2014), 
  #                    .(sample_date, 
  #                      wlObs = wl_initial_cm,
  #                      wlL1 = wl_l1_cm,
  #                      pet = pet_cm,
  #                      wa = normalized_wa_cm,
  #                      precip = best_precip_cm)][min(which(!is.na(wlObs))):.N]

# To Do:
# - Get water levels to go above cp (probably need to add piecewise component to p)
# - Add quickflow in/out from preceding day's precip (lagged precip term)
# - Try to vary esy for P and PET again. Maybe having P when wl > cp not depend
# on esy will help set a more accurate esy(s). Or try to work, previously they
# ended up negatively afffecting the esy fit
# - Consider changing distributions for priors (i.e. lognormal, or skew_normal for
# some things)
# - Try making G dependent on normalized_wa_cm (need to experiment for decent priors)
# - Adjust cp from constant prior based on observed percentiles to a strict 
# distributional prior
# - Try out detrending water availability across years. Potentially just using
# the previous year's final ytd_wa as the starting value
# - Include interception term
# - Split PET into wet and dry days or try Drooger & Allen 2002
# - PET does not seem to cause drawdown when ~ -pet/esy > -1, Had really good
# training data fit when I dropped Bg and and then did piecewise regression on
# PET with Bpet + step(pet - cpPET) * (Mpet * pet / esy). It overestimated
# drawdown for the validation data because PET rose above cpPET much earlier
# than drawdown occurred
# - Try two stage Sy relationship. In original Sy plots, it decreases again 
# after the critical water level
# - Try adding pan evap 

# daily_water_levels[site == "152",
#                    .(pet_drawdown = -pet_cm / (Besy + Mesy^(wl_initial_cm - cp)),
#                      p_rise = best_precip_cm / (Besy + Mesy^(wl_initial_cm - cp)),
#                      ds = Ds_cm,
#                      pet_cm = -pet_cm,
#                      best_precip_cm,
#                      sample_date)] %>% 
#   ggplot(aes(x = pet_drawdown,
#              y = ds - p_rise,
#              color = ifelse(p_rise > 0, "wet", "dry"))) +
#   geom_point()

wb_form <- 
  bf(
    # Generic Water Budget
    wlObs ~ wlL1 + p + et + q + g,
    
    # Precip response scaled by ecosystem specific yield
    nlf(p ~ step(wlL1 - cp) * (precip / esy) + step(cp - wlL1) * (Mp * precip)),
    
    # PET response scaled by ecosystem specific yield and 0 when annual water
    # availability (YTD P - YTD PET) is above a threshold value
    nlf(et ~ step(cpWA - wa) * (Mpet * pet / esy)),
    # nlf(et ~ step(pet - cpPET) * (Mpet * pet / esy)),
    
    # Streamflow and subsurface flow losses above a water level threshold
    nlf(q ~ step(wlL1 - cp) * (Bq + Mq * (wlL1 - cp)^2)),
    
    # Local subsurface inputs when annual water availability is above threshold
    nlf(g ~ step(cp - cpWA) * Bg),
    
    # Ecosystem Specific Yield
    nlf(esy ~ Besy + Mesy ^ (step(cp - wlL1) * (wlL1 - cp))),
    
    # All effects are population effects
    cpWA + cp + Mp + Mpet + Bg + Bq + Mq + Besy + Mesy ~ 1,
    nl = TRUE)

wb_priors <- 
  set_prior(paste0("normal(", get_mode(brm_test$wlObs), ", 1)"), nlpar = "cp") +
  prior(normal(-10, 1), nlpar = "cpWA", ub = 0) +
  # prior(normal(1.5, 1), nlpar = "Mp", lb = 1) +
  prior(constant(1), nlpar = "Mp", lb = 1) +
  prior(constant(-1), nlpar = "Mpet", ub = 0) +
  # prior(normal(1, 0.1), nlpar = "cpPET", lb = 0) +
  prior(normal(1, 10), nlpar = "Bq", ub = 0) +
  prior(normal(0.02, 1), nlpar = "Mq", ub = 0) +
  prior(normal(1, 5), nlpar = "Bg", lb = 0) +
  prior(normal(0.2, 1), nlpar = "Besy", lb = 0, ub = 0.5) +
  prior(normal(1.1, 1), nlpar = "Mesy", lb = 1, ub = 2) +
  prior(gamma(2, .1), class = nu)

wb_mod <- 
  brm(wb_form,
      data = brm_test,
      prior = wb_priors,
      family = student,
      iter = 500,
      warmup = 100,
      seed = 1234,
      chains = 1)

for(i in 1:nrow(val_dat)){
  
  if(i == 1){
    wb_coefs <- 
      fixef(wb_mod, robust = TRUE)
    
    wl_hat2 <- 
      numeric(nrow(val_dat))
    
    wl_hat2[1] <- 
      val_dat$wlObs[1]
    
    precip <- 
      val_dat$precip
    
    pet <- 
      val_dat$pet
    
    wa <- 
      val_dat$wa
    
    p <- et <- q <- g <- esy <- numeric(nrow(val_dat))
    
    cpWA <- wb_coefs["cpWA_Intercept", 1]
    cp <- wb_coefs["cp_Intercept", 1]
    Mpet <- wb_coefs["Mpet_Intercept", 1]
    Mp <- wb_coefs["Mp_Intercept", 1]
    Bq <- wb_coefs["Bq_Intercept", 1]
    Mq <- wb_coefs["Mq_Intercept", 1]
    Bg <- wb_coefs["Bg_Intercept", 1]
    Besy <- wb_coefs["Besy_Intercept", 1]
    Mesy <- wb_coefs["Mesy_Intercept", 1]
    # cpPET <- wb_coefs["cpPET_Intercept", 1]
    
    esy[1] <- 
      Besy + Mesy ^ ((cp > wl_hat2[1]) * (wl_hat2[1] - cp))
    
    next
  }
  
  p[i] <- 
    (cp > wl_hat2[i-1]) * (precip[i-1] / esy[i-1]) + (cp <= wl_hat2[i-1]) * (Mp * precip[i-1])
  
  et[i] <- 
    (cpWA > wa[i-1]) * (Mpet * pet[i-1] / esy[i-1])
  # ((pet[i-1] / esy[i-1]) < cpPET) * (Mpet * pet[i-1] / esy[i-1])  
  
  q[i] <- 
    (cp <= wl_hat2[i-1]) * (Bq + Mq * wl_hat2[i-1]^2)
  
  g[i] <-
    (cpWA <= wa[i-1]) * (Bg)
  
  wl_hat2[i] <- 
    wl_hat2[i-1] + p[i] + et[i] + q[i] + g[i]
  
  esy[i] <- 
    Besy + Mesy ^ ((cp > wl_hat2[i]) * (wl_hat2[i] - cp))
  
}

plot(wlObs ~ sample_date, data = val_dat, col = 'gray40', type = 'l'); lines(wl_hat2 ~ val_dat$sample_date)

# Simple Model ------------------------------------------------------------

get_mode <- 
  function(x){
    d <- density(na.omit(x))
    d$x[which.max(d$y)]
  }

kenton <- external_met[station_name == "kenton"]
kenton[, offset_pet_cm := pmax(0, pet_cm - mean_na(.SD[month(sample_date) %in% c(1:2, 11:12),
                                                       pet_cm]))]
kenton[, offset_wa_cm := cumsum(total_input_cm - offset_pet_cm)]
kenton[, ytd_offset_wa_cm := water_availability_cm - first(water_availability_cm), 
       by = .(sample_year)]
kenton[, l1_ytd_offset_wa_cm := shift(ytd_offset_wa_cm , 1, fill = 0), by = .(sample_year)]

water_budget[, wlL1 := shift(wl_initial_cm, 1), 
             by = .(site, sample_year)]

simple_tr <- 
  water_budget[site == "135" & sample_year %in% c(2012, 2013)
               ][kenton,
                 .(sample_date, 
                   wlObs = wl_initial_cm,
                   dObs = Ds_cm,
                   wlL1 = wlL1,
                   inputs = l1_total_input_cm,
                   pet = l1_pet_cm,
                   melt = l1_melt_cm,
                   rain = l1_rain_cm,
                   wa = i.l1_ytd_offset_wa_cm),
                 on = "sample_date"][!is.na(wlObs) & !is.na(wlL1)]

simple_val <- 
  kenton[sample_year %in% 2012:2014 & sample_date >= as.Date("2012-04-01"),
         .(sample_date,
           wa = ytd_offset_wa_cm,
           pet = pet_cm,
           melt = melt_cm,
           rain = rain_cm,
           inputs = total_input_cm)
         ][water_budget[site == "135" & sample_year %in% 2012:2014],
           `:=`(wlObs = i.wl_initial_cm,
                dObs = i.Ds_cm),
           on = "sample_date"
         ]

simple_form <- 
  bf(
    # Generic Water Budget
    wlObs ~ wlL1 + p - et - q,
    
    nlf(p ~ step(wa - cpWA) * (Mm * (melt + rain)) + step(cpWA - wa) * (rain / esy)),
    
    nlf(et ~ step(cpWA - wa) * (pet / esy)),
    
    # Ecosystem Specific Yield
    # nlf(esy ~ step(1 - (Besy + Mesy * wlL1 + M2esy * wlL1^2)) * step((Besy + Mesy * wlL1 + M2esy * wlL1^2) - 0.05) * (Besy + Mesy * wlL1 + M2esy * wlL1^2) +
    #       step(0.05 - (Besy + Mesy * wlL1 + M2esy * wlL1^2)) * 0.05 + 
    #       step((Besy + Mesy * wlL1 + M2esy * wlL1^2) - 1)),
    nlf(esy ~ Besy + Mesy ^ ((cp > wlL1) * (wlL1 - cp))),
    
    # Surface Flow
    nlf(q ~ step(wlL1 - cp) * (Mq * (wlL1 - cp)^2)),
    
    # All effects are population effects
    cpWA + cp + Besy + Mesy + Mq + Mm ~ 1,
    nl = TRUE)

simple_priors <- 
  set_prior(paste0("normal(", get_mode(simple_tr$wlObs), ", 1)"), nlpar = "cp") +
  prior(normal(10, 1), nlpar = "cpWA") + 
  prior(student_t(3, 0.4, 0.1), nlpar = "Besy", lb = 0) +
  prior(student_t(3, 1.3, 0.1), nlpar = "Mesy", lb = 1) +
  prior(gamma(2, .1), nlpar = "Mq", lb = 0) +
  prior(gamma(2, .1), nlpar = "Mm", lb = 0) +
  prior(gamma(2, .1), class = nu)

simple_mod <- 
  brm(simple_form,
      data = simple_tr,
      prior = simple_priors,
      family = student,
      iter = 500,
      warmup = 100,
      seed = 1234,
      chains = 1)

fit <- 
  fitted(simple_mod)

plot(wlObs ~ sample_date, data = simple_tr, col = 'gray40', type = 'l'); lines(fit[, "Estimate"] ~ simple_tr$sample_date)

?for(i in 1:(nrow(simple_val) - 1)){
  
  if(i == 1){
    simp_coefs <- 
      fixef(simple_mod, robust = TRUE)
    
    swl_hat <- 
      numeric(nrow(simple_val))
    
    swl_hat[1] <- 
      simple_val$wlObs[1]
    
    rain <- 
      simple_val$rain
    
    melt <- 
      simple_val$melt
    
    pet <- 
      simple_val$pet
    
    wa <- 
      simple_val$wa
    
    ds <- p <- et <- q <- esy <- numeric(nrow(simple_val))
    
    cpWA <- simp_coefs["cpWA_Intercept", 1]
    cp <- simp_coefs["cp_Intercept", 1]
    Besy <- simp_coefs["Besy_Intercept", 1]
    Mesy <- simp_coefs["Mesy_Intercept", 1]
    Mq <- simp_coefs["Mq_Intercept", 1]
    Mm <- simp_coefs["Mm_Intercept", 1]
    Bq <- simp_coefs["Bq_Intercept", 1]
    
    esy[1] <- 
      Besy + Mesy ^ ((cp > swl_hat[1]) * (swl_hat[1] - cp))
    
    rm(simp_coefs)
    
  }

  p[i] <- 
    (wa[i] > cpWA) * (Mm * (melt[i] + rain[i])) + 
    (wa[i] <= cpWA) * (rain[i] / esy[i])
  
  et[i] <- 
    (wa[i] <= cpWA) * (pet[i] / esy[i])
  
  q[i] <- 
    (swl_hat[i] >= cp) * (Bq + Mq * (swl_hat[i] - cp)^2)

  ds[i] <- 
    p[i] - et[i] - q[i]
      
  swl_hat[i + 1] <- 
    swl_hat[i] + ds[i]

  esy[i + 1] <- 
    Besy + Mesy ^ ((cp > swl_hat[i]) * (swl_hat[i] - cp))
  
}

plot(wlObs ~ sample_date, data = simple_val, col = 'gray40', type = 'l'); lines(swl_hat ~ simple_val$sample_date)


safe_brm <- 
  quietly(brm)

# 2012 and 2013 and 2012 & 2018 both worked well for 152 test models

update_cp <- 
  function(priors, new.cp){
    priors[priors$nlpar == "cp", ] <- 
      set_prior(paste0("constant(", new.cp, ")"), nlpar = "cp")
    
    priors
  }

wl_mods2012 <- 
  daily_water_levels[sample_year == 2012,
                     .(wlObs = wl_initial_cm,
                       wlL1 = wl_l1_cm,
                       precip = best_precip_cm,
                       pet = pet_cm,
                       wa = normalized_wa_cm,
                       site)][,
                              .(mod = list(safe_brm(wb_form,
                                                    data = .SD,
                                                    prior = update_cp(wb_priors, quantile(daily_water_levels[site == .BY[[1]], wl_initial_cm], 0.8, na.rm = TRUE)),
                                                    family = student,
                                                    iter = 500,
                                                    warmup = 100,
                                                    seed = 1234,
                                                    chains = 1))),
                              keyby = .(site)]

# saveRDS(wl_mods2012, "tmp/full_wetland_models.rds")


wl_mods2012[, fit := map(mod, "result")]

wl_mods2012_extras <- 
  wl_mods2012[,
              .(site,
                output = map(mod, "output"),
                warnings = map(mod, "warnings"),
                messages = map(mod, "messages"))]


wl_mods2012[, pop_eff := map(fit, ~fixef(.x, robust = TRUE)[, 1])]

daily_water_levels[sample_date == as.Date("2012-04-01"), 
                   .(site, sample_date, wl_initial_cm)]

wl_mods2012[daily_water_levels[sample_date == as.Date("2012-04-01")],
            initial.wl := i.wl_initial_cm,
            on = c("site")]

wl_mods2012[, met := list(kenton[sample_date >= as.Date("2012-04-01")])]

predict_water_levels <- 
  function(initial.wl,
           met,
           pop_eff){
    
    for(i in 1:nrow(met)){
      
      if(i == 1){
        wl_hat <- 
          numeric(nrow(met))
        
        wl_hat[1] <- 
          initial.wl
        
        precip <- 
          met$precip
        
        pet <- 
          met$pet
        
        wa <- 
          met$wa
        
        p <- et <- q <- g <- esy <- numeric(nrow(met))
        
        cpWA <- pop_eff[["cpWA_Intercept"]]
        cp <- pop_eff[["cp_Intercept"]]
        Mpet <- pop_eff[["Mpet_Intercept"]]
        Bq <- pop_eff[["Bq_Intercept"]]
        Mq <- pop_eff[["Mq_Intercept"]]
        Bg <- pop_eff[["Bg_Intercept"]]
        Besy <- pop_eff[["Besy_Intercept"]]
        Mesy <- pop_eff[["Mesy_Intercept"]]
        
        esy[1] <- 
          Besy + Mesy ^ ((cp > wl_hat[1]) * (wl_hat[1] - cp))
        
        next
      }
      
      p[i] <- 
        (cp > wl_hat[i-1]) * (precip[i-1] / esy[i-1]) + (cp <= wl_hat[i-1]) * (Mp * precip[i-1])
      
      et[i] <- 
        (cpWA > wa[i-1]) * (Mpet * pet[i-1] / esy[i-1])  
      
      q[i] <- 
        (cp <= wl_hat[i-1]) * (Bq + Mq * wl_hat[i-1]^2)
      
      g[i] <-
        (cpWA <= wa[i-1]) * (Bg)
      
      wl_hat[i] <- 
        wl_hat[i-1] + p[i] + et[i] + q[i] + g[i]
      
      esy[i] <- 
        Besy + Mesy ^ ((cp > wl_hat[i]) * (wl_hat[i] - cp))
      
    }
    
    data.table(sample_date = met$sample_date,
               wl_hat,
               p_hat = p,
               et_hat = et,
               q_hat = q,
               g_hat = g,
               esy_hat = esy)
  }


wl_mods2012[,
            out := pmap(.SD[, .(initial.wl, met, pop_eff)],
                        predict_water_levels)]

plot(wl_hat ~ sample_date, 
     data = wl_mods2012["135", out[[1]]],
     type = 'l')
lines(wl_initial_cm ~ sample_date,
       data = daily_water_levels[site == "135"],
       col = 'gray40')

daily_water_levels[wl_mods2012[, out[[1]], by = .(site)], 
                   wl_hat := i.wl_hat,
                   on = c("site", "sample_date")]

split(daily_water_levels[!is.na(wl_initial_cm) & site != "006"],
      by = "site") %>% 
  map(~as.data.table(t(gof(sim = .x$wl_hat, obs = .x$wl_initial_cm)))) %>% 
  rbindlist()

daily_water_levels[site == "152", 
                   as.data.table(t(gof(sim = .SD$wl_hat, obs = .SD$wl_initial_cm)))]




# Physical Model ----------------------------------------------------------

bank_water <- 
  function(x){
    
    b <- 0
    
    for(i in 1:length(x)){
      
      x[i] <- 
        x[i] + b/10
      
      if(x[i] > 0){
        b <- x[i]
        x[i] <- 0
      } else {
        b <- 0
      }
    }
    
    x
  }

kenton[sample_year != 2020, 
       banked_wa := bank_water(ytd_offset_wa_cm),
       by = .(sample_year)]

plot(ytd_offset_wa_cm ~ sample_date, data = kenton[year(sample_date) %in% 2012:2014], 
     type = 'l')


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

kenton[, new_ytd_water_availability_cm := (water_availability_cm - first(water_availability_cm)),
       by = .(station_name, sample_year)]

water_budget[kenton,
             new_ytd_water_availability_cm := i.new_ytd_water_availability_cm,
             on = "sample_date"]

# new_ytd_water_availability_cm shows more variation within a small range of 
# new_ytd_water_availability_cm

drying_curves <- 
  water_budget[site == "152", 
               .SD[
                 # Index of maximum water availability (melt period)
                 which.max(water_availability_cm):
                   # Index of max WA + Index of min WA that occurs after max WA
                   (which.max(water_availability_cm) + which.min(.SD[which.max(water_availability_cm):.N, water_availability_cm]))], 
               by = .(sample_year)]

ggplot(drying_curves,
       aes(x = sample_date,
           y = water_availability_cm)) +
  geom_line() +
  facet_wrap(~sample_year,
             scales = "free")

# The offset to align the curves lay in the selection of the start period for
# the drawdown. When I offset each curve to the water availability at the start
# of the drawdown, the curves sync much better. Still some work to identify the
# best start ands top points for the drawdown curve to get the best curve fit.
# It also seems like just using water avilability will work too
ggplot(drying_curves[drying_curves[, .SD[1, .(sample_date, water_availability_cm)], by = .(sample_year)], on = "sample_year", .(wl_initial_cm, sample_year,  water_availability_cm = water_availability_cm - i.water_availability_cm)], 
       aes(x = water_availability_cm, 
           y = wl_initial_cm)) + 
  geom_point(aes(color = factor(sample_year)))

drying_curves <- 
  drying_curves[drying_curves[, .SD[1, .(sample_date, water_availability_cm)], 
                              by = .(sample_year)], 
                on = "sample_year", 
                .(wl_initial_cm, sample_year, sample_date,
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





# The new drying curves seem to work pretty well too. Maybe they would be easier
# to get all the years aligned?
# ggplot(drying_curves[, 
#                      .(ytd_water_availability_cm = cumsum(pmin(total_input_cm - pet_cm, 0)),
#                        wl_initial_cm = cumsum(pmin(Ds_cm, 0))),
#                      by = .(sample_year)], 
#        aes(x = ytd_water_availability_cm, 
#            y = wl_initial_cm)) + 
#   geom_point(aes(color = factor(sample_year)))
# new_drying <-
#   drying_curves[,
#                 .(ytd_water_availability_cm = cumsum(pmin(total_input_cm - pet_cm, 0)),
#                   wl_initial_cm = cumsum(pmin(Ds_cm, 0))),
#                 by = .(sample_year)]
# 
# drying_mod <-
#   robustbase::nlrob(wl_initial_cm ~ b + m1 * ytd_water_availability_cm - m2**ytd_water_availability_cm,
#                     start = list(b = 20, m1 = 1, m2 = 0.86), na.action = na.exclude,
#                     data = new_drying,
#                     maxit = 50)
# 
# hyst_function <-
#   deriv(bquote(.(b) + .(m1)*water.availability - .(m2)**water.availability,
#                where = list(b = coef(drying_mod)[[1]],
#                             m1 = coef(drying_mod)[[2]],
#                             m2 = coef(drying_mod)[[3]])),
#         namevec = "water.availability",
#         func = TRUE)


# WA <-
#   water_budget[site == "157" & sample_year == 2012][min(which(!is.na(wl_initial_cm))):.N,
#     ytd_water_availability_cm]

WA <- kenton[between(sample_year, 2007, 2019), new_ytd_water_availability_cm]
D <- kenton[between(sample_year, 2007, 2019), sample_date]
PET <- kenton[between(sample_year, 2007, 2019), pet_cm]
P <- kenton[between(sample_year, 2007, 2019), total_input_cm]
R <- kenton[between(sample_year, 2007, 2019), rain_cm]
M <- kenton[between(sample_year, 2007, 2019), melt_cm]

WL <- 
  water_budget[site == "157" & sample_year == 2013][min(which(!is.na(wl_initial_cm))):.N, 
    wl_initial_cm]

# PET <- 
#   water_budget[site == "135"][min(which(!is.na(wl_initial_cm))):.N, 
#     pet_cm]

# P <- 
#   water_budget[site == "135"][min(which(!is.na(wl_initial_cm))):.N, 
#     total_input_cm]

# Adjust by Previous year's water deficit
kenton[,
       .(wd = .SD[!is.na(ytd_water_availability_cm), 
                  last(ytd_water_availability_cm) - first(ytd_water_availability_cm)]), 
       by = .(sample_year)][, .(sample_year, wd = shift(wd, 1))]

ggplot(drying_curves[kenton[,
                            .(wd = .SD[!is.na(ytd_water_availability_cm), 
                                       last(ytd_water_availability_cm) - first(ytd_water_availability_cm)]), 
                            by = .(sample_year)][, .(sample_year, wd = shift(pmin(5, pmax(wd, 0)), 1))], 
                     on = "sample_year", nomatch = NULL,
                     .(wl_initial_cm, sample_year, ytd_water_availability_cm = ytd_water_availability_cm + i.wd)], 
       aes(x = ytd_water_availability_cm, 
           y = wl_initial_cm)) + 
  geom_point(aes(color = factor(sample_year))) + 
  geom_function(fun = ~14.875 - 0.8455**.x)


for(t in 1:length(WA)){
  
  # Set up on first iteration
  if(t == 1){
   
    X <- 
      WA
    
    # Y is predicted water level
    # out is water level output from model
    # gradient is the derivative of the model slope 
    # dWA is the change in X
    Y <- out <- gradient <- dX <- 
      numeric(length(WL))
    
    # Start the predictions at water yield == 0
    Y[1] <- coef(drying_mod)[1]
    
    # Move onto t = 2
    next
  }
  
  # Calculate dX at step t
  dX[t] <- 
    X[t] - X[t-1]
  
  # Hystersis Loop
  if(dX[t] < 0){
    # If x is a drawdown, then 
    Y[t] <- 
      Y[t-1] + (hyst_function(X[t]) - hyst_function(X[t-1]))
  } else {
    # Y[t] <- Y[t-1] + (X[t] - X[t-1]) * gradient[t]
    # After the new water level is computed, I should recalculate a new water 
    # availability (based on the water level) and then use the maximum of the 
    # observed water availability or the offset water availability
    # Save slope of gradient at t
    gradient[t] <- 
      attr(hyst_function(X[t]), "gradient")[1, ]
    
    Y[t] <- 
      pmax(Y[t-1], hyst_function(X[t-1] + dX[t] * gradient[t]))
      
  }
  
}

ggplot(water_budget[site == "157"], 
       aes(y = wl_initial_cm, x = sample_date)) +
  geom_line(color = 'gray40') +
  # geom_line(data = data.table(sample_date = D, 
  #                             wl_initial_cm = out, 
  #                             sample_year = year(D))[sample_year %between% c(2012, 2019)]) +
  geom_line(data = data.table(sample_date = D, 
                              wl_initial_cm = Y,
                              sample_year = year(D))[sample_year %between% c(2012, 2019)], 
            color = "blue") + 
  facet_wrap(~sample_year, scales = "free")


plot(wl_initial_cm ~ sample_date, data = water_budget[site == "157" & sample_year == 2014], type = 'l', col = 'gray40'); lines(Y ~ D)

plot(WL, type = 'l', col = 'gray40'); lines(Y)
plot(WA, type = 'l')
plot(Y ~ D, type = 'l')



plot(-attr(deriv(bquote(.(c) + (.(d) - .(c))*(1-exp(-exp(.(b)*(water.availability - .(e))))),
                        where = as.list(coef(robustbase::nlrob(wl_initial_cm ~ c + (d-c)*(1-exp(-exp(b*(ytd_water_availability_cm - e)))),
                                                               drying_curves[sample_year == 2012],
                                                               na.action = na.exclude,
                                                               start = list(c = 20, d = -70, b = -0.25, e = -25))))),
                 namevec = "water.availability",
                 func = TRUE)(drying_curves[sample_year == 2012, ytd_water_availability_cm]), "gradient")[,1] * drying_curves[sample_year == 2012, pet_cm], type = 'h')

plot(attr(deriv(bquote(.(c) + (.(d) - .(c))*(1-exp(-exp(.(b)*(water.availability - .(e))))),
                       where = as.list(coef(robustbase::nlrob(wl_initial_cm ~ c + (d-c)*(1-exp(-exp(b*(ytd_water_availability_cm - e)))),
                                                              drying_curves[sample_year == 2012],
                                                              na.action = na.exclude,
                                                              start = list(c = 20, d = -70, b = -0.25, e = -25))))),
                namevec = "water.availability",
                func = TRUE)(drying_curves[sample_year == 2012, ytd_water_availability_cm]), "gradient")[,1] ~ drying_curves[sample_year == 2012, wl_initial_cm])

plot(-attr(hyst_function(drying_curves[sample_year == 2012, ytd_water_availability_cm]), "gradient")[,1] * drying_curves[sample_year == 2012, pet_cm], type = 'h')

plot(attr(deriv(bquote(.(c) + (.(d) - .(c))*(1-exp(-exp(.(b)*(water.availability - .(e))))),
                       where = as.list(coef(robustbase::nlrob(wl_initial_cm ~ c + (d-c)*(1-exp(-exp(b*(ytd_water_availability_cm - e)))),
                                                              drying_curves[sample_year == 2012],
                                                              na.action = na.exclude,
                                                              start = list(c = 20, d = -70, b = -0.25, e = -25))))),
                namevec = "water.availability",
                func = TRUE)(drying_curves[sample_year == 2012, ytd_water_availability_cm]), "gradient")[,1] * (cumsum(drying_curves[sample_year == 2012, total_input_cm] - drying_curves[sample_year == 2012, pet_cm])) + 20, type = 'l'); lines(drying_curves[sample_year == 2012, wl_initial_cm], col = 'red')



# This may be exatly what I want
plot((cumsum(attr(deriv(bquote(.(c) + (.(d) - .(c))*(1-exp(-exp(.(b)*(water.availability - .(e))))),
                               where = as.list(coef(robustbase::nlrob(wl_initial_cm ~ c + (d-c)*(1-exp(-exp(b*(ytd_water_availability_cm - e)))),
                                                                      drying_curves[sample_year == 2012],
                                                                      na.action = na.exclude,
                                                                      start = list(c = 20, d = -70, b = -0.25, e = -25))))),
                        namevec = "water.availability",
                        func = TRUE)(drying_curves[sample_year == 2012, ytd_water_availability_cm]), "gradient")[,1] * drying_curves[sample_year == 2012, total_input_cm]) -
        cumsum(attr(deriv(bquote(.(c) + (.(d) - .(c))*(1-exp(-exp(.(b)*(water.availability - .(e))))),
                                 where = as.list(coef(robustbase::nlrob(wl_initial_cm ~ c + (d-c)*(1-exp(-exp(b*(ytd_water_availability_cm - e)))),
                                                                        drying_curves[sample_year == 2012],
                                                                        na.action = na.exclude,
                                                                        start = list(c = 20, d = -70, b = -0.25, e = -25))))),
                          namevec = "water.availability",
                          func = TRUE)(drying_curves[sample_year == 2012, ytd_water_availability_cm]), "gradient")[,1] * drying_curves[sample_year == 2012, pet_cm])) + 20, type = 'l'); lines(drying_curves[sample_year == 2012, wl_initial_cm], col = 'red')


# Drawdown should be controlled by this function:
plot(x = drying_curves[sample_year != 2015, wl_initial_cm],
     y = attr(hyst_function(drying_curves[sample_year!= 2015, water_availability_cm]), "gradient")[,1],
     main = "Drying Curve Gradient",
     xlab = "Water Level (cm)",
     ylab = expression(delta/delta~X))

library(mcp)
gradient_mod <- 
  mcp(list(gradient ~ wl_initial_cm, ~ 0 + wl_initial_cm),
      data = drying_curves[sample_year != 2015 & !is.na(wl_initial_cm), 
                           .(wl_initial_cm,
                             gradient = attr(hyst_function(water_availability_cm), "gradient")[,1])])

grad_coefs <- 
  mcp::fixef(gradient_mod)

# gradient_function <- 
#   as.function(list(water.level = NULL,
#                    bquote((.(cp_1) >= water.level) * (.(int_1) + water.level * .(wl_initial_cm_1)) +
#                             (.(cp_1) < water.level) * (.(int_1) + .(cp_1) * .(wl_initial_cm_1) + (water.level - .(cp_1)) * .(wl_initial_cm_2)),
#                           where = set_names(as.list(grad_coefs$mean), grad_coefs$name))))

gradient_function <- 
  as.function(list(water.level = NULL,
                   bquote(.(int_1) + .(cp_1) * .(wl_initial_cm_1) + (water.level - .(cp_1)) * .(wl_initial_cm_2),
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
      Y[t] + pmin(0, (hyst_function(X[t]) - ifelse(t > 1, hyst_function(X[t-1]), 0)))
  } else {
    Y[t+1] <-
      Y[t] + gradient[t] * P[t]
  }
  
  # # Hystersis Loop
  # # Should the conditional take into account the expected precip vs pet response?
  # # hyst_function(X[t] + dX[t]) < hyst_function(X[t-1] + dX[t] * gradient[t])
  # if(dX[t] < 0){
  #   # If x is a drawdown, then 
  #   Y[t] <- 
  #     y_hat[t]
  #   
  #   X[t + 1] <- 
  #     X[t] + dX[t]
  #   
  # } else {
  #   # After the new water level is computed, I should recalculate a new water 
  #   # availability (based on the water level) and then use the maximum of the 
  #   # observed water availability or the offset water availability
  #   Y[t] <- 
  #     hyst_function(X[t-1] + dX[t] * gradient[t])
  #  
  #    X[t + 1] <- 
  #     # Set X[t+1] based on adding gradient of precip
  #     X[t] + (gradient[t]*(P[t]) - PET[t])
  #     # Set X[t+1] same as for drawdown days
  #      # X[t] + dX[t]
  #    # Set X[t+1] based on Y[t] 
  #    # ifelse(Y[t] > coef(drying_mod)[["d"]],
  #    #        e + log(-log(-((Y[t]-c)/(d-c) - 1)))/b,
  #    #        X[t] + P[t] - PET[t])
  # }
  # 
  # # if(coef(drying_mod)[["b"]] < Y[t]){
  # #   Y[t] <-
  # #     Y[t] - (Y[t] - coef(drying_mod)[["b"]])/2
  # # }
  # 
  
}


ggplot(kenton[sample_year %between% c(2012, 2019)], 
       aes(y = ytd_water_availability_cm, x = sample_date)) +
  geom_line(color = 'gray40') +
  geom_line(data = data.table(sample_date = D, 
                              ytd_water_availability_cm = X,
                              sample_year = year(D))[sample_year %between% c(2012, 2019)], 
            color = "blue",
            linetype = "dotted") + 
  facet_wrap(~sample_year, scales = "free_x")


ggplot(water_budget[site == "152"], 
       aes(y = wl_initial_cm, x = sample_date)) +
  geom_line(color = 'gray40') +
  geom_line(data = data.table(sample_date = D, 
                              wl_initial_cm = Y,
                              sample_year = year(D))[sample_year %between% c(2012, 2019)], 
            color = "blue",
            linetype = "dotted") + 
  facet_wrap(~sample_year, scales = "free_x")


ggplot(water_budget[site == "152", 
                    .(sample_date,
                      wl_initial_cm = as.numeric(scale(wl_initial_cm))),
                    by = .(sample_year)], 
       aes(y = wl_initial_cm, x = sample_date)) +
  geom_line(color = 'gray40') +
  geom_line(data = data.table(sample_date = D, 
                              wl_initial_cm = Y,
                              sample_year = year(D))[sample_year %between% c(2012, 2019)][,
                                                                                          .(sample_date,
                                                                                            wl_initial_cm = as.numeric(scale(wl_initial_cm))),
                                                                                          by = .(sample_year)], 
            color = "blue",
            linetype = "dotted") + 
  facet_wrap(~sample_year, scales = "free_x")



plot(wl_initial_cm ~ sample_date,
     data = water_budget[site == "157" & sample_year == 2014],
     col = 'gray40',
     type = 'l')
lines(Y ~ D)
lines(WA ~ D, col = palette()[4])

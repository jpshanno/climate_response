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

daily_water_levels[, wl_l1_cm := shift(wl_initial_cm, 1),
                   by = .(site)]

daily_water_levels[, threshold := coef(lm(wl_l1_cm ~ Ds_cumulative,
                                          data = .SD))[[1]],
                   by = .(site)]


# Big Bayesian Model ------------------------------------------------------

library(brms)
options(mc.cores = 4)

brm_test <- 
  daily_water_levels[site == "152" & sample_year %in% c(2012, 2013) & !is.na(wl_initial_cm), 
                     .(wlObs = wl_initial_cm,
                       wlL1 = wl_l1_cm,
                       pet = pet_cm,
                       wa = normalized_wa_cm,
                       precip = best_precip_cm)]

val_dat <- 
  daily_water_levels[site == "152" & sample_year %in% c(2012:2014), 
                     .(sample_date, 
                       wlObs = wl_initial_cm,
                       wlL1 = wl_l1_cm,
                       pet = pet_cm,
                       wa = normalized_wa_cm,
                       precip = best_precip_cm)][min(which(!is.na(wlObs))):.N]

# To Do:
# - Get water levels to go above cp (probably need to add piecewise component to p)
# - Add quickflow in/out from preceding day's precip
# - Try to vary esy for P and PET again. Maybe having P when wl > cp not depend
# on esy will help set a more accurate esy(s). Or try to work, previously they
# ended up negatively afffecting the esy fit
# - Consider changing distributions for priors (i.e. lognormal, or skew_normal for
# some things)
# - Try making G dependent on normalized_wa_cm (need to experiment for decent priors)
# - Adjust cp from constant prior based on observed percentiles to a strict 
# distributional prior

wb_form <- 
  bf(
    # Generic Water Budget
    wlObs ~ wlL1 + p + et + q + g,
    
    # Precip response scaled by ecosystem specific yield
    nlf(p ~ step(wlL1 - cp) * (precip / esy) + step(cp - wlL1) * (Mp * precip)),
    
    # PET response scaled by ecosystem specific yield and 0 when annual water
    # availability (YTD P - YTD PET) is above a threshold value
    nlf(et ~ step(cpWA - wa) * (Mpet * pet / esy)),
    
    # Streamflow and subsurface flow losses above a water level threshold
    nlf(q ~ step(wlL1 - cp) * (Bq + Mq * wlL1^2)),
    
    # Local subsurface inputs when annual water availability is above threshold
    nlf(g ~ step(cp - cpWA) * Bg),
    
    # Ecosystem Specific Yield
    nlf(esy ~ Besy + Mesy ^ (step(cp - wlL1) * (wlL1 - cp))),
    
    # All effects are population effects
    cpWA + cp + Mp + Mpet + Bg + Bq + Mq + Besy + Mesy ~ 1,
    nl = TRUE)

wb_priors <- 
  # Not sure how to set priors on the intermediate values, got an error that
  # they are not parameters in the model. It would be nice to be able to control
  # them directly rather than through multiple parameters
  # prior(gamma(0.16, 1.54), nlpar = "p", lb = 0) + # determined via lmomco for brm_test$precip_cm
  # prior(normal(5, 10), nlpar = "et", ub = 0) +
  # prior(gamma(0.5, 1), nlpar = "q", ub = 0) +
  # prior(normal(5, 10), nlpar = "g", lb = 0) +
  # prior(normal(1, 1), nlpar = "Mp", lb = 0) +
  # prior(normal(0.5, 1), nlpar = "esy", lb = 0) +
  set_prior(paste0("constant(", quantile(brm_test$wlObs, 0.8), ")"), nlpar = "cp") +
  prior(normal(-15, 1), nlpar = "cpWA", ub = 0) +
  prior(normal(1.5, 1), nlpar = "Mp", lb = 1) +
  prior(constant(-1), nlpar = "Mpet", ub = 0) +
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
      fixef(wb_mod)
    
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
    
    esy[1] <- 
      Besy + Mesy ^ ((cp > wl_hat2[1]) * (wl_hat2[1] - cp))
    
    next
  }
  
  p[i] <- 
    (cp > wl_hat2[i-1]) * (precip[i-1] / esy[i-1]) + (cp <= wl_hat2[i-1]) * (Mp * precip[i-1])
  
  et[i] <- 
    (cpWA > wa[i-1]) * (Mpet * pet[i-1] / esy[i-1])  
  
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
                                                    prior = update_cp(wb_priors, quantile(.SD$wlObs, 0.8, na.rm = TRUE)),
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

wl_mods2012[, met := list(kenton[sample_date >= as.Date("2012-04-01"),
                                 .(sample_date, 
                                   pet = pet_cm,
                                   wa = normalized_wa_cm,
                                   precip = precip_cm)])]

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
     data = wl_mods2012["152", out[[1]]],
     type = 'l')
lines(wl_initial_cm ~ sample_date,
       data = daily_water_levels[site == "152"],
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




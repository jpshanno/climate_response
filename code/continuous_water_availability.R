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
  # external_met[format(sample_date, "%m%d") %in% c("0401", "1231"),
  #            .(water_year, wa_trend = ((shift(water_availability_cm, -1) - shift(water_availability_cm, 1)) / (365*3))),
  #            by = .(station_name,
  #                   sample_year)][, .(station_name, water_year, wa_trend)]

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


ggplot(external_met[water_year %in% 2012:2020], 
       aes(x = sample_date, y = normalized_wa_cm)) + 
  geom_line() + 
  geom_hline(aes(yintercept = 0)) +
  facet_grid(station_name ~ sample_year, scales = "free") + 
  scale_x_date(date_breaks = "1 month", date_labels = "%m", expand = expansion()) + 
  theme_bw()


ggplot(external_met, 
       aes(x = sample_date, y = detrend_wa_cm)) + 
  geom_line() + 
  geom_point(data = external_met[format(sample_date, "%m%d") == "1101"], 
             color = "red") + 
  facet_wrap(~station_name, ncol = 1) + 
  scale_x_date(breaks = as.Date(paste0(2000:2020, "-11-01")))


daily_water_levels[, Ds_continuous := diff_na(wl_initial_cm),
                   by = .(site)]

daily_water_levels[, Ds_cumulative := ytd_sum(Ds_continuous),
                   by = .(site)]

ggplot(daily_water_levels,
       aes(x = wl_initial_cm,
           y = normalized_wa_cm/Ds_cumulative)) +
  geom_point(shape = 20,
             alpha = 0.2) +
  coord_cartesian(ylim = c(-1, 3)) +
  facet_wrap(~site, 
             scales = "free")


ggplot(daily_water_levels,
       aes(x = sample_date,
           y = Ds_cumulative)) +
  geom_line() +
  facet_grid(site ~ water_year, 
             scales = "free")

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


ggplot(daily_water_levels,
       aes(x = wl_initial_cm,
           y = Ds_cumulative/normalized_wa_cm)) +
  geom_point() +
  coord_cartesian(ylim = c(-20, 20)) +
  facet_grid(site ~ sample_year, 
             scales = "free")

ggplot(daily_water_levels,
       aes(x = wl_initial_cm,
           y = Ds_cumulative/normalized_wa_cm,
           color = factor(sample_year))) +
  geom_point() +
  facet_wrap(~ site, 
             scales = "free")

ggplot(daily_water_levels,
       aes(x = wl_initial_cm,
           y = Ds_cumulative/normalized_wa_cm)) +
  geom_point(aes(color = factor(sample_year))) +
  geom_smooth(method = lmrob,
              formula = y ~ x,
              method.args = list(setting = "KS2014")) +
  coord_cartesian(ylim = c(-20, 20)) +
  facet_wrap(~ site, 
             scales = "free")

daily_water_levels[, esy_wa_emp := Ds_cumulative / normalized_wa_cm]

daily_water_levels[, wl_l1_cm := shift(wl_initial_cm, 1),
                   by = .(site)]

daily_water_levels[, threshold := coef(lm(wl_l1_cm ~ Ds_cumulative,
                                          data = .SD))[[1]],
                   by = .(site)]


ggplot(daily_water_levels[wl_l1_cm < threshold],
       aes(x = wl_l1_cm,
           y = esy_wa_emp,
           color = factor(sample_year))) +
  geom_point() +
  facet_wrap(~ site, 
             scales = "free")

esy_wa_mods <-
  daily_water_levels[is.finite(esy_wa_emp) & wl_l1_cm < threshold,
                     .(mod_lm = list(lmrob(esy_wa_emp ~ wl_l1_cm,
                                           data = .SD,
                                           setting = "KS2014"))
                       , mod_quad = list(lmrob(esy_wa_emp ~ wl_l1_cm + I(wl_l1_cm^2),
                                               data = .SD,
                                               setting = "KS2014"))
                       , b = min_na(esy_wa_emp)
                       , mod_nls = list(nls(esy_wa_emp ~ m**(wl_l1_cm - o),
                                            start = list(m = 1.1, o = unique(.SD$threshold)),
                                            control = nls.control(warnOnly = TRUE),
                                            data = .SD))
                       # , mod_mcp = list(mcp(list(esy_wa_emp ~ wl_initial_cm,
                       #                           ~ 0 + wl_initial_cm),
                       #                      prior = list(cp_1 = "dnorm(0, 1)",
                       #                                   int_1 = "dnorm(1, 1)"),
                       #                      cores = 1,
                       #                      data = .SD))
                     ),
                     keyby = .(site)]

esy_wa_mods[,
            f_nls := map2(mod_nls, b,
                          ~{
                            coefs <-
                              coef(.x)
                            as.function(list(water.level = NULL,
                                             bquote(.(b) + .(m)**(water.level - .(o)),
                                                    where = list(b = .y,
                                                                 m = coefs[[1]],
                                                                 o = coefs[[2]]))))
                          })]



esy_wa_mods[,
            f_lm := map(mod_lm,
                        ~{
                          coefs <-
                            coef(.x)

                          as.function(list(water.level = NULL,
                                           bquote(.(b) + .(m) * water.level,
                                                  where = list(b = coefs[[1]],
                                                               m = coefs[[2]]))))
                        })]

esy_wa_mods[,
            f_quad := map(mod_quad,
                        ~{
                          coefs <-
                            coef(.x)

                          as.function(list(water.level = NULL,
                                           bquote(.(b) + .(m) * water.level + .(m2) * water.level^2,
                                                  where = list(b = coefs[[1]],
                                                               m = coefs[[2]],
                                                               m2 = coefs[[3]]))))
                        })]

# esy_wa_mods[,
#             f_mcp := map(mod_mcp,
#                          ~{
#                            coefs <-
#                              mcp::fixef(.x)
# 
#                            as.function(list(water.level = NULL,
#                                             bquote(.(b) +
#                                                      (water.level <= .(cp)) * (water.level * .(m1)) +
#                                                      (water.level > .(cp)) * (.(cp) * .(m1) + .(m2) * (water.level - .(cp))),
#                                                    where = list(cp = coefs[1, "mean"],
#                                                                 b = coefs[2, "mean"],
#                                                                 m1 = coefs[4, "mean"],
#                                                                 m2 = coefs[5, "mean"]))))
#                          })]

daily_water_levels[, esy_wa_lm := esy_wa_mods[.BY[[1]], f_lm[[1]]](wl_l1_cm),
                   by = .(site)]
daily_water_levels[, Ds_hat_lm := normalized_wa_cm * esy_wa_lm]

daily_water_levels[, esy_wa_quad := esy_wa_mods[.BY[[1]], f_quad[[1]]](wl_l1_cm),
                   by = .(site)]
daily_water_levels[, Ds_hat_quad := normalized_wa_cm * esy_wa_quad,
                   by = .(site)]

daily_water_levels[, esy_wa_nls := esy_wa_mods[.BY[[1]], f_predict[[1]]](wl_initial_cm),
                   by = .(site)]
daily_water_levels[, Ds_hat_nls := normalized_wa_cm / esy_wa_nls]

# daily_water_levels[, esy_wa_mcp := esy_wa_mods[.BY[[1]], f_mcp[[1]]](wl_l1_cm),
#                    by = .(site)]
# daily_water_levels[, Ds_hat_mcp := normalized_wa_cm * esy_wa_mcp,
#                    by = .(site)]

ggplot(daily_water_levels[wl_initial_cm < (threshold - 1)],
       aes(x = wl_initial_cm,
           y = esy_wa_emp)) +
  geom_point() +
  geom_line(aes(y = esy_wa_nls), 
            color = "red") +
  # geom_line(aes(y = esy_wa_quad),
  #           color = "blue") +
  # geom_line(aes(y = esy_wa_mcp), 
  #           color = "blue") +
  facet_wrap(~ site, 
             scales = "free")

ggplot(daily_water_levels,
       aes(x = sample_date,
           y = Ds_cumulative)) +
  geom_line(color = "gray40") +
  geom_line(aes(y = Ds_hat_lm),
            color = "blue") +
  facet_grid(site ~ sample_year, 
             scales = "free")

ggplot(daily_water_levels[wl_initial_cm < (threshold - 1)],
       aes(x = Ds_cumulative,
           y = Ds_hat_nls)) +
  geom_point() +
  geom_abline(color = "red") +
  facet_grid(site ~ sample_year, 
             scales = "free")

test <- 
  daily_water_levels[site == "152" & sample_year == 2013][min(which(!is.na(wl_initial_cm))):.N]

f_esy <- 
  esy_wa_mods["152", f_lm[[1]]]

WL <- 
  test$wl_initial_cm

threshold <- 
  unique(test$threshold)

WA <- 
  test$normalized_wa_cm

P <- 
  test$best_precip_cm

PET <- 
  test[external_met[station_name == "kenton"], pet_cm, nomatch = NULL, on = "sample_date"]

SY <- 
  test$esy_wa_emp

Dobs <- 
  test$Ds_cumulative

Dcm <- 
  test$Ds_cm

d_hat <- 
  numeric(nrow(test))

d_hat[1] <- 
  0

wl_hat <- 
  numeric(nrow(test))

wl_hat[1] <- 
  WL[1]

sy_hat <- 
  numeric(nrow(test))

sy_hat[1] <- 
  f_esy(WL[1])

for(i in 2:length(d_hat)){

  # if(wl_hat[i-1] < threshold-0.6){
  #   sy_hat[i] <- 
  #     0.1932 + 1.1532**(wl_hat[i-1] - 17.2069)
  #   # daily_water_levels[, sy_emp := (best_precip_cm + melt_cm - pet_cm) / Ds_cm]
  #   # nls(sy_emp ~ b + m**(wl_l1_cm-c),
  #   #     data = daily_water_levels[site == "152" & wl_l1_cm < threshold & between(sy_emp, 0, 1)],
  #   #     start = list(b = 0.1, m = 1.1, c = 20),
  #   #     na.action = na.exclude)
  #   
  # } else {
  #   sy_hat[i] <- 
  #     10
  # }
  # 
  # d_hat[i] <- 
  #   (P[i]-PET[i])/sy_hat[i]
  
  # sy_hat[i] <- 
  #   0.1328917 + 1.1109099**(wl_hat[i-1] - 31.4566861)
  
  # if(P[i] < PET[i]){
    sy_hat[i] <-
      0.306868 + 0.002693 * wl_hat[i-1]
      # ifelse(wl_hat[i-1] < threshold,
      #        0.2047 + 1.1045**(wl_hat[i-1] - 21.3802),
      #        # Test as below
      #        # nlrob({(-pet_cm)/Ds_cm} ~ b + m **(wl_l1_cm - c), data = test[best_precip_cm < 0.2 & wl_l1_cm < threshold],
      #        #       na.action = na.exclude,
      #        #       start = c(b = 0.2, m = 1.1, c = 20))
      #        2)
      # pmin(1,
      #      0.5 + 1.2437181**(wl_hat[i-1] - 12.8420301))

      d_hat[i] <-
        (-PET[i]) / sy_hat[i]

      # if(wl_hat[i-1] > threshold) {
      #
      #   d_hat[i] <-
      #     d_hat[i] + ((threshold - wl_hat[i-1]) / 2)
      #
      # }

  # } else {

    sy_hat[i] <-
      0.216396 + 0.002312 * wl_hat[i - 1]
      # ifelse(wl_hat[i-1] < threshold,
      #        0.1432 + 1.1785**(wl_hat[i-1] - 15.0037),
      #        # test <- daily_water_levels[site == "152"]
      #        # nlrob({(best_precip_cm - 0.2)/Ds_cm} ~ b + m**(wl_l1_cm - c), data = test[wl_l1_cm < threshold & !is.na(best_precip_cm) & !is.na(Ds_cm)],
      #        # start = c(b = 0.05, m = 1.1, c = 50))
      #        1)
      # pmin(1,
      #      0.4 + 1.2437181**(wl_hat[i-1] - 12.8420301))

    d_hat[i] <-
     (P[i]) / sy_hat[i] + d_hat[i]
  # }
  
  wl_hat[i] <- 
    wl_hat[i-1] + d_hat[i]
  
  if(wl_hat[i] > (threshold)){
    wl_hat[i] <-
      wl_hat[i] + (threshold - wl_hat[i]) / 2
  }

}

plot(WL, type = 'l', col = "gray40"); lines(wl_hat)
plot(SY, type = 'l', col = "gray40"); lines(sy_hat)
plot(d_hat, type = 'h')
plot(wl_hat, type = 'l')
plot(sy_hat, type = 'l')
plot(WL[1] + cumsum(d_hat), type = 'l'); lines(WL, col = 'gray40')
lines(Dobs, lty = 3)
lines(WL[1] + as.numeric(filter(WA, 0.0198785393130021, "rec", sides = 1)), type = 'l')
plot(Dobs)

library(brms)
options(mc.cores = 4)

brm_test <- 
  daily_water_levels[site == "152" & sample_year %in% c(2012, 2018) & !is.na(wl_initial_cm), 
       .(wlObs = wl_initial_cm,
           wlL1 = wl_l1_cm,
           pet = pet_cm,
           precip = best_precip_cm)]

# Original Formula
brm_form <- 
  brmsformula(wlObs ~ wlL1 + 
                (precip / (step(cp - wlL1) * (bW + mW * wlL1) + step(wlL1 - cp) * (bW2))) -
                (pet / (step(cp - wlL1) * (bD + mD * wlL1) + step(wlL1 - cp) * (bD2)) + 
                   step(wlL1 - cp) * (bQ + mQ * wlL1)),
              cp + bW + mW + bW2 + bD + mD + bD2 + bQ + mQ ~ 1,
              nl = TRUE)

brm_prior <- 
  prior(normal(0, 10), nlpar = "cp", lb = 0) +
  prior(normal(0.1, 1), nlpar = "bW", lb = 0) +
  prior(normal(0.002, 0.1), nlpar = "mW", lb = 0) +
  prior(normal(0.1, 1), nlpar = "bW2") +
  prior(normal(0.1, 1), nlpar = "bD") +
  prior(normal(0.002, 0.1), nlpar = "mD", lb = 0) +
  prior(normal(0.1, 1), nlpar = "bD2") +
  prior(normal(-10, 10), nlpar = "bQ", ub = 0) +
  prior(normal(-0.1, 1), nlpar = "mQ", ub = 0) +
  prior(gamma(2, .1), class = nu)

# The scale does not match. The trends are good, but the scales don't match,
# weird. It may be good to switch to water_availabilty rather than P & PET
# because it would work well between years to capture when drawdown begins. Or
# if not between years, it may provide a signal for when drawdown starts

wl_mod <- 
  brm(brm_form,
      data = brm_test,
      prior = brm_prior,
      family = student,
      iter = 500,
      warmup = 100,
      seed = 1234,
      chains = 1)


# wl ~ wlL1 + d_in - d_out
# d_in ~ precip / sy_in
# sy_in ~ ((wlL1 < cp) * (bW + mW * wlL1) + (wlL1 >= cp) * (mW2 * wlL1))
# 
# d_out ~ (pet / sy_out) + Q
# sy_out ~ (wlL1 < cp) * (bD + mD * wlL1) + (wlL1 >= cp) * (mD2 * wlL1) 
# Q ~ (wlL1 >= cp) * (bQ + bM * wlL1)

val_dat <- 
  daily_water_levels[site == "152" & sample_year %in% c(2012:2013), 
                     .(sample_date, 
                       wlObs = wl_initial_cm,
                       wlL1 = wl_l1_cm,
                       pet = pet_cm,
                       precip = best_precip_cm)][min(which(!is.na(wlObs))):.N]

wl_hat2 <- 
  numeric(nrow(val_dat))

wl_hat2[1] <- 
  val_dat$wlObs[1]

precip <- 
  val_dat$precip

pet <- 
  val_dat$pet

q <- sy_in <- sy_out <- d_in <- d_out <- numeric(length(wl_hat2))

wl_coefs <- 
  fixef(wl_mod)

cp <- wl_coefs[1, 1]
bW <- wl_coefs[2, 1]
mW <- wl_coefs[3, 1]
bW2 <- wl_coefs[4, 1]
bD <- wl_coefs[5, 1]
mD <- wl_coefs[6, 1]
bD2 <- wl_coefs[7, 1]
bQ <- wl_coefs[8, 1]
mQ <- wl_coefs[9, 1]


for(i in 2:length(wl_hat2)){
  
  sy_in[i] <- 
    ((wl_hat2[i-1] < cp) * (bW + mW * wl_hat2[i-1]) + (wl_hat2[i-1] >= cp) * (bW2))
  
  d_in[i] <- 
    (precip[i-1] / sy_in[i])
  
  sy_out[i] <- 
    ((wl_hat2[i-1] < cp) * (bD + mD * wl_hat2[i-1]) + (wl_hat2[i-1] >= cp) * (bD2)) 
  
  q[i] <-  
    (wl_hat2[i-1] >= cp) * (bQ + mQ * wl_hat2[i-1])
  
  d_out[i] <- 
    (pet[i-1] / sy_out[i]) - q[i]
  
  wl_hat2[i] <- wl_hat2[i-1] + d_in[i] + d_out[i]

  if(is.na(wl_hat2[i])){
    stop("wl_hat2 is NA at i = ", i)
  }
  
}

plot(val_dat$wlObs, col = 'gray40', type = 'l'); lines(wl_hat2)
plot(scale(val_dat$wlObs[1:244]), col = 'gray40', type = 'l'); lines(scale(wl_hat2[1:244]))

# Big Bayesian Model ------------------------------------------------------


val_dat <- 
  daily_water_levels[site == "152" & sample_year %in% c(2012:2014), 
                     .(sample_date, 
                       wlObs = wl_initial_cm,
                       wlL1 = wl_l1_cm,
                       pet = pet_cm,
                       precip = best_precip_cm)][min(which(!is.na(wlObs))):.N]

wb_form <- 
  bf(wlObs ~ wlL1 + p + et + q + g,
     nlf(p ~ precip / esyP),
     nlf(et ~ step(cp - wlL1) * (Mpet * pet / esyET)),
     # This could have a below threshold component that response to previous day's precip
     nlf(q ~ step(wlL1 - cp) * (Bq + Mq * wlL1^2)),
     nlf(g ~ step(cp - wlL1) * Bg),
     nlf(esy ~ Besy + Mesy ^ (step(cp - wlL1) * (wlL1 - cp))),
     cp + Mpet + Bg + Bq + Mq + Besy + Mesy ~ 1,
     nl = TRUE)

wb_priors <- 
  # Not sure how to set priors on the intermediate values, got an error that
  # they are not parameters in the model
  # prior(gamma(0.16, 1.54), nlpar = "p", lb = 0) + # determined via lmomco for brm_test$precip_cm
  # prior(normal(5, 10), nlpar = "et", ub = 0) +
  # prior(gamma(0.5, 1), nlpar = "q", ub = 0) +
  # prior(normal(5, 10), nlpar = "g", lb = 0) +
  # prior(normal(1, 1), nlpar = "Mp", lb = 0) +
  # prior(normal(0.5, 1), nlpar = "esy", lb = 0) +
  prior(normal(5, 1), nlpar = "cp", lb = 0) +
  prior(constant(-1), nlpar = "Mpet", ub = 0) +
  prior(normal(1, 10), nlpar = "Bq", ub = 0) +
  prior(normal(0.02, 1), nlpar = "Mq", ub = 0) +
  prior(normal(1, 5), nlpar = "Bg", lb = 0) +
  # prior(normal(0.02, 1), nlpar = "Mg", lb = 0) +
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
    
    p <- et <- q <- g <- esy <- numeric(nrow(val_dat))
    
    # Mp <- wb_coefs["Mp_Intercept", 1]
    cp <- wb_coefs["cp_Intercept", 1]
    Mpet <- wb_coefs["Mpet_Intercept", 1]
    Bq <- wb_coefs["Bq_Intercept", 1]
    Mq <- wb_coefs["Mq_Intercept", 1]
    Bg <- wb_coefs["Bg_Intercept", 1]
    # Mg <- wb_coefs["Mg_Intercept", 1]
    Besy <- wb_coefs["Besy_Intercept", 1]
    Mesy <- wb_coefs["Mesy_Intercept", 1]
    
    esy[1] <- 
      Besy + Mesy ^ ((cp > wl_hat2[1]) * (wl_hat2[1] - cp))
    
    next
  }
  
  p[i] <- 
    precip[i-1] / esy[i-1]
  
  et[i] <- 
    (cp > wl_hat2[i-1]) * (Mpet * pet[i-1] / esy[i-1])  
  
  q[i] <- 
    (cp <= wl_hat2[i-1]) * (Bq + Mq * wl_hat2[i-1]^2)
  
  g[i] <-
    (cp <= wl_hat2[i-1]) * (Bg) # + Mg * wl_hat2[i-1])
  
  wl_hat2[i] <- 
    wl_hat2[i-1] + p[i] + et[i] + q[i] + g[i]

  esy[i] <- 
    Besy + Mesy ^ ((cp > wl_hat2[i]) * (wl_hat2[i] - cp))
  
}

plot(wlObs ~ sample_date, data = val_dat, col = 'gray40', type = 'l'); lines(wl_hat2 ~ val_dat$sample_date)


predict_water_levels <- 
  function(wl.obs,
           water.availability,
           esy.model){
    
    n <- 
      length(water.availability) + 1
    
    d_hat <- 
      numeric(n)
    
    init <- 
      1
    
    if(length(wl.obs) > 1){
      init <- 
        min(which(!is.na(wl.obs)))
    }

    d_hat[init] <- 
      0
    
    for(i in 2:n){
      d_hat[i] <- 
        (f_esy(d_hat[i-1]) * WA[i-1])
    }
    
    wl.obs[init] + d_hat[-c(1)]
  }
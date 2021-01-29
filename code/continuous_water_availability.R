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
                   eos = as.Date(paste0(.BY[[2]], "-11-01")),
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

external_met[,
             normalized_wa_cm := (water_availability_cm - first(water_availability_cm)) - max(water_availability_cm - first(water_availability_cm)),
             by = .(station_name, water_year)]


ggplot(external_met[sample_year %in% 2012:2020], 
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
                   by = .(site, water_year)]

daily_water_levels[, Ds_cumulative := ytd_sum(Ds_continuous),
                   by = .(site, water_year)]


ggplot(daily_water_levels,
       aes(x = sample_date,
           y = Ds_cumulative)) +
  geom_line() +
  facet_grid(site ~ water_year, 
             scales = "free")

daily_water_levels[external_met[station_name == "kenton"],
                   `:=`(detrend_wa_cm = i.detrend_wa_cm,
                        normalized_wa_cm = i.normalized_wa_cm),
                   on = c("dowy", "water_year")]


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

daily_water_levels[, normalized_Ds := Ds_cumulative - max_na(Ds_cumulative),
                   by = .(site, water_year)]
daily_water_levels[, esy_wa_emp := Ds_cumulative/normalized_wa_cm]
daily_water_levels[, wl_l1_cm := shift(wl_initial_cm, -1),
                   by = .(site, water_year)]

# MCP DID HORRIBLY (with esy_wa_emp made using detrend_wa_cm & wl_initial_cm as X)
# Ds_cumulative as X does great, but that's because then it's in X & Y

esy_wa_mods <-
  daily_water_levels[is.finite(esy_wa_emp), 
                     .(mod_lm = list(lmrob(esy_wa_emp ~ wl_l1_cm,
                                           data = .SD,
                                           setting = "KS2014"))
                       , mod_quad = list(lmrob(esy_wa_emp ~ wl_l1_cm + I(wl_l1_cm^2),
                                               data = .SD,
                                               setting = "KS2014"))
                       # , mod_mcp = list(mcp(list(esy_wa_emp ~ wl_initial_cm,
                       #                         ~ 0 + wl_initial_cm),
                       #                    data = .SD))
                       ),
                     by = .(site)]


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
daily_water_levels[, Ds_hat_lm := normalized_wa_cm * esy_wa_lm,
                   by = .(site)]

daily_water_levels[, esy_wa_quad := esy_wa_mods[.BY[[1]], f_quad[[1]]](wl_l1_cm),
                   by = .(site)]
daily_water_levels[, Ds_hat_quad := normalized_wa_cm * esy_wa_quad,
                   by = .(site)]

# daily_water_levels[, esy_wa_mcp := esy_wa_mods[.BY[[1]], f_mcp[[1]]](wl_l1_cm),
#                    by = .(site)]
# daily_water_levels[, Ds_hat_mcp := detrend_wa_cm * esy_wa_mcp,
#                    by = .(site)]

ggplot(daily_water_levels,
       aes(x = wl_l1_cm,
           y = esy_wa_emp)) +
  geom_point() +
  geom_line(aes(y = esy_wa_lm), 
            color = "red") +
  geom_line(aes(y = esy_wa_quad),
            color = "blue") +
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

ggplot(daily_water_levels,
       aes(x = Ds_cumulative,
           y = Ds_hat_lm)) +
  geom_point() +
  geom_abline(color = "red") +
  facet_grid(site ~ sample_year, 
             scales = "free")

test <- 
  daily_water_levels[site == "135" & water_year == 2017 & !is.na(wl_initial_cm)]

f_esy <- 
  esy_wa_mods["135", f_lm[[1]]]

m <- 
  esy_wa_mods["135", coef(mod_lm[[1]])[[2]]]
  
WL <- 
  test$wl_initial_cm

WA <- 
  test$normalized_wa_cm

SY <- 
  test$esy_wa_emp

d_hat <- 
  numeric(nrow(test))

Dobs <- 
  test$Ds_cumulative

d_hat[1] <- 
  0

wl_hat <- 
  numeric(nrow(test))

wl_hat[1] <- 
  0

for(i in 2:length(d_hat)){
  d_hat[i] <- 
    (f_esy(wl_hat[i-1]) * WA[i-1])
  
  wl_hat[i] <- 
    wl_hat[1] + d_hat[i]
}

plot(WL, type = 'l', col = "gray40")
lines(WL[1] + wl_hat)
plot(wl_hat, type = 'l')
plot(d_hat, type = 'l')
lines(Dobs, lty = 3)
lines(WL[1] + as.numeric(filter(WA, 0.0198785393130021, "rec", sides = 1)), type = 'l')
plot(Dobs)

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
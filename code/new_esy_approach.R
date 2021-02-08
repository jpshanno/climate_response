# New ESy
tar_load(external_met)
tar_load(daily_water_levels)

external_met[, gdd_10 := (tmax_c + tmin_c) / 2 - 10]

# Drop leap-day so all doys match across year
regional_met <- 
  external_met[format(sample_date, "%m%d") != "0229",  
               lapply(.SD, median, na.rm = TRUE),
               by = .(sample_date),
               .SDcols = c("pet_cm", "precip_cm", "gdd_10")]

# Using DIN-4049 Water year (Germany & Austria)
# ggplot(tr_dat, aes(x = doy, y = wl_initial_cm)) +
#   geom_line(stat = "summary", fun = median, color = "blue") +
#   geom_line(data = regional_met, aes(x = dowy, y = -20 * roll_pet_cm)) +
#   facet_wrap(~sample_year)
regional_met[, sample_year := year(sample_date)]
regional_met[, water_year := ifelse(sample_date >= as.Date(paste0(.BY[[1]], "-11-01")),
                                    .BY[[1]] + 1,
                                    .BY[[1]]),
             by = .(sample_year)]

setorder(regional_met, "sample_date")

regional_met[,
             doy := 1:.N,
             by = .(sample_year)]
regional_met[sample_year == 2011,
             doy := yday(sample_date)]
regional_met[,
             dowy := 1:.N,
             by = .(water_year)]

# Cumulative Values by water year
regional_met[order(water_year, dowy), 
             `:=`(wytd_gdd_10 = cumsum(gdd_10),
                  wytd_precip_cm = cumsum(precip_cm),
                  wytd_pet_cm = cumsum(pet_cm)),
             by = .(water_year)]

regional_met[order(sample_year, doy), 
             `:=`(ytd_gdd_10 = cumsum(gdd_10),
                  ytd_precip_cm = cumsum(precip_cm),
                  ytd_pet_cm = cumsum(pet_cm)),
             by = .(sample_year)]

regional_met[, smooth_pet_cm := smooth.spline(pet_cm)$y]

regional_met[,
             `:=`(wy_water_availability_cm = wytd_precip_cm - wytd_pet_cm)]

regional_met[,
             `:=`(water_availability_cm = ytd_precip_cm - ytd_pet_cm)]

regional_met[, smooth_wa_cm := smooth.spline(water_availability_cm)$y]

regional_met[, 
             `:=`(rot_wa_cm = water_availability_cm - doy * last(water_availability_cm) / 365,
                  rot_pet_cm = ytd_pet_cm - doy * last(ytd_pet_cm) / 365),
             by = .(sample_year)]

daily_water_levels[regional_met, 
                   `:=`(regional_pet_cm = i.pet_cm,
                        regional_p_cm = i.precip_cm,
                        ytd_precip_cm = i.ytd_precip_cm,
                        wytd_precip_cm = i.wytd_precip_cm,
                        ytd_pet_cm = i.ytd_pet_cm,
                        wytd_pet_cm = i.wytd_pet_cm,
                        wy_water_availability_cm = i.wy_water_availability_cm,
                        water_availability_cm = i.water_availability_cm,
                        rot_wa_cm = i.rot_wa_cm),
                   on = "sample_date"]

# Just linear interception
# Using gross precip shows problems at lower water levels with much higher 
# Esy values. But previous interception estimates seemed to overestimate 
# interception and wreck the exponential shape of the Esy ~ WL curve
new_i_lm <-
  daily_water_levels[obs_precip_cm > 0,
                     .(mod = list(lmrob(Dl_signed_cm ~ obs_precip_cm,
                                        data = .SD,
                                        setting = "KS2014"))),
                     keyby = .(site_status)]

new_i_lm[, i_lm_cm := map_dbl(mod, ~{-coef(.x)[[1]] / coef(.x)[[2]]})]

daily_water_levels[new_i_lm[, .(site_status, i_lm_cm)], i_lm_cm := i.i_lm_cm, on = "site_status"]

# Old Interception estimates
tar_load(mod_interception)
daily_water_levels[, i_old_cm := mod_interception[.BY[[1]], f_predict[[1]]](doy),
                   by = .(site_status)]

daily_water_levels[, net_precip_cm := pmax(0, best_precip_cm - i_lm_cm)]

daily_water_levels[, Ds_resid := Ds_cm - net_precip_cm + regional_pet_cm - melt_cm]

daily_water_levels[, wl_bin := cut_interval(wl_initial_cm,
                                            n = uniqueN(round(wl_initial_cm, 0))/2),
                   by = .(site)]

reduced_dat <-
  daily_water_levels[Ds_resid < 5 & wl_initial_cm > -10,
                     .(Ds_resid = quantile(Ds_resid, 0.25, na.rm = TRUE),
                       wl_initial_cm = median_na(wl_initial_cm)),
                     by = .(site, wl_bin)]

ggplot(daily_water_levels[wl_initial_cm > -10], 
       aes(x = wl_initial_cm, y = Ds_resid)) + 
  geom_point() + 
  geom_point(data = reduced_dat,
             color = "red") +
  facet_wrap(~site, scales = "free")


flow_mods <- 
  reduced_dat[, .(mod = list(mcp(list(Ds_resid ~ 1, ~0 + I(wl_initial_cm^2)), 
                                 data = .SD))), 
              keyby = .(site)]

# map(flow_mods$mod, 
#     plot,
#     q_predict = TRUE) %>% 
#   reduce(`+`)

flow_mods[, coefs := map(mod, mcp::fixef)]
flow_mods[, f_predict := map(coefs,
                             ~{as.function(list(water.level = NULL,
                                                bquote({
                                                  ifelse(water.level < .(cp),
                                                         0,
                                                         .(b) + .(m) * (water.level - .(cp))^2)
                                                },
                                                where = list(cp = .x[1, "mean"],
                                                             b = .x[2, "mean"],
                                                             m = .x[4, "mean"]))))})]

daily_water_levels[, q_cm := flow_mods[.BY[[1]], f_predict[[1]]](wl_initial_cm),
                   by = .(site)]

ggplot(daily_water_levels[wl_initial_cm > -10], 
       aes(x = wl_initial_cm, y = Ds_resid)) + 
  geom_point() + 
  geom_point(data = reduced_dat,
             color = "red") +
  geom_line(aes(y = q_cm),
             color = 'blue') +
  facet_wrap(~site, scales = "free")

daily_water_levels[, esy_emp := (best_precip_cm + melt_cm - q_cm) / Dl_signed_cm]

daily_water_levels[, wl_ub := quantile(.SD[best_precip_cm > i_lm_cm & between(esy_emp, 0, 2.5), 
                                           wl_min_cm],
                                       0.975, 
                                       na.rm = TRUE), 
                   by = .(site)]

daily_water_levels[, min_esy := min_na(.SD[best_precip_cm > i_lm_cm & between(esy_emp, 0, 2.5) & wl_initial_cm < wl_ub, esy_emp]),
                   by = .(site)]

esy_mods <- 
  daily_water_levels[best_precip_cm > i_lm_cm & between(esy_emp, 0, 2.5) & wl_initial_cm < wl_ub, 
                     .(mod = list(nls(esy_emp ~ min_esy + m^(wl_min_cm - c),
                                      data = .SD,
                                      start = list(m = 1.5, c = -50),
                                      control = nls.control(warnOnly = TRUE, maxiter = 200)))),
                     keyby = .(site)]

# This would need a limit to avoid esy < 0
esy_mods[daily_water_levels[best_precip_cm > i_lm_cm & between(esy_emp, 0, 2.5) & wl_initial_cm < wl_ub, 
                            .(mod = list(mcp(list(esy_emp ~ 1, ~ 0 + I(wl_min_cm ^ 2)),
                                             # prior = list(int_1 = .SD[1, min_esy]),
                                             data = .SD))),
                            keyby = .(site)], 
         mcp_mod := i.mod]

esy_mods[, mcp_predict := map(mcp_mod,
                              ~{
                                coefs <- 
                                  mcp::fixef(.x)
                                
                                as.function(list(wl = NULL,
                                                 bquote(.(b) + 
                                                          (wl > .(cp)) * (.(m2) * (wl - .(cp))^2),
                                                        where = list(cp = coefs[1, "mean"],
                                                                     b = coefs[2, "mean"],
                                                                     # m1 = coefs[4, "mean"],
                                                                     m2 = coefs[4, "mean"]))))
                              })]

esy_mods[daily_water_levels[best_precip_cm > i_lm_cm & between(esy_emp, 0, 2.5) & wl_initial_cm < wl_ub & !is.na(wl_min_cm), 
                            .(mod = list(nls(esy_emp ~ b + m*exp(d*wl_min_cm),
                                             data = .SD,
                                             start = list(b = .SD[1, min_esy], m = 0.2, d = 0.1),
                                             algorithm = "port",
                                             lower = list(b = .SD[1, min_esy], m = 0, d = 0),
                                             control = nls.control(warnOnly = TRUE, maxiter = 200)))),
                            keyby = .(site)], 
         exp_mod := i.mod]

daily_water_levels[,
                   esy_mcp := esy_mods[.BY[[1]], mcp_predict[[1]]](wl_min_cm), 
                   by = .(site)]

daily_water_levels[,
                   esy_hat := predict(esy_mods[.BY[[1]], mod[[1]]], newdata = .SD),
                   by = .(site)]

daily_water_levels[,
                   esy_exp:= predict(esy_mods[.BY[[1]], exp_mod[[1]]], newdata = .SD),
                   by = .(site)]

ggplot(daily_water_levels[best_precip_cm > i_lm_cm & between(esy_emp, 0, 2.5) & wl_min_cm < wl_ub], 
       aes(x = wl_min_cm, y = esy_emp)) +
  geom_point() + 
  geom_line(aes(y = esy_mcp),
            color = "blue") +
  geom_line(aes(y = esy_hat),
            color = "red") +
  geom_line(aes(y = esy_exp),
            color = "green", linetype = "dashed") +
  facet_wrap(~site, scales = "free")


ggplot(daily_water_levels,
       aes(x = water_availability_cm / esy_hat,
           y = wl_initial_cm)) + 
  geom_point(aes(color = factor(sample_year))) +
  facet_wrap(~site, scales = "free")


ggplot(daily_water_levels[sample_year == 2013],
       aes(x = water_availability_cm,
           y = wl_initial_cm)) + 
  geom_path(aes(color = doy)) +
  facet_wrap(~site, scales = "free")

daily_water_levels[!is.na(wl_initial_cm), 
                   ytd_precip_sy := cumsum((ytd_precip_cm - shift(ytd_precip_cm, 1, 0)) / esy_brm),
                   by = .(site, sample_year)]

daily_water_levels[!is.na(wl_initial_cm), 
                   ytd_pet_sy := cumsum((ytd_pet_cm - shift(ytd_pet_cm, 1, 0)) / esy_brm),
                   by = .(site, sample_year)]

daily_water_levels[, wl_offset := first(na.omit(wl_initial_cm - ytd_precip_sy + ytd_pet_sy)),
                   by = .(site, sample_year)]

# Esy pooled is created from running nls() model on var_dat below
# var_dat <-
#   daily_water_levels[best_precip_cm == 0 & !(melt_cm > 0),
#                      .(esy_bounds = qnorm(0.975) * sd((-regional_pet_cm) / Dl_signed_cm, na.rm = TRUE),
#                        wl_initial_cm = median_na(wl_initial_cm)),
#                      by = .(site, wl_bin = cut_interval(wl_initial_cm,
#                                                         n = uniqueN(round(wl_initial_cm, 0))/5))]
# daily_water_levels[, esy_pooled := pmin(1, 0.14199 + 0.39584*(exp(0.03978*wl_initial_cm)))]

ggplot(daily_water_levels[sample_year == 2012],
       aes(x = doy,
           y = wl_initial_cm)) + 
  geom_line() +
  geom_line(color = "blue",
            aes(y = wl_offset + ytd_precip_sy - ytd_pet_sy)) +
  facet_wrap(~site, scales = "free")




# New PET ESY -------------------------------------------------------------

ggplot(daily_water_levels[daily_water_levels[, .(sample_date, mdp = frollsum(best_precip_cm, 3, align = "right")), by = .(site, sample_year)], on = c("sample_date", "site", "sample_year")][mdp == 0 & wl_max_cm > wl_final_cm, .(doy, site, wl_max_cm, regional_pet_cm, pet_response = wl_max_cm - wl_final_cm)], aes(x = wl_max_cm, y = regional_pet_cm / pet_response)) + geom_point() + facet_wrap(~site, scales = "free_x") + coord_cartesian(ylim = c(-0, 2.5))

pet_esy_dat <- 
  daily_water_levels[
    daily_water_levels[, 
                       .(sample_date, mdp = frollsum(regional_p_cm, 3, align = "right")), 
                       by = .(site, sample_year)], 
    on = c("sample_date", "site", "sample_year")
    ][mdp == 0 & wl_max_cm > wl_final_cm, 
      .(doy, site, wl_initial_cm, regional_pet_cm, 
        pet_response = wl_max_cm - wl_final_cm,
        pet_sy = regional_pet_cm / (wl_max_cm - wl_final_cm))
      ][pet_sy <= 1.5]

library(brms)
options(mc.cores = 4)

brm_form <- 
  brmsformula(nl = TRUE,
              pet_sy ~ (wl_initial_cm <= cp) * (b + m*exp(d*wl_initial_cm)) + (wl_initial_cm > cp) * b0,
              cp + b + m + d + b0~ 1)

brm_priors <- 
  prior(nlpar = "cp", normal(5, 1)) +
  prior(nlpar = "b", normal(0, 1), lb = 0.02) +
  prior(nlpar = "m", normal(0.2, 1), lb = 0) +
  prior(nlpar = "d", normal(0.1, 1), lb = 0) + 
  prior(nlpar = "b0", normal(0.5, 1), lb = 0)

brm_pet <- 
  pet_esy_dat[, 
              .(mod = list(brm(brm_form,
                               prior = brm_priors, 
                               data = .SD))),
              keyby = .(site)]

brm_pet[, f_predict := map(mod,
                           ~{
                             coefs <- 
                               fixef(.x)[, "Estimate"]
                             
                             as.function(
                               list(water.level = NULL,
                                    bquote(
                                      (water.level <= .(cp)) * (.(b) + .(m)*exp(.(d)*water.level)) + (water.level > .(cp)) * .(b0),
                                      where = 
                                        list(cp = coefs[[1]],
                                             b = coefs[[2]],
                                             m = coefs[[3]],
                                             d = coefs[[4]],
                                             b0 = coefs[[5]])
                                    ))
                             )
                           })]

saveRDS(brm_pet, "tmp/brms_pet_esy.rds")

brm_pet[, f_predict := map(mod,
                           ~{
                             coefs <- 
                               fixef(.x)[, "Estimate"]
                             
                             as.function(
                               list(water.level = NULL,
                                    bquote(
                                      (water.level <= .(cp)) * (.(b) + .(m)*exp(.(d)*water.level)) + (water.level > .(cp)) * .(b0),
                                      where = 
                                        list(cp = coefs[[1]],
                                             b = coefs[[2]],
                                             m = coefs[[3]],
                                             d = coefs[[4]],
                                             b0 = coefs[[5]])
                                    ))
                             )
                           })]

brm_pet[, f_predict := map(mod,
                           ~{
                             coefs <- 
                               fixef(.x)[["Estimate"]]
                             
                             as.function(
                               list(water.level = NULL,
                                    bquote(
                                      (water.level <= .(cp)) * (.(b) + .(m)*exp(.(d)*water.level)) + (water.level > .(cp)) * .(b0),
                                      where = 
                                        list(cp = coefs[1],
                                             b = coefs[2],
                                             m = coefs[3],
                                             d = coefs[4],
                                             b0 = coefs[5])
                                    ))
                             )
                           })]

pet_esy_dat[, esy_brm := brm_pet[.BY[[1]], f_predict[[1]]](wl_initial_cm),
                   by = .(site)]

daily_water_levels[, esy_brm := brm_pet[.BY[[1]], f_predict[[1]]](wl_initial_cm),
                   by = .(site)]

# New Precip ESY ----------------------------------------------------------

ggplot(daily_water_levels[, .(wl_bin, esy_pooled, wl_min_cm, doy, best_precip_cm = frollsum(best_precip_cm, n = 3, align = "right"), response = pmax(wl_max_cm, shift(wl_max_cm, 1), shift(wl_max_cm, 2) ) - wl_min_cm, sample_year), by = .(site)][best_precip_cm > 0.04], aes(x = wl_min_cm, y = best_precip_cm / response)) + geom_point() + geom_smooth(method = "glm", method.args = list(family = Gamma(log))) + facet_wrap(~site, scales = "free")

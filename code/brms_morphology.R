# Changepoint Models in brms from 
# https://discourse.mc-stan.org/t/piecewise-linear-mixed-models-with-a-random-change-point/5306
library(brms)
options(mc.cores = 3)

# Filtering by esy < 5 barely changes the medians
# water_budget[, esy := (-pet_cm - Q? net_precip_cm + melt_cm) / (wl_final_cm - wl_initial_cm)]


tar_load(water_budget)

water_budget[, `:=`(Ds_cm = Ds_cm / sy,
                    Dl_signed_cm = Dl_signed_cm / sy,
                    sy = NULL)]

water_budget[, net_flow_cm := Ds_cm + pet_cm - net_precip_cm - melt_cm]

water_budget[, wl_bin := cut_interval(wl_initial_cm,
                                      n = uniqueN(round(wl_initial_cm, 0))/3), 
             by = .(site)]

# Reduced data to median flow for each centimeter of water level to remove 
# noise
reduced_dat <- 
  water_budget[net_flow_cm < 3, 
               .(net_out = quantile(net_flow_cm,
                                    probs = 0.1,
                                    na.rm = TRUE),
                 net_in = quantile(net_flow_cm,
                                   probs = 0.9,
                                   na.rm = TRUE),
                 wl_initial_cm = median_na(wl_initial_cm)), 
               by = .(site, wl_bin)]

# Remove wl_bin from original data
# water_budget[, wl_bin := NULL]

reduced_dat <- 
  reduced_dat[!is.na(net_in + net_out) & !is.na(wl_initial_cm)]

ggplot(water_budget,
       aes(x = wl_initial_cm,
           y = net_flow_cm)) + 
  geom_point() +
  geom_point(data = reduced_dat,
             aes(y = net_in),
             color = "blue",
             shape = 20) +
  geom_point(data = reduced_dat,
             aes(y = net_out),
             color = "brown",
             shape = 20) +
  facet_wrap(~site,
             scales = "free")

# The model
q_form <- bf(
  # Section 1 of piece-wise regression (not strictly necessary to specify)
  net_flow_cm ~ 0 + 
    # Section 2 of piece-wise regression
    (slope2 * (wl_initial_cm - change)^2) * step(wl_initial_cm - change),
  # Allow intercept, slope1, slope2, and change point to vary with site
  slope2 + change ~ 1 + (1|site),
  # Specify non-linear formula (ie `*` means multiply not crossed-effects and
  # exponents can be used directly)
  nl = TRUE
)

# Priors (should have a limit on slope to make sure it's negative. Just specifying
# ub = 0 doesn't ensure that the group-level effects are negative
q_prior <- 
  prior(normal(0, 1), nlpar = "slope2") +
  prior(normal(5, 2), nlpar = "change")

q_mod <- 
  brm(q_form, 
      data = reduced_dat, 
      prior = q_prior, 
      chains = 3)

q_coefs <- 
  water_budget[, setNames(as.list(fixef(q_mod)[, "Estimate"] + ranef(q_mod)$site[.BY[[1]], "Estimate", ]),
                          c("m", "cp")),
               keyby = .(site)]

q_coefs[, f_predict := map2(m, cp,
                            ~as.function(list(water.level = NULL, 
                                              bquote(.(quadm) * (water.level - .(change))^2 * as.numeric(.(change) < water.level), 
                                                     where = list(quadm = .x, change = .y)))))]

water_budget[, q_hat := q_coefs[.BY[[1]], f_predict[[1]]](wl_initial_cm), 
             by = .(site)]

water_budget[, meteo_flow_cm := net_flow_cm - q_hat]


ggplot(water_budget,
       aes(x = wl_initial_cm,
           y = net_flow_cm)) + 
  geom_point() +
  geom_line(aes(y = q_hat),
             color = "red",
             shape = 20) +
  facet_wrap(~site,
             scales = "free")

water_budget[, esy := (-pet_cm + net_precip_cm + melt_cm) / (wl_max_cm - wl_min_cm)]

reduced_sy <- 
  water_budget[(pet_cm) < (net_precip_cm + melt_cm) & !is.na(pet_cm + net_precip_cm + melt_cm), 
               .(esy = median_na(esy)), 
               by = .(site, wl_min_cm = round(wl_min_cm, 0))]

reduced_sy <- reduced_sy[!is.na(esy)]

sy_mods <- 
  water_budget[(pet_cm) < (net_precip_cm + melt_cm) & !is.na(pet_cm + net_precip_cm + melt_cm + esy), 
             .(mod = list(nlrob(esy ~ b + m ** (wl_min_cm - c),
                                data = .SD,
                                start = list(b = min(.SD$esy, na.rm = TRUE), m = 1.25, c = -50),
                                algorithm = "port",
                                maxit = 100,
                                lower = list(b = min(.SD$esy, na.rm = TRUE), m = 0, c =-100000),
                                control = nls.control(warnOnly = TRUE, maxiter = 100)))),
             keyby = .(site)]

sy_mods[, coefs := map(mod, coef)]
sy_mods[, f_predict := map(coefs,
                           ~as.function(list(water.level = NULL,
                                             bquote(.(B) + .(M)**(water.level-.(C)),
                                                    where = list(B = .x[["b"]],
                                                                 M = .x[["m"]],
                                                                 C = .x[["c"]])))))]

water_budget[, esy_hat := sy_mods[.BY[[1]], f_predict[[1]]](wl_initial_cm),
             by = .(site)]

ggplot(water_budget[(pet_cm) < (net_precip_cm + melt_cm)],
       aes(x = wl_initial_cm,
           y = esy)) + 
  geom_point() +
  geom_line(aes(y = esy_hat),
            color = "red") +
  geom_line(aes(y = pmin(sy_orig, 1)), 
            color = "blue",
            linetype = "dashed") +
  facet_wrap(~site,
             scales = "free")


water_budget[, Ds_esy := Ds_cm * esy_hat]

water_budget[, net_flow_cm := Ds_esy + pet_cm - net_precip_cm - melt_cm]



ggplot(water_budget,
       aes(x = wl_initial_cm,
           y = Ds_esy + pet_cm - net_precip_cm - melt_cm - slow_flow_cm)) + 
  geom_point() +
  facet_wrap(~site,
             scales = "free")

ggplot(water_budget,
       aes(x = pet_cm,
           y = Ds_cm)) + 
  geom_point() +
  geom_point(aes(y = Ds_cm * esy_hat),
             color = "red") +
  facet_wrap(~site,
             scales = "free") +
  coord_cartesian(ylim = c(-3, 3))


water_budget[abs((-pet_cm + net_precip_cm + melt_cm) / (wl_final_cm - wl_initial_cm)) < 25, 
             .(esy = median_na(abs((-pet_cm + net_precip_cm + melt_cm) / (wl_final_cm - wl_initial_cm))),
               wl_initial_cm = median_na(wl_initial_cm)), 
             by = .(site,
                    wl_bin = cut_interval(wl_initial_cm,
                                          n = uniqueN(round(wl_initial_cm, 0))/5))] %>% 
  ggplot(aes(x = wl_initial_cm,
             y = esy)) + 
  geom_point() +
  geom_smooth(method = lm,
              formula = y ~ x + I(x^2),
              se = FALSE) +
  facet_wrap(~site,
             scales = "free")





water_budget[, meteo_flow_cm := net_flow_cm - q_hat]
water_budget[, et_scaled := as.numeric(scale(et_cm)), by = .(site)]


ggplot(water_budget[between(meteo_flow_cm, -2, 4)],
       aes(x = pet_cm,
           y = meteo_flow_cm)) + 
  geom_point() +
  facet_wrap(~site,
             scales = "free")

meteo_form <- 
  bf(
    meteo_flow_cm ~ Intercept + petslope * pet_cm,
    # Allow intercept, slope1, slope2, and change point to vary with site
    Intercept + petslope ~ 1 + (1|site),
    # Specify non-linear formula (ie `*` means multiply not crossed-effects and
    # exponents can be used directly)
    nl = TRUE
  )

# Priors
meteo_priors <- 
  prior(normal(0, 5), nlpar = "Intercept") +
  prior(normal(0, 2), nlpar = "petslope")

mod_meteo <- 
  water_budget[, .(mod = list(brm(meteo_form, 
                                  data = .SD, 
                                  prior = meteo_priors, 
                                  chains = 3))),
               keyby = .(site_status)]
  
options(mc.cores = 1)

meteo_coefs <- 
  water_budget[, setNames(as.list(fixef(mod_meteo[.BY[[2]], mod[[1]]])[, "Estimate"] + ranef(mod_meteo[.BY[[2]], mod[[1]]])$site[.BY[[1]], "Estimate", ]),
                          c("intercept", "m_pet")),
               keyby = .(site, site_status)]

meteo_coefs[, f_predict := pmap(list(intercept, m_pet), 
                             function(intercept, m_pet){
                               as.function(list(pet = NULL, 
                                                bquote(.(b) + .(m) * pet, 
                                                       where = list(b = intercept, 
                                                                    m = m_pet))))})]


water_budget[, meteo_hat := meteo_coefs[CJ(.BY[[1]], .BY[[2]]), f_predict[[1]]](pet_cm),
             by = .(site, site_status)]

ggplot(water_budget,
       aes(x = meteo_flow_cm,
           y = meteo_hat)) + 
  geom_point() +
  geom_abline(color = "red") +
  facet_wrap(~site,
             scales = "free")

water_budget[, net_flow_hat := q_hat + meteo_hat]

ggplot(water_budget,
       aes(x = net_flow_cm,
           y = net_flow_hat)) + 
  geom_point() +
  geom_abline(color = "red") +
  facet_wrap(~site,
             scales = "free")


water_budget[, et_cm := (q_hat + meteo_hat) - Ds_cm - net_precip_cm - melt_cm]

ggplot(water_budget,
       aes(x = pet_cm,
           y = et_cm)) + 
  geom_point() +
  geom_abline(color = "red") +
  facet_wrap(~site,
             scales = "free")

et_form <- 
  bf(
    et_cm ~ Intercept + petslope * pet_cm,
    # Allow intercept, slope1, slope2, and change point to vary with site
    Intercept + petslope ~ 1 + (1|site),
    # Specify non-linear formula (ie `*` means multiply not crossed-effects and
    # exponents can be used directly)
    nl = TRUE
  )

# Priors
et_priors <- 
  prior(normal(0, 5), nlpar = "Intercept") +
  prior(normal(0, 2), nlpar = "petslope")

mod_et <- 
  water_budget[, .(mod = list(brm(et_form, 
                                  data = .SD, 
                                  prior = et_priors, 
                                  chains = 3))),
               keyby = .(site_status)]

options(mc.cores = 1)

et_coefs <- 
  water_budget[, setNames(as.list(fixef(mod_et[.BY[[2]], mod[[1]]])[, "Estimate"] + ranef(mod_et[.BY[[2]], mod[[1]]])$site[.BY[[1]], "Estimate", ]),
                          c("intercept", "m_pet")),
               keyby = .(site, site_status)]

et_coefs[, f_predict := pmap(list(intercept, m_pet), 
                             function(intercept, m_pet){
                               as.function(list(pet = NULL, 
                                                bquote(.(b) + .(m) * pet, 
                                                       where = list(b = intercept, 
                                                                    m = m_pet))))})]


water_budget[, et_hat := et_coefs[CJ(.BY[[1]], .BY[[2]]), f_predict[[1]]](pet_cm),
             by = .(site, site_status)]


ggplot(water_budget,
       aes(x = pet_cm,
           y = et_cm)) + 
  geom_point() +
  geom_line(aes(y = et_hat),
            color = "red") +
  facet_wrap(~site,
             scales = "free")


water_budget[, Ds_hat := (q_hat + meteo_hat) - et_hat + net_precip_cm + melt_cm]


ggplot(water_budget,
       aes(x = Ds_cm,
           y = Ds_hat)) + 
  geom_point() +
  geom_point(aes(y = (q_hat + meteo_hat) - et_hat + net_precip_cm + melt_cm),
             color = "blue",
             alpha = 0.25,
             shape = 20) +
  geom_abline(color = "red") +
  facet_wrap(~site,
             scales = "free")


dat <- 
  water_budget[net_flow_cm < 5 & site == "139" & site_status == "Treated" & !is.na(net_flow_cm + wl_initial_cm)]

full_form <- bf(
  # Section 1 of piece-wise regression
  Ds_cm ~ Intercept + precipslope * net_precip_cm + petslope * pet_cm + gddslope * gdd_10 + slope1 * wl_initial_cm * step(change - wl_initial_cm) +
    # Section 2 of piece-wise regression
    (slope1 * change + slope2 * (wl_initial_cm - change)^2) * step(wl_initial_cm - change),
  # Allow intercept, slope1, slope2, and change point to vary with site
  Intercept + slope1 + slope2 + change + petslope + gddslope + precipslope ~ 1,
  # Specify non-linear formula (ie `*` means multiply not crossed-effects and
  # exponents can be used directly)
  nl = TRUE
)

# Priors
full_prior <- 
  prior(normal(0, 5), nlpar = "Intercept") +
  prior(normal(0, 2), nlpar = "slope1") +
  prior(normal(0, 2), nlpar = "slope2") + 
  prior(normal(0, 2), nlpar = "petslope") +
  prior(normal(0, 2), nlpar = "gddslope") +
  prior(normal(0, 2), nlpar = "precipslope") +
  prior(normal(5, 3), nlpar = "change")  # Within observed range

full_fit <- 
  brm(full_form, 
      data = dat, 
      prior = full_prior, 
      chains = 3)

dat[, ds_hat := predict(full_fit, newdata = .SD)[, "Estimate"]]

ggplot(dat,
       aes(x = Ds_cm,
           y = ds_hat)) + 
  geom_point() +
  geom_abline(color = "red") +
  coord_cartesian(xlim = c(-1.5, 1.5),
                  ylim = c(-1.5, 1.5))

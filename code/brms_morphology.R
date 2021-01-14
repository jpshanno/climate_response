library(brms)
options(mc.cores = 3)

tar_load(water_budget)

# Filtering by esy < 5 barely changes the medians


water_budget[, wl_bin := cut_interval(wl_initial_cm,
                                      n = uniqueN(round(wl_initial_cm, 0))/3), 
             by = .(site)]

# Reduced data to median flow for each centimeter of water level to remove 
# noise
reduced_dat <- 
  water_budget[Ds_cm < 5, 
               .(Ds_cm = median_na(Ds_cm + pet_cm - net_precip_cm - melt_cm),
                 wl_initial_cm = median_na(wl_initial_cm)), 
               by = .(site, wl_bin)]

# Remove wl_bin from original data
water_budget[, wl_bin := NULL]

# Remove NAs for modeling & flows > 5 cm as these all seem to occur at high
# water levels and likely represent melt/rain on snow or missed precip 
# events. Removing data from before May removes many of these too, but then 
# the most extreme high water levels are missed

reduced_dat <- 
  subset(x = reduced_dat,
         subset = 
           !is.na(Ds_cm) & 
           !is.na(wl_initial_cm))

ggplot(water_budget,
       aes(x = wl_initial_cm,
           y = Ds_cm)) + 
  geom_point() +
  geom_point(data = reduced_dat,
             color = "red",
             shape = 20) +
  facet_wrap(~site,
             scales = "free")

# The model
morpho_form <- bf(
  # Section 1 of piece-wise regression
  Ds_cm ~ Intercept + slope1 * wl_initial_cm * step(change - wl_initial_cm) +
    # Section 2 of piece-wise regression
    (slope1 * change + slope2 * (wl_initial_cm - change)^2) * step(wl_initial_cm - change),
  # Allow intercept, slope1, slope2, and change point to vary with site
  Intercept + slope1 + slope2 + change ~ 1 + (1|site),
  # Specify non-linear formula (ie `*` means multiply not crossed-effects and
  # exponents can be used directly)
  nl = TRUE
)

# Priors
morpho_prior <- 
  prior(normal(0, 5), nlpar = "Intercept") +
  prior(normal(0, 2), nlpar = "slope1") +
  prior(normal(0, 2), nlpar = "slope2") +
  prior(normal(5, 3), nlpar = "change")

morpho_mod <- 
  brm(morpho_form, 
      data = reduced_dat, 
      prior = morpho_prior, 
      chains = 3)

morpho_coefs <- 
  water_budget[, setNames(as.list(fixef(morpho_mod)[, "Estimate"] + ranef(morpho_mod)$site[.BY[[1]], "Estimate", ]),
                          c("intercept", "m", "m2", "cp")),
               keyby = .(site)]

morpho_coefs[, f_predict := pmap(list(intercept, m, m2, cp), 
                                 function(intercept, m, m2, cp){
                                   as.function(list(water.level = NULL, 
                                                    bquote(.(b) + .(m1) * water.level * as.numeric(.(change) >= water.level) + .(m1) * .(change) + .(quadm) * (water.level - .(change))^2 * as.numeric(.(change) < water.level), 
                                                           where = list(b = intercept, m1 = m, quadm = m2, change = cp))))})]

plot(conditional_effects(morpho_mod,
                         prob = 0,
                         re_formula = NULL,
                         conditions = data.frame(site = unique(reduced_dat$site))),
     points = TRUE)[[1]] + 
  facet_wrap(~ site, scales = "free")

water_budget[, q_hat := morpho_coefs[.BY[[1]], f_predict[[1]]](wl_initial_cm), 
             by = .(site)]

ggplot(water_budget,
       aes(x = wl_initial_cm,
           y = Ds_cm)) + 
  geom_point() +
  geom_point(aes(y = q_hat),
             color = "red",
             shape = 20) +
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

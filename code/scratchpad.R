# Drop unused scenarios on input
swg_simulations_loca <- simplify_scenarios(swg_simulations_loca)[scenario != "scenario"]

future.climate <- simplify_scenarios(swg_simulations_loca)[scenario != "historical"]
obs.climate <- tar_read(swg_simulations_observed)

obs.climate[j = `:=`(gcm = "observed", scenario = "current")]

climate <- rbind(
  obs.climate[j = .(gcm, scenario, simulation_date, simulation_id, pet_cm, precip_cm, rain_cm, snow_cm, swe_cm, melt_cm, total_input_cm)],
  future.climate[j = .(gcm, scenario, simulation_date, simulation_id, pet_cm, precip_cm, rain_cm, snow_cm, swe_cm, melt_cm, total_input_cm)]
)

climate[
  j = `:=`(
    pet_cm = filled_rolling_mean(pet_cm),
    # daily_precip_cm = filled_rolling_mean(precip_cm),
    # monthly_precip_cm = filled_rolling_sum(precip_cm),
    daily_rain_cm = filled_rolling_mean(rain_cm),
    rain_days = filled_rolling_mean(rain_cm > 0)
  ),
  by = .(gcm, scenario, simulation_id)
]

drivers <- data.table::melt(
  data = climate,
  measure.vars = c("pet_cm", "daily_rain_cm", "rain_days", "melt_cm"),
  id.vars = c("gcm", "scenario", "simulation_date", "simulation_id")
)

ribbon_dat <- drivers[
  j = data.table::as.data.table(ggdist::mean_hdci(value, .width = 0.67, na.rm = TRUE)),
  by = .(gcm, scenario, simulation_date, variable)
]

ribbon_dat[
  j = climate := fcase(
    scenario == "current", "Current",
    scenario == "rcp45", "Less Sensitive",
    scenario == "rcp85", "More Sensitive"
  )
]

ggplot(ribbon_dat) +
  aes(x = simulation_date,
      y = abs(y),
      ymin = abs(ymin),
      ymax = abs(ymax),
      color = climate,
      fill = climate) +
  geom_line() +
  geom_ribbon(color = NA, alpha = 0.2) +
  facet_wrap(~variable, ncol = 1, scales = "free")


filled_rolling_sum <- function(x, na.rm = TRUE) {
  
  # Mirror first and last 15 rows to get rid of NA values from rolling mean
  x <- c(
    head(x, 15),
    x,
    tail(x, 15)
  )
  
  # Get 30 day rolling mean of precip variables
  x <- data.table::frollsum(x = x, n = 30, align = "center", na.rm = na.rm)
  
  # Drop mirrored values
  x <- head(x, -15)
  x <- tail(x, -15)
  
}

filled_rolling_mean <- function(x) {
  
  # Mirror first and last 15 rows to get rid of NA values from rolling mean
  x <- c(
    head(x, 15),
    x,
    tail(x, 15)
  )
  
  # Get 30 day rolling mean of precip variables
  x <- data.table::frollmean(x = x, n = 30, align = "center", na.rm = TRUE)
  
  # Drop mirrored values
  x <- head(x, -15)
  x <- tail(x, -15)
  
}

params <- tar_read(model_params)[site %in% tar_read(treatment_sites)]
params <- params[, .(params = list(as.data.table(head(params[[1]], -1)))), by = .(site, site_status, future.forest.change)]
# The control and treated parameters are not in the same order. Have to do 
# rbindlist rather than just convert the list to a datatable and let rbind
params <- params[, .(site, site_status, future.forest.change, rbindlist(params, use.names = TRUE))]
params[j = `:=`(
  MPET = fcoalesce(future.forest.change * MPET, MPET),
  future.forest.change = NULL
)]

params <- melt(params, id.vars = c("site", "site_status"))
ggplot(params[site %in% tar_read(treatment_sites) & variable %in% c("MPET", "MP")]) +
  aes(x = site,
      y = value,
      color = site_status) +
  geom_boxplot() +
  facet_wrap(~variable, scales = "free")

dcast(params[variable == "MPET"], site ~ site_status, value.var = "value")


rbindlist(list(tar_read(test_data_fits), tar_read(test_data_future_forest_fits)))[j = pet_hat := filled_rolling_mean(pet_hat), by = .(site, site_status)][j = .(pet_hat = mean(pet_hat, na.rm = TRUE)), by = .(site, site_status, doy = as.dowy(sample_date, 11))] %>% ggplot() + aes(x = doy, y = pet_hat, color = site_status) + geom_line() + facet_wrap(~site, scales = "free")
data <- 
  tar_read(training_data)


ggplot(data,
       aes(x = wl_initial_cm,
           y = Ds_cm)) +
  geom_point() +
  facet_wrap(~site,
             scales = "free")

data[, wl_bin := cut_interval(wl_initial_cm,
                              n = uniqueN(round(wl_initial_cm, 0))/1), 
     by = .(site)]

# Reduced data to median flow for each centimeter of water level to remove 
# noise
reduced_dat <- 
  data[Ds_cm < 5, 
       .(Ds_cm = median_na(Ds_cm),
         wl_initial_cm = median_na(wl_initial_cm)), 
       by = .(site, wl_bin)]

ggplot(data,
       aes(x = wl_initial_cm,
           y = Ds_cm)) +
  geom_point() +
  geom_line(data = reduced_dat, color = "red") +
  facet_wrap(~site,
             scales = "free")

# Remove wl_bin from original data
data[, wl_bin := NULL]

# Remove NAs for modeling & flows > 5 cm as these all seem to occur at high
# water levels and likely represent melt/rain on snow or missed precip 
# events. Removing data from before May removes many of these too, but then 
# the most extreme high water levels are missed

reduced_dat <- 
  subset(x = reduced_dat,
         subset = 
           !is.na(Ds_cm) & 
           !is.na(wl_initial_cm))

outflow_mods <- 
  reduced_dat[, .(mod = list(mcp(list(Ds_cm ~ 0,
                                      ~0 + wl_initial_cm + I(wl_initial_cm^2)),
                                 data = .SD))),
              keyby = .(site)]

outflow_mods <- 
  outflow_mods[, cbind(.SD,
                   map_dfr(mod,
                           ~dcast(setDT(fixef(.x)[, c("name", "mean")]), 
                                  . ~ name, 
                                  value.var = "mean")[, -c(".")]))]

# Drop model sigma estimate
outflow_mods[, sigma_1 := NULL]

setnames(outflow_mods,
         c("cp_1", "wl_initial_cm_2", "wl_initial_cm_2_E2"),
         c("changepoint", "linear_slope_2", "quadratic_slope_2"))

predict_functions <- 
  split(outflow_mods,
        by = "site") %>%
  map(~as.function(list(water.level = NULL,
                        substitute({
                          q <- 
                            ifelse(water.level <= changepoint,
                                   0,
                                   linear_slope_2 * (water.level - changepoint) + quadratic_slope_2 * (water.level - changepoint)^2)
                          q},
                          env = .x))))


outflow_mods[, f_predict := predict_functions]

data[, streamflow_cm := -1 * outflow_mods[.BY[[1]], f_predict[[1]]](wl_initial_cm),
     by = .(site)]


ggplot(data,
       aes(x = wl_initial_cm,
           y = Ds_cm)) +
  geom_point() +
  geom_line(aes(y = -streamflow_cm), 
            color = "red") +
  facet_wrap(~site,
             scales = "free")


data[, g_cm := Ds_cm - net_precip_cm + pet_cm + streamflow_cm]

ggplot(data,
       aes(x = water_availability_cm,
           y = g_cm)) +
  geom_point() +
  geom_smooth(method = lmrob,
              formula = y ~ x + I(x^2),
              method.args = list(setting = "KS2014"),
              se = FALSE) +
  facet_wrap(~site,
             scales = "free") +
  coord_cartesian(ylim = c(-2, 4))

new_flow_mods <-
  data[,
       .(mod = list(lmrob(net_flow_cm ~
                             water_availability_cm + I(water_availability_cm^2),
                           data = .SD,
                           setting = "KS2014"))),
       keyby = .(site, site_status)]




# options(future.fork.enable = TRUE)
# new_flow_mods <- 
#   data[, .(mod = list(mcp(list(net_flow_cm ~ water_availability_cm,
#                                       ~0 + water_availability_cm + I(water_availability_cm^2)),
#                           cores = 3,
#                           data = .SD))),
#               keyby = .(site)]
# 
# new_flow_mods <- 
#   new_flow_mods[, cbind(.SD,
#                        map_dfr(mod,
#                                ~dcast(setDT(fixef(.x)[, c("name", "mean")]), 
#                                       . ~ name, 
#                                       value.var = "mean")[, -c(".")]))]
# 
# # Drop model sigma estimate
# new_flow_mods[, sigma_1 := NULL]
# 
# setnames(new_flow_mods,
#          c("cp_1", "int_1", "water_availability_cm_1", "water_availability_cm_2", "water_availability_cm_2_E2"),
#          c("changepoint", "intercept", "linear_slope_1", "linear_slope_2", 
#            "quadratic_slope_2"))
# 
# predict_functions <- 
#   split(new_flow_mods,
#         by = "site") %>%
#   map(~as.function(list(water_availability_cm = NULL,
#                         substitute({
#                           q <- 
#                             ifelse(water_availability_cm <= changepoint,
#                                    intercept + linear_slope_1 * water_availability_cm,
#                                    intercept + linear_slope_1 * (water_availability_cm - changepoint) + linear_slope_2 * (water_availability_cm - changepoint) + quadratic_slope_2 * (water_availability_cm - changepoint)^2)
#                           q},
#                           env = .x))))
# 
# 
# new_flow_mods[, f_predict := predict_functions]


data[, pred_net_flow_cm := predict(new_flow_mods[CJ(.BY[[1]], .BY[[2]]), 
                                                 mod[[1]]],
                                   newdata = .SD),
     by = .(site, site_status)]

data[, drawdown_cm := pred_net_flow_cm - Ds_cm + net_precip_cm - streamflow_cm]

ggplot(data,
       aes(x = pet_cm,
           y = drawdown_cm)) +
  geom_point() +
  geom_smooth(method = lmrob,
              formula = y ~ 0 + x + I(x^0.5),
              method.args = list(setting = "KS2014"),
              se = FALSE) +
  facet_wrap(~site,
             scales = "free")

new_drawdown_mods <- 
  data[, 
       .(mod = list(lmrob(drawdown_cm ~ 
                             0 + pet_cm + I(pet_cm^0.5), 
                           data = .SD,
                           setting = "KS2014"))),
       keyby = .(site, site_status)]

data[, pred_drawdown_cm := predict(new_drawdown_mods[CJ(.BY[[1]], .BY[[2]]), 
                                                 mod[[1]]],
                                   newdata = .SD),
     by = .(site, site_status)]


data[, pred_ds := net_precip_cm + pred_net_flow_cm - pred_drawdown_cm - streamflow_cm]

ggplot(data = data,
       aes(x = Ds_cm,
           y = pred_ds)) +
  geom_point() +
  geom_abline(color = "red") +
  facet_wrap(~site,
             scales = "free") +
  coord_cartesian(ylim = c(-2, 2),
                  xlim = c(-2, 2))


f_interception <- 
  intercept.models[site.status, f_predict[[1]]]
f_esy <- 
  esy.models[site.id, f_predict[[1]]]
f_rise <- 
  precip.models[site.id, f_predict[[1]]]



if(any(is.na(weather$precip_cm))){
  message("Some missing precip values were filled in as 0.")
  
  weather[, precip_cm := nafill(precip_cm,
                                type = "const",
                                fill = 0)]
}


# Calculate Solar Radiation & PET -----------------------------------------

# bc() needs sample date but just converts it back to doy, so year doesn't
# matter
# calculate_solar_radiation uses station_name (just for grouping) and lat
weather[, `:=`(lat = lat,
               station_name = site.id,
               sample_date = as.Date(doy, origin = "2020-01-01"))]

# calculate_solar_radiation() calculate_mean_temp() and
# calculate_hargreaves_pet() were designed to work in a pipe so they return
# modified data as a side effect. This may change in the future
calculate_solar_radiation(weather, solrad.coefs)
calculate_mean_temp(weather)
calculate_hargreaves_pet(weather, lambda.MJ.kg = 2.45)

# Calculate Interception ---------------------------------------------------

weather[, interception_cm := f_interception(doy)]
weather[, net_precip_cm := pmax(0, precip_cm - interception_cm)]

wl <- Ds <- esy <- net_flow <- water_availability <- drawdown <- streamflow <- rain_rise <- 
  numeric(nrow(weather))

wl[1] <- 
  initial.wl

water_availability[1] <- 
  initial.water.availability

for(i in 1:(nrow(weather)-1)){
  
  esy[i] <- 
    f_esy(wl[i])
  
  streamflow[i] <- 
    -1 * outflow_mods[site.id, f_predict[[1]]](wl[i])
  
  net_flow[i] <- 
    predict(new_flow_mods[CJ(site.id, site.status),
                          mod[[1]]],
            newdata = data.frame(water_availability_cm = water_availability[i]))
  
  drawdown[i] <-
    predict(new_drawdown_mods[CJ(site.id, site.status),
                              mod[[1]]],
            newdata = data.frame(pet_cm = weather[i, pet_cm]))
  
  rain_rise[i] <-
    f_rise(weather[i, net_precip_cm])
  
  # if(rain_rise[i] > 0){
  #   drawdown[i] <-
  #     0.5 * drawdown[i]*0.5
  # }
  
  # HAVE TO EXCLUDE MELT for now
  Ds[i] <- 
    rain_rise[i] + net_flow[i] - drawdown[i] - streamflow[i]
  
  water_availability[i+1] <- 
    water_availability[i] + (weather[i, precip_cm] - weather[i, pet_cm])
  
  wl[i+1] <- wl[i] + Ds[i] / esy[i]
}

pred_budget <- 
  data.table(doy = weather$doy, 
             g_precip = weather$precip_cm,
             pred_i = weather$interception_cm,
             pred_wl = wl, 
             pred_Ds = Ds, 
             pred_esy = esy, 
             pred_netflow = net_flow,
             drawdown, 
             streamflow = streamflow,
             pred_water_availability = water_availability)

test_dat[pred_budget,
         `:=`(pred_wl = i.pred_wl,
              g_precip = i.g_precip,
              pred_water_availability = i.pred_water_availability),
         on = "doy"]

test_dat[, `:=`(cum_precip_cm = cumsum(nafill(best_precip_cm, "const", 0)),
                cum_g_precip = cumsum(g_precip))]

test_dat[, wa_diff_cm := cum_precip_cm - cum_g_precip]

ggplot(test_dat,
       aes(x = doy)) +
  geom_line(aes(y = wl_initial_cm)) +
  geom_line(aes(y = pred_wl),
            color = "blue",
            linetype = "dotted") +
  geom_line(aes(y = (pred_wl + wa_diff_cm)),
            color = "blue",
            linetype = "dashed") +
  theme_bw()

ggplot(test_dat,
       aes(x = doy)) +
  geom_line(aes(y = Ds_cm)) +
  geom_line(data = pred_budget,
            aes(y = pred_Ds),
            color = "red",
            linetype = "dashed") +
  theme_bw()

ggplot(test_dat,
       aes(x = doy)) +
  geom_line(aes(y = net_flow_cm)) +
  geom_line(data = pred_budget,
            aes(y = pred_netflow - streamflow),
            color = "red",
            linetype = "dashed") +
  theme_bw()

ggplot(test_dat,
       aes(x = doy)) +
  geom_line(aes(y = pet_cm)) +
  geom_line(data = pred_budget[drawdown != 0],
            aes(y = drawdown),
            color = "red",
            linetype = "dashed") +
  theme_bw()

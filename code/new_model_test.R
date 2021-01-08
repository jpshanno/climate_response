# Model G
# Model ET
# Model P
# Calculate pred Ds
# Model Q

# MAY BE BETTER TO TRY DOING BELOW APPROACH AND THEN MODEL STREAMFLOW AS RESIDUAL
# OF UNEXPLAINED DS
# Good performance with net_flow_cm ~ wl_initial_cm*(water_availability_m + I(water_availability_m^2))
# with no separate streamflow model

# Load Data ---------------------------------------------------------------

data <- 
  tar_read(training_data)


# Model G ----------------------------------------------------------

ggplot(data,
       aes(x = water_availability_cm,
           y = net_flow_cm)) +
  geom_point() +
  geom_smooth(method = lmrob,
              formula = y ~ x + I(x^2),
              method.args = list(setting = "KS2014"),
              se = FALSE) +
  facet_wrap(~site,
             scales = "free") +
  coord_cartesian(ylim = c(-2, 4))

mod_g <-
  data[,
       .(mod = list(lmrob(net_flow_cm ~
                            water_availability_cm + I(water_availability_cm^2),
                          data = .SD,
                          setting = "KS2014"))),
       keyby = .(site, site_status)]

mod_g[, c("intercept", "slope", "quad_slope") := map_dfr(mod, coef)]

func_g <-
  split(mod_g,
        by = c("site", "site_status")) %>%
  map(~as.function(list(water.availability = NULL,
                        substitute({
                         intercept + slope * water.availability + quad_slope * water.availability^2},
                          env = .x))))


mod_g[, f_predict := func_g]
  

# data[, wa_bin := cut_interval(water_availability_cm,
#                               n = uniqueN(round(water_availability_cm, 0))/5),
#      by = .(site)]
# 
# # Reduced data to median flow for each centimeter of water level to remove
# # noise
# reduced_dat_wa <-
#   data[abs(net_flow_cm) < 5,
#        .(net_flow_cm = median_na(net_flow_cm),
#          water_availability_cm = median_na(water_availability_cm)),
#        by = .(site, wa_bin)]
# 
# ggplot(data[abs(net_flow_cm) < 5],
#        aes(x = water_availability_cm,
#            y = net_flow_cm)) +
#   geom_point() +
#   geom_point(data = reduced_dat_wa, 
#              shape = 19,
#              color = "red") +
#   facet_wrap(~site,
#              scales = "free")
# 
# mod_g <-
#   reduced_dat_wa[, 
#                  .(mod = list(mcp(list(net_flow_cm ~ water_availability_cm,
#                                        ~0 + water_availability_cm),
#                                   data = .SD))),
#                  keyby = .(site)]
# 
# mod_g <-
#   mod_g[, cbind(.SD,
#                 map_dfr(mod,
#                         ~dcast(setDT(fixef(.x)[, c("name", "mean")]),
#                                . ~ name,
#                                value.var = "mean")[, -c(".")]))]
# 
# # Drop model sigma estimate
# mod_g[, sigma_1 := NULL]
# 
# setnames(mod_g,
#          c("cp_1", "int_1", "water_availability_cm_1", "water_availability_cm_2"),
#          c("changepoint", "intercept", "linear_slope_1", "linear_slope_2"))
# 
# func_g <-
#   split(mod_g,
#         by = "site") %>%
#   map(~as.function(list(water_availability_cm = NULL,
#                         substitute({
#                           q <-
#                             ifelse(water_availability_cm <= changepoint,
#                                    intercept + linear_slope_1 * water_availability_cm,
#                                    intercept + linear_slope_1 * changepoint + linear_slope_2 * (water_availability_cm - changepoint))
#                           q},
#                           env = .x))))
# 
# 
# mod_g[, f_predict := func_g]

data[, g_hat := mod_g[CJ(.BY[[1]]),
                      f_predict[[1]]](water_availability_cm),
     by = .(site)]

# Model ET ----------------------------------------------------------------

data[, et_cm := g_hat - Ds_cm + net_precip_cm]

ggplot(data,
       aes(x = pet_cm,
           y = et_cm)) +
  geom_point() +
  geom_smooth(method = lmrob,
              formula = y ~ x + I(x^0.5),
              method.args = list(setting = "KS2014"),
              se = FALSE) +
  facet_wrap(~site,
             scales = "free")

plan(multisession,
     workers = 2)

mod_et <- 
  split(data,
        by = c("site_status")) %>% 
  future_map(~.x[,.(site_status = first(site_status),
                    mod = list(rlmerRcpp(et_cm ~ 
                                           pet_cm + I(pet_cm^0.5) + (1 + pet_cm + I(pet_cm^0.5) || site), 
                                         data = .SD)))]) %>% 
  rbindlist() %>% 
  setkey("site_status")

# data[, 
#      .(mod = list(rlmerRcpp(et_cm ~ 
#                               pet_cm + I(pet_cm^0.5) + (1 + pet_cm + I(pet_cm^0.5) || site), 
#                             data = .SD))),
#      keyby = .(site_status)]


# Fixed effects predictions:
mod_et <-
  mod_et[, c(.SD, map_dfr(mod, lme4::fixef))]

# Clean up coefficient names
setnames(mod_et,
         c("site_status", "mod", "intercept", "slope1", "slope2"))

# Created prediction functions
func_et_fixed <- 
  split(mod_et,
        by = "site_status") %>%
  map(~as.function(list(pet = NULL,
                        substitute({
                          
                          intercept + 
                            slope1 * pet + 
                            slope2 * sqrt(pet)
                          },
                          env = .x))))

# Add prediction functions to data.table
mod_et[, f_predict_fixed := func_et_fixed]

# Calculate predicted fixed-effect values with predict()
true_predict <- 
  data[, predict(mod_et[.BY[[1]], mod[[1]]],
                 re.form = NA,
                 newdata = .SD),
       by = .(site_status)]$V1

# Calculate predicted fixed-effect values with function
test_predict <- 
  data[, mod_et[.BY[[1]], f_predict_fixed[[1]]](pet_cm),
       by = .(site_status)]$V1

# Check that manual and predict() predictions match
if(isFALSE(all.equal(true_predict,test_predict))){
  stop("The fixed effects prediction function does not match the output from predict()")
}

# Random effects predictions:
# Expand mod_et to all site/status combinations
mod_et <- 
  mod_et[unique(data[, .(site_status, site)])]

# Extract the random effects from the site_status models
mod_et[, random_effects := map2(mod, site, ~lme4::ranef(.x)$site[.y, ])]

# Adjust the fixed-effects coefficients for random effects
mod_et[, c("intercept", "slope1", "slope2") := 
         map2(.SD[, .(intercept, slope1, slope2)], 
              random_effects[[1]], 
              `+`), 
       by = .(site)]

# Set key to ensure ordering and look-up for prediction
setkey(mod_et, "site", "site_status")

# Create prediction functions
func_et_random <- 
  split(mod_et,
        by = c("site", "site_status")) %>%
  map(~as.function(list(pet = NULL,
                        substitute({
                          
                          intercept + 
                            slope1 * pet + 
                            slope2 * sqrt(pet)
                        },
                        env = .x))))

# Add prediction functions to data.table
mod_et[, f_predict_random := func_et_random]

# Calculate predicted random-effect values with predict()
true_rand_predict <- 
  data[, predict(mod_et[CJ(.BY[[1]], .BY[[2]]), mod[[1]]],
                 newdata = cbind(.SD, site = .BY[[1]])),
       by = .(site, site_status)]$V1

# Calculate predicted random-effect values with function
test_rand_predict <- 
  data[, mod_et[CJ(.BY[[1]], .BY[[2]]), f_predict_random[[1]]](pet_cm),
       by = .(site, site_status)]$V1

# Check that manual and predict() predictions match for random predictions
all.equal(true_rand_predict,test_rand_predict)

data[, et_hat := mod_et[CJ(.BY[[1]], .BY[[2]]),
                        f_predict_random[[1]]](pet_cm),
     by = .(site, site_status)]


# Model P -----------------------------------------------------------------

data[, p_cm := et_hat - g_hat + Ds_cm]

ggplot(data[net_precip_cm > 0],
       aes(x = net_precip_cm,
           y = p_cm)) + 
  geom_point() +
  geom_smooth(method = lmrob,
              formula = y ~ 0 + x + I(x^0.5),
              se = FALSE,
              method.args = list(setting = "KS2014")) +
  facet_wrap(~site,
             scales = "free")

mod_p <- 
  data[net_precip_cm > 0, 
       .(mod = list(lmrob(p_cm ~ 0 + net_precip_cm + I(net_precip_cm^0.5),
                          data = .SD,
                          setting = "KS2014"))),
       keyby = .(site)]

data[, p_hat := predict(mod_p[.BY[[1]], mod[[1]]],
                        newdata = .SD),
     by = .(site)]


# Model Q -----------------------------------------------------------------

data[, q_cm := p_hat - et_hat + g_hat - Ds_cm]

ggplot(data,
       aes(x = wl_initial_cm,
           y = q_cm)) + 
  geom_point() +
  facet_wrap(~site,
             scales = "free")


data[, wl_bin := cut_interval(wl_initial_cm,
                              n = uniqueN(round(wl_initial_cm, 0))/2), 
     by = .(site)]

# Reduced data to median flow for each centimeter of water level to remove 
# noise
q_dat <- 
  data[q_cm > -2, 
       .(q_cm = quantile(q_cm, probs = 0.5, na.rm = TRUE),
         wl_initial_cm = median_na(wl_initial_cm)), 
       by = .(site, wl_bin)]

ggplot(data,
       aes(x = wl_initial_cm,
           y = q_cm)) +
  geom_point() +
  geom_point(data = q_dat, 
             color = "red") +
  facet_wrap(~site,
             scales = "free")


mod_q <-
  q_dat[, 
        .(mod = list(mcp(list(q_cm ~ wl_initial_cm,
                              ~0 + wl_initial_cm + I(wl_initial_cm^2)),
                         data = .SD))),
        keyby = .(site)]

mod_q <-
  mod_q[, cbind(.SD,
                map_dfr(mod,
                        ~dcast(setDT(fixef(.x)[, c("name", "mean")]),
                               . ~ name,
                               value.var = "mean")[, -c(".")]))]

# Drop model sigma estimate
mod_q[, sigma_1 := NULL]

setnames(mod_q,
         c("cp_1", "int_1", "wl_initial_cm_1", "wl_initial_cm_2", "wl_initial_cm_2_E2"),
         c("changepoint", "intercept", "linear_slope_1", "linear_slope_2", "quad_slope_2"))

func_q <-
  split(mod_q,
        by = "site") %>%
  map(~as.function(list(water.level = NULL,
                        substitute({
                          q <-
                            ifelse(water.level <= changepoint,
                                   intercept + linear_slope_1 * water.level,
                                   intercept + linear_slope_1 * changepoint + linear_slope_2 * (water.level - changepoint) +  + quad_slope_2 * (water.level - changepoint)^2)
                          q},
                          env = .x))))


mod_q[, f_predict := func_q]

data[, q_hat := mod_q[CJ(.BY[[1]]),
                      f_predict[[1]]](wl_initial_cm),
     by = .(site)]


# Calculate Ds_hat --------------------------------------------------------

data[, Ds_hat := p_hat - et_hat - q_hat + g_hat]

ggplot(data,
       aes(x = Ds_cm,
           y = Ds_hat,
           color = site_status)) +
  geom_point() + 
  geom_abline() +
  facet_wrap(~site,
             scales = "free")

ggplot(data,
       aes(x = Ds_cm,
           y = Ds_hat,
           color = site_status)) +
  geom_point() + 
  geom_abline() +
  facet_wrap(~site) +
  lims(x = c(-1.5, 1.5),
       y = c(-1.5, 1.5))

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
  reduced_dat[, .(mod = list(mcp(list(Ds_cm ~ wl_initial_cm,
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
         c("cp_1", "int_1", "wl_initial_cm_1", "wl_initial_cm_2", "wl_initial_cm_2_E2"),
         c("changepoint", "intercept", "linear_slope_1", "linear_slope_2",
           "quadratic_slope_2"))


predict_functions <- 
  split(outflow_mods,
        by = "site") %>%
  map(~as.function(list(water.level = NULL,
                        substitute({
                          q <- 
                            ifelse(water.level <= changepoint,
                                   intercept + linear_slope_1 * water.level,
                                   intercept + linear_slope_1 * (water.level - changepoint) + linear_slope_2 * (water.level - changepoint) + quadratic_slope_2 * (water.level - changepoint)^2)
                          
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



ggplot(data,
       aes(x = water_availability_cm,
           y = g_cm)) +
  geom_point() +
  geom_line(aes(y = pred_net_flow_cm),
            color = "red") +
  facet_wrap(~site,
             scales = "free")





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



# Simulation --------------------------------------------------------------

wl <- Ds <- esy <- g <- water_availability <- et <- q <- p <- 
  numeric(nrow(weather))

wl[1] <- 
  initial.wl

water_availability[1] <- 
  initial.water.availability

for(i in 1:(nrow(weather)-1)){
  
  esy[i] <- 
    f_esy(wl[i])
  
  q[i] <- 
    mod_q[site.id, f_predict[[1]]](wl[i])
  
  # if(wl[i] < mod_q[site.id, changepoint]){
  #   streamflow[i] <- 0
  # }
  
  # net_flow[i] <- 
  #   predict(new_flow_mods[CJ(site.id, site.status),
  #                         mod[[1]]],
  #           newdata = data.frame(water_availability_cm = water_availability[i]))
  
  g[i] <- 
    mod_g[CJ(site.id, site.status),
          f_predict[[1]]](water_availability[i])
  
  et[i] <-
    mod_et[CJ(site.id, site.status), f_predict_random[[1]]](weather[i, pet_cm])
  
  p[i] <-
    predict(mod_p[site.id, mod[[1]]], 
            newdata = data.frame(net_precip_cm = weather[i, net_precip_cm]))
  
  # if(rain_rise[i] > 0){
  #   drawdown[i] <-
  #     0.5 * drawdown[i]*0.5
  # }
  
  # HAVE TO EXCLUDE MELT for now
  Ds[i] <- 
    p[i] + g[i] - et[i]
  
  water_availability[i+1] <- 
    water_availability[i] + (weather[i, precip_cm] - weather[i, pet_cm])
  
  wl[i+1] <- wl[i] + Ds[i]
}

pred_budget <- 
  data.table(doy = weather$doy, 
             g_precip = weather$precip_cm,
             pred_i = weather$interception_cm,
             pred_wl = wl, 
             pred_Ds = Ds, 
             pred_esy = esy, 
             pred_g = g,
             pred_et = et, 
             pred_q = q,
             pred_p = p,
             pred_water_availability = water_availability)

test_dat[pred_budget,
         `:=`(pred_wl = i.pred_wl,
              pred_Ds = i.pred_Ds,
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
  labs(caption = paste(site.id, site.status)) +
  theme_bw()

ggplot(test_dat,
       aes(x = doy)) +
  geom_col(aes(y = Ds_cm),
           alpha = 0.5) +
  geom_col(data = pred_budget,
            aes(y = pred_Ds),
            fill = "red",
           alpha = 0.5) +
  theme_bw()

ggplot(test_dat,
       aes(x = doy)) +
  geom_col(aes(y = net_flow_cm),
               alpha = 0.5) +
  geom_col(data = pred_budget,
           aes(y = pred_g),
           alpha = 0.5,
           fill = "red") +
  theme_bw()

ggplot(test_dat,
       aes(x = doy)) +
  geom_line(aes(y = pet_cm)) +
  geom_line(data = pred_budget[pred_et != 0],
            aes(y = pred_et),
            color = "red",
            linetype = "dashed") +
  theme_bw()


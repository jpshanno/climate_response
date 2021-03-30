source("code/load_project.R")
library(magrittr, include.only = "%T>%")
tar_load(solrad_coefs)


present_simulations <- 
  fread("output/tabular/weather_sequences_gfdl-cm3_1980-2009.csv") %>% 
  .[, precip := 0.1 * precip] %>% 
  setnames(c("tmin", "tmax", "precip"),
           c("tmin_c", "tmax_c", "precip_cm")) %>% 
  set(j = "station_name",
      value = "PIEM4") %>% 
  calculate_mean_temp()

present_simulations[, solrad_MJ_m2 := calculate_solar_radiation(.SD, 
                                                                coefs = solrad_coefs[station_name == "PIEM4"],
                                                                stop.on.error = FALSE,
                                                                return.vector = TRUE),
                    by = .(seq_ID)]

present_simulations[, pet_cm := -pmax(0, 0.1 * 0.0135 * solrad_MJ_m2 / 2.45 * (tmean_c + 17.8))]

present_sims <- 
  present_simulations %>%
  set(j = "seq_set", value = present_simulations$seq_ID %/% 30) %>% 
  set(j = "water_year", value = as.water_year(present_simulations$sample_date, 11)) %>% 
  .[!(water_year %in% c(1980, 2010))] %>% 
  split(by = "seq_set") %>% 
  map(~possibly(calculate_snowmelt, otherwise = NULL)(.x)) %>% 
  rbindlist() %>% 
  .[water_year != 2009 & water_year != 1981] %>% 
  .[format(sample_date, "%m%d") != "0229"]

future_simulations <- 
  fread("output/tabular/weather_sequences_gfdl-cm3_2070-2099.csv") %>% 
  .[, precip := 0.1 * precip] %>% 
  setnames(c("tmin", "tmax", "precip"),
           c("tmin_c", "tmax_c", "precip_cm")) %>% 
  set(j = "station_name",
      value = "PIEM4") %>% 
  calculate_mean_temp()

future_simulations[, solrad_MJ_m2 := calculate_solar_radiation(.SD, 
                                                                coefs = solrad_coefs[station_name == "PIEM4"],
                                                                stop.on.error = FALSE,
                                                                return.vector = TRUE),
                    by = .(seq_ID)]

future_simulations[, pet_cm := -pmax(0, 0.1 * 0.0135 * solrad_MJ_m2 / 2.45 * (tmean_c + 17.8))]

future_sims <- 
  future_simulations %>%
  set(j = "seq_set", value = future_simulations$seq_ID %/% 30) %>% 
  set(j = "water_year", value = as.water_year(future_simulations$sample_date, 11)) %>% 
  .[!(water_year %in% c(2070, 2100))] %>% 
  split(by = "seq_set") %>% 
  map(~possibly(calculate_snowmelt, otherwise = NULL)(.x)) %>% 
  rbindlist() %>% 
  .[water_year != 2099 & water_year != 2071] %>% 
  .[format(sample_date, "%m%d") != "0229"]

site_info <- 
  fread(text = 
          "site,lon,lat,study,treatment,treatment_date,morphology
             006,-89.19036,46.32392,planting,Girdle,2013-10-31,flowthrough
             009,-89.71207,46.41624,eco,Ash Cut,2013-10-31,threshold
             053,-89.53043,46.61028,pws,Ash Cut,2014-10-31,flowthrough
             077,-88.86684,46.71778,eco,Ash Cut,2013-10-31,threshold
             111,-89.71369,46.57636,planting,Control,2013-10-31,threshold
             113,-89.80990,46.35090,pws,Control,2014-10-31,flowthrough
             119,-89.59433,46.28794,eco,Girdle,2013-10-31,closed
             135,-88.87838,46.75382,eco,Control,2013-10-31,flowthrough
             139,-88.95577,46.59619,planting,Ash Cut,2013-10-31,threshold
             140,-89.61962,46.43279,eco,Girdle,2013-10-31,flowthrough
             151,-88.89502,46.67307,eco,Girdle,2013-10-31,threshold
             152,-89.71468,46.41277,eco,Control,2013-10-31,flowthrough
             156,-89.58431,46.28327,eco,Ash Cut,2013-10-31,closed
             157,-89.58479,46.28173,eco,Control,2013-10-31,threshold",
        select = c(site = "character",
                   lon = "numeric",
                   lat = "numeric",
                   study = "character",
                   treatment = "character",
                   treatment_date = "Date",
                   morphology = "character"),
        key = "site")


full_dat <- 
  water_budget[CJ(site = unique(water_budget[, site]), 
                  sample_date = kenton[between(water_year, 2012, 2019), sample_date])
  ][kenton[between(water_year, 2012, 2019)],
    .(sample_date, 
      site, 
      wl_initial_cm, 
      idw_precip_cm,
      water_year = i.water_year,
      pet_cm = -i.pet_cm,
      rain_cm = i.rain_cm,
      melt_cm = i.melt_cm,
      total_input_cm = i.total_input_cm,
      tmax_c = i.tmax_c,
      tmin_c = i.tmin_c,
      tmean_c = i.tmean_c),
    on = "sample_date"]

setkey(full_dat, site, sample_date)

# Set Site status
full_dat[site_info,
         site_status := fifelse(sample_date > i.treatment_date & i.treatment != "Control",
                                "Treated",
                                "Control"),
         on = "site"]

# Drop Leap year
full_dat <- 
  full_dat[format(sample_date, "%m%d") != "0229"]

# Drop anomalous melt
full_dat[between(melt_cm, 0, 0.1) & month(sample_date) %in% 6:9, 
         `:=`(total_input_cm = total_input_cm - melt_cm,
              melt_cm = 0)]

# Drop anomalous pet
full_dat[tmax_c <= 0,
         pet_cm := 0]

# Calculate YTD after dropping above points
full_dat[, 
         `:=`(ytd_pet_cm = cumsum(pet_cm),
              ytd_p_cm = cumsum(rain_cm),
              ytd_m_cm = cumsum(melt_cm),
              ytd_input_cm = cumsum(rain_cm + melt_cm)),
         by = .(water_year)]

full_dat[, 
         `:=`(ytd_water_balance_cm = ytd_pet_cm + ytd_p_cm + ytd_m_cm,
              water_balance_cm = cumsum(pet_cm + rain_cm + melt_cm))]

full_dat[,
         max_wl_cm := possibly(function(x){density(x, na.rm = TRUE)$x[which.max(density(x, na.rm = TRUE)$y)]}, 
                               NA_real_)(wl_initial_cm),
         by = .(site)]

mod_dat <- 
  full_dat[, .(max_wl_cm = first(max_wl_cm),
               training_dat = list(.SD[water_year == ifelse(.BY[[2]] == "Control", 2012, 2015)]),
               valid_dat = list(.SD[water_year != ifelse(.BY[[2]] == "Control", 2012, 2015)])),
           by = .(site, site_status)
  ][map_dbl(training_dat, ~sum(!is.na(.x$wl_initial_cm))) > 0]

mod_dat[, params := map2(training_dat, max_wl_cm,
                         ~possibly(quad_opt, 
                                   otherwise = NULL)(wobs = .x$wl_initial_cm, 
                                                     met = .x[, .(pet_cm, rain_cm, melt_cm)],
                                                     max.wl = .y,
                                                     control = list(fnscale = -1)))]

sim_dat <- 
  mod_dat[, .(site, site_status, max_wl_cm, params)]

sim_dat[, `:=`(present_dat = map2(params, max_wl_cm,
                                  ~quad_curve_wl(PET = present_sims$pet_cm,
                                                 P = present_sims$rain_cm, 
                                                 M = present_sims$melt_cm, 
                                                 max.wl = first(.y),
                                                 b0 = .x$par[["b0"]],
                                                 b1 = .x$par[["b1"]],
                                                 b2 = .x$par[["b2"]])),
               future_dat = map2(params, max_wl_cm,
                                 ~quad_curve_wl(PET = future_sims$pet_cm,
                                                P = future_sims$rain_cm, 
                                                M = future_sims$melt_cm, 
                                                max.wl = first(.y),
                                                b0 = .x$par[["b0"]],
                                                b1 = .x$par[["b1"]],
                                                b2 = .x$par[["b2"]])))]

summary(sim_dat$future_dat[[17]][, wl_hat])
plot(sim_dat$future_dat[[17]][1:(365*10), wl_hat], type = 'l'); for(i in 365*(1:10)){abline(v = i)}
plot(sim_dat$future_dat[[17]][1:(365*10), wa_hat], type = 'l'); for(i in 365*(1:10)){abline(v = i)}

sim_dat[, present_summary := 
          map(present_dat,
              ~.x[, 
                  .(wl_hat_mode = rethinking::chainmode(wl_hat),
                    wl_hat_mean = mean(wl_hat),
                    wl_hat_med = median(wl_hat),
                    wl_hat_10 = quantile(wl_hat, 0.1),
                    wl_hat_025 = quantile(wl_hat, 0.025),
                    wl_hat_975 = quantile(wl_hat, 0.975),
                    wl_hpdi_lower = rethinking::HPDI(wl_hat, 0.66)[[1]],
                    wl_hpdi_upper = rethinking::HPDI(wl_hat, 0.66)[[2]]), 
                  by = .(doy = rep(1:365, times = nrow(.x)/365))])]

sim_dat[, future_summary := 
          map(future_dat,
              ~.x[, 
                  .(wl_hat_mode = rethinking::chainmode(wl_hat),
                    wl_hat_mean = mean(wl_hat),
                    wl_hat_med = median(wl_hat),
                    wl_hat_10 = quantile(wl_hat, 0.1),
                    wl_hat_025 = quantile(wl_hat, 0.025),
                    wl_hat_975 = quantile(wl_hat, 0.975),
                    wl_hpdi_lower = rethinking::HPDI(wl_hat, 0.66)[[1]],
                    wl_hpdi_upper = rethinking::HPDI(wl_hat, 0.66)[[2]]), 
                  by = .(doy = rep(1:365, times = nrow(.x)/365))])]

ggplot(sim_dat$present_summary[[17]],
       aes(x = doy,
           y = wl_hat_10)) +
  geom_ribbon(aes(ymin = wl_hpdi_lower,
                  ymax = wl_hpdi_upper),
              alpha = 0.6) +
  geom_line() +
  geom_line(color = "red",
            aes(y = wl_hat_cm),
            data = water_budget[site == "152", .(wl_hat_cm = mean_na(wl_initial_cm)), by = .(doy = dowy)])


simulations_wb <- 
  present_sims[, 
               .(dowy, 
                 water_balance_cm = cumsum(melt_cm + rain_cm + pet_cm) - max(cumsum(melt_cm + rain_cm + pet_cm))), 
               by = .(seq_ID = rep(1:(nrow(present_sims)/365), each = 365))
  ][, 
    .(wb_hat_cm = mean(water_balance_cm),
      wb_hat_025 = quantile(water_balance_cm, 0.025),
      wb_hat_975 = quantile(water_balance_cm, 0.975),
      wb_hpdi_lower = rethinking::HPDI(water_balance_cm, 0.66)[[1]],
      wb_hpdi_upper = rethinking::HPDI(water_balance_cm, 0.66)[[2]]), 
    by = .(dowy)]

observed_wb <- 
  kenton[water_year %in% 2012:2019, 
         .(wb_hat_cm = mean(ytd_water_availability_cm),
           wb_hat_025 = quantile(ytd_water_availability_cm, 0.025),
           wb_hat_975 = quantile(ytd_water_availability_cm, 0.975),
           wb_hpdi_lower = rethinking::HPDI(ytd_water_availability_cm, 0.66)[[1]],
           wb_hpdi_upper = rethinking::HPDI(ytd_water_availability_cm, 0.66)[[2]]), 
         by = .(dowy)]

ggplot(simulations_wb,
       aes(x = dowy,
           y = wb_hat_cm)) +
  geom_ribbon(aes(ymin = wb_hat_025,
                  ymax = wb_hat_975),
              alpha = 0.6) +
  geom_line() +
  geom_line(color = "red",
            data = observed_wb)



# Compare Treatment Effect ------------------------------------------------

present_dat <- 
  sim_dat[site %in% c("009", "077", "119", "140", "151", "139", "156"), 
          present_summary[[1]], 
          by = .(site, site_status)]

future_dat <- 
  sim_dat[site %in% c("009", "077", "119", "140", "151", "139", "156"), 
          future_summary[[1]], 
          by = .(site, site_status)]

simulations <- 
  rbindlist(list(Present = present_dat,
                 Future = future_dat),
            idcol = "period")

ggplot(present_dat,
       aes(x = doy, 
           y = wl_hat_mode, 
           color = site_status,
           fill = site_status,
           linetype = site_status)) +
  geom_line() +
  geom_ribbon(aes(ymin = wl_hpdi_lower,
                  ymax = wl_hpdi_upper),
              alpha = 0.6,
              color = NA) +
  facet_wrap(~site,
             scales = "free")


ggplot(future_dat,
       aes(x = doy, 
           y = wl_hat_mode, 
           color = site_status,
           fill = site_status,
           linetype = site_status)) +
  geom_line() +
  geom_ribbon(aes(ymin = wl_hpdi_lower,
                  ymax = wl_hpdi_upper),
              alpha = 0.6,
              color = NA) +
  facet_wrap(~site,
             scales = "free")


{ggplot(simulations[site_status == "Control"],
       aes(x = doy,
           y = wl_hat_med,
           color = period,
           fill = period,
           linetype = period)) +
  # geom_ribbon(aes(ymin = wl_hpdi_lower,
  #                 ymax = wl_hpdi_upper),
  #             alpha = 0.6,
  #             color = NA) +
  geom_line() +
  facet_wrap(~site) +
  colorblindr::scale_color_OkabeIto() +
  labs(x = "Day of Water Year",
       y = "Water Level (cm)",
       title = "Control Conditions") +
  ggthemes::theme_few(base_size = 16) +
  guides(color = guide_legend(title = NULL),
         linetype = guide_legend(title = NULL)) +
  theme(legend.position = c(0.65, 0.15),
        legend.text = element_text(size = rel(1.5)))} %T>%
  ggsave(filename = "tmp/climate_figures/control_conditions.png",
         device = "png",
         width = 11.75,
         height = 7.25,
         units = "in")

{ggplot(simulations[site_status == "Treated"],
       aes(x = doy,
           y = wl_hat_med,
           color = period,
           fill = period,
           linetype = period)) +
  # geom_ribbon(aes(ymin = wl_hpdi_lower,
  #                 ymax = wl_hpdi_upper),
  #             alpha = 0.6,
  #             color = NA) +
  geom_line() +
  facet_wrap(~site) +
  colorblindr::scale_color_OkabeIto() +
  labs(x = "Day of Water Year",
       y = "Water Level (cm)",
       title = "Treated Conditions") +
  ggthemes::theme_few(base_size = 16) +
  guides(color = guide_legend(title = NULL),
         linetype = guide_legend(title = NULL)) +
  theme(legend.position = c(0.65, 0.15),
        legend.text = element_text(size = rel(1.5)))} %T>%
  ggsave(filename = "tmp/climate_figures/treated_conditions.png",
         device = "png",
         width = 11.75,
         height = 7.25,
         units = "in")


{ggplot(simulations[(site_status == "Control" & period == "Present") | (site_status == "Treated" & period == "Future"),
                   .(scenario = factor(paste(site_status, period, sep = " - ")),
                     wl_hat_med,
                     doy,
                     site)],
       aes(x = doy,
           y = wl_hat_med,
           color = scenario,
           fill = scenario,
           linetype = scenario)) +
  geom_line() +
  facet_wrap(~site) +
  colorblindr::scale_color_OkabeIto() +
  labs(x = "Day of Water Year",
       y = "Water Level (cm)",
       title = "Compare EAB & Climate Impacts") +
  ggthemes::theme_few(base_size = 16) +
  guides(color = guide_legend(title = NULL),
         linetype = guide_legend(title = NULL)) +
  theme(legend.position = c(0.65, 0.15),
        legend.text = element_text(size = rel(1.5)))} %T>%
  ggsave(filename = "tmp/climate_figures/scenarios.png",
         device = "png",
         width = 11.75,
         height = 7.25,
         units = "in")

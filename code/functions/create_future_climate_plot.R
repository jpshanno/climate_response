create_future_climate_plot <- function(observed.data, gcm.data, solrad.coefs, output.file, ...) {
  # observed.data <- tar_read(swg_data)[station_name == "bergland_dam"]
  # gcm.data <- simplify_scenarios(tar_read(loca_simulations)[station_name == "bergland_dam"])
  # solrad.coefs <- tar_read(solrad_coefs)[station_name == "PIEM4"]

  darkblue <- "#3085c7"
  orange <- "#d55e00"
  
  gcm.data <- gcm.data[scenario != "historical"]
  gcm.data[, Climate := fifelse(scenario == "rcp85" & gcm == "gfdl-cm3", "More Sensitive Future Scenario", "Less Sensitive Future Scenario")]
  
  gcm.data <- gcm.data %>%
    calculate_mean_temp() %>% 
    calculate_solar_radiation(coefs = solrad.coefs,
                              stop.on.error = FALSE,
                              return.vector = FALSE) %>% 
    calculate_hargreaves_pet(lambda.MJ.kg = 2.45)
  
  observed.data <- observed.data %>%
    calculate_mean_temp() %>% 
    calculate_solar_radiation(coefs = solrad.coefs,
                              stop.on.error = FALSE,
                              return.vector = FALSE) %>% 
    calculate_hargreaves_pet(lambda.MJ.kg = 2.45)

  dat <- 
    rbind(observed.data[, .(Climate = "Observed Climate", sample_year, sample_date = as.Date(sample_date), pet_cm = abs(pet_cm), precip_cm)],
          gcm.data[, .(Climate, sample_year, sample_date, pet_cm = abs(pet_cm), precip_cm)])
    
  dat[, sample_season := as.climate_season(sample_date, FALSE)]
  
  dat <- dat[
    .SDcols = c("precip_cm", "pet_cm"),
    j = c(list(sample_season = sample_season), lapply(.SD, filled_rolling_sum)),
    by = .(Climate, sample_year)
  ]
  
  dat[
    j = `:=`(
      Climate = factor(
        Climate,
        levels = c("More Sensitive Future Scenario", "Less Sensitive Future Scenario", "Observed Climate"),
        ordered = TRUE
      ),
      sample_season = factor(
        sample_season,
        levels = c("djf", "mam", "jja", "son"),
        labels = c("DJF", "MAM", "JJA", "SON"),
        ordered = TRUE
      )
    )
  ]
  
  setnames(
    x = dat, 
    old = c("pet_cm", "precip_cm"),
    new = c("30-day Potential Evapotranspiration (cm)", "30-day Precipitation (cm)")
  )

  
  base_plot <- ggplot(dat) +
    aes(
      fill = Climate,
      y = Climate
    ) +
    scale_fill_manual(
      values = c(
        "Observed Climate" = "gray30",
        "Less Sensitive Future Scenario" = darkblue,
        "More Sensitive Future Scenario" = orange),
      breaks = c(
        c(
          "Observed Climate",
          "Less Sensitive Future Scenario",
          "More Sensitive Future Scenario")
      )
    ) +
    facet_wrap(~ sample_season, nrow = 1)
  
  fig <- 
    {base_plot + geom_density_ridges(scale = 5, alpha = 0.8, color = "gray70", aes(x = `30-day Potential Evapotranspiration (cm)`))} /
    {base_plot + geom_density_ridges(scale = 5, alpha = 0.8, color = "gray70", aes(x = `30-day Precipitation (cm)`))} +
    plot_layout(guides = "collect") &
    theme_minimal(base_size = 12) &
    theme(
      legend.position = "bottom",
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      legend.title = element_blank()
    )
  
  ggsave(plot = fig, 
         filename = output.file,
         ...)
  
  output.file
}

create_future_climate_table <- function(observed.data, gcm.data, solrad.coefs) {
  # observed.data <- tar_read(swg_data)[station_name == "bergland_dam"]
  # gcm.data <- simplify_scenarios(tar_read(loca_simulations)[station_name == "bergland_dam"])
  # solrad.coefs <- tar_read(solrad_coefs)[station_name == "PIEM4"]
  
  darkblue <- "#3085c7"
  orange <- "#d55e00"
  
  gcm.data <- gcm.data[scenario != "historical"]
  gcm.data[, Climate := fifelse(scenario == "rcp85" & gcm == "gfdl-cm3", "More Sensitive Future Scenario", "Less Sensitive Future Scenario")]
  
  gcm.data <- gcm.data %>%
    calculate_mean_temp() %>% 
    calculate_solar_radiation(coefs = solrad.coefs,
                              stop.on.error = FALSE,
                              return.vector = FALSE) %>% 
    calculate_hargreaves_pet(lambda.MJ.kg = 2.45)
  
  observed.data <- observed.data %>%
    calculate_mean_temp() %>% 
    calculate_solar_radiation(coefs = solrad.coefs,
                              stop.on.error = FALSE,
                              return.vector = FALSE) %>% 
    calculate_hargreaves_pet(lambda.MJ.kg = 2.45)
  
  dat <- 
    rbind(observed.data[, .(Climate = "Observed Climate", sample_year, sample_date = as.Date(sample_date), pet_cm = abs(pet_cm), precip_cm)],
          gcm.data[, .(Climate, sample_year, sample_date, pet_cm = abs(pet_cm), precip_cm)])
  
  dat[, sample_season := as.climate_season(sample_date, FALSE)]
  
  dat[
    j = `:=`(
      Climate = factor(
        Climate,
        levels = c("More Sensitive Future Scenario", "Less Sensitive Future Scenario", "Observed Climate"),
        ordered = TRUE
      ),
      sample_season = factor(
        sample_season,
        levels = c("djf", "mam", "jja", "son"),
        labels = c("DJF", "MAM", "JJA", "SON"),
        ordered = TRUE
      )
    )
  ]
  
  dat <- dat[
    i = sample_season == "JJA",
    j = lapply(.SD[,c("precip_cm", "pet_cm")], sum),
    by = .(Climate, sample_year)]
  
  dat[, deficit := precip_cm - pet_cm]
  
  tab <- dat[
    j = ggdist::median_hdci(.SD, precip_cm, pet_cm, deficit, .width = 0.67),
    by = .(Climate)]
  
  tab <- tab[
    j = .(Climate, 
          `Precipitation (cm)` = glue::glue_data(.SD, "{pretty_round(precip_cm, 2)} ({pretty_round(precip_cm.lower, 2)}, {pretty_round(precip_cm.upper, 2)})"), 
          `Potential Evapotranspiration (cm)` = glue::glue_data(.SD, "{pretty_round(pet_cm, 2)} ({pretty_round(pet_cm.lower, 2)}, {pretty_round(pet_cm.upper, 2)})"), 
          `Water Defecit (cm)` = glue::glue_data(.SD, "{pretty_round(deficit, 2)} ({pretty_round(deficit.lower, 2)}, {pretty_round(deficit.upper, 2)})"))]
  
  flextable::flextable(tab)
}

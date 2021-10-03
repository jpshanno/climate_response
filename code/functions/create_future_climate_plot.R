create_future_climate_plot <- function(observed.data, gcm.data, output.file, ...) {
  # observed.data <- tar_read(external_met)
  # gcm.data <- simplify_scenarios(tar_read(loca_simulations)[station_name == "bergland_dam"])

  darkblue <- "#3085c7"
  orange <- "#d55e00"
  
  gcm.data <- gcm.data[scenario != "historical"]
  gcm.data[, Climate := fifelse(scenario == "rcp85" & gcm == "gfdl-cm3", "More Sensitive Future Scenario", "Less Sensitive Future Scenario")]
  
  dat <- 
    rbind(observed.data[, .(Climate = "Observed Climate", sample_date = as.Date(sample_date), tmin_c, tmax_c, precip_cm)],
          gcm.data[, .(Climate, sample_date, tmin_c, tmax_c, precip_cm)])
    
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
  
  setnames(
    x = dat, 
    old = c("tmin_c", "tmax_c", "precip_cm"),
    new = c("Maximum Temperature (ºC)", "Minimum Temperature (ºC)", "Daily Precipitation (cm)")
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
    facet_wrap(~ sample_season, nrow = 1, scales = "free")
  
  fig <- 
    {base_plot + geom_density_ridges(scale = 5, alpha = 0.8, color = "gray70", aes(x = `Maximum Temperature (ºC)`))} /
    {base_plot + geom_density_ridges(scale = 5, alpha = 0.8, color = "gray70", aes(x = `Minimum Temperature (ºC)`))} /
    {base_plot + geom_density_ridges(scale = 5, alpha = 0.8, color = "gray70", aes(x = `Daily Precipitation (cm)`)) + scale_x_sqrt()} +
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
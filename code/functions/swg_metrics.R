create_gcm_check_plot <- function(observed.data, gcm.data, output.file, ...) {
  # obs <- tar_read(swg_data)[station_name == "bergland_dam"]
  # loca <- tar_read(loca_simulations)[station_name == "bergland_dam" & scenario == "historical"]
  
  green <- as.character(palette.colors()[4])
  blue <- as.character(palette.colors()[6])
  orange <- as.character(palette.colors()[7])
  
  dat <- 
    rbind(observed.data[, .(type = 'Observed', gcm = "Observed", gcm.wrap = "GFDL-CM3", sample_date, sample_month, sample_year, tmin_c, tmax_c, precip_cm)],
          observed.data[, .(type = 'Observed', gcm = "Observed", gcm.wrap = "CCSM4", sample_date, sample_month, sample_year, tmin_c, tmax_c, precip_cm)],
          gcm.data[, .(type = 'LOCA', gcm = toupper(gcm), gcm.wrap = toupper(gcm), sample_date, sample_month, sample_year, tmin_c, tmax_c, precip_cm)])
  
  dat[, sample_season := as.climate_season(sample_date, TRUE)]
  
  base_plot <- 
    ggplot(dat) +
    aes(y = sample_season,
        fill = gcm) +
    scale_y_discrete(labels = toupper) +
    facet_wrap(~gcm.wrap, ncol = 1) +
    theme_minimal() +
    theme(strip.text = element_blank())
  
  
  fig <- 
    {
    {base_plot + geom_density_ridges(aes(x = tmin_c), color = "gray60", alpha = 0.5, scale = 1)} +
    {base_plot + geom_density_ridges(aes(x = tmax_c), color = "gray60", alpha = 0.5, scale = 1)} +
    {base_plot + geom_density_ridges(aes(x = precip_cm), color = "gray60", alpha = 0.5, scale = 1) + scale_x_sqrt()}
  } +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom",
          legend.title = element_blank())
  
  ggsave(plot = fig, 
         filename = output.file,
         ...)

  output.file
}

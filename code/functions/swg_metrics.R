create_gcm_check_plot <- function(observed.data, gcm.data, output.file, ...) {
  # observed.data <- tar_read(swg_data)[station_name == "bergland_dam"]
  # gcm.data <- tar_read(loca_simulations)[station_name == "bergland_dam" & scenario == "historical"]
  
  green <- as.character(palette.colors()[4])
  blue <- as.character(palette.colors()[6])
  orange <- as.character(palette.colors()[7])
  
  dat <- 
    rbind(observed.data[, .(type = 'Observed', gcm = "Observed", gcm.wrap = "GFDL-CM3", sample_date, sample_month, sample_year, tmin_c, tmax_c, precip_cm)],
          observed.data[, .(type = 'Observed', gcm = "Observed", gcm.wrap = "CCSM4", sample_date, sample_month, sample_year, tmin_c, tmax_c, precip_cm)],
          gcm.data[, .(type = 'LOCA, Current Climate', gcm = toupper(gcm), gcm.wrap = toupper(gcm), sample_date, sample_month, sample_year, tmin_c, tmax_c, precip_cm)])
  
  dat[, sample_season := as.climate_season(sample_date, TRUE)]
  
  base_plot <- function(data){
    ggplot(data) +
    aes(y = sample_season,
        fill = type) +
    scale_y_discrete(labels = toupper) +
    scale_fill_manual(values = c("LOCA, Current Climate" = green, "Observed" = "gray70")) +
    # facet_wrap(~gcm.wrap, ncol = 1) +
    theme_minimal()} #+
    # theme(strip.text = element_blank())
  
  fig <- 
    {
      {
      {base_plot(dat[gcm.wrap == "CCSM4"]) + geom_density_ridges(aes(x = tmin_c), color = "gray60", alpha = 0.5, scale = 1) + labs(title = "A. CCSM4", x = "", y = "Season")} +
      {base_plot(dat[gcm.wrap == "CCSM4"]) + geom_density_ridges(aes(x = tmax_c), color = "gray60", alpha = 0.5, scale = 1) + labs(x = "", y = "")} +
      {base_plot(dat[gcm.wrap == "CCSM4"]) + geom_density_ridges(aes(x = precip_cm), color = "gray60", alpha = 0.5, scale = 1) + scale_x_sqrt() + labs(x = "", y = "")}
      } / {
      {base_plot(dat[gcm.wrap == "GFDL-CM3"]) + geom_density_ridges(aes(x = tmin_c), color = "gray60", alpha = 0.5, scale = 1) + labs(title = "B. GFDL-CM3", x = "Minimum Temperature (C)", y = "Season")} +
      {base_plot(dat[gcm.wrap == "GFDL-CM3"]) + geom_density_ridges(aes(x = tmax_c), color = "gray60", alpha = 0.5, scale = 1) + labs(x = "Maximum Temperature (C)", y = "")} +
      {base_plot(dat[gcm.wrap == "GFDL-CM3"]) + geom_density_ridges(aes(x = precip_cm), color = "gray60", alpha = 0.5, scale = 1) + scale_x_sqrt() + labs(x = "Daily Precipitatation (cm)", y = "")}
      }
    } +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom",
          legend.title = element_blank())
  
  ggsave(plot = fig, 
         filename = output.file,
         ...)

  output.file
}

create_swg_table <- function(loca.data, swg.data) {
  # swg.data <- tar_read(swg_simulations_loca)[scenario == "historical"]
  # loca.data <- tar_read(loca_simulations)[station_name == "bergland_dam" & scenario == "historical"]
  
  dat <- 
    rbind(swg.data[, .(type = "SWG", gcm = toupper(gcm), doy = yday(simulation_date), sample_year = simulation_id, sample_month = simulation_month, season = simulation_season, tmin_c, tmax_c, precip_cm)],
          loca.data[, .(type = "LOCA", gcm = toupper(gcm), doy = yday(sample_date), sample_year, sample_month, season = as.climate_season(sample_month, TRUE), tmin_c, tmax_c, precip_cm)])
  
  # Tables for 
  tabs <- 
    melt(dat, measure.vars = c("tmin_c", "tmax_c", "precip_cm"))[
      j = threshold := qnorm(0.975) * sd(value[.SD[["type"]] == "LOCA"]),
      by = .(gcm, season, variable)
    ][
      by = .(type, gcm, season, variable),
      j = .(val_summary = glue::glue_data(.SD,
                                          "{pretty_round(mean(value), 2)} ({pretty_round(mean(value) - qnorm(0.975) * sd(value), 2)}, {pretty_round(mean(value) + qnorm(0.975) * sd(value), 2)}); {pretty_round(100 * sum(abs(value - mean(value)) > threshold) / .N, 2)}% outliers")
      )] %>%
    split(f = .$variable) %>%
    map_dfr(~dcast(.x, gcm + season ~ type,
                   value.var = "val_summary"),
            .id = "Variable")
  
  
  tabs %>% 
    transform(Variable = 
                str_replace_all(
                  Variable,
                  pattern = c(
                    "tmin_c" = "Minimum Temperature",
                    "tmax_c" = "Maximum Temperature",
                    "precip_cm" = "Precipitation"
                  ))) %>%
    transform(season = toupper(season)) %>%
    setnames(old = c("gcm", "season"), new = c("GCM", "Season")) %>%
    flextable::flextable() %>%
    flextable::merge_v(j = 1:2) %>%
    flextable::valign(j = 1:2, valign = "top") %>%
    flextable::autofit()
}

create_swg_plot <- function(loca.data, swg.data, output.file, ...) {
  # swg.data <- tar_read(swg_simulations_loca)[scenario == "historical"]
  # loca.data <- tar_read(loca_simulations)[station_name == "bergland_dam" & scenario == "historical"]
  
  dat <- 
    rbind(swg.data[, .(type = "SWG", gcm = toupper(gcm), doy = yday(simulation_date), sample_year = simulation_id, sample_month = simulation_month, season = simulation_season, tmin_c, tmax_c, precip_cm)],
          loca.data[, .(type = "LOCA", gcm = toupper(gcm), doy = yday(sample_date), sample_year, sample_month, season = as.climate_season(sample_month, TRUE), tmin_c, tmax_c, precip_cm)])

  dat[, season := toupper(season)]
  
  green <- as.character(palette.colors()[4])
  blue <- as.character(palette.colors()[6])
  orange <- as.character(palette.colors()[7])
  
  set.seed(1234)
  
  simulation_samples <- 
    replicate(n = 30,
              expr = sample(1:10000, size = 30, replace = FALSE),
              simplify = FALSE)
  
  base_plot <- 
    ggplot() +
    aes(y = ..scaled..) +
    facet_grid(.~season)
  
  
  for(i in seq_len(30)){
    if(i == 1){
      tmin_plot <- 
        base_plot + aes(x = tmin_c)
    }
    
    tmin_plot <- 
      tmin_plot +
      geom_density(data = dat[type == "SWG" & sample_year %in% simulation_samples[[i]]],
                   alpha = 0.9,
                   size = 0.05,
                   color = 'gray30')
    
    if(i == 30){
      tmin_plot <- 
        tmin_plot +
        geom_density(data = dat[type == "LOCA"],
                     color = green)
    }
  }
  
  
  for(i in seq_len(30)){
    
    if(i == 1){
      tmax_plot <- 
        base_plot + aes(x = tmax_c)
    }
    
    tmax_plot <- 
      tmax_plot +
      geom_density(data = dat[type == "SWG" & sample_year %in% simulation_samples[[i]]],
                   alpha = 0.9,
                   size = 0.05,
                   color = 'gray30')
    
    if(i == 30){
      tmax_plot <- 
        tmax_plot +
        geom_density(data = dat[type == "LOCA"],
                     color = green)
    }
  }
  
  for(i in seq_len(30)){
    
    if(i == 1){
      t_plot <- 
        base_plot
    }
    
    t_plot <- 
      t_plot +
      geom_density(data = dat[type == "SWG" & sample_year %in% simulation_samples[[i]]],
                   aes(x = tmin_c),
                   alpha = 0.9,
                   size = 0.05,
                   color = blue) +
      geom_density(data = dat[type == "SWG" & sample_year %in% simulation_samples[[i]]],
                   aes(x = tmax_c),
                   alpha = 0.9,
                   size = 0.05,
                   color = orange)
    
    if(i == 30){
      t_plot <- 
        t_plot +
        geom_density(data = dat[type == "LOCA"],
                     aes(x = tmin_c),
                     color = blue)+
        geom_density(data = dat[type == "LOCA"],
                     aes(x = tmax_c),
                     color = orange)
    }
  }
  
  for(i in seq_len(30)){
    
    if(i == 1){
      precip_plot <-
        base_plot + aes(x = precip_cm)
    }
    
    precip_plot <- 
      precip_plot +
      geom_density(data = dat[type == "SWG" & sample_year %in% simulation_samples[[i]], 
                              .(precip_cm = sum(precip_cm)),
                              by = .(sample_year, season)],
                   alpha = 0.9,
                   size = 0.05,
                   color = 'gray30')
    
    if(i == 30){
      precip_plot <- 
        precip_plot +
        geom_density(data = dat[type == "LOCA", 
                                .(precip_cm = sum(precip_cm)),
                                by = .(sample_year, season)],
                     color = green)
    }
  }
  
  
  fig <- {t_plot + 
      labs(y = NULL,
           x = "Daily Minimum and Maximum Temperature (C)")} / 
    {precip_plot +  
        labs(y = NULL, 
             x = "Monthly Precipitation (cm)")} &
    plot_annotation(tag_levels = "A") & 
    theme_minimal(base_size = 13)
 
  ggsave(plot = fig, filename = output.file, ...)
  
  output.file
}
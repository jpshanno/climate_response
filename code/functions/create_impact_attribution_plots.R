create_total_impact_plot <- function(proportions, output.file, ...){
  
  # proportions <- tar_read(analysis_simulations)$proportions
  green <- "#2b730b"
  brown <- "#756953" 
  blue <- "#3096bd"
  
  orange <- "#d55e00"
  darkblue <- "#3085c7"
  gold <- "#e69f00"
  red <- "#dc2b2bff"
  
  ci_width <- 0.67
  ci_function <- function(...) {hdci(...)}
  point_function <- function(...) {median(...)}

  baseline <- 
    proportions[site_status == "Control" & scenario == "historical" & month(simulation_date) %in% 5:10, 
                .(Connectivity = prop_within_5cm_max_wl,
                  Drawdown = (1-prop_above_neg_50),
                  Inundation = prop_above_neg_10,
                  simulation_month = month(simulation_date, label = TRUE, abbr = TRUE))] %>% 
    melt(id.vars = "simulation_month",
         value.name = "Probability")
  
  zeroing <- baseline[
    j = .(med_prob = median(Probability)),
    by = .(simulation_month, variable)
  ]
  
  baseline[
    i = zeroing,
    on = c("simulation_month", "variable"),
    j = Probability := Probability - med_prob
  ]
  
  summary_dat <- 
    proportions[site_status != "Control" & scenario != "historical" & month(simulation_date) %in% 5:10,
                .(connectivity = point_function(prop_within_5cm_max_wl),
                  connectivity_min = ci_function(prop_within_5cm_max_wl, .width = ci_width)[[1]],
                  connectivity_max = ci_function(prop_within_5cm_max_wl, .width = ci_width)[[2]],
                  drawdown = point_function((1-prop_above_neg_50)),
                  drawdown_min = ci_function((1-prop_above_neg_50), .width = ci_width)[[1]],
                  drawdown_max = ci_function((1-prop_above_neg_50), .width = ci_width)[[2]],
                  inundation = point_function(prop_above_neg_10),
                  inundation_min = ci_function(prop_above_neg_10, .width = ci_width)[[1]],
                  inundation_max = ci_function(prop_above_neg_10, .width = ci_width)[[2]]),
                by = .(site_status, 
                       scenario = fcase(scenario == "rcp45", "Warm & Dry", 
                                        scenario == "rcp85", "Hot & Wet"),
                       simulation_month = month(simulation_date, label = TRUE, abbr = TRUE))] %>% 
    transform(scenario = factor(scenario, levels = c("Warm & Dry", "Hot & Wet"), ordered = TRUE)) %>% 
    melt(id.vars = c("site_status", "scenario", "simulation_month"),
         measure.vars = list(Probability = c("connectivity", "inundation", "drawdown"),
                             ymin = c("connectivity_min", "inundation_min", "drawdown_min"),
                             ymax = c("connectivity_max", "inundation_max", "drawdown_max"))) %>% 
    transform(variable = fct_relabel(variable, ~fcase(.x=="1", "Connectivity", .x=="2", "Inundation", .x=="3", "Drawdown"))) %>%
    transform(site_status = str_replace_all(site_status, "Treated", "Non-Forested")) %>%
    transform(site_status = factor(site_status, levels = c("Non-Forested", "Future Forested"), ordered = TRUE)) %>%
    .[
      i = zeroing,
      on = c("simulation_month", "variable"),
      j = `:=`(
        Probability = Probability - med_prob,
        ymin = ymin - med_prob,
        ymax = ymax - med_prob
      )
    ]
  
  
  
  fig <-
    ggplot(summary_dat) + 
    ylab("Change in Probability of Occurrence") +
    stat_gradientinterval(data = baseline,
                          color = NA,
                          aes(x = simulation_month,
                              y = Probability,
                              slab_alpha = stat(-pmin(ci_width, pmax(abs(1-2*cdf), 0))))) + # Shading covers 67% of data, starts to fade at 10%
    geom_crossbar(data = baseline,
                  stat = "summary",
                  fun = median,
                  show.legend = FALSE,
                  aes(x = simulation_month,
                      y = Probability,
                      color = "Baseline: Modeled Black Ash, Current Climate")
    ) +
    geom_text(
      data = summary_dat[j = .(site_status = unique(site_status), max_prob = max(ymax)), by = .(variable)][order(variable, site_status), .(variable, site_status, label = LETTERS[1:6], x = 0.75, y = 0.75 * max_prob)],
      aes(x = x, y = y, label = label),
      size = 5,
      vjust = 0
    ) +
    geom_pointrange(aes(x = simulation_month,
                        y = Probability,
                        ymin = ymin,
                        ymax = ymax,
                        color = scenario),
                    # show.legend = FALSE,
                    position = position_dodge(width = 0.6)) +
    scale_color_manual(values = c(`Baseline: Modeled Black Ash, Current Climate` = "gray30",
                                  `Warm & Dry` = darkblue,
                                  `Hot & Wet` = orange),
                       breaks = c("Warm & Dry", "Hot & Wet", "Baseline: Modeled Black Ash, Current Climate")) +
    facet_grid(variable ~ site_status,
               scales = "free") +
    guides(
      color = guide_legend(
        title = NULL)
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom",
          axis.title.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.border = element_rect(fill = NA))
  
  ggsave(plot= fig,
         filename = output.file,
         ...)
  
  output.file

}


create_eab_impact_plot <- function(proportions, output.file, ...){
  
  # proportions <- tar_read(analysis_simulations)
  green <- "#2b730b"
  brown <- "#756953" 
  blue <- "#3096bd"
  
  orange <- "#d55e00"
  darkblue <- "#3085c7"
  gold <- "#e69f00"
  red <- "#dc2b2bff"
  
  ci_width <- 0.67
  ci_function <- function(...) {hdci(...)}
  point_function <- function(...) {median(...)}
  
  # Create baseline of black ash forest under reach climate condition
  eab_baseline <- proportions[
    i = month(simulation_date) %in% 5:10 & site_status == "Control",
    j = .(
      simulation_month = month(simulation_date, label = TRUE, abbr = TRUE),
      scenario = fcase(scenario == "historical", "Current Climate",
                       scenario == "rcp45", "Warm & Dry", 
                       scenario == "rcp85", "Hot & Wet"),
      Connectivity = prop_within_5cm_max_wl,
      Drawdown = 1-prop_above_neg_50,
      Inundation = prop_above_neg_10
    )
  ] %>% 
    melt(id.vars = c("scenario", "simulation_month"),
         value.name = "Probability")
  
  eab_zeroing <- eab_baseline[
    j = .(med_prob = median(Probability)),
    by = .(scenario, simulation_month, variable)
  ]
  
  eab_baseline[
    i = eab_zeroing,
    on = c("scenario", "simulation_month", "variable"),
    j = Probability := Probability - med_prob
  ]
  
  eab_summary_dat <- 
    proportions[month(simulation_date) %in% 5:10 & site_status != "Control",
                .(connectivity = point_function(prop_within_5cm_max_wl),
                  connectivity_min = ci_function(prop_within_5cm_max_wl, .width = ci_width)[[1]],
                  connectivity_max = ci_function(prop_within_5cm_max_wl, .width = ci_width)[[2]],
                  drawdown = point_function((1-prop_above_neg_50)),
                  drawdown_min = ci_function((1-prop_above_neg_50), .width = ci_width)[[1]],
                  drawdown_max = ci_function((1-prop_above_neg_50), .width = ci_width)[[2]],
                  inundation = point_function(prop_above_neg_10),
                  inundation_min = ci_function(prop_above_neg_10, .width = ci_width)[[1]],
                  inundation_max = ci_function(prop_above_neg_10, .width = ci_width)[[2]]),
                by = .(site_status, 
                       scenario = fcase(scenario == "historical", "Current Climate",
                                        scenario == "rcp45", "Warm & Dry", 
                                        scenario == "rcp85", "Hot & Wet"),
                       simulation_month = month(simulation_date, label = TRUE, abbr = TRUE))] %>% 
    transform(scenario = factor(scenario, levels = c("Current Climate", "Warm & Dry", "Hot & Wet"), ordered = TRUE)) %>% 
    melt(id.vars = c("site_status", "scenario", "simulation_month"),
         measure.vars = list(Probability = c("connectivity", "inundation", "drawdown"),
                             ymin = c("connectivity_min", "inundation_min", "drawdown_min"),
                             ymax = c("connectivity_max", "inundation_max", "drawdown_max"))) %>% 
    transform(variable = fct_relabel(variable, ~fcase(.x=="1", "Connectivity", .x=="2", "Inundation", .x=="3", "Drawdown"))) %>%
    transform(site_status = factor(site_status,
                                   levels = c("Treated", "Future Forested"),
                                   labels = c("Non-Forested", "Future Forested"),
                                   ordered = TRUE)) %>%
    .[
      i = eab_zeroing,
      on = c("scenario", "simulation_month", "variable"),
      j = `:=`(
        Probability = Probability - med_prob,
        ymin = ymin - med_prob,
        ymax = ymax - med_prob
      )
    ]
  
  fig <-
    ggplot(eab_summary_dat) + 
    ylab("Change in Probability of Occurrence") +
    stat_gradientinterval(data = eab_baseline,
                          color = NA,
                          aes(x = simulation_month,
                              y = Probability,
                              slab_alpha = stat(-pmin(ci_width, pmax(abs(1-2*cdf), 0))))) + # Shading covers 67% of data, starts to fade at 10%
    geom_crossbar(data = eab_baseline,
                  stat = "summary",
                  fun = median,
                  show.legend = FALSE,
                  aes(x = simulation_month,
                      y = Probability,
                      color = "Baseline: Modeled Black Ash")
    ) +
    geom_pointrange(aes(x = simulation_month,
                        y = Probability,
                        ymin = ymin,
                        ymax = ymax,
                        color = site_status),
                    position = position_dodge(width = 0.6)) +
    scale_color_manual(values = c(`Baseline: Modeled Black Ash` = "gray30",
                                  `Non-Forested` = gold,
                                  `Future Forested` = green),
                       breaks = c("Non-Forested", "Future Forested", "Baseline: Modeled Black Ash")) +
    facet_grid(variable ~ scenario,
               scales = "free") +
    guides(
      color = guide_legend(
        title = NULL)
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom",
          axis.title.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.border = element_rect(fill = NA))
  
  ggsave(plot= fig,
         filename = output.file,
         ...)
  
  output.file
  
}

create_climate_impact_plot <- function(proportions, output.file, ...){
  
  # proportions <- tar_read(analysis_simulations)
  green <- "#2b730b"
  brown <- "#756953" 
  blue <- "#3096bd"
  
  orange <- "#d55e00"
  darkblue <- "#3085c7"
  gold <- "#e69f00"
  red <- "#dc2b2bff"
  
  ci_width <- 0.67
  ci_function <- function(...) {hdci(...)}
  point_function <- function(...) {median(...)}
  
  # Create baseline of each alternative vegetative condition under current climate conditions 
  climate_baseline <- proportions[
    i = month(simulation_date) %in% 5:10 & scenario == "historical",
    j = .(
      simulation_month = month(simulation_date, label = TRUE, abbr = TRUE),
      site_status = fcase(site_status == "Control", "Black Ash",
                          site_status == "Treated", "Non-Forested",
                          site_status == "Future Forested", "Future Forested"),
      Connectivity = prop_within_5cm_max_wl,
      Drawdown = 1-prop_above_neg_50,
      Inundation = prop_above_neg_10
    )
  ] %>% 
    melt(id.vars = c("site_status", "simulation_month"),
         value.name = "Probability")
  
  climate_zeroing <- climate_baseline[
    j = .(med_prob = median(Probability)),
    by = .(site_status, simulation_month, variable)
  ]
  
  climate_baseline[
    i = climate_zeroing,
    on = c("site_status", "simulation_month", "variable"),
    j = Probability := Probability - med_prob
  ]
  
  climate_summary_dat <- 
    proportions[month(simulation_date) %in% 5:10 & scenario != "historical",
                .(connectivity = point_function(prop_within_5cm_max_wl),
                  connectivity_min = ci_function(prop_within_5cm_max_wl, .width = ci_width)[[1]],
                  connectivity_max = ci_function(prop_within_5cm_max_wl, .width = ci_width)[[2]],
                  drawdown = point_function((1-prop_above_neg_50)),
                  drawdown_min = ci_function((1-prop_above_neg_50), .width = ci_width)[[1]],
                  drawdown_max = ci_function((1-prop_above_neg_50), .width = ci_width)[[2]],
                  inundation = point_function(prop_above_neg_10),
                  inundation_min = ci_function(prop_above_neg_10, .width = ci_width)[[1]],
                  inundation_max = ci_function(prop_above_neg_10, .width = ci_width)[[2]]),
                by = .(site_status, 
                       scenario = fcase(scenario == "rcp45", "Warm & Dry", 
                                        scenario == "rcp85", "Hot & Wet"),
                       simulation_month = month(simulation_date, label = TRUE, abbr = TRUE))] %>% 
    transform(scenario = factor(scenario, levels = c("Warm & Dry", "Hot & Wet"), ordered = TRUE)) %>% 
    melt(id.vars = c("site_status", "scenario", "simulation_month"),
         measure.vars = list(Probability = c("connectivity", "inundation", "drawdown"),
                             ymin = c("connectivity_min", "inundation_min", "drawdown_min"),
                             ymax = c("connectivity_max", "inundation_max", "drawdown_max"))) %>% 
    transform(variable = fct_relabel(variable, ~fcase(.x=="1", "Connectivity", .x=="2", "Inundation", .x=="3", "Drawdown"))) %>%
    transform(site_status = factor(site_status,
                                   levels = c("Control", "Treated", "Future Forested"),
                                   labels = c("Black Ash", "Non-Forested", "Future Forested"),
                                   ordered = TRUE)) %>%
    .[
      i = climate_zeroing,
      on = c("site_status", "simulation_month", "variable"),
      j = `:=`(
        Probability = Probability - med_prob,
        ymin = ymin - med_prob,
        ymax = ymax - med_prob
      )
    ]
  
  
  fig <- ggplot(climate_summary_dat) + 
    ylab("Change in Probability of Occurrence") +
    stat_gradientinterval(data = climate_baseline,
                          color = NA,
                          aes(x = simulation_month,
                              y = Probability,
                              slab_alpha = stat(-pmin(ci_width, pmax(abs(1-2*cdf), 0))))) + # Shading covers 67% of data, starts to fade at 10%
    geom_crossbar(data = climate_baseline,
                  stat = "summary",
                  fun = median,
                  show.legend = FALSE,
                  aes(x = simulation_month,
                      y = Probability,
                      color = "Baseline: Current Climate")
    ) +
    geom_pointrange(aes(x = simulation_month,
                        y = Probability,
                        ymin = ymin,
                        ymax = ymax,
                        color = scenario),
                    position = position_dodge(width = 0.6)) +
    scale_color_manual(values = c(`Baseline: Current Climate` = "gray30",
                                  `Warm & Dry` = darkblue,
                                  `Hot & Wet` = orange),
                       breaks = c("Warm & Dry", "Hot & Wet", "Baseline: Current Climate")) +
    facet_grid(variable ~ site_status,
               scales = "free") +
    guides(
      color = guide_legend(
        title = NULL)
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom",
          axis.title.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.border = element_rect(fill = NA))
  
  ggsave(plot= fig,
         filename = output.file,
         ...)
  
  output.file
  
}

create_impact_attribution_plot <- function(proportions, output.file, ...) {
  
  green <- "#2b730b"
  brown <- "#756953" 
  blue <- "#3096bd"
  
  orange <- "#d55e00"
  darkblue <- "#3085c7"
  gold <- "#e69f00"
  red <- "#dc2b2b"
  
  baseline_proportions <- proportions[
    i = month(simulation_date) %in% 6:8 & scenario == "historical" & site_status == "Control",
    j = .(
      site,
      simulation_date,
      gcm,
      connectivity = prop_within_5cm_max_wl,
      drawdown = 1-prop_above_neg_50,
      inundation = prop_above_neg_10
    )
  ] %>%
    melt(measure.vars = c("connectivity", "drawdown", "inundation"),
         value.name = "Control")
  
  # Create baseline of each vegetation type under current climate conditions
  current_climate_proportions <- proportions[
    i = month(simulation_date) %in% 6:8 & scenario == "historical" & site_status != "Control",
    j = .(
      site,
      site_status,
      simulation_date,
      gcm,
      connectivity = prop_within_5cm_max_wl,
      drawdown = 1-prop_above_neg_50,
      inundation = prop_above_neg_10
    )
  ] %>%
    melt(measure.vars = c("connectivity", "drawdown", "inundation"),
         value.name = "historical")
  
  # Calculate impact of EAB as difference relative to black ash stands under
  # the current climate
  eab_impact <- 
    proportions[
      i = month(simulation_date) %in% 6:8 & scenario == "historical" & site_status != "Control",
      j = .(
        site,
        site_status,
        simulation_date,
        gcm,
        connectivity = prop_within_5cm_max_wl,
        drawdown = 1-prop_above_neg_50,
        inundation = prop_above_neg_10
      )
    ] %>%
    melt(measure.vars = c("connectivity", "drawdown", "inundation"),
         value.name = "Alternative") %>%
    .[baseline_proportions,
      on = c("site", "simulation_date", "gcm", "variable")
    ] %>%
    .[
      j = `Vegetation\nChange Impact` := Alternative - Control
    ] %>%
    .[
      j = median_hdci(.SD, `Vegetation\nChange Impact`, .width = 0.67),
      by = .(variable, site_status),
    ] %>% 
    melt(
      measure.vars = "Vegetation\nChange Impact",
      variable.name = "impact_source"
    )
  
  # Calculate impact of Warm & Dry climate scenario relative to the same
  # vegetation conditions under the current climate
  cc_l_impact <- 
    proportions[
      i = month(simulation_date) %in% 6:8 & scenario == "rcp45" & gcm == "ccsm4" & site_status != "Control",
      j = .(
        site,
        site_status,
        simulation_date,
        gcm,
        connectivity = prop_within_5cm_max_wl,
        drawdown = 1-prop_above_neg_50,
        inundation = prop_above_neg_10
      )
    ] %>%
    melt(measure.vars = c("connectivity", "drawdown", "inundation"),
         value.name = "future") %>%
    .[current_climate_proportions,
      on = c("site", "site_status", "simulation_date", "gcm", "variable")
    ] %>%
    .[!is.na(future)] %>% 
    .[
      j = `Climate Impact\n(Warm & Dry)` := future - historical
    ] %>%
    .[
      j = median_hdci(.SD, `Climate Impact\n(Warm & Dry)`, .width = 0.67),
      by = .(variable, site_status),
    ] %>% 
    melt(
      measure.vars = "Climate Impact\n(Warm & Dry)",
      variable.name = "impact_source"
    )
  
  # Calculate impact of Hot & Wet climate scenario relative to the same
  # vegetation conditions under the current climate
  cc_h_impact <- 
    proportions[
      i = month(simulation_date) %in% 6:8 & scenario == "rcp85" & gcm == "gfdl-cm3" & site_status != "Control",
      j = .(
        site,
        site_status,
        simulation_date,
        gcm,
        connectivity = prop_within_5cm_max_wl,
        drawdown = 1-prop_above_neg_50,
        inundation = prop_above_neg_10
      )
    ] %>%
    melt(measure.vars = c("connectivity", "drawdown", "inundation"),
         value.name = "future") %>%
    .[current_climate_proportions,
      on = c("site", "site_status", "simulation_date", "gcm", "variable")
    ] %>%
    .[!is.na(future)] %>% 
    .[
      j = `Climate Impact\n(Hot & Wet)` := future - historical
    ] %>%
    .[
      j = median_hdci(.SD, `Climate Impact\n(Hot & Wet)`, .width = 0.67),
      by = .(variable, site_status),
    ] %>%
    melt(
      measure.vars = "Climate Impact\n(Hot & Wet)",
      variable.name = "impact_source",
    )
  
  tab_dat <- rbindlist(
    list(
      # eab_impact[, c(.SD, list(scenario = "Warm & Dry"))],
      # eab_impact[, c(.SD, list(scenario = "Hot & Wet"))],
      eab_impact,
      cc_l_impact,
      cc_h_impact
    )
  )
  
  # Format Output Table
  tab_dat[, site_status := factor(site_status,
                                  levels = c("Treated", "Future Forested"),
                                  labels = c("Non-Forested", "Future Forested"),
                                  ordered = TRUE)]
  tab_dat[, variable := factor(variable,
                               levels = c("connectivity", "inundation", "drawdown"),
                               labels = c("Connectivity", "Inundation", "Drawdown"),
                               ordered = TRUE)]
  
  attribution_fig <-
    ggplot(tab_dat) +
    aes(
      x = impact_source,
      y = value,
      ymin = .lower,
      ymax = .upper,
      color = impact_source,
      shape = impact_source
    ) +
    geom_hline(
      aes(
        yintercept = 0
      ),
      color = "gray10",
      linetype = "dashed"
    ) +
    geom_pointrange(
      fill = "white",
      position = position_dodge(width = 0.5)
    ) +
    ylab("Change in Probability of Occurrence") +
    facet_grid(
      variable ~ site_status,
      scales = "free"
    ) +
    theme_minimal(base_size = 14) +
    scale_color_manual(
      name = NULL,
      values = c(
        "Vegetation\nChange Impact" = brown,
        "Climate Impact\n(Warm & Dry)" = darkblue,
        "Climate Impact\n(Hot & Wet)" = orange
      )
    ) +
    scale_shape_manual(
      name = NULL,
      values = c(
        "Vegetation\nChange Impact" = 21,
        "Climate Impact\n(Warm & Dry)" = 22,
        "Climate Impact\n(Hot & Wet)" = 25
      )
    ) +
    theme(legend.position = "bottom",
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.border = element_rect(fill = NA))
  
  ggsave(
    output.file,
    attribution_fig,
    ...
  )
}


create_impact_attribution_table <- function(proportions) {
  
  # Create baseline of black ash forest under current climate conditions
  baseline_proportions <- proportions[
    i = month(simulation_date) %in% 7:9 & scenario == "historical" & site_status == "Control",
    j = .(
      site,
      simulation_date,
      gcm,
      connectivity = prop_within_5cm_max_wl,
      drawdown = 1-prop_above_neg_50,
      inundation = prop_above_neg_10
    )
  ] %>%
  melt(measure.vars = c("connectivity", "drawdown", "inundation"),
         value.name = "Control")
  
  # Create baseline of each vegetation type under current climate conditions
  current_climate_proportions <- proportions[
    i = month(simulation_date) %in% 7:9 & scenario == "historical" & site_status != "Control",
    j = .(
      site,
      site_status,
      simulation_date,
      gcm,
      connectivity = prop_within_5cm_max_wl,
      drawdown = 1-prop_above_neg_50,
      inundation = prop_above_neg_10
    )
  ] %>%
    melt(measure.vars = c("connectivity", "drawdown", "inundation"),
         value.name = "historical")
  
  # Calculate impact of EAB as difference relative to black ash stands under
  # the current climate
  eab_impact <- 
    proportions[
      i = month(simulation_date) %in% 7:9 & scenario == "historical" & site_status != "Control",
      j = .(
        site,
        site_status,
        simulation_date,
        gcm,
        connectivity = prop_within_5cm_max_wl,
        drawdown = 1-prop_above_neg_50,
        inundation = prop_above_neg_10
        )
      ] %>%
    melt(measure.vars = c("connectivity", "drawdown", "inundation"),
         value.name = "Alternative") %>%
    .[baseline_proportions,
      on = c("site", "simulation_date", "gcm", "variable")
    ] %>%
    .[
      j = `EAB Impact` := Alternative - Control
    ] %>%
    .[
      j = .(`EAB Impact` = pretty_hdci(100 * `EAB Impact`)),
      by = .(variable, site_status),
    ]
  
  # Calculate impact of less sensitive climate scenario relative to the same
  # vegetation conditions under the current climate
  cc_l_impact <- 
    proportions[
      i = month(simulation_date) %in% 7:9 & scenario == "rcp45" & gcm == "ccsm4" & site_status != "Control",
      j = .(
        site,
        site_status,
        simulation_date,
        gcm,
        connectivity = prop_within_5cm_max_wl,
        drawdown = 1-prop_above_neg_50,
        inundation = prop_above_neg_10
      )
    ] %>%
      melt(measure.vars = c("connectivity", "drawdown", "inundation"),
           value.name = "future") %>%
      .[current_climate_proportions,
        on = c("site", "site_status", "simulation_date", "gcm", "variable")
       ] %>%
      .[!is.na(future)] %>% 
      .[
        j = `Climate Impact (Less Sensitive)` := future - historical
      ] %>%
      .[
        j = .(`Climate Impact\n(Less Sensitive)` = pretty_hdci(100 * `Climate Impact (Less Sensitive)`)),
        by = .(variable, site_status),
      ] %>% 
      .[]
  
  # Calculate impact of more sensitive climate scenario relative to the same
  # vegetation conditions under the current climate
  cc_h_impact <- 
    proportions[
      i = month(simulation_date) %in% 7:9 & scenario == "rcp85" & gcm == "gfdl-cm3" & site_status != "Control",
      j = .(
        site,
        site_status,
        simulation_date,
        gcm,
        connectivity = prop_within_5cm_max_wl,
        drawdown = 1-prop_above_neg_50,
        inundation = prop_above_neg_10
      )
    ] %>%
    melt(measure.vars = c("connectivity", "drawdown", "inundation"),
         value.name = "future") %>%
    .[current_climate_proportions,
      on = c("site", "site_status", "simulation_date", "gcm", "variable")
    ] %>%
    .[!is.na(future)] %>% 
    .[
      j = `Climate Impact (More Sensitive)` := future - historical
    ] %>%
    .[
      j = .(`Climate Impact\n(More Sensitive)` = pretty_hdci(100 * `Climate Impact (More Sensitive)`)),
      by = .(variable, site_status),
    ] %>% 
    .[]
  
  tab_dat <- reduce(
    list(eab_impact, cc_l_impact, cc_h_impact),
    merge,
    by = c("variable", "site_status")
  )
  
  # Format Output Table
  tab_dat[, site_status := factor(site_status,
                                    levels = c("Treated", "Future Forested"),
                                        labels = c("Non-Forested", "Future Forested"),
                                        ordered = TRUE)]
  tab_dat[, variable := factor(variable,
                               levels = c("connectivity", "inundation", "drawdown"),
                               labels = c("Connectivity", "Inundation", "Drawdown"),
                               ordered = TRUE)]
  setorder(tab_dat, variable, site_status)
  setnames(tab_dat, c("variable", "site_status"), c("Variable", "Vegetation Condition"))
  
  flextable::flextable(tab_dat) %>% 
    flextable::merge_v(j = 1) %>%
    flextable::valign(j = 1, valign = "top") %>%
    flextable::align(align = "center", part = c("all")) %>% 
    flextable::autofit()
  
}

pretty_hdci <- function(x, unit = "%") {
  
  int <- median_hdci(x, .width = 0.67)
  
  glue::glue("{pretty_round(int$y, 2)}{unit}; ({pretty_round(int$ymin, 2)}, {pretty_round(int$ymax, 2)})")
  
}

simplify_scenarios <- function(data) {
  
  data[!(gcm == "gfdl-cm3" & scenario == "rcp45") & !(gcm == "ccsm4" & scenario == "rcp85")]
  
}

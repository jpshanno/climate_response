source("code/load_project.R")
proportions <- tar_read(analysis_simulations)$proportions
green <- "#2b730b"
brown <- "#756953" 
blue <- "#3096bd"

orange <- "#d55e00"
darkblue <- "#3085c7"
gold <- "#e69f00"
red <- "#dc2b2b"

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
                     scenario = fcase(scenario == "rcp45", "Less Sensitive", 
                                      scenario == "rcp85", "More Sensitive"),
                     simulation_month = month(simulation_date, label = TRUE, abbr = TRUE))] %>% 
  transform(scenario = factor(scenario, levels = c("Less Sensitive", "More Sensitive"), ordered = TRUE)) %>% 
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
  geom_pointrange(aes(x = simulation_month,
                      y = Probability,
                      ymin = ymin,
                      ymax = ymax,
                      color = scenario),
                  # show.legend = FALSE,
                  position = position_dodge(width = 0.6)) +
  scale_color_manual(values = c(`Baseline: Modeled Black Ash, Current Climate` = "gray30",
                                `Less Sensitive` = darkblue,
                                `More Sensitive` = orange),
                     breaks = c("Less Sensitive", "More Sensitive", "Baseline: Modeled Black Ash, Current Climate")) +
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

ggsave(
  "~/phd/defense_figures/total_impact_zeroed.png",
  fig,
  width = 10,
  height = 5,
  units = "in"
)



# EAB Impact --------------------------------------------------------------

# Create baseline of black ash forest under reach climate condition
eab_baseline <- proportions[
  i = month(simulation_date) %in% 5:10 & site_status == "Control",
  j = .(
    simulation_month = month(simulation_date, label = TRUE, abbr = TRUE),
    scenario = fcase(scenario == "historical", "Current Climate",
                     scenario == "rcp45", "Less Sensitive", 
                     scenario == "rcp85", "More Sensitive"),
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
                                      scenario == "rcp45", "Less Sensitive", 
                                      scenario == "rcp85", "More Sensitive"),
                     simulation_month = month(simulation_date, label = TRUE, abbr = TRUE))] %>% 
  transform(scenario = factor(scenario, levels = c("Current Climate", "Less Sensitive", "More Sensitive"), ordered = TRUE)) %>% 
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

eab_fig <-
  ggplot(eab_summary_dat) + 
  ylab("Probability of Occurrence") +
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

ggsave(
  "~/phd/defense_figures/eab_impact_zeroed.png",
  eab_fig,
  width = 10,
  height = 5,
  units = "in"
)

# Climate Impacts ---------------------------------------------------------

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
                     scenario = fcase(scenario == "rcp45", "Less Sensitive", 
                                      scenario == "rcp85", "More Sensitive"),
                     simulation_month = month(simulation_date, label = TRUE, abbr = TRUE))] %>% 
  transform(scenario = factor(scenario, levels = c("Less Sensitive", "More Sensitive"), ordered = TRUE)) %>% 
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


climate_fig <- ggplot(climate_summary_dat) + 
  ylab("Probability of Occurrence") +
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
                                `Less Sensitive` = darkblue,
                                `More Sensitive` = orange),
                     breaks = c("Less Sensitive", "More Sensitive", "Baseline: Current Climate")) +
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

ggsave(
  "~/phd/defense_figures/climate_impact_zeroed.png",
  climate_fig,
  width = 10,
  height = 5,
  units = "in"
)

# Impact comparison table -------------------------------------------------

# Create baseline of black ash forest under current climate conditions
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
    j = `EAB Impact` := Alternative - Control
  ] %>%
  .[
    j = median_hdci(.SD, `EAB Impact`, .width = 0.67),
    by = .(variable, site_status),
  ] %>% 
  melt(
    measure.vars = "EAB Impact",
    variable.name = "impact_source"
  )

# Calculate impact of less sensitive climate scenario relative to the same
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
    j = `Climate Impact (Less Sensitive)` := future - historical
  ] %>%
  .[
    j = median_hdci(.SD, `Climate Impact (Less Sensitive)`, .width = 0.67),
    by = .(variable, site_status),
  ] %>% 
  melt(
    measure.vars = "Climate Impact (Less Sensitive)",
    variable.name = "impact_source"
  )

# Calculate impact of more sensitive climate scenario relative to the same
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
    j = `Climate Impact (More Sensitive)` := future - historical
  ] %>%
  .[
    j = median_hdci(.SD, `Climate Impact (More Sensitive)`, .width = 0.67),
    by = .(variable, site_status),
  ] %>%
  melt(
    measure.vars = "Climate Impact (More Sensitive)",
    variable.name = "impact_source",
  )

tab_dat <- rbindlist(
  list(
    # eab_impact[, c(.SD, list(scenario = "Less Sensitive"))],
    # eab_impact[, c(.SD, list(scenario = "More Sensitive"))],
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
      "EAB Impact" = brown,
      "Climate Impact (Less Sensitive)" = darkblue,
      "Climate Impact (More Sensitive)" = orange
    )
  ) +
  scale_shape_manual(
    name = NULL,
    values = c(
      "EAB Impact" = 21,
      "Climate Impact (Less Sensitive)" = 22,
      "Climate Impact (More Sensitive)" = 25
    )
  ) +
  theme(legend.position = "bottom",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_rect(fill = NA))

ggsave(
  "~/phd/defense_figures/impact_attribution.png",
  attribution_fig,
  dpi = 300,
  width = 7,
  height = 8,
  units = "in"
)

# Summer Changes ----------------------------------------------------------

observed.data <- tar_read(swg_data)[station_name == "bergland_dam"]
gcm.data <- simplify_scenarios(tar_read(loca_simulations)[station_name == "bergland_dam"])
solrad.coefs <- tar_read(solrad_coefs)[station_name == "PIEM4"]

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
  rbind(observed.data[as.climate_season(sample_date, FALSE) == "jja", .(Climate = "Observed Climate", sample_year, sample_date = as.Date(sample_date), pet_cm = abs(pet_cm), precip_cm)],
        gcm.data[as.climate_season(sample_date, FALSE) == "jja", .(Climate, sample_year, sample_date, pet_cm = abs(pet_cm), precip_cm)])

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
  .SDcols = c("precip_cm", "pet_cm"),
  j = lapply(.SD, sum),
  by = .(Climate, sample_year)
]

dat[, deficit_cm := precip_cm - pet_cm]

dat[
  j = ggdist::median_hdci(.SD, precip_cm, pet_cm, deficit_cm, .width = 0.67),
  by = .(Climate)][
    j = .(Climate, 
          `Precipitation (cm)` = glue::glue_data(.SD, "{pretty_round(precip_cm, 2)} ({pretty_round(precip_cm.lower, 2)}, {pretty_round(precip_cm.upper, 2)})"), 
          `Potential Evapotranspiration (cm)` = glue::glue_data(.SD, "{pretty_round(pet_cm, 2)} ({pretty_round(pet_cm.lower, 2)}, {pretty_round(pet_cm.upper, 2)})"), 
          `Water Balance (cm)` = glue::glue_data(.SD, "{pretty_round(deficit_cm, 2)} ({pretty_round(deficit_cm.lower, 2)}, {pretty_round(deficit_cm.upper, 2)})"))
    ]


setnames(
  dat,
  c("precip_cm", "pet_cm", "deficit_cm"),
  c("JJA Precipitation (cm)", "JJA Potential Evapotranspiration (cm)",
    "JJA Water Balance (cm)")
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
  )

drivers_plots <- 
  {base_plot + geom_density_ridges(scale = 5, alpha = 0.8, color = "gray70", aes(x = `JJA Precipitation (cm)`))} /
  {base_plot + geom_density_ridges(scale = 5, alpha = 0.8, color = "gray70", aes(x = `JJA Potential Evapotranspiration (cm)`))} /
  {base_plot + geom_density_ridges(scale = 5, alpha = 0.8, color = "gray70", aes(x = `JJA Water Balance (cm)`))} +
  plot_layout(guides = "collect") &
  theme_minimal(base_size = 12) &
  theme(
    legend.position = "bottom",
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    legend.title = element_blank()
  )

ggsave(
  "~/phd/defense_figures/drivers_of_change.png",
  drivers_plots,
  width = 6,
  height = 8,
  units = "in",
  dpi = 300
)

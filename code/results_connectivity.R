source("code/load_project.R")
library(ggtext)

tar_load(wetland_simulation_summaries)
tar_load(model_params)

scenarios <- 
  data.table(scenario = c("historical", "rcp45", "rcp85"),
             gcm = c("ccsm4", "ccsm4", "gfdl-cm3"))

green <- "#2b730b"
brown <- "#756953" 
blue <- "#3096bd"

orange <- "#d55e00"
darkblue <- "#3085c7"
gold <- "#e69f00"
red <- "#dc2b2bff"


model_params[site_status == "Control", .(maxWL = map_dbl(params, pluck, "maxWL")), by = .(site)]



proportions <- 
  wetland_simulation_summaries[["proportions"]]

connectivity <- 
  proportions[scenarios,
              on = c("gcm", "scenario"),
              ggdist::median_qi(prop_within_5cm_max_wl, .width = 0.5),
              by = .(site_status = factor(fcase(site_status == "Control", "Ash Forest", 
                                                site_status == "Future Forested", "Replacement Forest",
                                                site_status == "Treated", "Non-Forested"),
                                          levels = c("Ash Forest",
                                                     "Replacement Forest",
                                                     "Non-Forested"),
                                          ordered = TRUE), 
                     scenario = fcase(scenario == "historical", "Current Climate", 
                                      scenario == "rcp45", "Less Sensitive", 
                                      scenario == "rcp85", "More Sensitive"), 
                     simulation_date,
                     simulation_season = as.climate_season(simulation_date, TRUE)),
              nomatch = NULL] 


connectivity_baseline <- 
  rbind(connectivity[site_status == "Ash Forest" & scenario == "Current Climate",
                     .(site_status = "Non-Forested", simulation_date, scenario = "control", y, ymin, ymax)],
        connectivity[site_status == "Ash Forest" & scenario == "Current Climate",
                     .(site_status = "Replacement Forest", simulation_date, scenario = "control", y, ymin, ymax)]) %>% 
  transform(site_status = factor(site_status,
                                 levels = c("Replacement Forest",
                                            "Non-Forested"),
                                 ordered = TRUE))


# Minimum median daily probabilties:
dcast(connectivity[, .(prop = min(y)), by = .(site_status, scenario)],
      scenario ~ site_status,
      value.var = "prop")

connectivity_panel <- 
  ggplot(connectivity[site_status != "Ash Forest"]) +
  aes(x = simulation_date, 
      y = y, 
      color = scenario, 
      linetype = scenario, 
      fill = scenario) +
  geom_rect(data = data.frame(season = c("mam", "son"),
                              xmin = as.Date(c("2009-03-01","2009-09-01")),
                              xmax = as.Date(c("2009-05-31","2009-11-30")),
                              ymin = c(0, 0),
                              ymax = c(1, 1),
                              y = 0,
                              scenario = c("Replacement Forest"),
                              simulation_date = as.Date("2009-01-01")),
            aes(xmin = xmin,
                ymin = ymin,
                xmax = xmax,
                ymax = ymax),
            color = NA,
            fill = 'gray80',
            linetype = 'solid') +
  geom_ribbon(data = connectivity_baseline,
              aes(ymin = ymin, 
                  ymax = ymax),
              color = NA, 
              alpha = 0.2,
              show.legend = FALSE) +
  geom_line(data = connectivity_baseline,
            show.legend = FALSE) +
  geom_ribbon(aes(ymin = ymin, 
                  ymax = ymax),
              color = NA, 
              alpha = 0.3,
              show.legend = FALSE) + 
  geom_line(show.legend = FALSE) +
  scale_fill_manual(values = c(`Current Climate` = green,
                               `Less Sensitive` = darkblue,
                               `More Sensitive` = orange,
                               control = 'gray20')) +
  scale_color_manual(values = c(`Current Climate` = green,
                                `Less Sensitive` = darkblue,
                                `More Sensitive` = orange,
                                control = 'black')) +
  scale_linetype_manual(values = c(`Current Climate` = "longdash",
                                   `Less Sensitive` = "dotdash",
                                   `More Sensitive` = "twodash",
                                   control = "solid")) +
  scale_x_date(name = NULL, 
               expand = c(0, 0), 
               # date_labels = "%b",
               breaks = function(lims){seq(lims[1], lims[2], by = "month")},
               labels = function(brks){labs <- format(brks, "%b"); labs[1 + 2*(0:5)] <- ""; labs}) +
  ylab("Probability of Connectivity") +
  geom_richtext(data = connectivity[site_status != "Ash Forest", .SD[1], by = .(site_status)],
                aes(x = as.Date("2009-01-01"), 
                    y = 0.1,
                    label = site_status),
                fill = NA,
                label.color = NA,
                color = 'black',
                vjust = 0.3,
                size = 16*5/14,
                hjust = 0) +
  facet_wrap(~site_status) +
  theme_minimal(base_size = 18) +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.text = element_blank(),
        axis.ticks.x = element_line(color = "grey60",
                                    size = rel(0.5)))


{connectivity_panel / 
  {{water_levels_panel / 
  connectivity_panel} *
  scale_x_date(name = NULL, 
               expand = c(0, 0), 
               date_labels = "%b") * 
      theme_minimal()}}  +
  plot_layout(design = 
                "AAA
                 BBB
                 BBB")

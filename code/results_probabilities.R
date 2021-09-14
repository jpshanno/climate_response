source("code/load_project.R")
library(ggtext)
library(ggdist)


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

month_breaks <- 
  function(lims){seq(lims[1], lims[2], by = "month")}

month_lims <- 
  function(brks){labs <- format(brks, "%b"); labs[1 + 2*(0:5)] <- ""; labs}


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

drawdown <- 
  proportions[scenarios,
              on = c("gcm", "scenario"),
              ggdist::median_qi((1-prop_above_neg_50), .width = 0.5),
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
                     simulation_date),
              nomatch = NULL] 


drawdown_baseline <- 
  rbind(drawdown[site_status == "Ash Forest" & scenario == "Current Climate",
                 .(site_status = "Non-Forested", simulation_date, scenario = "control", y, ymin, ymax)],
        drawdown[site_status == "Ash Forest" & scenario == "Current Climate",
                 .(site_status = "Replacement Forest", simulation_date, scenario = "control", y, ymin, ymax)]) %>% 
  transform(site_status = factor(site_status,
                                 levels = c("Replacement Forest",
                                            "Non-Forested"),
                                 ordered = TRUE))


inundation <- 
  proportions[scenarios,
              on = c("gcm", "scenario"),
              ggdist::median_qi(prop_above_neg_10, .width = 0.5),
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
                     simulation_date),
              nomatch = NULL] 


inundation_baseline <- 
  rbind(inundation[site_status == "Ash Forest" & scenario == "Current Climate",
                   .(site_status = "Non-Forested", simulation_date, scenario = "control", y, ymin, ymax)],
        inundation[site_status == "Ash Forest" & scenario == "Current Climate",
                   .(site_status = "Replacement Forest", simulation_date, scenario = "control", y, ymin, ymax)]) %>% 
  transform(site_status = factor(site_status,
                                 levels = c("Replacement Forest",
                                            "Non-Forested"),
                                 ordered = TRUE))


# Table -------------------------------------------------------------------

# Minimum median daily probabilties:
dcast(inundation[, .(prop = min(y)), by = .(site_status, scenario)],
      scenario ~ site_status,
      value.var = "prop")

# Minimum median daily probabilties:
dcast(drawdown[, .(prop = max(y)), by = .(site_status, scenario)],
      scenario ~ site_status,
      value.var = "prop")


# Minimum median daily probabilties:
dcast(connectivity[, .(prop = min(y)), by = .(site_status, scenario)],
      scenario ~ site_status,
      value.var = "prop")

# Color or denote cells (maybe up/down arrows) based on current conditions (do
# current conditions fall within, above, or below the QI)

table_baseline <- 
  proportions[site_status == "Control" & scenario == "historical" & month(simulation_date) %in% 5:10, 
              .(connectivity_base = median(prop_within_5cm_max_wl),
                drawdown_base = median((1-prop_above_neg_50)),
                inundation_base = median(prop_above_neg_10)),
              by = .(simulation_season = as.climate_season(simulation_date, TRUE))]

table_dat <- 
  proportions[site_status != "Control" & month(simulation_date) %in% 5:10, 
              .(connectivity = median(prop_within_5cm_max_wl),
                connectivity_min = ggdist::qi(prop_within_5cm_max_wl, .width = 0.5)[[1]],
                connectivity_max = ggdist::qi(prop_within_5cm_max_wl, .width = 0.5)[[2]],
                drawdown = median((1-prop_above_neg_50)),
                drawdown_min = ggdist::qi((1-prop_above_neg_50), .width = 0.5)[[1]],
                drawdown_max = ggdist::qi((1-prop_above_neg_50), .width = 0.5)[[2]],
                inundation = median(prop_above_neg_10),
                inundation_min = ggdist::qi(prop_above_neg_10, .width = 0.5)[[1]],
                inundation_max = ggdist::qi(prop_above_neg_10, .width = 0.5)[[2]]),
              by = .(site_status, 
                     scenario = fcase(scenario == "historical", "Current Climate", 
                                      scenario == "rcp45", "CCSM4 (4.5)", 
                                      scenario == "rcp85", "GFDL-CM3 (8.5)"),
                     simulation_season = as.climate_season(simulation_date, TRUE))] %>% 
  .[table_baseline, 
    `:=`(connectivity_status = fcase(between(connectivity_base, connectivity_min, connectivity_max), "⟷", 
                                         connectivity_base > connectivity_max, "↓",
                                         connectivity_base < connectivity_min, "↑"),
         drawdown_status = fcase(between(drawdown_base, drawdown_min, drawdown_max), "⟷",
                                     drawdown_base > drawdown_max, "↓",
                                     drawdown_base < drawdown_min, "↑"),
         inundation_status = fcase(between(inundation_base, inundation_min, inundation_max), "⟷",
                                       inundation_base > inundation_max, "↓",
                                       inundation_base < inundation_min, "↑")),
    on = "simulation_season"] %>% 
  .[, .(site_status, scenario, simulation_season,
        connectivity = sprintf(paste0("%0.3f\n(%0.2f-%0.2f)\n", connectivity_status), connectivity, connectivity_min, connectivity_max),
        drawdown = sprintf(paste0("%0.3f\n(%0.2f-%0.2f)\n", drawdown_status), drawdown, drawdown_min, drawdown_max),
        inundation = sprintf(paste0("%0.3f\n(%0.2f-%0.2f)\n",inundation_status), inundation, inundation_min, inundation_max))] %>% 
  melt(id.vars = c("site_status", "scenario", "simulation_season")) %>%
  split(f = .$site_status) %>% 
  imap(~{rbind(dcast(.x, variable ~ simulation_season + scenario, value.var = "scenario")[1,],
         dcast(.x, variable ~ simulation_season + scenario))} %>% 
        set_names(~substr(.x, 1, 3)) %>%
        .[1, var:=.y]) %>% 
  rbindlist() %>%
  setnames("var", "")

table_header <- 
  names(table_dat)

table_dat %>% 
  set_names(paste0("X", 1:ncol(.))) %>% 
  flextable::flextable() %>% 
  flextable::set_header_labels(values = set_names(as.list(table_header), paste0("X", 1:ncol(table_dat)))) %>% 
  flextable::merge_h(part = "header") %>% 
  flextable::merge_v()

# Panels ------------------------------------------------------------------

inundation_panel <- 
  ggplot(inundation[site_status != "Ash Forest"]) +
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
  geom_ribbon(data = inundation_baseline,
              aes(ymin = ymin, 
                  ymax = ymax),
              color = NA, 
              alpha = 0.2,
              show.legend = FALSE) +
  geom_line(data = inundation_baseline,
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
               breaks = month_breaks,
               labels = month_lims) +
  ylab("Probability of Inundation") +
  geom_richtext(data = inundation[site_status != "Ash Forest", .SD[1], by = .(site_status)],
                aes(x = as.Date("2009-01-01"), 
                    y = 0.1,
                    label = site_status),
                fill = NA,
                label.color = NA,
                color = 'black',
                vjust = 0.3,
                size = 16*5/14,
                hjust = 0) +
  facet_wrap(~site_status) 

drawdown_panel <- 
  ggplot(drawdown[site_status != "Ash Forest"]) +
  aes(x = simulation_date, 
      y = y, 
      color = scenario, 
      linetype = scenario, 
      fill = scenario) +
  geom_rect(data = data.frame(season = c("mam", "son"),
                              xmin = as.Date(c("2009-03-01","2009-09-01")),
                              xmax = as.Date(c("2009-05-31","2009-11-30")),
                              ymin = c(0, 0),
                              ymax = c(0.5, 0.5),
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
  geom_ribbon(data = drawdown_baseline,
              aes(ymin = ymin, 
                  ymax = ymax),
              color = NA, 
              alpha = 0.2,
              show.legend = FALSE) +
  geom_line(data = drawdown_baseline,
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
               breaks = month_breaks,
               labels = month_lims) +
  ylab("Probability of Drawdown") +
  geom_richtext(data = drawdown[site_status != "Ash Forest", .SD[1], by = .(site_status)],
                aes(x = as.Date("2009-01-01"), 
                    y = 0.1,
                    label = site_status),
                fill = NA,
                label.color = NA,
                color = 'black',
                vjust = 0.3,
                size = 16*5/14,
                hjust = 0) +
  facet_wrap(~site_status)


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
               breaks = month_breaks,
               labels = month_lims) +
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
  facet_wrap(~site_status) 


{connectivity_panel /
  inundation_panel /
  drawdown_panel} *
  theme_minimal(base_size = 18) *
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.text = element_blank(),
        axis.ticks.x = element_line(color = "grey60",
                                    size = rel(0.5)))




# Summary Panel -----------------------------------------------------------
library(forcats)
ci_width <- 0.67
ci_function <- function(...) {ggdist::hdci(...)}
point_function <- function(...) {median(...)}


tar_load(c(training_data, testing_data, control_optimization))

obs_cols <- c(
  "site", "sample_date", "sample_year", "wl_initial_cm", "treatment",
  "site_status")

observations <- rbindlist(list(
  testing_data[site_status == "Control" & month(sample_date) %in% 5:10, ..obs_cols],
  training_data$control[month(sample_date) %in% 5:10, ..obs_cols]
  ))

observations[, simulation_month := month(sample_date, label = TRUE, abbr = TRUE)]

observations[control_optimization[, .(site, max_wl = map_dbl(params, "maxWL"))],
             on = c("site"),
             max_wl := i.max_wl]

obs_probs <- 
  observations[!is.na(wl_initial_cm),
             .(Connectivity = sum(wl_initial_cm >= (max_wl - 5), na.rm = TRUE) / .N,
               Inundation = sum(wl_initial_cm >= -10, na.rm = TRUE) / .N,
               Drawdown = sum(wl_initial_cm <= -50, na.rm = TRUE) / .N),
             by = .(site, simulation_month)] %>%
  melt(id.vars = c("site", "simulation_month"), value.name = "Probability") %>%
  {rbind(.[, c(.SD, site_status = "Future Forested")],
        .[, c(.SD, site_status = "Treated")])}


summary_baseline <- 
  proportions[site_status == "Control" & scenario == "historical" & month(simulation_date) %in% 5:10, 
              .(Connectivity = point_function(prop_within_5cm_max_wl),
                Drawdown = point_function((1-prop_above_neg_50)),
                Inundation = point_function(prop_above_neg_10)),
              by = .(simulation_month = month(simulation_date, label = TRUE, abbr = TRUE))] %>% 
  melt(id.vars = "simulation_month",
       value.name = "Probability")

# Using qi() not HDCI because of how skewed the probabilities are
summary_dat <- 
  proportions[site_status != "Control" & month(simulation_date) %in% 5:10,
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
  transform(variable = fct_relabel(variable, ~fcase(.x=="1", "Connectivity", .x=="2", "Inundation", .x=="3", "Drawdown")))


test_dat <- 
  proportions[site_status != "Control" & month(simulation_date) %in% 5:10, 
              .(Connectivity = prop_within_5cm_max_wl,
                Drawdown = (1-prop_above_neg_50),
                Inundation = prop_above_neg_10,
                simulation_month = month(simulation_date, label = TRUE, abbr = TRUE),
                scenario = factor(scenario, 
                                  levels = c("historical", "rcp45", "rcp85"),
                                  labels = c("Current Climate", "Less Sensitive", "More Sensitive"), 
                                  ordered = TRUE),
                site_status = factor(site_status,
                                     levels = c("Future Forested", 
                                                "Treated"),
                                     labels = c("Replacement Forest",
                                                "Non-Forested"),
                                     ordered = TRUE))] %>% 
  melt(measure.vars = c("Connectivity", "Inundation", "Drawdown"),
       value.name = "Probability",
       variable.factor = TRUE)

test <- 
  proportions[site_status == "Control" & scenario == "historical" & month(simulation_date) %in% 5:10, 
              .(Connectivity = prop_within_5cm_max_wl,
                Drawdown = (1-prop_above_neg_50),
                Inundation = prop_above_neg_10,
                simulation_month = month(simulation_date, label = TRUE, abbr = TRUE))] %>% 
  melt(id.vars = "simulation_month",
       value.name = "Probability")


ggplot(summary_dat) + 
  stat_gradientinterval(data = test,
                        color = NA,
                        aes(x = simulation_month,
                            y = Probability,
                            slab_alpha = stat(-pmin(ci_width, pmax(abs(1-2*cdf), 0))))) + # Shading covers 67% of data, starts to fade at 10%
  geom_boxplot(data = summary_baseline,
               aes(x = simulation_month,
                   y = Probability),
               color = "gray30") +
  geom_pointrange(aes(x = simulation_month,
                      y = Probability,
                      ymin = ymin,
                      ymax = ymax,
                      color = scenario),
                  position = position_dodge(width = 0.6)) +
  scale_fill_manual(values = c(`Current Climate` = green,
                               `Less Sensitive` = darkblue,
                               `More Sensitive` = orange,
                               control = 'gray20')) +
  scale_color_manual(values = c(`Current Climate` = green,
                                `Less Sensitive` = darkblue,
                                `More Sensitive` = orange,
                                control = 'black')) +
  facet_grid(variable ~ site_status,
             scales = "free") +
  labs(caption = str_wrap("All ranges represent 67% of the simulations as HDCI. Points represent median. Gray shaded area represents model simulations for control conditions under the current climate.", width = 100)) +
  guides(color = guide_legend(title = NULL)) +
  theme_minimal(base_size = 18) +
  theme(legend.position = c(0.01, 0.3),
        legend.justification = c(0, 1),
        panel.grid.minor.y = element_blank(),
        panel.border = element_rect(fill = NA))

ggplot(rbindlist(list(simulation = test, obs = obs_probs), idcol = "type", fill = TRUE)) + 
  geom_boxplot(aes(x = simulation_month,
                   y = Probability,
                   fill = type)) +
  facet_grid(~ variable,
             scales = "free") 
  
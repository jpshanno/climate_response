create_probabilities_plot <- function(proportions, output.file, ...){
  
  # proportions <- tar_read(wetland_simulation_summaries)[["proportions"]]
  green <- "#2b730b"
  brown <- "#756953" 
  blue <- "#3096bd"
  
  orange <- "#d55e00"
  darkblue <- "#3085c7"
  gold <- "#e69f00"
  red <- "#dc2b2bff"
  
  ci_width <- 0.67
  ci_function <- function(...) {ggdist::hdci(...)}
  point_function <- function(...) {median(...)}

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
    transform(site_status = factor(site_status, levels = c("Non-Forested", "Future Forested"), ordered = TRUE))
  
  baseline <- 
    proportions[site_status == "Control" & scenario == "historical" & month(simulation_date) %in% 5:10, 
                .(Connectivity = prop_within_5cm_max_wl,
                  Drawdown = (1-prop_above_neg_50),
                  Inundation = prop_above_neg_10,
                  simulation_month = month(simulation_date, label = TRUE, abbr = TRUE))] %>% 
    melt(id.vars = "simulation_month",
         value.name = "Probability")
  
  fig <- ggplot(summary_dat) + 
    ggdist::stat_gradientinterval(data = baseline,
                          color = NA,
                          aes(x = simulation_month,
                              y = Probability,
                              slab_alpha = stat(-pmin(ci_width, pmax(abs(1-2*cdf), 0))))) + # Shading covers 67% of data, starts to fade at 10%
    geom_crossbar(data = baseline,
                  stat = "summary",
                  fun = median,
                 aes(x = simulation_month,
                     y = Probability,
                     color = "Modeled Control")
                 ) +
    geom_pointrange(aes(x = simulation_month,
                        y = Probability,
                        ymin = ymin,
                        ymax = ymax,
                        color = scenario),
                    position = position_dodge(width = 0.6)) +
    geom_text(
      data = summary_dat[j = .(site_status = unique(site_status), max_prob = max(ymax)), by = .(variable)][order(variable, site_status), .(variable, site_status, label = LETTERS[1:6], x = 0.75, y = 0.75 * max_prob)],
      aes(x = x, y = y, label = label),
      size = 5,
      vjust = 0
    ) +
    scale_color_manual(values = c(`Modeled Control` = "gray30",
                                  `Less Sensitive` = darkblue,
                                  `More Sensitive` = orange),
                       breaks = c("Less Sensitive", "More Sensitive", "Modeled Control")) +
    facet_grid(variable ~ site_status,
               scales = "free") +
    # labs(caption = str_wrap("All ranges represent 67% of the simulations as HDCI. Points represent median. Gray shaded area represents model simulations for control conditions under the current climate.", width = 100)) +
    guides(color = guide_legend(title = NULL)) +
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

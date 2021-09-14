wetland_model_metrics <- function(data) {
  
  data[!is.na(wl_initial_cm + wl_hat),
        .(r2 = hydroGOF::rPearson(wl_hat, wl_initial_cm),
          mae = hydroGOF::mae(wl_hat, wl_initial_cm)),
        by = .(site, water_year, site_status)] %>% 
    melt(id.vars = c("site", "water_year", "site_status"),
         variable.name = "metric")
  
}
  

create_wetland_model_metrics_plot <- function(metrics, outlier.sites, output.file, ...) {
  
  metrics <- copy(metrics)
  
  metrics[, site_class := fifelse(site %in% outlier.sites, "Outlier Sites", "Analysis Sites")]
  metrics[, metric := str_replace_all(metric, c("r2" = "R<sup>2</sup>", "med_err" = "Median Error (cm)", "rmse" = "RMSE"))]
  
  fig <- metrics %>%
    split(f = interaction(.$site_class, .$metric)) %>%
    map(~{
      ggplot(.x) +
        aes(x = site, y = value) +
        geom_boxplot() +
        labs(y = .x[["metric"]][1]) +
        facet_grid(~ site_class, scales = "free", space = "free")
      }) %>%
    reduce(`+`)  +
    plot_layout(design = "AAAAB\nCCCCD\nEEEEF", guides = "collect") &
    theme_minimal(base_size = 12) +
    theme(axis.title.x = element_blank(),
          axis.title.y = ggtext::element_markdown(),
          legend.title = element_blank())
  
  ggsave(plot = fig, filename = output.file, ...)
  
  output.file
  
}

calculate_predicted_probabilities <- function(data, max_wl_data) {
  
  data[max_wl_data, max_wl := max_wl, on = "site"]
  
  inundation_probs <- 
    data[!is.na(wl_initial_cm),
        .(Observed = prob_inundated(wl_initial_cm),
          Predicted = prob_inundated(wl_hat)),
        by = .(site, water_year, site_status)] %>% 
    melt(id.vars = c("site", "water_year", "site_status"), variable.name = "type", value.name = "Inundation")
  
  connectivity_probs <- 
    data[!is.na(wl_initial_cm),
        .(Observed = prob_connected(wl_initial_cm, max_wl),
          Predicted = prob_connected(wl_hat, max_wl)),
        by = .(site, water_year, site_status)] %>% 
    melt(id.vars = c("site", "water_year", "site_status"), variable.name = "type", value.name = "Connectivity")
  
  drawdown_probs <- 
    data[!is.na(wl_initial_cm),
        .(Observed = prob_drawdown(wl_initial_cm),
          Predicted = prob_drawdown(wl_hat)),
        by = .(site, water_year, site_status)] %>% 
    melt(id.vars = c("site", "water_year", "site_status"), variable.name = "type", value.name = "Drawdown")
  
  reduce(list(inundation_probs, connectivity_probs, drawdown_probs), merge, by = c("site", "water_year", "site_status", "type")) %>%
    melt(id.vars = c("site", "water_year", "site_status", "type"), value.name = "Probability")
  
}
    
prob_inundated <- function(water_level) {
  
  assertthat::assert_that(is.numeric(water_level))
  assertthat::assert_that(any(!is.na(water_level)))
  
  sum(water_level >= -10) / length(water_level)
  
}

prob_connected <- function(water_level, spill_level) {
  
  assertthat::assert_that(is.numeric(water_level))
  assertthat::assert_that(any(!is.na(water_level)))
  
  sum(water_level >= (spill_level - 5)) / length(water_level)
  
}

prob_drawdown <- function(water_level) {
  
  assertthat::assert_that(is.numeric(water_level))
  assertthat::assert_that(any(!is.na(water_level)))
  
  sum(water_level < -50) / length(water_level)
  
}

check_probability_differences <- function(data) {
  
  mod <- lme4::lmer(Probability ~ 0 + type*site_status + (0 + type*site_status | site), data = data)
  
  emm <- emmeans::emmeans(mod, specs = ~ type | site_status)
  
  out <- broom::tidy(multcomp::cld(emm))
  
  out
  
}

test_predicted_probabilities <- function(data) {
  
  purrr::map_dfr(
    .x = split(data, f = data$variable),
    .f = check_probability_differences,
    .id = "variable"
    ) %>%
  as.data.table()
  
}

create_wetland_model_probability_table <- function(tests) {
  
  simp <- tests[,
    .(variable, site_status, type,
      estimate = glue::glue_data(tests, "{pretty_round(estimate, 2)} ({pretty_round(conf.low, 2)}, {pretty_round(conf.high, 2)})"),
      within_variable_group = letters[as.integer(.group)])] %>%
  dcast(variable + site_status + within_variable_group ~ type, value.var = "estimate") %>%
  .[, .(variable, site_status, Observed, Predicted, within_variable_group)]
  
  simp[, variable := factor(variable, levels = c("Connectivity", "Inundation", "Drawdown"), labels = c("Connectivity", "Inundation", "Drawdown"), ordered = TRUE)]
  simp[order(variable)]
  
  flextable::flextable(simp) %>%
    flextable::merge_v(j = 1:2) %>%
    flextable::valign(j = 1:2, valign = "top") %>%
    flextable::autofit()
  
}

create_wetland_model_probability_plot <- function(data, tests, output.file, ...) {
  
  fig <-
    ggplot(data) +
    aes(x = site_status, y = Probability) +
    geom_point(aes(color = type), position = position_jitterdodge()) +
    geom_pointrange(data = tests,
                    aes(y = estimate, ymin = conf.low, ymax = conf.high, group = type),
                    position = position_dodge(width = 0.75)) +
    facet_wrap(. ~ variable, scales = "free") +
    theme_minimal(base_size = 16) +
    theme(legend.position = "bottom",
          axis.title.x = element_blank(),
          legend.title = element_blank())
  
  ggsave(filename = output.file, plot = fig, ...)
  
  output.file
  
}

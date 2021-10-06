calculate_wetland_model_metrics <- function(data, max_wl_data) {
  
  data[max_wl_data, max_wl := max_wl, on = "site"]
  # data[, wl_range := median_na(.SD[!is.na(wl_initial_cm + wl_hat), .(wl_range = diff(range(wl_initial_cm, na.rm = TRUE))), by = .(water_year)][["wl_range"]]), by = .(site)]
  data[!is.na(wl_initial_cm + wl_hat),
        .(r2 = cor(wl_hat, wl_initial_cm, use = "pairwise.complete.obs"),
          med_err = median(wl_hat - wl_initial_cm, na.rm = TRUE),
          rmedse = rmedse(wl_hat, wl_initial_cm),
          rmedse_range = rmedse(wl_hat, wl_initial_cm) /  diff(range(wl_initial_cm, na.rm = TRUE))),
        by = .(site, water_year, site_status)] %>% 
    melt(id.vars = c("site", "water_year", "site_status"),
         variable.name = "metric")
  
}

rmedse <- function(sim, obs){
  sqrt(median((sim - obs)^2, na.rm = TRUE))
}

rrmedse <- function(sim, obs){
  rmedse(sim, obs) / diff(range(obs, na.rm = TRUE))
}

create_wetland_model_metrics_plot <- function(data, output.file, metrics, ...) {
  
  data <- copy(data[metric %in% metrics])
  
  data[, metric := str_replace_all(metric, c("r2" = "R<sup>2</sup>", "med_err" = "Median Error (cm)", "^rmedse$" = "RMedSE", "rmedse_range" = "rRMedSE"))]
  
  # Keeping split + combine to label y-axis with metric. Could melt the data and
  # then facet by variable. But this was set up from when there were outlier 
  # sites removed
  fig <- data %>%
    split(f = .$metric) %>%
    map(~{
      ggplot(.x) +
        aes(x = site, y = value) +
        geom_boxplot() +
        labs(y = .x[["metric"]][1])
      }) %>%
    reduce(`+`)  +
    plot_layout(ncol = 1, guides = "collect") &
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
    flextable::valign(j = 1:2, valign = "top")
  
}

create_wetland_model_probability_plot <- function(data, tests, output.file, ...) {
  
  data <- copy(data)
  tests <- copy(tests)
  
  data[, variable := factor(variable, levels = c("Connectivity", "Inundation", "Drawdown"), labels = c("Connectivity", "Inundation", "Drawdown"), ordered = TRUE)]
  tests[, variable := factor(variable, levels = c("Connectivity", "Inundation", "Drawdown"), labels = c("Connectivity", "Inundation", "Drawdown"), ordered = TRUE)]
  
  fig <-
    ggplot(data) +
    aes(x = site_status, y = Probability) +
    ylab("Probability of Occurrence") +
    geom_point(aes(color = type, shape = type), position = position_jitterdodge()) +
    geom_pointrange(data = tests,
                    aes(y = estimate, ymin = conf.low, ymax = conf.high, group = type),
                    position = position_dodge(width = 0.75)) +
    facet_wrap(. ~ variable, scales = "free") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom",
          axis.title.x = element_blank(),
          legend.title = element_blank())
  
  ggsave(filename = output.file, plot = fig, ...)
  
  output.file
  
}


create_wetland_model_metrics_table <- function(data) {
  
  tab <- data[,
    .(`Minimum\nObserved` = pretty_round(min(value), 2),
      `Maximum\nObserved` = pretty_round(max(value), 2),
      Mean = pretty_round(mean(value), 2),
      `Median` = pretty_round(median(value), 2)),
    by = .(`Model\nMetric` = metric, `Site\nCondition` = site_status)]
  
  tab[, 
      `Model\nMetric` := str_replace_all(`Model\nMetric`, c("r2" = "R~2~", "med_err" = "Median Error (cm)", "^rmedse$" = "RMedSE", "rmedse_range" = "rRMedSE"))]
  
  flextable::flextable(tab) %>%
    flextable::merge_v(j = 1:2) %>%
    flextable::valign(j = 1, valign = "top") %>%
    flextable::align(align = "right") %>%
    flextable::align(j = 1, align = "left") %>%
    flextable::align(j = 2, align = "center") %>%
    flextable::align(part = "header", align = "center")
  
}
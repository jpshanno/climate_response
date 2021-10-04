#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param data
#' @param ...
#' @return
#' @author Joe Shannon
#' @export
#' 
#' 
create_pet_impact_plot <- function(data, output.file, ...) {

  # data <- tar_read(hydrological_components)
  
  green <- "#2b730b"
  brown <- "#756953" 
  blue <- "#3096bd"
  
  orange <- "#d55e00"
  darkblue <- "#3085c7"
  gold <- "#e69f00"
  red <- "#dc2b2bff"
  
  data[
    j = `:=`(
      Climate = fcase(
        scenario == "historical", "Current Climate",
        scenario == "rcp45", "Less Sensitive",
        scenario == "rcp85", "More Sensitive"
      ),
      site_status = fcase(
        site_status == "Control", "Modeled Black Ash",
        site_status == "Treated", "Non-Forested",
        site_status == "Future Forested", "Future Forested"
      )
    )
  ]

  data <- data[
    .SDcols = patterns("^aet_"),
    j = c(list(simulation_date = simulation_date),
          lapply(.SD, filled_rolling_mean)),
    by = .(site_status, Climate)]
  
  data <- data[
    month(simulation_date) %in% 4:10]
  
  # Combine the two gcms current climate conditions
  # TODO: come up with a better solution
  data[
    i = Climate == "Current Climate",
    j = `:=`(
      aet_hat = mean(aet_hat),
      aet_hat.lower = mean(aet_hat.lower),
      aet_hat.upper = mean(aet_hat.upper),
      aet_effect = mean(aet_effect),
      aet_effect.lower = mean(aet_effect.lower),
      aet_effect.upper = mean(aet_effect.upper)
    ),
    by = .(site_status, Climate, simulation_date)]
  
  fig <- ggplot(data) +
    aes(x = simulation_date) +
    ggpattern::geom_ribbon_pattern(aes(ymin = aet_hat.lower, ymax = aet_hat.upper, fill = "Est. AET"), color = NA, alpha = 0.2, pattern_fill = "#ffffffff", pattern_colour = "#ffffffff") +
    geom_ribbon(aes(ymin = aet_effect.lower, ymax = aet_effect.upper, fill = "Effective AET"), color = NA, alpha = 0.2) +
    geom_line(aes(y = aet_hat, color = "Est. AET", linetype = "Est. AET")) +
    geom_line(aes(y = aet_effect, color = "Effective AET", linetype = "Effective AET")) +
    facet_grid(site_status~Climate) +
    theme_minimal()
  
  # fixed_climate_fig <- ggplot(data) +
  #   aes(
  #     x = simulation_date,
  #     y = y,
  #     ymin = ymin,
  #     ymax = ymax,
  #     fill = site_status,
  #     color = site_status
  #   ) +
  #   ylab("30-day Mean PET (cm)") +
  #   xlab(NULL) +
  #   ggtitle("B.") +
  #   geom_ribbon(alpha = 0.3, color = NA) +
  #   geom_line() +
  #   scale_color_manual(
  #     name = "Vegetation Conditions",
  #     values = c(
  #       "Modeled Black Ash" = brown,
  #       "Future Forested" = gold,
  #       "Non-Forested" = green 
  #     ),
  #     breaks = c(
  #       "Modeled Black Ash",
  #       "Non-Forested",
  #       "Future Forested"
  #     )
  #     ) +
  #   scale_fill_manual(
  #     name = "Vegetation Conditions",
  #     values = c(
  #       "Modeled Black Ash" = brown,
  #       "Future Forested" = gold,
  #       "Non-Forested" = green 
  #     ),
  #     breaks = c(
  #       "Modeled Black Ash",
  #       "Non-Forested",
  #       "Future Forested"
  #     )
  #   ) +
  #   facet_wrap(~Climate)
  # 
  # fixed_veg_fig <- ggplot(data) +
  #   aes(
  #     x = simulation_date,
  #     y = y,
  #     ymin = ymin,
  #     ymax = ymax,
  #     fill = Climate,
  #     color = Climate
  #   ) +
  #   ylab("30-day Mean PET (cm)") +
  #   xlab(NULL) +
  #   ggtitle("B.") +
  #   geom_ribbon(alpha = 0.3, color = NA) +
  #   geom_line() +
  #   scale_color_manual(
  #     name = "Vegetation Conditions",
  #     values = c(
  #       "Current Climate" = "gray30",
  #       "Less Sensitive" = blue,
  #       "More Sensitive" = orange 
  #     )
  #   ) +
  #   scale_fill_manual(
  #     name = "Vegetation Conditions",
  #     values = c(
  #       "Current Climate" = "gray30",
  #       "Less Sensitive" = blue,
  #       "More Sensitive" = orange 
  #     )
  #   ) +
  #   facet_wrap(~site_status)
  # 
  # fig <- 
  #   {fixed_climate_fig / fixed_veg_fig} +
  #   plot_layout(guides = "collect") &
  #   theme_minimal(base_size = 9) &
  #   theme(
  #     legend.position = "bottom"
  #   )

  ggsave(
    plot = fig,
    filename = output.file,
    ...)
  
  output.file
  
}

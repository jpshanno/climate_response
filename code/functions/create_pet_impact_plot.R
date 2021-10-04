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

  data[
    j = `:=`(
      y = abs(y),
      ymin = abs(ymax),
      ymax = abs(ymin)
    )
  ]
    
  # data <- data[
  #   .SDcols = c("y", "ymin", "ymax"),
  #   j = c(list(simulation_date = simulation_date),
  #         lapply(.SD, filled_rolling_mean)),
  #   by = .(site_status, Climate)]
  
  data <- data[
    month(simulation_date) %in% 4:10]
  
  fixed_climate_fig <- ggplot(data) +
    aes(
      x = simulation_date,
      y = y,
      ymin = ymin,
      ymax = ymax,
      fill = site_status,
      color = site_status
    ) +
    ylab("30-day Mean PET (cm)") +
    xlab(NULL) +
    ggtitle("B.") +
    geom_ribbon(alpha = 0.3, color = NA) +
    geom_line() +
    scale_color_manual(
      name = "Vegetation Conditions",
      values = c(
        "Modeled Black Ash" = brown,
        "Future Forested" = gold,
        "Non-Forested" = green 
      ),
      breaks = c(
        "Modeled Black Ash",
        "Non-Forested",
        "Future Forested"
      )
      ) +
    scale_fill_manual(
      name = "Vegetation Conditions",
      values = c(
        "Modeled Black Ash" = brown,
        "Future Forested" = gold,
        "Non-Forested" = green 
      ),
      breaks = c(
        "Modeled Black Ash",
        "Non-Forested",
        "Future Forested"
      )
    ) +
    facet_wrap(~Climate)
  
  fixed_veg_fig <- ggplot(data) +
    aes(
      x = simulation_date,
      y = y,
      ymin = ymin,
      ymax = ymax,
      fill = Climate,
      color = Climate
    ) +
    ylab("30-day Mean PET (cm)") +
    xlab(NULL) +
    ggtitle("B.") +
    geom_ribbon(alpha = 0.3, color = NA) +
    geom_line() +
    scale_color_manual(
      name = "Vegetation Conditions",
      values = c(
        "Current Climate" = "gray30",
        "Less Sensitive" = blue,
        "More Sensitive" = orange 
      )
    ) +
    scale_fill_manual(
      name = "Vegetation Conditions",
      values = c(
        "Current Climate" = "gray30",
        "Less Sensitive" = blue,
        "More Sensitive" = orange 
      )
    ) +
    facet_wrap(~site_status)
  
  fig <- 
    {fixed_climate_fig / fixed_veg_fig} +
    plot_layout(guides = "collect") &
    theme_minimal(base_size = 9) &
    theme(
      legend.position = "bottom"
    )

  ggsave(
    plot = fig,
    filename = output.file,
    ...)
  
  output.file
  
}

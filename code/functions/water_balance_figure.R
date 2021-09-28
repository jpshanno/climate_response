create_water_balance_plot <- function(
  data,
  output.file,
  ...
) {
  dat <- 
    data[j = .(water_availability_cm = mean(total_input_cm + pet_cm)),
         by = .(sample_date)
         ][
          j = .(sample_date,
                water_availability_cm = cumsum(water_availability_cm))]
  
  fig <- 
    ggplot(dat) +
    geom_line(aes(x = sample_date, y = water_availability_cm), color = "gray30") +
    labs(x = NULL, y = "Cumulative Water Availability (cm)") +
    scale_x_date(minor_breaks = "1 year") +
    theme_minimal(base_size = 12)
  
  ggsave(plot = fig, filename = output.file, ...)
  
  output.file
}


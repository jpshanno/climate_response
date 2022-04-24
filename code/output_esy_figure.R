source("code/load_project.R")
tar_load(control_optimization)
data <- tar_read(training_data)$control


drawdown <- 
  data[, 
       .SD[which(format(sample_date, "%m%d") == "0401"):which.min(wl_initial_cm)],
       by = .(site)]

drawdown[, ytd_water_balance := cumsum(rain_cm + pet_cm + melt_cm),
         by = .(site, water_year)]

drawdown[, 
         c("drawdown_emp", "esy_emp",
           "esy_x_intercept") := {
             
             mod <- 
               robustbase::nlrob(wl_initial_cm ~ quad(ytd_water_balance, 
                                                      b0 = b0, b1 = b1, b2 = b2),
                                 data = .SD, 
                                 na.action = na.exclude,
                                 maxit = 50,
                                 start = list(b0 = 8, b1 = 1, b2 = -1))
             
             list(drawdown_emp = predict(mod, newdata = .SD),
                  esy_emp = quad_prime(mod = mod, wa = .SD[, ytd_water_balance]),
                  esy_x_intercept = quad_x_intercept(mod))
             
           },
         by = .(site)]

esy_form <- 
  bf(esy_emp ~ a - (a - b) * exp (c * wl_initial_cm),
     a + b + c ~ 0 + site,
     nl = TRUE)

dat <- 
  drawdown[site == "135"]

pub_fig <-
  {ggplot(dat,
          aes(x = ytd_water_balance,
              y = wl_initial_cm)) + 
      geom_path() +
      geom_line(aes(y = drawdown_emp),
                color = 'red') + 
      scale_x_continuous(name = "Year-to-Date Water Balance (cm)",
                         position = "top") +
      labs(y = "Water Level Relative to\nGround Surface (cm)")} /
  {ggplot(dat,
          aes(x = wl_initial_cm,
              y = esy_emp)) + 
      geom_point() +
      geom_function(fun = ~control_optimization[site == "135", params[[1]]$funESY](.x, control_optimization[site == "135", params[[1]]$esy_min]),
                    color = 'red') + 
      labs(x = "Water Level Relative to Ground Surface (cm)",
           y = expression(Empiricial~E[Sy]))} * theme_bw(base_size = 16) +
  plot_annotation(tag_levels = "A")

ggsave(plot = pub_fig,
       filename = "output/figures/scatter_plot_and_line_empirical_esy.png",
       width = 12,
       height = 8,
       dpi = 300)


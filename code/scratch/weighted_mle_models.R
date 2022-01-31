source("code/load_project.R")
tar_load(training_data)
tar_load(testing_data)
tar_load(treatment_sites)
tar_load(esy_functions)

optimize_params <- 
  function(data, par, fixed = NULL, ...){
    opt <- optim(par = par,
                 fixed = fixed,
                 ...,
                 control = list(fnscale = 1,
                                maxit = 2000),
                 fn = 
                   function(params, fixed = NULL){

                    params <- 
                      c(params, fixed)

                    wl_hat <- 
                      wetland_model(data = data, params)$wl_hat

                    resids <-
                      (wl_hat - data$wl_initial_cm)

                    # Initial weights are equal
                    init_weight <- 1 / sum(is.na(data$wl_initial_cm))

                    # Weights increase asymmetrically as water levels drop
                    wghts <- pmax(init_weight, init_weight * (params$maxWL - data$wl_initial_cm))

                    # Weights are squared
                    wghts <- wghts^2
                    wghts[is.na(wghts)] <- init_weight^2

                    # Weight residuals
                    weighted_resids <- resids * wghts
                    weighted_resids <- weighted_resids[!is.na(weighted_resids)]
                    obs_se <-
                      sd(diff(data$wl_initial_cm), na.rm = TRUE) / sum(!is.na(data$wl_initial_cm))

                    -sum(
                      wghts[!is.na(resids)] * dnorm(resids[!is.na(resids)],
                            mean = 0,
                            sd = obs_se,
                            log = TRUE)
                      )
                   })

    opt$par <- c(opt$par, unlist(fixed))

    opt
  }

control_optimization <-
    fit_models(training_data[["control"]],
               esy_functions,
               par = list(MPET = 1,
                          MP = 1.5,
                          MM = 1,
                          MQ = 0.5,
                          minESY = 1,
                          phiM = 0.9,
                          phiP = 0.5))

training_data[["control"]][esy_functions, `:=`(max_wl = i.max_wl, funESY = i.pred_fun)]

control_population_optim <-
  optim(
    par = list(
      MPET = 1,
      MP = 1.5,
      MM = 1,
      MQ = 0.5,
      minESY = 1,
      phiM = 0.9,
      phiP = 0.5
    ),
    control = list(
      fnscale = 1,
      maxit = 2000
    ),
    fn =
      function(params) {

        split_data <- split(
          training_data[["control"]],
          f = training_data[["control"]]$site
        )

        sum(
          purrr::map_dbl(
          .x = split_data,
          .f = ~{
            site_params <- 
              c(params, list(maxWL = unique(.x$max_wl), funESY = .x$funESY[[1]]))

            wl_hat <-
              wetland_model(data = .x, site_params)$wl_hat

            resids <-
              (wl_hat - .x$wl_initial_cm)

            # Initial weights are equal
            init_weight <- 1 / sum(is.na(.x$wl_initial_cm))

            # Weights increase asymmetrically as water levels drop
            wghts <- pmax(init_weight, init_weight * (site_params$maxWL - .x$wl_initial_cm))

            # Weights are squared
            wghts <- wghts^2
            wghts[is.na(wghts)] <- init_weight^2

            obs_se <-
              sd(diff(.x$wl_initial_cm), na.rm = TRUE) / sum(!is.na(.x$wl_initial_cm))

            -sum(
              wghts[!is.na(resids)] * dnorm(resids[!is.na(resids)],
                mean = 0,
                sd = obs_se,
                log = TRUE
              )
              )
          }
        )
        )

      }
)

training_data[["control"]][,
  population_wl_hat :=
    wetland_model(
      .SD,
      as.list(c(control_population_optim$par, maxWL = max_wl[[1]], funESY = funESY[[1]]))
    ),
    by = .(site)
  ]

treated_optimization <-
    refit_model(training_data[["treated"]],
                control_optimization,
                refit = list(MPET = 1,
                             MP = 1))

model_params <-
    concatenate_model_params(Control = control_optimization, 
                             Treated = treated_optimization)

train_data_fits <- 
    predict_train_period(data = training_data, model_params)

test_data_fits <-
    predict_test_period(data = testing_data[site %in% treatment_sites], 
                        model.params = model_params)

wetland_model_metrics <-
    calculate_wetland_model_metrics(data = test_data_fits, max_wl_data = esy_functions)

create_wetland_model_metrics_table(wetland_model_metrics)

create_wetland_model_metrics_plot(
      data = wetland_model_metrics,
      metrics = c("r2", "med_err", "rmedse", "rmedse_range"),
      output.file = "/tmp/boxplot_wetland_model_metrics.pdf",
      width = 4,
      height = 6
      )

wetland_model_predicted_probabilities <-
    calculate_predicted_probabilities(test_data_fits, max_wl_data = esy_functions)

wetland_model_predicted_probs_tests <-
    test_predicted_probabilities(wetland_model_predicted_probabilities)

wetland_model_predicted_probs_table <-
    create_wetland_model_probability_table(wetland_model_predicted_probs_tests)

create_wetland_model_probability_plot(
      data = wetland_model_predicted_probabilities,
      tests = wetland_model_predicted_probs_tests,
      output.file = "/tmp/jitterpoint_and_pointrange_wetland_model_predicted_probability.pdf",
      width = 7.5,
      height = 3.5)


control_fits <- control_population_optim

# Population level MLE:
ggplot(training_data[["control"]]) + aes(x = as.dowy(sample_date, 11)) + geom_line(aes(y = wl_initial_cm)) + geom_line(aes(y = population_wl_hat), color = 'red', linetype = 'dashed') + facet_wrap(~site + water_year, scales = "free_y") + ggtitle("Population-Level MLE")

# MLE by site:
ggplot(train_data_fits[["control"]]) + aes(x = as.dowy(sample_date, 11)) + geom_line(aes(y = wl_initial_cm)) + geom_line(aes(y = wl_hat), color = 'red', linetype = 'dashed') + facet_wrap(~site + water_year, scales = "free_y") + ggtitle("Group-Level MLE")

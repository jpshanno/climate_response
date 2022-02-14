####### WORK WITH ONLY TREATED SITES TO AVOID UNEVEN GROUP SIZES

### Likelihood should probably go back to dnorm with obs sd
### Need to think about matching scale of ll for likelihood and prior
### Need to get proposal function better tuned

source("code/load_project.R")
tar_load(training_data)
tar_load(testing_data)
tar_load(treatment_sites)
tar_load(esy_functions)
training_data[["control"]][esy_functions, `:=`(max_wl = i.max_wl, funESY = i.pred_fun)]
con_dat <- training_data[["control"]][site %in% treatment_sites]

# https://kevintshoemaker.github.io/NRES-746/LECTURE7.html#Myxomatosis_revisited_(again)
ll_function <- function(params) {

        split_data <- split(
          con_dat,
          f = con_dat$site
        )

        mean(
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
            wghts <- wghts / sum(wghts, na.rm = TRUE)

            obs_sd <-
              sd(diff(.x$wl_initial_cm), na.rm = TRUE)

            sum(
              wghts[!is.na(resids)] *
              dt(
                resids[!is.na(resids)],
                df = sum(!is.na(.x$wl_initial_cm)),
                log = TRUE
              )
            )
          }
        )
        )
      }

start_values <- c(
  MPET = 1,
  MP = 1.5,
  MM = 1,
  MQ = 0.5,
  minESY = 1,
  phiM = 0.9,
  phiP = 0.5
)

prior_ll <- function(params) {
    mpet_prior <- dt(params[["MPET"]], df = sum(!is.na(con_dat$wl_initial_cm)), log = TRUE)
    mp_prior <- dt(params[["MP"]], df = sum(!is.na(con_dat$wl_initial_cm)), log = TRUE)
    mm_prior <- dt(params[["MM"]], df = sum(!is.na(con_dat$wl_initial_cm)), log = TRUE)
    mq_prior <- dt(params[["MQ"]], df = sum(!is.na(con_dat$wl_initial_cm)), log = TRUE)
    esy_prior <- dt(params[["minESY"]], df = sum(!is.na(con_dat$wl_initial_cm)), log = TRUE)
    phim_prior <- dt(params[["phiM"]], df = sum(!is.na(con_dat$wl_initial_cm)), log = TRUE)
    phip_prior <- dt(params[["phiP"]], df = sum(!is.na(con_dat$wl_initial_cm)), log = TRUE)

    length(unique(con_dat$site)) *
      (mpet_prior +
      mp_prior +
      mm_prior +
      mq_prior +
      esy_prior +
      phim_prior +
      phip_prior)
}

PosteriorRatio2 <- function(old_guess, new_guess) {
  oldLogLik <- ll_function(old_guess)   # compute likelihood and prior density at old guess
  oldLogPrior <- prior_ll(old_guess)
  newLogLik <- ll_function(new_guess)             # compute likelihood and prior density at new guess
  newLogPrior <- prior_ll(new_guess)
  return(exp((newLogLik+newLogPrior)-(oldLogLik+oldLogPrior)))          # compute ratio of weighted likelihoods
}

proposal_function <- function(oldguess) {
  # This function is built from Florian's post on implementing MCMC
  res <- oldguess + rep(-100, length(oldguess))

  # Ensure that we only return positive new guesses
  while(any(res < .Machine$double.neg.eps)) {
    res <- oldguess + rnorm(length(oldguess), mean = 0, sd = 0.01)
  }
  res
}

sample_mcmc <- function(chain_length, start_values, warmup) {
  old_guess <- start_values
  guesses <- matrix(0, nrow = chain_length, ncol = length(start_values))
  colnames(guesses) <- names(start_values)
  guesses[1, ] <- start_values
  counter <- 2
  while(counter <= chain_length) {
    new_guess <- proposal_function(old_guess)
    post.rat <- PosteriorRatio2(old_guess, new_guess)
    prob.accept <- min(1, post.rat)
    rand <- runif(1)
    if(rand <= prob.accept) {
      old_guess <- new_guess
      guesses[counter, ] <- new_guess
      counter = counter + 1
    }
  }
  guesses[(warmup+1):chain_length,]
}

chain <- sample_mcmc(10000, start_values, 1000)

colMeans(chain[-c(1:500),])
apply(chain[-c(1:500), ], 2, median)

for(i in seq_len(ncol(chain))) {
  plot(chain[, i], type = 'l', main = colnames(chain)[i])
}
for(i in seq_len(ncol(chain))) {
  plot(density(chain[, i]), main = colnames(chain)[i])
}






library(BayesianTools)
point_priors <- c(
  MPET = 1,
  MP = 1.5,
  MM = 1,
  MQ = 0.5,
  minESY = 1,
  phiM = 0.9,
  phiP = 0.5
)
b_setup <- createBayesianSetup(
  ll_function,
  prior = createUniformPrior(
    lower = c(
      MPET = 0, 
      MP = 0, 
      MM = 0, 
      MQ = 0, 
      minESY = 0, 
      phiM = 0, 
      phiP = 0
    ),
    upper = c(
      MPET = 10,
      MP = 10,
      MM = 10,
      MQ = 1-.Machine$double.neg.eps,
      minESY = 1,
      phiM = 1-.Machine$double.neg.eps,
      phiP = 1-.Machine$double.neg.eps
    ),
    best = point_priors
  ),
  names = names(point_priors)
)

runMCMC(b_setup, sampler = "Metropolis")



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

        # future::plan(future::multicore, workers = 4)
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
        # future::plan(future::sequential)

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

####### WORK WITH ONLY TREATED SITES TO AVOID UNEVEN GROUP SIZES

### Likelihood should probably go back to dnorm with obs sd
### Need to think about matching scale of ll for likelihood and prior
### Need to get proposal function better tuned

source("code/load_project.R")
tar_load(training_data)
tar_load(testing_data)
tar_load(treatment_sites)
tar_load(esy_functions)
esy_functions[, min_esy := pred_fun[[1]](max_wl, -9999), by = "site"]
training_data[["control"]][esy_functions, `:=`(max_wl = i.max_wl, funESY = i.pred_fun, minESY = i.min_esy)]
con_dat <- training_data[["control"]][site %in% treatment_sites]

# https://kevintshoemaker.github.io/NRES-746/LECTURE7.html#Myxomatosis_revisited_(again)
ll_function <- function(params) {

        split_data <- split(
          con_dat,
          f = con_dat$site
        )

        sum(
          purrr::map_dbl(
          .x = split_data,
          .f = ~{
            site_params <-
              c(params,
                list(
                  maxWL = unique(.x$max_wl),
                  funESY = .x$funESY[[1]]
                )
              )

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
              dnorm(
                resids[!is.na(resids)],
                mean = 0,
                sd = obs_sd,
                log = TRUE
              )
            )
          }
        )
        )
      }

prior_ll <- function(params, start_values) {
  mvtnorm::dmvnorm(
    x = params,
    mean = start_values,
    sigma = 0.5  * start_values * diag(length(params)),
    log = TRUE
  )
}

PosteriorRatio2 <- function(old_guess, new_guess, start_values) {
  oldLogLik <- ll_function(old_guess)   # compute likelihood and prior density at old guess
  oldLogPrior <- prior_ll(old_guess, start_values)
  newLogLik <- ll_function(new_guess)             # compute likelihood and prior density at new guess
  newLogPrior <- prior_ll(new_guess, start_values)
  exp((newLogLik+newLogPrior)-(oldLogLik+oldLogPrior))          # compute ratio of weighted likelihoods
}

proposal_function <- function(oldguess, sigma_scale) {
  # This function is built from Florian's post on implementing MCMC
  res <- oldguess + rep(-100, length(oldguess))

  # Ensure that we only return positive new guesses
  while(any(res < .Machine$double.neg.eps)) {
    jumps <- mvtnorm::rmvnorm(
      n = 1,
      mean = rep(0, length(oldguess)),
      sigma = sigma_scale * diag(length(oldguess))
    )
    res <- oldguess + jumps
  }
  res <- res[1,]
  names(res) <- names(oldguess)
  res
}

sample_mcmc <- function(chain_length, start_values, warmup, adapt_rate) {

  # Set-up
  guesses <- matrix(0, nrow = chain_length, ncol = length(start_values))
  sigma_scale <- acceptance_rate <- numeric(length(chain_length))
  colnames(guesses) <- names(start_values)

  # Define Initial Conditions
  old_guess <- start_values
  guesses[1, ] <- start_values
  counter <- 2
  pb <- txtProgressBar(max = chain_length)
  proposals <- 1
  sigma_scale[c(1, 2)] <- c(NA_real_, 0.0001)
  acceptance_rate[1] <- NA_real_

  # Loop through chain length saving accepted proposals
  while(counter <= chain_length) {
    # Get a new guess from a proposed jump
    new_guess <- proposal_function(old_guess, sigma_scale[counter])
    # Get the ration between liklihood and prior likelihood
    post.rat <- PosteriorRatio2(old_guess, new_guess, start_values)
    # Get the probability of acceptance (max of 1)
    prob.accept <- min(1, post.rat)
    # Get a random probability between 0 and 1
    rand <- runif(1)
    # Increment number of proposals
    proposals <- proposals + 1
    # Check if probability of acceptance is greater than the random number
    if(rand <= prob.accept) {
      # Store new guess as old guess
      old_guess <- new_guess
      # Store new guess in accepted guesses (draws)
      guesses[counter, ] <- new_guess
      # Get acceptance rate of proposals
      acceptance_rate[counter] <- (counter - 1) / proposals
      # Tune towards 30% acceptance rate by adjusting jump size
      sigma_scale[counter + 1] <- dplyr::case_when(
        acceptance_rate[counter] > 0.4 ~ sigma_scale[counter] * (1 + adapt_rate),
        acceptance_rate[counter] < 0.2 ~ sigma_scale[counter] * (1 - adapt_rate),
        TRUE ~ sigma_scale[counter]
      )
      # Increment counter
      counter <- counter + 1
      setTxtProgressBar(pb, counter)
    }
  }
  # Return results
  list(
    warmup = ifelse(warmup == 0, numeric(), guesses[1:warmup,]),
    draws = guesses[(warmup+1):chain_length,],
    acceptance_rate = acceptance_rate,
    sigma_scale = sigma_scale
  )
}

init_values <- c(
  MPET = 1,
  MP = 1.5,
  MM = 1,
  MQ = 0.5,
  minESY = 1,
  phiM = 0.9,
  phiP = 0.5
)

chain <- sample_mcmc(200, init_values, 0, 0.001)

colMeans(chain$draws)
apply(chain$draws, 2, median)

par(mfrow = c(4, 2))
for(i in seq_len(ncol(chain$draws))) {
  plot(chain$draws[, i], type = 'l', main = colnames(chain$draws)[i])
}

for(i in seq_len(ncol(chain$draws))) {
  plot(density(chain$draws[, i]), main = colnames(chain$draws)[i])
}
par(mfrow = c(1, 1))
pairs(chain$draws)

plot(chain$acceptance_rate, type = 'l')
plot(chain$sigma_scale, type = 'l')

future::plan(future::multicore, workers = 3)
chains <- furrr::future_map(
  list(
    chain_1 = init_values, # set_names(runif(length(start_values)), names(start_values)),
    chain_2 = init_values, # set_names(runif(length(start_values)), names(start_values)),
    chain_3 = init_values # set_names(runif(length(start_values)), names(start_values))
  ),
  ~ sample_mcmc(
      chain_length = 500,
      start_values = .x,
      warmup = 100,
      adapt_rate = 0.001
    ),
  .options = furrr::furrr_options(seed = TRUE)
)
future::plan(future::sequential)

saveRDS(chains, "~/phd/climate_response/tmp/chains.rds")

map(chains, ~colMeans(.x$draws))
map(chains, ~apply(.x$draws, 2, median))
draws <- map_dfr(chains, ~as.data.table(.x$draws)[, .draw := 1:nrow(.SD)], .id = "chain") %>%
  tidyr::pivot_longer(cols = -c(chain, .draw)) %>%
  as.data.table()

ggplot(draws) +
  aes(x = .draw, y = value, color = chain) +
  geom_line() +
  facet_wrap(~name, scales = "free")

# 7998 is the largest number less than 8000 and divisible by 3
sample_frame <-
  data.table(
    chain = rep(c("chain_1", "chain_2", "chain_3"), each = 7998/3),
    .draw = sample(8000, size = 7998)
  )

thinned <- draws[sample_frame, on = c("chain", ".draw")]
ggplot(thinned) +
  aes(x = .draw, y = value, color = chain) +
  geom_line() +
  facet_wrap(~name, scales = "free")



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

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
con_dat <- con_dat[(!(month(sample_date) == 2 & mday(sample_date) == 29)),]

site <- "151"
dat <- con_dat[site == "151"]

# Weights increase asymmetrically as water levels drop
create_weights <- function(x, scale) {
  init_weight <- 1 / sum(!is.na(x))
  wghts <- pmax(init_weight, init_weight * (scale - x))
  # Weights are squared
  wghts <- wghts^2
  wghts[is.na(x)] <- 0
  wghts <- wghts / sum(wghts)
}

create_data_matrix <- function(data, var) {
  lapply(
    X = sort(unique(data[["site"]])),
    FUN = \(x) matrix(data[site == x][[var]], nrow = 1)
  ) %>%
  purrr::reduce(rbind)
}

con_dat[, wghts := create_weights(wl_initial_cm, max_wl), by = .(site)]
con_dat[, filled_wl := nafill(wl_initial_cm ,"const", -9999)]
esy_functions[sort(unique(con_dat$site))]$pred_fun
esy_params <- matrix(
  c(
    0.9440797, 9.68796681237216, 0.963832093292945, 0.011490387703557,
    1.2087133, 9.6930777345254, 2.10532835889272, 0.0157893963100476,
    1.1956961, 9.35343568266112, 1.62678535681258, 0.00984811084065675,
    1.1942750, 8.85796577611301, 2.0742188380092, 0.0148425995465377,
    1.3711603, 9.79628207399244, 1.31799862570157, 0.00835140512795191,
    1.4753461, 9.92020216902697, 2.40233114598982, 0.0103268508715977,
    1.0764391, 9.56079966326987, 1.17368355521875, 0.00762101076778987,
    0.2738756, 9.64372934079638, 2.0767145808409, 0.0127361070535151),
  byrow = TRUE,
  ncol = 4
)

obs_sds <- con_dat[, .(v = sd(diff(wl_initial_cm), na.rm = TRUE)), by = .(site)][["v"]]

stan_data <- list(
  D = 365L,
  K = length(unique(con_dat$site)),
  pet = create_data_matrix(con_dat, "pet_cm"),
  rain = create_data_matrix(con_dat, "rain_cm"),
  melt = create_data_matrix(con_dat, "melt_cm"),
  wghts = create_data_matrix(con_dat, "wghts"),
  maxWL = create_data_matrix(con_dat, "max_wl")[, 1],
  y = create_data_matrix(con_dat, "filled_wl"),
  esyParams = esy_params,
  ySD = obs_sds
)


library(cmdstanr)
library(posterior)
library(bayesplot)
mod <- cmdstan_model(stan_file = "code/scratch/wetland_model.stan")
# On the last fit with all 4 parameters it look like 1 chain went bad when bQ
# was included. Probably have to fine tune the priors or the specification in
# the model
fit <- mod$sample(
  data = stan_data,
  seed = 1234567,
  chains = 4,
  parallel_chains = 4,
  refresh = 500,
  adapt_delta = 0.9,
  max_treedepth = 10
)

fit$summary()
mcmc_hist(fit$draws(c("bPET", "bRain", "bMelt", "bQ")))

fit_mcmc <- as_mcmc.list(fit)

# color_scheme_set("mix-blue-pink")
mcmc_trace(fit_mcmc,  pars = c("bPET", "bRain", "bMelt", "bQ"), n_warmup = 1000,
                facet_args = list(nrow = 2, labeller = label_parsed))


new_params <- list(
  MPET = 1.36,
  MP = 1.48,
  MM = 1.04,
  MQ = 0.516,
  minESY = esy_functions[site][["min_esy"]],
  phiM = 0,
  phiP =0,
  maxWL = unique(dat$max_wl),
  funESY = esy_functions[site, pred_fun]
)
test <- wetland_model(dat, new_params)
test
plot(dat$wl_initial_cm, type = "l")
lines(test$wl_hat, col = 'red', lty = "dashed")

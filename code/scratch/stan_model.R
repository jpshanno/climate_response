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

init_weight <- 1 / 365

# Weights increase asymmetrically as water levels drop
# create_weights <- function()
wghts <- pmax(init_weight, init_weight * (dat$max_wl - dat$wl_initial_cm))

# Weights are squared
wghts <- wghts^2
wghts[is.na(dat$wl_initial)] <- 0

stan_data <- list(
  D = 365L,
  K = 1L,
  pet = matrix(dat$pet_cm, nrow = 1),
  rain = matrix(dat$rain_cm, nrow = 1),
  melt = matrix(dat$melt_cm, nrow = 1),
  wghts = matrix(wghts, nrow = 1),
  maxWL = as.array(unique(dat$max_wl)), # wrapping in as.array found as recommended fix for dims 'declared=(1); dims found=()' error
  y = matrix(nafill(dat$wl_initial_cm, "const", -9999), nrow = 1),
  esyParams = matrix(c(1.0764391, 9.56079966326987, 1.17368355521875, 0.00762101076778987), nrow = 1),
  ySD = as.array(sd(diff(dat$wl_initial_cm), na.rm = TRUE))
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

####### WORK WITH ONLY TREATED SITES TO AVOID UNEVEN GROUP SIZES

# It looks like the step size is being set much smaller on the chains that end
# up not mixing. I think the issue is that the warm-up period is leading to
# poor step sizes for some of the chains. Need to determine how to either set
# some minimum, or explore how init values relate to final step size. It may be
# that init values to close to the final parameters are leading to small step
# sizes? I'd have to check to be sure

# Getting multi-modal estimates of invididual sigma values. Likely have to
# contrain sigma (not desirable) or update the weights (better path). It's still
# too high of a likelihood to fit a mean model and just have terrible errors any
# time the water level is not hovering around max

source("code/load_project.R")
tar_load(training_data)
tar_load(testing_data)
tar_load(treatment_sites)
tar_load(esy_functions)
esy_functions[, min_esy := pred_fun[[1]](max_wl, -9999), by = "site"]
training_data[["control"]][esy_functions, `:=`(max_wl = i.max_wl, funESY = i.pred_fun, minESY = i.min_esy)]
con_dat <- training_data[["control"]][site %in% treatment_sites]
con_dat <- con_dat[(!(month(sample_date) == 2 & mday(sample_date) == 29)),]

library(cmdstanr)
library(posterior)
library(bayesplot)

# Weights increase asymmetrically as water levels drop
create_weights <- function(x, scale) {
  init_weight <- 1 / 365
  wghts <- pmax(init_weight, init_weight * (scale - x))
  # Weights are squared
  wghts <- wghts^2
  wghts[is.na(x)] <- 0
  wghts <- wghts / sum(wghts)
  assertthat::assert_that(assertthat::noNA(wghts))
  wghts
}

create_data_matrix <- function(data, var) {
  res <- lapply(
    X = sort(unique(data[["site"]])),
    FUN = \(x) matrix(data[site == x][[var]], nrow = 1)
  ) %>%
  purrr::reduce(rbind)

  assertthat::assert_that(assertthat::noNA(res))
  res
}

init_values <- function(pooling) {
  set.seed(8675309)
  generate_values <- function(pooling) {
    x <- list(
      bPET = runif(1, 0.9, 1.1),
      bRain = runif(1, 0.9, 1.1),
      bMelt = runif(1, 0.9, 1.1),
      bQ = runif(1, 0.4, 0.7)
    )
    # x$sigma <- switch(pooling,
    #   "full" = 0.1,
    #   "partial" = rep(0.01, 8),
    #   "none" = rep(0.01, 8)
    # )
    x
  }
  replicate(generate_values(pooling), n = 4, simplify = FALSE)
}

con_dat[, wghts := create_weights(wl_initial_cm, max_wl), by = .(site)]
con_dat[, filled_wl := nafill(wl_initial_cm ,"const", -9999)]
# esy_functions[sort(unique(con_dat$site))]$pred_fun
# esy_functions[sort(unique(con_dat$site))]$min_esy
esy_params <- matrix(
  c(
      0.9440797, 9.68796681237216, 0.963832093292945, 0.011490387703557
    , 1.2087133, 9.6930777345254, 2.10532835889272, 0.0157893963100476
    , 1.1956961, 9.35343568266112, 1.62678535681258, 0.00984811084065675
    , 1.1942750, 8.85796577611301, 2.0742188380092, 0.0148425995465377
    , 1.3711603, 9.79628207399244, 1.31799862570157, 0.00835140512795191
    , 1.4753461, 9.92020216902697, 2.40233114598982, 0.0103268508715977
    , 1.0764391, 9.56079966326987, 1.17368355521875, 0.00762101076778987
    , 0.2738756, 9.64372934079638, 2.0767145808409, 0.0127361070535151
  ),
  byrow = TRUE,
  ncol = 4
)

obs_sigma <- con_dat[order(site), .(v = median(abs(diff(wl_initial_cm)), na.rm = TRUE)), by = .(site)][["v"]]

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
  obs_sigma = obs_sigma
)

mod <- cmdstan_model(
  stan_file = "code/scratch/wetland_model.stan",
  dir = "/tmp",
  include_paths = file.path(getwd(), "code/scratch/"),
  force_recompile = TRUE
  )

fit <- mod$sample(
  data = stan_data,
  seed = 1234567,
  chains = 4,
  parallel_chains = 4,
  adapt_delta = 0.70,
  iter_warmup = 500,
  iter_sampling = 100,
  refresh = 50,
  save_warmup = TRUE,
  init = 0
)

# saveRDS(fit, "/Users/jpshanno/phd/climate_response/tmp/fit_20220304.rds")
# list.files(
#     "/var/folders/h_/8p62bvmn4b73006tnwkgfrlc0000gn/T/",
#     recursive = TRUE,
#     pattern = "wetland_model-202203032325",
#     full.names=TRUE
#   ) %>%
#   file.copy(file.path("/Users/jpshanno/phd/climate_response/tmp", basename(.)))

fit$summary()
mcmc_hist(fit$draws(c("bPET", "bRain", "bMelt", "bQ", "sigma"), inc_warmup = TRUE))

fit_mcmc <- as_mcmc.list(fit)
# color_scheme_set("mix-blue-pink")
mcmc_trace(
  fit_mcmc,
  pars = c("bPET", "bRain", "bMelt", "bQ"),
  n_warmup = 500,
  facet_args = list(nrow = 2, labeller = label_parsed)
  )

draws <- map(
  seq_len(4),
  ~list.files(
    "/var/folders/h_/8p62bvmn4b73006tnwkgfrlc0000gn/T/",
    recursive = TRUE,
    pattern = "wetland_model-202203072326",
    full.names=TRUE
  )[.x] %>%
  readr::read_csv(skip = 45) %>%
  dplyr::filter(dplyr::if_all(.fn = ~!is.na(.x)))
)

par(mfrow = c(2,2)); iwalk(draws, ~plot(log10(.x$stepsize__), type = "l", main = .y)); par(mfrow = c(1,1))

test_site <- "009"
dat <- con_dat[site == test_site]
new_params <- list(
  MPET = 0.958,
  MP = 0.854,
  MM = 0.832,
  MQ = 0.644,
  minESY = esy_functions[test_site][["min_esy"]],
  phiM = 0,
  phiP =0,
  maxWL = unique(dat$max_wl),
  funESY = esy_functions[test_site, pred_fun]
)
test <- wetland_model(dat, new_params)
plot(dat$wl_initial_cm, type = "l")
lines(test$wl_hat, col = 'red', lty = "dashed")


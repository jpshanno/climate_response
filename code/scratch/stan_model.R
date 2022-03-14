####### WORK WITH ONLY TREATED SITES TO AVOID UNEVEN GROUP SIZES

# TODO: adjust weights to weight more heavily for days above max_wl
# TODO: set lower values to 1 not 0

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

generate_init_params <- function(){
     c(
        runif(1, 0.9, 1.1),
        runif(1, 0.9, 1.1),
        runif(1, 0.9, 1.1),
        runif(1, 0.4, 0.7)
      )
  }

generate_values <- function() {
    list(
      bPop = generate_init_params(),
      bGroup = t(replicate(n = 8, generate_init_params())),
      sigma = 1,
      tau = runif(4, 0.1, 0.3)
    )
  }

init_values <- function(chain_id) {
  set.seed(8675309 + chain_id)
  generate_values()
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

# Priors -----------------------------------------------------------------------
dist <- distributional::generate(distributional::dist_gamma(10, 20), 10000)[[1]]
plot(density(dist))
ggdist::mean_hdci(dist, .width = 0.9)
ggdist::median_hdci(dist, .width = 0.9)


gdist <- distributional::generate(distributional::dist_normal(0, 0.05), 10000)[[1]]
plot(density(gdist))
ggdist::mean_hdci(gdist, .width = 0.9)
ggdist::median_hdci(gdist, .width = 0.9)


mod <- cmdstan_model(
  stan_file = "code/scratch/wetland_model.stan",
  dir = "/tmp",
  force_recompile = TRUE
  )

fit <- mod$sample(
  data = stan_data,
  seed = 1234567,
  chains = 4,
  parallel_chains = 4,
  adapt_delta = 0.8,
  # max_treedepth = 11,
  iter_warmup = 1000,
  iter_sampling = 1000,
  refresh = 50,
  save_warmup = TRUE,
  init = 0
)
  
# saveRDS(fit, "/Users/jpshanno/phd/climate_response/tmp/fit_20220304.rds")
# list.files(
#     "/var/folders/h_/8p62bvmn4b73006tnwkgfrlc0000gn/T/",
#     recursive = TRUE,
#     pattern = "wetland_model-202203132138",
#     full.names=TRUE
#   ) %>%
#   file.copy(file.path("/Users/jpshanno/phd/climate_response/tmp", basename(.)))

fit$summary() %>% dplyr::select(variable, mean, median, rhat) %>% print(n = nrow(.))
mcmc_hist(fit$draws(c("bGroup[3,1]", "bGroup[3,2]", "bGroup[3,3]", "bGroup[3,4]")), inc_warmup = TRUE)

fit_mcmc <- as_mcmc.list(fit)
color_scheme_set("mix-blue-pink")
mcmc_trace(
  fit_mcmc,
  pars = glue::glue("bPop[{p}]", p = 1:4)
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


plot_test_site <- function(site_id) {
  idx <- which(c("009", "053", "077", "119", "139", "140", "151", "156") == site_id)
  test_params <- fit$summary() %>% dplyr::filter(grepl(glue::glue("\\[{idx},"), .$variable)) %>% dplyr::pull(median)
  print(test_params)
  test_site <- unique(con_dat$site)[idx]
  dat <- testing_data[site == test_site][water_year == min(water_year)]
  new_params <- list(
    MPET = test_params[1],
    MP = test_params[2],
    MM = test_params[3],
    MQ = test_params[4],
    minESY = esy_functions[test_site][["min_esy"]],
    phiM = 0,
    phiP =0,
    maxWL = esy_functions[test_site, max_wl],
    funESY = esy_functions[test_site, pred_fun]
  )
  test <- wetland_model(dat, new_params)
  plot(dat$wl_initial_cm, type = "l")
  lines(test$wl_hat, col = 'red', lty = "dashed")
}
plot_test_site("151")
plot_test_site("077")

plot_train_site <- function(site_id) {
  idx <- which(c("009", "053", "077", "119", "139", "140", "151", "156") == site_id)
  test_params <- fit$summary() %>% dplyr::filter(grepl(glue::glue("\\[{idx},"), .$variable)) %>% dplyr::pull(median)
  print(test_params)
  test_site <- unique(con_dat$site)[idx]
  dat <- con_dat[site == test_site][water_year == min(water_year)]
  new_params <- list(
    MPET = test_params[1],
    MP = test_params[2],
    MM = test_params[3],
    MQ = test_params[4],
    minESY = esy_functions[test_site][["min_esy"]],
    phiM = 0,
    phiP =0,
    maxWL = esy_functions[test_site, max_wl],
    funESY = esy_functions[test_site, pred_fun]
  )
  test <- wetland_model(dat, new_params)
  plot(dat$wl_initial_cm, type = "l")
  lines(test$wl_hat, col = 'red', lty = "dashed")
}
plot_train_site("151")

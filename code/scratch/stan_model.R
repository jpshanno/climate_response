####### WORK WITH ONLY TREATED SITES TO AVOID UNEVEN GROUP SIZES


# TODO: Loosen priors
# TODO: get init working
# TODO: Consider heavier weighting (cubic)
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
  init_weight <- 1 / !is.na(x)
  wghts <- pmax(init_weight, init_weight * abs(scale - x))
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

generate_init_params <- function() {
  matrix(
    c(
      runif(3, 0.9, 1.1),
      runif(1, 0.4, 0.7)
    ),
    nrow = 1
  )
  }

generate_values <- function() {
    init_params <- generate_init_params()
    init_tau <- matrix(runif(4, 0.01, 0.15), ncol = 4)
    list(
      bPop = init_params,
      bGroup = replicate(8, generate_init_params()),
      tau = init_tau,
      sigma = 1
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
      0.9440797, 9.68796681237216, 0.963832093292945, 0.011490387703557, 1.78, -0.0660
    , 1.2087133, 9.6930777345254, 2.10532835889272, 0.0157893963100476, 0.653, -0.0915
    , 1.1956961, 9.35343568266112, 1.62678535681258, 0.00984811084065675, 2.29, -0.0724
    , 1.1942750, 8.85796577611301, 2.0742188380092, 0.0148425995465377, 2.32, -0.0427
    , 1.3711603, 9.79628207399244, 1.31799862570157, 0.00835140512795191, 1.93, -0.0635
    , 1.4753461, 9.92020216902697, 2.40233114598982, 0.0103268508715977, 2.32, -0.0518
    , 1.0764391, 9.56079966326987, 1.17368355521875, 0.00762101076778987, 2.04, -0.0871
    , 0.2738756, 9.64372934079638, 2.0767145808409, 0.0127361070535151, 3.03, -0.0720
  ),
  byrow = TRUE,
  ncol = 6
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
dist <- distributional::generate(distributional::dist_gamma(1, 4), 10000)[[1]]
ex_dist <- distributional::generate(distributional::dist_exponential(5), 10000)[[1]]
ggdist::mean_hdci(dist, .width = 0.9)
ggdist::median_hdci(dist, .width = 0.9)
plot(density(dist))
lines(density(ex_dist), col = "red")


gdist <- distributional::generate(distributional::dist_normal(0, 0.05), 10000)[[1]]
plot(density(gdist))
ggdist::mean_hdci(gdist, .width = 0.9)
ggdist::median_hdci(gdist, .width = 0.9)


mod <- cmdstan_model(
  stan_file = "code/scratch/wetland_model.stan",
  dir = "/tmp",
  force_recompile = TRUE
  )

# Tightening priors and increasing adapt_delta reduces divergences (see last
# commit)
fit <- mod$sample(
  data = stan_data,
  seed = 1234567,
  chains = 4,
  parallel_chains = 4,
  adapt_delta = 0.99,
  # max_treedepth = 11,
  iter_warmup = 200,
  iter_sampling = 200,
  refresh = 50,
  save_warmup = TRUE,
  init = 0
)

fit$summary() %>% dplyr::select(variable, mean, median, rhat) %>% print(n = nrow(.))
mcmc_dens(fit$draws(c("bPop[1]", "bPop[2]", "bPop[3]"))) +
  geom_density(
    color = 'red',
    linetype = 'dashed',
    aes(x = distributional::generate(distributional::dist_gamma(10, 10), 12000)[[1]]))

mcmc_dens(fit$draws(c("bPop[4]"))) +
  geom_density(
    color = 'red',
    linetype = 'dashed',
    aes(x = distributional::generate(distributional::dist_gamma(1, 4), 4000)[[1]]))

mcmc_dens(fit$draws(c("bPop[5]", "bPop[6]"))) +
  geom_density(
    color = 'red',
    linetype = 'dashed',
    aes(x = distributional::generate(distributional::dist_exponential(10), 8000)[[1]]))

# draws <- map(
#   seq_len(4),
#   ~list.files(
#     "/var/folders/h_/8p62bvmn4b73006tnwkgfrlc0000gn/T/",
#     recursive = TRUE,
#     pattern = "wetland_model-202203072326",
#     full.names=TRUE
#   )[.x] %>%
#   readr::read_csv(skip = 45) %>%
#   dplyr::filter(dplyr::if_all(.fn = ~!is.na(.x)))
# )

# par(mfrow = c(2,2)); iwalk(draws, ~plot(log10(.x$stepsize__), type = "l", main = .y)); par(mfrow = c(1,1))


plot_test_site <- function(site_id, pop_level = FALSE) {
  idx <- which(c("009", "053", "077", "119", "139", "140", "151", "156") == site_id)
  if(pop_level) {
    test_params <- fit$summary() %>% dplyr::filter(grepl("^[bp]", .$variable)) %>% dplyr::pull(mean)
  } else {
    test_params <- fit$summary() %>% dplyr::filter(grepl(glue::glue("\\[{idx}\\]"), .$variable)) %>% dplyr::pull(mean)
  }
  print(c(site_id, test_params))
  dat <- testing_data[site == site_id & site_status == "Control"][, c(.SD, list(n_records = !is.na(wl_initial_cm))), by = .(water_year)][n_records > 0]
  if(nrow(dat) == 0) return(NULL)
  dat <- dat[water_year == sample(water_year, 1)]
  new_params <- list(
    MPET = test_params[1],
    MP = test_params[2],
    MM = test_params[3],
    MQ = test_params[4],
    minESY = esy_functions[site_id][["min_esy"]],
    phiP = test_params[5],
    phiM =test_params[6],
    maxWL = esy_functions[site_id, max_wl],
    funESY = esy_functions[site_id, pred_fun]
  )
  test <- wetland_model(dat, new_params)
  ylim <- c(min(c(test$wl_hat, dat$wl_initial_cm)), max(c(test$wl_hat, dat$wl_initial_cm)))
  plot(dat$wl_initial_cm, type = "l", main = site_id, ylim = ylim)
  lines(test$wl_hat, col = 'red', lty = "dashed")
}
sites_to_test <- intersect(testing_data[site_status == "Control"]$site, con_dat$site)
par(mfrow = c(2, 4)); purrr::walk(sites_to_test, ~plot_test_site(.x)); par(mfrow = c(1,1))
par(mfrow = c(2, 4)); purrr::walk(sites_to_test, ~plot_test_site(.x, TRUE)); par(mfrow = c(1,1))

plot_train_site <- function(site_id, pop_level = FALSE) {
  idx <- which(c("009", "053", "077", "119", "139", "140", "151", "156") == site_id)
  if(pop_level) {
    test_params <- fit$summary() %>% dplyr::filter(grepl("^[bp]", .$variable)) %>% dplyr::pull(mean)
  } else {
    test_params <- fit$summary() %>% dplyr::filter(grepl(glue::glue("\\[{idx}\\]"), .$variable)) %>% dplyr::pull(mean)
  }
  test_site <- unique(con_dat$site)[idx]
  dat <- con_dat[site == test_site][water_year == min(water_year)]
  new_params <- list(
    MPET = test_params[1],
    MP = test_params[2],
    MM = test_params[3],
    MQ = test_params[4],
    minESY = esy_functions[test_site][["min_esy"]],
    phiP = test_params[5],
    phiM = test_params[6],
    maxWL = esy_functions[test_site, max_wl],
    funESY = esy_functions[test_site, pred_fun]
  )
  test <- wetland_model(dat, new_params)
  ylim <- c(
    min(c(test$wl_hat, dat$wl_initial_cm), na.rm = TRUE),
    max(c(test$wl_hat, dat$wl_initial_cm), na.rm = TRUE)
  )
  plot(dat$wl_initial_cm, type = "l", ylim = ylim)
  lines(test$wl_hat, col = 'red', lty = "dashed")
}
par(mfrow = c(2, 4)); purrr::walk(unique(con_dat$site), ~plot_train_site(.x)); par(mfrow = c(1,1))
par(mfrow = c(2, 4)); purrr::walk(unique(con_dat$site), ~plot_train_site(.x, TRUE)); par(mfrow = c(1,1))

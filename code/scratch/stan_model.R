####### WORK WITH ONLY TREATED SITES TO AVOID UNEVEN GROUP SIZES


# TODO: Limit rain inputs to only periods when rain is greater than PET
# TODO: compare symmetric and asymmetric weights
# TODO: limit number of params by setting melt params equal to rain params
# TODO: Try same starting params for all chains
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
  # wghts <- wghts / sum(wghts)
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

generate_values <- function() {
    list(
      bPET = runif(1, 0.9, 1.1),
      bRain = runif(1, 0.9, 1.1),
      # bMelt = runif(1, 0.9, 1.1),
      bQ = runif(1, 0.4, 0.6),
      # bphiRain = runif(1, 0.1, 0.2),
      # bphiMelt = runif(1, 0.1, 0.2),
      bEsyInt = runif(1, 1, 2),
      bEsySlope = runif(1, 0, 0.1),
      bEsyMin = runif(1, 0, 1),
      # taubPET = runif(1, 0.01, 0.03),
      # taubRain = runif(1, 0.01, 0.03),
      # taubMelt = runif(1, 0.01, 0.03),
      # taubQ = runif(1, 0.01, 0.03),
      # tauphiRain = runif(1, 0.01, 0.03),
      # tauphiMelt = runif(1, 0.01, 0.03),
      # taubEsyInt = runif(1, 0.01, 0.03),
      # taubEsySlope = runif(1, 0.01, 0.03),
      # taubEsyMin = runif(1, 0.01, 0.03),
      # gPET = runif(8, 0.9, 1.1),
      # gRain = runif(8, 0.9, 1.1),
      # gMelt = runif(8, 0.9, 1.1),
      # gQ = runif(8, 0.9, 1.1),
      # gphiRain = runif(8, 0.9, 1.1),
      # gphiMelt = runif(8, 0.9, 1.1),
      # gEsyInt = runif(8, 1, 2),
      # gEsySlope = runif(8, 0, 0.1),
      # gEsyMin = runif(8, 0, 1),
      omega = runif(1, 0, 1),
      alpha = runif(1, -1, 0)
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
      0.9440797, 9.68796681237216, 0.963832093292945, 0.011490387703557, 1.78, 0.0660
    , 1.2087133, 9.6930777345254, 2.10532835889272, 0.0157893963100476, 0.653, 0.0915
    , 1.1956961, 9.35343568266112, 1.62678535681258, 0.00984811084065675, 2.29, 0.0724
    , 1.1942750, 8.85796577611301, 2.0742188380092, 0.0148425995465377, 2.32, 0.0427
    , 1.3711603, 9.79628207399244, 1.31799862570157, 0.00835140512795191, 1.93, 0.0635
    , 1.4753461, 9.92020216902697, 2.40233114598982, 0.0103268508715977, 2.32, 0.0518
    , 1.0764391, 9.56079966326987, 1.17368355521875, 0.00762101076778987, 2.04, 0.0871
    , 0.2738756, 9.64372934079638, 2.0767145808409, 0.0127361070535151, 3.03, 0.0720
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
  y = create_data_matrix(con_dat, "filled_wl"),
  esyParams = esy_params,
  maxWL = create_data_matrix(con_dat, "max_wl")[, 1]
)

wl_model <- function(
    # bPET,
    # bRain,
    # bMelt,
    # bQ,
    # phiRain,
    # phiMelt,
    pet,
    rain,
    melt,
    D,
    maxWL,
    esyint,
    esyslope,
    esymin
){
  wlHat <- numeric(D)
  Rain <- numeric(D)
  Melt <- numeric(D)
  gradient <- numeric(D)
  qHat <- numeric(D)
  mHat <- numeric(D)
  pHat <- numeric(D)
  petHat <- numeric(D)
      
      # Initialize model at full water level
      wlHat[1] <- maxWL
      
      # Loop through weather data
      for(t in 2:D){
            # Rain[t] <- rain[t] + rain[t-1] * phiRain
            # Melt[t] <- melt[t] + melt[t-1] * phiMelt
            wlHat[t] <- wlHat[t - 1]
            # Esy
            # Calculate gradient of drawdown 
            gradient[t] <- pmax(esyint - esyslope * wlHat[t], esymin)
            
            # PET or P times Esy
            # Use net input to determine if water level increases or decreases
            # Assuming AET is negligible on days where P >= PET
            # Water level drawdown = PET2, if P2 <= PET2 then it can be assumed to be
            # less than interception (not necessarily true, but works as a
            # simplifying assumption)
            petHat[t] <- pet[t] * gradient[t]
                  # pet_fun(wl_hat[t-1], maxWL, future.forest.change) * (MPET * PET[t]) * gradient[t]
            wlHat[t] <- wlHat[t] + petHat[t]
            pHat[t] <- rain[t] * gradient[t]
            wlHat[t] <- wlHat[t] + pHat[t]
            

            # Snowmelt
            mHat[t] <- Melt[t] * gradient[t]
            wlHat[t] <- wlHat[t] + mHat[t]

            # Streamflow
            # If WL is above spill point threshold then lose some to streamflow. 
            # This could probably be improved using the morphology models to determine
            # streamflow
            # if(wlHat[t-1] > maxWL){
            # qHat[t] <- bQ * (wlHat[t-1] - maxWL)
            # if(qHat[t] > wlHat[t] - maxWL){
            #       qHat[t] <- wlHat[t] - maxWL
            # }
            # wlHat[t] <- wlHat[t] - qHat[t]
            # }

      }
    wlHat
   }

# Priors -----------------------------------------------------------------------
dist <- distributional::generate(distributional::dist_gamma(1, 1), 10000)[[1]]
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
# vals <- lapply(1:4, \(x){init_values(x)})
fit <- mod$sample(
  data = stan_data,
  seed = 1234567,
  parallel_chains = 4,
  adapt_delta = 0.99,
  iter_warmup = 1000,
  iter_sampling = 200,
  refresh = 250,
  save_warmup = TRUE,
  init = init_values
)

fit$summary()
mcmc_trace(fit$draws(inc_warmup = FALSE))
mcmc_dens(fit$draws(inc_warmup = FALSE))
mcmc_trace(fit$draws(inc_warmup = FALSE))
mcmc_trace(fit$draws(c("bTilde"), TRUE))
mcmc_trace(fit$draws(c("bTilde[1]", "bTilde[2]"), TRUE))
mcmc_trace(fit$draws(c("z[1]", "z[2]"), TRUE))
mcmc_trace(fit$draws(c("tau[1]", "tau[2]"), TRUE))
mcmc_pairs(fit$draws(c("bParams[1]", "bParams[2]"), TRUE))
mcmc_pairs(fit$draws(c("bTilde[1]", "bTilde[2]"), TRUE))
mcmc_pairs(fit$draws(c("z[1]", "z[2]"), TRUE))
mcmc_pairs(fit$draws(c("bPET", "bEsyMin")))
mcmc_dens(fit$draws(c("bPET", "bRain", "bMelt")))
mcmc_dens(fit$draws(c("bParams[1]", "bParams[2]", "bTilde", "z[1]", "z[2]")))
mcmc_hist(fit$draws(c("bParams[1]", "bParams[2]", "bTilde", "z[1]", "z[2]")))

mcmc_dens(fit$draws(c("bQ"))) +
  geom_density(
    color = 'red',
    linetype = 'dashed',
    aes(x = distributional::generate(distributional::dist_gamma(1, 4), 800)[[1]]))

mcmc_dens(fit$draws(c("bphiRain", "bphiMelt"))) +
  geom_density(
    color = 'red',
    linetype = 'dashed',
    aes(x = distributional::generate(distributional::dist_exponential(10), 1600)[[1]]))

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

plot_train_site <- function(site_id, pop_level = FALSE) {
  idx <- which(c("009", "053", "077", "119", "139", "140", "151", "156") == site_id)
  if(pop_level) {
    test_params <- fit$summary() %>% dplyr::filter(grepl("^[bp]", .$variable)) %>% {set_names(.$mean, .$variable)}
  } else {
    test_params <- fit$summary() %>% dplyr::filter(grepl(glue::glue("\\[{idx}\\]"), .$variable)) %>% {set_names(.$mean, stringr::str_replace(.$variable, "(g([a-zA-Z]{1,})\\[[0-9]{1,}\\])", "b\\2"))}
  }
  test_site <- unique(con_dat$site)[idx]
  dat <- con_dat[site == test_site][water_year == min(water_year)]
  wl_hat <- wl_model(
    # bPET = test_params[["bPET"]],
    # bRain = test_params[["bRain"]],
    # bMelt = test_params[["bMelt"]],
    # bQ = test_params[["bQ"]],
    # phiRain = test_params[["bphiRain"]],
    # phiMelt = test_params[["bphiMelt"]],
    pet = dat$pet_cm,
    rain = dat$rain_cm,
    melt = dat$melt_cm,
    D = nrow(dat),
    maxWL = esy_functions[site_id, max_wl],
    esyint = test_params[["bEsyInt"]],
    esyslope = test_params[["bEsySlope"]],
    esymin = test_params[["bEsyMin"]]
  )
  ylim <- c(min(c(wl_hat, dat$wl_initial_cm), na.rm = TRUE), max(c(wl_hat, dat$wl_initial_cm), na.rm = TRUE))
  plot(dat$wl_initial_cm, type = "l", main = site_id, ylim = ylim)
  lines(wl_hat, col = 'red', lty = "dashed")
}

plot_test_site <- function(site_id, pop_level = FALSE) {
  idx <- which(c("009", "053", "077", "119", "139", "140", "151", "156") == site_id)
  if(pop_level) {
    test_params <- fit$summary() %>% dplyr::filter(grepl("^[bp]", .$variable)) %>% {set_names(.$mean, .$variable)}
  } else {
    test_params <- fit$summary() %>% dplyr::filter(grepl(glue::glue("\\[{idx}\\]"), .$variable)) %>% {set_names(.$mean, stringr::str_replace(.$variable, "(g([a-zA-Z]{1,})\\[[0-9]{1,}\\])", "b\\2"))}
  }
  print(c(site_id, test_params))
  dat <- testing_data[site == site_id & site_status == "Control"][, c(.SD, list(n_records = !is.na(wl_initial_cm))), by = .(water_year)][n_records > 0]
  if(nrow(dat) == 0) return(NULL)
  dat <- dat[water_year == sample(water_year, 1)]
  wl_hat <- wl_model(
    # bPET = test_params[["bPET"]],
    # bRain = test_params[["bRain"]],
    # bMelt = test_params[["bMelt"]],
    # bQ = test_params[["bQ"]],
    # phiRain = test_params[["bphiRain"]],
    # phiMelt = test_params[["bphiMelt"]],
    pet = dat$pet_cm,
    rain = dat$rain_cm,
    melt = dat$melt_cm,
    D = nrow(dat),
    maxWL = esy_functions[site_id, max_wl],
    esyint = test_params[["bEsyInt"]],
    esyslope = test_params[["bEsySlope"]],
    esymin = test_params[["bEsyMin"]]
  )
  ylim <- c(min(c(wl_hat, dat$wl_initial_cm)), max(c(wl_hat, dat$wl_initial_cm)))
  plot(dat$wl_initial_cm, type = "l", main = site_id, ylim = ylim)
  lines(wl_hat, col = 'red', lty = "dashed")
}

par(mfrow = c(2, 4)); purrr::walk(unique(con_dat$site), ~plot_train_site(.x)); par(mfrow = c(1,1))
par(mfrow = c(2, 4)); purrr::walk(unique(con_dat$site), ~plot_train_site(.x, TRUE)); par(mfrow = c(1,1))
sites_to_test <- intersect(testing_data[site_status == "Control"]$site, con_dat$site)
par(mfrow = c(2, 4)); purrr::walk(sites_to_test, ~plot_test_site(.x)); par(mfrow = c(1,1))
par(mfrow = c(2, 4)); purrr::walk(sites_to_test, ~plot_test_site(.x, TRUE)); par(mfrow = c(1,1))

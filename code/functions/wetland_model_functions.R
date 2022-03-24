# Functions to fit wetland model

##' @title Quadratic curve

##' @return
##' @author Joe Shannon
##' @export
##' 
quad <- 
  function(wa, b0, b1, b2){
    b0 + b1*wa + b2*wa^2
  }

##' @title Derivative of a quadratic curve model

##' @return
##' @author Joe Shannon
##' @export
##' 
##' Supplies the empirical estimate of ESy

quad_prime <- 
  function(mod = NULL, wa, b1, b2){
    if(!is.null(mod)){
      b1 <- coef(mod)[["b1"]]
      b2 <- coef(mod)[["b2"]]
    }
    b1 + b2 * 2 * wa
  }

##' @title X-intercept for a quadratic model

##' @return
##' @author Joe Shannon
##' @export
##' 
##' One method of determining wetland max water level

quad_x_intercept <- 
  function(mod = NULL, b1, b2){
    if(!is.null(mod)){
      b1 <- coef(mod)[["b1"]]
      b2 <- coef(mod)[["b2"]]
    }
    
    -b1 / (b2 * 2)
  }

##' @title Calculate the point with highest density

##' @return
##' @author Joe Shannon
##' @export
##' 
##' One method of determining wetland max water level

numeric_mode <- 
  function(x, ...){
    den <- density(x, ...)
    den$x[which.max(den$y)]
  }


##' @title build_esy_functions

##' @return
##' @author Joe Shannon
##' @export
##' 
build_esy_functions <- 
  function(data){
    
    drawdown <- 
      data[, 
           .SD[which(format(sample_date, "%m%d") == "0401"):which.min(wl_initial_cm)],
           by = .(site)]
    
    drawdown[, ytd_water_balance := cumsum(rain_cm + pet_cm + melt_cm),
             by = .(site, water_year)]
    # New potential Esy model:
    # TODO: figure out why 053 water balance is so different from the other
    # sites. The shape is right, but there's a huge shift into positive water
    # balances
    # mod <- brm(bf(wl_initial_cm ~ a - (a - b) * exp(-c * ytd_water_balance), a + b + c ~ 1 + (1 | site), nl = TRUE), data = drawdown[site!="053"], cores = 4, warmup = 500, iter = 600)
    # deriv = (a - b) * (exp(-c * x) * c)
    # Alternatively could so some sort of mixture of two models (a curve and 
    # a linear model), or more simply two linear models. A steep model for low
    # water levels and a shallow model for high water levels. Then have a mixing
    # parameter alpha that ranges from 0-1 around some inflection point.
    # drawdown ~ alpha * high_esy + (1 - alpha) * low_esy
    # emperical_esy = deriv(drawdown_mod)
    brm
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
    # split(drawdown, f = drawdown$site) %>% purrr::map_dfr(~coef(lm((34.77 - 12.894) * (exp(-0.059 * ytd_water_balance) * 0.059) ~ wl_initial_cm, data = .x)), .id = "site") 
    esy_form <- 
      bf(esy_emp ~ a - (a - b) * exp (c * wl_initial_cm),
         a + b + c ~ 0 + site,
         nl = TRUE)
    
    esy_priors <- 
      prior(nlpar = "a", normal(10, 0.5), lb = 0) +
      prior(nlpar = "b", normal(2, 0.25), lb = 0) +
      prior(nlpar = "c", normal(0.01, 0.002), lb = 0) +
      prior(gamma(2, .1), class = nu)
    
    brm_esy <- 
      brm(esy_form,
          prior = esy_priors,
          family = student,
          cores = 4,
          save_model = "output/models/esy_model.stan", 
          data = drawdown[!is.na(wl_initial_cm)])
    
    esy_coefs <- 
      fixef(brm_esy)
    
    esy_a <- 
      esy_coefs[str_detect(rownames(esy_coefs), "a_"), "Estimate"] %>% 
      set_names(~str_extract(., "[0-9]{3}"))
    
    esy_b <- 
      esy_coefs[str_detect(rownames(esy_coefs), "b_"), "Estimate"] %>% 
      set_names(~str_extract(., "[0-9]{3}"))
    
    esy_c <- 
      esy_coefs[str_detect(rownames(esy_coefs), "c_"), "Estimate"] %>% 
      set_names(~str_extract(., "[0-9]{3}"))
    
    drawdown[, `:=`(a = esy_a[[.BY[[1]]]],
                    b = esy_b[[.BY[[1]]]],
                    c = esy_c[[.BY[[1]]]]), 
             by = .(site)]
    
    drawdown[, esy_hat := a - (a - b) * exp (c * wl_initial_cm)]
    
    esy_functions <- 
      drawdown[, 
               .(pred_fun = list(as.function(list(wl = NULL, 
                                                  min.esy = NULL,
                                                  bquote(pmax(min.esy,
                                                              .(a) - (.(a) - .(b)) * exp (.(c) * wl)),
                                                         where = .SD[1, .(a, b, c)])))),
                 max_wl = numeric_mode(wl_initial_cm, na.rm = TRUE)),
            keyby = .(site)]
    
    esy_functions
  }


##' @title Wetland model

##' PET adjustment has pmin(1, ...) because the adjustment regression
##' could show an increase in PET for extermely high water levels, but
##' I don't think that's realistic

##' @return
##' @author Joe Shannon
##' @export
##' 
wetland_model <- 
  function(data, params, future.forest.change = NA){

    # Percentage of Basal Area As Ash
    # 1 009   0.752355 - 2012 Veg Surveys
    #   053   0.756437 - 2014 Veg Surveys
    # 2 077   0.440609 - 2012 Veg Surveys
    #   111   0.48
    #   113   0.663741 - 2015 Veg Surveys
    # 3 119   0.838836 - 2012 Veg Surveys
    # 4 135   0.261113 - 2012 Veg Surveys
    #   139   0.38
    # 5 140   0.703470 - 2012 Veg Surveys
    # 6 151   0.535355 - 2012 Veg Surveys
    # 7 152   0.663909 - 2012 Veg Surveys
    # 8 156   0.850605 - 2012 Veg Surveys
    # 9 157   0.821684 - 2012 Veg Surveys
    
    MPET <-
      params[["MPET"]]
    
    if(is.na(future.forest.change)){
      pet_fun <- 
        function(wl, max.wl, future.forest.change){ 1 }
    } else {
      pet_fun <- 
        function(wl, max.wl, future.forest.change){
          (1 - future.forest.change) + future.forest.change * pmin(1, abs(1/(1.45077 - 0.05869 * (wl - max.wl))))
        }
    }
    
    MP <- 
      params[["MP"]]
    
    MM <-
      params[["MM"]]
    
    MQ <-
      params[["MQ"]]
    
    minESY <- 
      params[["minESY"]]
    
    # MG <- 
    #   params[["MG"]]
    
    phiM <- 
      params[["phiM"]]
    
    phiP <- 
      params[["phiP"]]
    
    maxWL <- 
      params[["maxWL"]]
    
    if(inherits(params[["funESY"]], "list")){
      esy_fun <- 
        params[["funESY"]][[1]]
      
    } else {
      esy_fun <- 
        params[["funESY"]]
    }
    
    PET <- data$pet_cm
    P <- as.numeric(filter(data$rain_cm, filter = phiP, method = "rec", sides = 1))
    M <- as.numeric(filter(data$melt_cm, filter = phiM, method = "rec", sides = 1))
    WL <- data$wl_initial_cm
    
    n <- 
      length(PET)
    
    # Create empty vectors
    wl_hat <- 
      gradient <- 
      q_hat <- 
      m_hat <- 
      p_hat <- 
      pet_hat <- 
      g_hat <- 
      numeric(n)
    
    # Initialize model at full water level
    wl_hat[1] <- 
      maxWL
    
    # Loop through weather data
    for(t in 2:n){

      ######## Esy
      # Calculate gradient2 of drawdown 
      gradient[t] <- 
        esy_fun(wl_hat[t-1], minESY)

      ######## PET or P times Esy
      # Use net input to determine if water level increases or decreases
      # Assuming AET is negligible on days where P >= PET
      if((P[t] + PET[t]) <= 0) {
        # Water level drawdown = PET2, if P2 <= PET2 then it can be assumed to be
        # less than interception (not necessarily true, but works as a
        # simplifying assumption)
        pet_hat[t] <- 
          pet_fun(wl_hat[t-1], maxWL, future.forest.change) * (MPET * PET[t]) * gradient[t]
        wl_hat[t] <-
          wl_hat[t-1] + pet_hat[t]
      } else {
        p_hat[t] <- 
          (MP * P[t]) * gradient[t]
        wl_hat[t] <-
          wl_hat[t-1] + p_hat[t]
      }

      ######## Snowmelt
      # Directly add melt to water level. This should probably have some sort of
      # multiplier
      m_hat[t] <- 
        MM * M[t] * gradient[t]
      wl_hat[t] <-
        wl_hat[t] + m_hat[t]

      ####### Streamflow
      # If WL is above spill point threshold then lose some to streamflow. 
      # This could probably be improved using the morphology models to determine
      # streamflow
      if(wl_hat[t-1] > maxWL){
        q_hat[t] <-
          MQ * (wl_hat[t-1] - maxWL)
        q_hat[t] <- 
          pmin(q_hat[t], wl_hat[t] - maxWL)
        wl_hat[t] <-
          wl_hat[t] - q_hat[t]
      }

    }
    data.table::data.table(wl_hat, q_hat, m_hat, p_hat, pet_hat, gradient)
  }

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
                     
                     # obs_weights <- 
                     #   as.data.table(density(data[, wl_initial_cm], 
                     #                         na.rm = TRUE)[c("x", "y")]
                     #                 )[, .(x, weight = -log(y))]
                     # 
                     # weights_index <- 
                     #   map_int(wl_hat, 
                     #           ~if(is.na(.x)){
                     #             return(NA_integer_)
                     #           } else {
                     #             which.min(abs(.x - obs_weights$x))
                     #           })
                     # 
                     # resid_weights <- 
                     #   obs_weights[weights_index, weight]
                     
                     # -sum(
                     #   # resid_weights[!is.na(data$wl_initial_cm)] * 
                     #     dnorm(resids[!is.na(data$wl_initial_cm)],
                     #           mean = 0,
                     #           sd = sd(diff(data$wl_initial_cm), na.rm = TRUE) / sum(!is.na(data$wl_initial_cm)),
                     #           log = TRUE))
                     
                     # hydroGOF::rsr(wl_hat[!is.na(data$wl_initial_cm)], data$wl_initial_cm[!is.na(data$wl_initial_cm)])
                     
                     # Initial weights are equal
                     init_weight <- 1 / sum(is.na(data$wl_initial_cm))
                     
                     # Weights increase asymmetrically as water levels drop
                     wghts <- pmax(init_weight, init_weight * (params$maxWL - data$wl_initial_cm))
                     
                     # Weights are squared
                     sqrt(weighted.mean(resids^2, w = wghts^2, na.rm = TRUE))
                     
                     
                   })
    
    opt$par <- c(opt$par, unlist(fixed))
    
    opt
  }



fit_models <- 
  function(data, 
           esy_mods,
           par = NULL,
           fixed = NULL,
           ...){
    
    optim_res <- 
      data[, .(res = 
                  list(optimize_params(
                    data = .SD, 
                    par = par,
                    fixed = c(fixed,
                              maxWL = esy_mods[.BY[[1]], max_wl],
                              funESY = esy_mods[.BY[[1]], pred_fun[[1]]]),
                    ...,
                    method = "L-BFGS-B",
                    lower = list(MPET = 0, 
                                 MP = 0, 
                                 MM = 0, 
                                 MQ = 0, 
                                 minESY = 0, 
                                 phiM = 0, 
                                 phiP = 0),
                    upper = list(MPET = Inf, 
                                 MP = Inf, 
                                 MM = Inf,  
                                 MQ = 1-.Machine$double.neg.eps,
                                 minESY = Inf, 
                                 phiM = 1-.Machine$double.neg.eps, 
                                 phiP = 1-.Machine$double.neg.eps)))),
            keyby = .(site)]
    
    
    optim_res[, `:=`(params = list(pluck(res, 1, "par")),
                     optimize_value = pluck(res, 1, "value"),
                     iterations = pluck(res, 1, "counts", 1),
                     convergence = pluck(res, 1, "convergence"),
                     message = pluck(res, 1, "message")),
              by = .(site)]
    
    optim_res
  }

refit_model <- 
  function(data, 
           optimized.models,
           refit = NULL,
           ...){
    
    if(is.null(refit)) stop("Specify at least one parameter to refit")
    
    existing_params <- 
      names(optimized.models[["params"]][[1]])
    
    refit_params <- 
      names(refit)

    if(any(!(refit_params %in% existing_params))){
      stop("The following refit parameters are not in the original models",
           refit[which(!(refit_params %in% existing_params))])
    }
        
    fixed <- 
      optimized.models[, 
                       map_dfr(params, as.data.table), 
                       by = .(site)]
    
    lower_bound <- 
      fixed[, .SD, .SDcols = c("site", refit_params)]
    
    upper_bound <- 
      fixed[, .SD, .SDcols = c("site", refit_params)]
      
    fixed <- 
      fixed[, c(refit_params) :=NULL ]
    
    # Set upper and lower bounds based on expected increase or decrease due to
    # treatment response (or based on physical parameters (0, -Inf, Inf, etc))
    # variables not modified here will use the control optimization values for
    # the upper/lower bound
    
    lower_bound <- 
      modify_at(lower_bound, c("MPET", "MM", "MQ", "minESY", "phiM", "phiP"), ~0)
    
    upper_bound <- 
      modify_at(upper_bound, c("MP", "MM", "minESY", "maxWL"), ~Inf) %>% 
      modify_at(c("MQ", "phiM", "phiP"), ~1-.Machine$double.neg.eps)
    
    refit_res <- 
      data[, .(res = 
                 list(optimize_params(
                   data = .SD, 
                   par = refit,
                   fixed = fixed[.BY[[1]], -c("site")],
                   ...,
                   method = "L-BFGS-B",
                   lower = lower_bound[.BY[[1]], -c("site")],
                   upper = upper_bound[.BY[[1]], -c("site")]))),
           keyby = .(site)]
    
    refit_res[, `:=`(params = list(pluck(res, 1, "par")),
                     optimize_value = pluck(res, 1, "value"),
                     iterations = pluck(res, 1, "counts", 1),
                     convergence = pluck(res, 1, "convergence"),
                     message = pluck(res, 1, "message")),
              by = .(site)]
    
    refit_res
    
  }
# Need to get IDW data for all external met sites

# Derived ESY Method

# 1. Fit asymptotic regression model Esy ~ a - (a - b) * exp (- c * wl) with 
#    dlst errors using optim. This supersedes having to fit a drawdown model and
#    then fit a second model to the gradient of that function to estimate ESy
# 2. Optimize MPET, MP, MM, and MQ for the training year
# 3. Validate optimized model against test datasets

quad <- 
  function(wa, b0, b1, b2){
    b0 + b1*wa + b2*wa^2
  }

# No convergence for 113, 119, 139. Initial looks make it seem like a better fit
# for the drawdown curve, so it should be used
quad_cp <- 
  function(wa, b0, cp, b2){
    b0 + (wa <= cp) * b2*(wa + cp)^2
  }

quad_prime_cp <- 
  function(mod = NULL, wa, cp, b2){
    if(!is.null(mod)){
      cp <- coef(mod)[["cp"]]
      b2 <- coef(mod)[["b2"]]
    }
    
    (wa <= cp) * b2 * 2 * (wa + cp)
  }

quad_x_intercept_cp <- 
  function(mod = NULL, cp, b2){
    
  }

quad_prime <- 
  function(mod = NULL, wa, b1, b2){
    if(!is.null(mod)){
      b1 <- coef(mod)[["b1"]]
      b2 <- coef(mod)[["b2"]]
      }
    
    grads <- b1 + b2 * 2 * wa
    
    attr(grads, "intercept") <- 
      b1
    
    grads
  }

quad_x_intercept <- 
  function(mod = NULL, b1, b2){
    if(!is.null(mod)){
      b1 <- coef(mod)[["b1"]]
      b2 <- coef(mod)[["b2"]]
    }
    
    -b1 / (b2 * 2)
  }

wetland_model <- 
  function(data, params){
    MPET <-
      params[["MPET"]]
    
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
    
    esy_fun <- 
      params[["minESY"]]
    
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
      
      # Calculate gradient2 of drawdown 
      gradient[t] <- 
        esy_fun(wl_hat[t-1], minESY)
      
      # Use net input to determine if water level increases or decreases
      if((P[t] + PET[t]) <= 0){
        
        # Water level drawdown = PET2, if P2 <= PET2 then it can be assumed to be
        # less than interception (not necessarily true, but works as a 
        # simplifying assumption)
        pet_hat[t] <- 
          (MPET * PET[t]) * gradient[t]
        
        wl_hat[t] <-
          wl_hat[t-1] + pet_hat[t]
        
      } else {
        
        # Water rise is P2 - PET2, which fits better, but need to rethink my
        # justification for this
        p_hat[t] <- 
          (MP * P[t]) * gradient[t]
        
        wl_hat[t] <-
          wl_hat[t-1] + p_hat[t]
        
      }
      
      # # Add in G
      # if(P[t-1] + PET[t-1] > 0){
      #   wl_hat[t] <- 
      #     wl_hat[t] + MG * P[t-1]
      # }
      
      
      # Directly add melt to water level. This should probably have some sort of
      # multiplier
      
      m_hat[t] <- 
        MM * M[t] * gradient[t]
      
      wl_hat[t] <-
        wl_hat[t] + m_hat
      
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
    data.table::data.table(wl_hat, q_hat, m_hat, p_hat, pet_hat)
  }

optimize_params <- 
  function(data, par, fixed = NULL, ...){
    opt <- optim(par = par,
                 fixed = fixed,
                 ...,
                 control = list(fnscale = -1, 
                                maxit = 2000),
                 fn = 
                   function(params, fixed = NULL, data = data){
                     
                     params <- 
                       c(params, fixed)
                     
                     wl_hat <- 
                       wetland_model(data, params)$wl_hat
                     
                     resids <- 
                       (wl_hat - data$wl_initial_cm)[!is.na(data$wl_initial_cm)]
                     
                     # -sum(dnorm(resids, 
                     #            mean = 0,
                     #            sd = sd(resids),
                     #            log = TRUE))
                     
                     hydroGOF::md(wl_hat[!is.na(data$wl_initial_cm)], data$wl_initial_cm[!is.na(data$wl_initial_cm)])
                   })
    
    opt$par <- c(opt$par, unlist(fixed))
    
    opt
  }


# Prep Data ---------------------------------------------------------------


source("code/load_project.R")
library(stringr)
library(extraDistr, 
        include.only = "dlst")
tar_load(external_met)
tar_load(water_budget)


kenton <- 
  external_met[station_name == "kenton"]


train <- 
  kenton[water_year %in% 2012:2019 & format(sample_date, "%m%d") != "0229",
         .(sample_date, 
           water_year,
           pet_cm = -pet_cm,
           rain_cm,
           melt_cm,
           total_input_cm,
           tmax_c,
           tmin_c,
           tmean_c)
  ]

train[between(melt_cm, 0, 0.1) & month(sample_date) %in% 6:9, 
      `:=`(total_input_cm = total_input_cm - melt_cm,
           melt_cm = 0)]

# Drop anomalous pet
train[tmax_c <= 0,
      pet_cm := 0]

# Calculate YTD after dropping above points
train[, 
      `:=`(ytd_pet_cm = cumsum(pet_cm),
           ytd_p_cm = cumsum(rain_cm),
           ytd_m_cm = cumsum(melt_cm),
           ytd_input_cm = cumsum(rain_cm + melt_cm)),
      by = .(water_year)]

train[, ytd_water_balance := ytd_pet_cm + ytd_p_cm + ytd_m_cm]
train[, dWB := c(0, diff(ytd_water_balance))]

train <-
  train[water_year == 2012][
    water_budget[, .(site, sample_date, wl_initial_cm)],
    on = "sample_date",
    nomatch = NULL
  ]

train[, 
      c("esy_emp",
        "esy_x_intercept") := {
          
          mod <- 
            robustbase::nlrob(wl_initial_cm ~ quad(ytd_water_balance, 
                                                   b0 = b0, b1 = b1, b2 = b2),
                              data = .SD[1:which.min(wl_initial_cm)], 
                              na.action = na.exclude,
                              maxit = 50,
                              start = list(b0 = 8, b1 = 1, b2 = -1))
            
          list(esy_emp = quad_prime(mod = mod, wa = .SD[, ytd_water_balance]),
               esy_x_intercept = quad_x_intercept(mod))
          
          },
      by = .(site)]

train[, hydro_period := rep(c("drawdown", "rebound"), 
                            times = c(which.min(wl_initial_cm), .N - which.min(wl_initial_cm))),
      by = .(site)]

ggplot(train[hydro_period == "drawdown"],
       aes(x = wl_initial_cm,
           y = esy_emp,
           color = site)) +
  geom_point()




# Calculate ESy Functions -------------------------------------------------
# Could potentially be improved by also looking at water balance. Something like
# a 'ghost' water level measurement that doesn't recover until water
# availability increases

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
      data = train[hydro_period == "drawdown"])

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

train[, `:=`(a = esy_a[[.BY[[1]]]],
             b = esy_b[[.BY[[1]]]],
             c = esy_c[[.BY[[1]]]]), 
      by = .(site)]

train[, esy_hat := a - (a - b) * exp (c * wl_initial_cm)]

esy_functions <- 
  train[, .(pred_fun = list(as.function(list(wl = NULL, min.esy = NULL,
                                        bquote(pmin(min.esy,
                                                    .(a) - (.(a) - .(b)) * exp (.(c) * wl)),
                                               where = .SD[1, .(a, b, c)])))),
            max_wl = density(na.omit(wl_initial_cm))$x[which.max(density(na.omit(wl_initial_cm))$y)]), 
        by = .(site)]

# 119 isn't great (has an early peak in spring). Could also work to identify 
# better local min(wl_initial_)
ggplot(train[hydro_period == "drawdown"],
       aes(x = wl_initial_cm,
           y = esy_emp)) +
  geom_point(aes(color = site)) +
  geom_line(aes(y = esy_hat)) +
  geom_vline(data = esy_fun,
             aes(xintercept = max_wl)) +
  facet_wrap(~site)







SITE <- 
  "135"



dat <-
  kenton[water_year %in% 2012:2019 & format(sample_date, "%m%d") != "0229",
         .(sample_date,
           water_year,
           pet_cm = -pet_cm,
           rain_cm,
           melt_cm,
           total_input_cm,
           tmax_c,
           tmin_c,
           tmean_c)
  ][water_budget[site == SITE],
    `:=`(site = i.site,
         best_precip_cm = i.best_precip_cm,
         idw_precip_cm = i.idw_precip_cm,
         wl_initial_cm = i.wl_initial_cm),
    on = "sample_date"]

dat[, site := unique(na.omit(site))]

dat[, coal_precip_cm := pmax(rain_cm, best_precip_cm, na.rm = TRUE)]

# Drop anomalous melt
dat[between(melt_cm, 0, 0.1) & month(sample_date) %in% 6:9, 
    `:=`(total_input_cm = total_input_cm - melt_cm,
         melt_cm = 0)]

# Drop anomalous pet
dat[tmax_c <= 0,
    pet_cm := 0]

# Calculate YTD after dropping above points
dat[, 
    `:=`(ytd_pet_cm = cumsum(pet_cm),
         ytd_p_cm = cumsum(rain_cm),
         ytd_m_cm = cumsum(melt_cm),
         ytd_input_cm = cumsum(rain_cm + melt_cm)),
    by = .(water_year)]

dat[, ytd_water_balance := ytd_pet_cm + ytd_p_cm + ytd_m_cm]
# dat[, ytd_water_balance := ytd_water_balance - max(ytd_water_balance),
#     by = .(water_year)]

train <- 
  dat[water_year == 2012]

# max.wl <- 
  # density(water_budget[site == SITE, na.omit(wl_initial_cm)])$x[which.max(density(water_budget[site == SITE, na.omit(wl_initial_cm)])$y)]

# Optimize ESy -------------------------------------------------------------

# Need to still fit the drawdown model then use quad_prime to get reasonable
# starting values for esy optimization

esy_fun <- 
  as.function(list(wl = NULL,
                   min.esy = 0,
                   max.wl = FALSE,
                   bquote(
                     if(max.wl) {
                       return(log((.(a)/(.(a)-.(b)))) / .(c))
                     } else {
                       pmax(min.esy, .(a) - (.(a) - .(b)) * exp (.(c) * wl))
                     },
                     where = as.list(optim(par = list(a = 10, b = 1, c = 0.02, sigma = 0.5, nu = 3),
                                           # method = "L-BFGS-B",
                                           # lower = list(a = .Machine$double.eps,
                                           #              b = .Machine$double.eps,
                                           #              c = .Machine$double.eps,
                                           #              sigma = .Machine$double.eps,
                                           #              nu = 1),
                                           # upper = list(a = 12,
                                           #              b = 20,
                                           #              c = 1 - .Machine$double.eps, 
                                           #              sigma = 2 * sqrt(var(train[1:which.min(wl_initial_cm)][!is.na(wl_initial_cm), wl_initial_cm])/length(train[1:which.min(wl_initial_cm)][!is.na(wl_initial_cm), wl_initial_cm])),
                                           #              nu = length(train[1:which.min(wl_initial_cm)][!is.na(wl_initial_cm), wl_initial_cm]) - 1),
                                           control = list(fnscale = 1, 
                                                          maxit = 2000),
                                           fn = 
                                             function(params){
                                               
                                               a <- params[["a"]]
                                               b <- params[["b"]]
                                               c <- params[["c"]]
                                               sigma <- params[["sigma"]]
                                               nu <- params[["nu"]]
                                               
                                               WL <- train[1:which.min(wl_initial_cm)][!is.na(wl_initial_cm), wl_initial_cm]
                                               dWA <- train[1:which.min(wl_initial_cm)][!is.na(wl_initial_cm), c(0, diff(ytd_water_balance))]
                                               
                                               n <- 
                                                 length(dWA)
                                               
                                               wl_hat <- 
                                                 gradient <- 
                                                 numeric(n)
                                               
                                               # Initialize model at full water level
                                               wl_hat[1] <- 
                                                 WL[1]
                                               
                                               # Loop through weather data
                                               for(t in 2:n){
                                                 
                                                 # Calculate gradient of drawdown 
                                                 gradient[t] <- 
                                                   a - (a - b) * exp (c * wl_hat[t-1])
                                                 
                                                 wl_hat[t] <- 
                                                   dWA[t-1] * gradient[t] + wl_hat[t-1]
                                               }
                                               
                                               resids <- 
                                                 (wl_hat - WL)[!is.na(WL)]
                                               
                                               # -sum(dnorm(resids, log = TRUE))
                                               -sum(dlst(resids,
                                                         df = nu,
                                                         sigma = sigma,
                                                         log = TRUE))
                                               
                                             })$par))))


max.wl <- 
  esy_fun(max.wl = TRUE)

# Parametize Water Balance Coefs ------------------------------------------

ggplot(dat,
       aes(x = sample_date,
           y = wl_initial_cm)) +
  geom_line() +
  facet_wrap(~water_year,
             scales = "free_x")

param_train <- 
  dat[water_year %in% c(2013)]

(mod_opt <- 
    optimize_params(par = list(MP = 1.5, MM = 0.1, ESY = 1, MPET = 1, phiM = 0.9, MQ = 0.5, phiP = 0.5),
                    method = "L-BFGS-B",
                    lower = list(MP = 0, MM = 0, ESY = 0, MPET = 0, phiM = 0, MQ = 0, phiP = 0),
                    upper = list(MP = Inf, MM = Inf, ESY = Inf, MPET = Inf, phiM = 0.9, MQ = 0.9, phiP = 0.9)))

# Test Fit Model ----------------------------------------------------------

test <- 
  dat[water_year %in% 2017]

{PET2 <- test$pet_cm
P2 <- test$rain_cm
P2 <- M2 <- as.numeric(filter(P2, filter = mod_opt$par[["phiP"]], method = "rec", sides = 1))
M2 <- test$melt_cm
M2 <- as.numeric(filter(M2, filter = mod_opt$par[["phiM"]], method = "rec", sides = 1))
D2 <- test$sample_date
WL2 <- test$wl_initial_cm

# Length of weather vector
n <- 
  length(PET2)

# Create empty vectors
wl_hat2 <- 
  gradient2 <- 
  q_hat2 <- 
  g_hat2 <- 
  m_hat2 <- 
  numeric(n)

# Initialize model at full water level
wl_hat2[1] <- 
  max.wl

# Loop through weather data
for(t in 2:n){
  
  
  # Calculate gradient2 of drawdown 
  gradient2[t] <- 
    esy_fun(wl_hat2[t-1], mod_opt$par[["ESY"]])
  
  # Use net input to determine if water level increases or decreases
  if((P2[t] + PET2[t]) <= 0){
    
    # Water level drawdown = PET2, if P2 <= PET2 then it can be assumed to be
    # less than interception (not necessarily true, but works as a 
    # simplifying assumption)
    wl_hat2[t] <-
      wl_hat2[t-1] + (mod_opt$par[["MPET"]] * PET2[t]) * gradient2[t]
    
  } else {
    
    # Water rise is P2 - PET2, which fits better, but need to rethink my
    # justification for this
    wl_hat2[t] <-
      wl_hat2[t-1] + (mod_opt$par[["MP"]]*P2[t]) * gradient2[t]
    
  }
  # Add in G
  # if(P2[t-1] + PET2[t-1] > 0){
  #   wl_hat2[t] <- 
  #     wl_hat2[t] + mod_opt$par[["MG"]] * P2[t-1]
  # }
  # Directly add melt to water level. This should probably have some sort of
  # multiplier
  
  m_hat2[t] <- 
    mod_opt$par[["MM"]] * M2[t] * gradient2[t]
  
  wl_hat2[t] <-
    wl_hat2[t] + m_hat2[t]
  
  # If WL is above spill point threshold then lose some to streamflow. 
  # This could probably be improved using the morphology models to determine
  # streamflow
  
  if(wl_hat2[t] > max.wl){

    q_hat2[t] <-
      mod_opt$par[["MQ"]] * (wl_hat2[t] - max.wl)

    wl_hat2[t] <-
      wl_hat2[t] - q_hat2[t]
  }
  
}
lplot(wl_hat2 ~ D2, ylim = c(min(c(wl_hat2, WL2), na.rm = TRUE), max(c(wl_hat2, WL2), na.rm = TRUE))); lines(WL2 ~ D2, col = 'gray40')
hydroGOF::md(wl_hat2[!is.na(WL2)], WL2[!is.na(WL2)])
}




################################################################################
# MODEL PERFROMANCE SEEMED TO DEGRADE SINCE INITIAL COMMMIT. MAY NEED TO REWIND
# AND TEST AGAIN. LOOKS LIKE ESY MAY BE GETTING OVER ESTIMATED FOR SOME SITES
# PLACE BOUNDS ON PARAMS
################################################################################

# Derived ESY Method

# 1. Fit drawdown model of wl_initial_cm ~ quad(ytd_water_balance) for selected
#    year (largest drawdown year is best)
# 2. Use b1, b2 to calculate quad_prime() as an approximation of Esy
# 3. Fit asymptotic regression model Esy ~ a - (a - b) * exp (- c * wl)
# 4. Optimize MPET, MP, MM, and MQ for the training year
# 5. Validate optimized model agains test datasets

quad <- 
  function(wa, b0, b1, b2){
    b0 + b1*wa + b2*wa^2
  }

quad_prime <- 
  function(wa, b1, b2){
    b1 + b2 * 2 * wa
  }

optimize_params <- 
  function(par, fixed = NULL, ...){
    opt <- optim(par = par,
                 fixed = fixed,
                 ...,
                 control = list(fnscale = -1, 
                                maxit = 2000),
                 fn = 
                   function(params, fixed = NULL){
                     
                     PET <- param_train$pet_cm
                     P <- param_train$rain_cm
                     M <- param_train$melt_cm
                     WL <- param_train$wl_initial_cm
                     
                     params <- 
                       c(params, fixed)
                     
                     MPET <-
                       params[["MPET"]]
                     
                     MP <- 
                       params[["MP"]]
                     
                     # MM <- 
                     #   params[["MM"]]
                     
                     MQ <-
                       params[["MQ"]]
                     
                     ESY <- 
                       params[["ESY"]]
                     
                     # MG <- 
                     #   params[["MG"]]
                     
                     n <- 
                       length(PET)
                     
                     # Create empty vectors
                     wl_hat <- 
                       gradient <- 
                       q_hat <- 
                       g_hat <- 
                       numeric(n)
                     
                     # Initialize model at full water level
                     wl_hat[1] <- 
                       max.wl
                     
                     # Loop through weather data
                     for(t in 2:n){
                       
                       # Calculate gradient2 of drawdown 
                       gradient[t] <- 
                         esy_fun(wl_hat[t-1], ESY)
                       
                       # Use net input to determine if water level increases or decreases
                       if((P[t] + PET[t]) <= 0){
                         
                         # Water level drawdown = PET2, if P2 <= PET2 then it can be assumed to be
                         # less than interception (not necessarily true, but works as a 
                         # simplifying assumption)
                         wl_hat[t] <-
                           wl_hat[t-1] + (MPET * PET[t]) * gradient[t]
                         
                       } else {
                         
                         # Water rise is P2 - PET2, which fits better, but need to rethink my
                         # justification for this
                         wl_hat[t] <-
                           wl_hat[t-1] + (MP * P[t]) * gradient[t]
                         
                       }
                       
                       # # Add in G
                       # if(P[t-1] + PET[t-1] > 0){
                       #   wl_hat[t] <- 
                       #     wl_hat[t] + MG * P[t-1]
                       # }
                       
                       
                       # Directly add melt to water level. This should probably have some sort of
                       # multiplier
                       wl_hat[t] <-
                         wl_hat[t] + MP * M[t] * gradient[t]
                       
                       # If WL is above spill point threshold then lose some to streamflow. 
                       # This could probably be improved using the morphology models to determine
                       # streamflow
                       
                       if(wl_hat[t] > max.wl){
                         
                         q_hat[t] <-
                           MQ * (wl_hat[t] - max.wl)
                         
                         wl_hat[t] <-
                           wl_hat[t] - q_hat[t]
                       }
                       
                     }
                     
                     resids <- 
                       (wl_hat - WL)[!is.na(WL)]
                     
                     # -sum(dnorm(resids, 
                     #            mean = 0,
                     #            sd = sd(resids),
                     #            log = TRUE))
                     
                     hydroGOF::md(wl_hat[!is.na(WL)], WL[!is.na(WL)])
                   })
    
    opt$par <- c(opt$par, unlist(fixed))
    
    opt
  }


# Prep Data ---------------------------------------------------------------


source("code/load_project.R")
tar_load(external_met)
tar_load(water_budget)


SITE <- 
  "135"

kenton <- 
  external_met[station_name == "kenton"]

dat <- 
  kenton[water_year %in% 2012:2019 & format(sample_date, "%m%d") != "0229",
         .(sample_date, 
           water_year,
           pet_cm = -pet_cm,
           # ytd_pet_cm = cumsum(-pet_cm),
           rain_cm,
           # ytd_p_cm = cumsum(rain_cm),
           melt_cm,
           # ytd_m_cm = cumsum(melt_cm),
           total_input_cm,
           # ytd_input_cm = cumsum(total_input_cm),
           tmax_c,
           tmin_c,
           tmean_c)
  ][water_budget[site == SITE],
    `:=`(site = i.site,
         wl_initial_cm = i.wl_initial_cm),
    on = "sample_date"]

dat[, site := unique(na.omit(site))]

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

# Calculate Drawdown ------------------------------------------------------

max.wl <- 
  density(water_budget[site == SITE, na.omit(wl_initial_cm)])$x[which.max(density(water_budget[site == SITE, na.omit(wl_initial_cm)])$y)]

(init_mod <- 
    robustbase::nlrob(wl_initial_cm ~ quad(ytd_water_balance, 
                                           b0 = b0, b1 = b1, b2 = b2),
                      data = train[1:which.min(wl_initial_cm)],
                      na.action = na.exclude,
                      maxit = 50,
                      start = list(b0 = 8, b1 = 1, b2 = -1)))

ggplot(data = train[1:which.min(wl_initial_cm), 
                  .(wl_initial_cm = wl_initial_cm, 
                    ytd_water_balance = ytd_water_balance)],
       aes(x = ytd_water_balance,
           y = wl_initial_cm)) +
  geom_path() +
  geom_function(fun = ~predict(init_mod, newdata = data.frame(ytd_water_balance = .x)),
                color = palette()[4]) 

b1 <- coef(init_mod)[["b1"]]
b2 <- coef(init_mod)[["b2"]]

train[, esy := quad_prime(ytd_water_balance, b1, b2)]

ggplot(data = train[1:which.min(wl_initial_cm)]) + 
  aes(x = wl_initial_cm, 
      y = esy) +
  geom_point() 

esy_mod <- 
  robustbase::nlrob(esy ~ (a - (a - b) * exp (- c * wl_initial_cm)),
                    start = c(a = 10, b = 0, c = -0.3), na.action = na.exclude,
                    data = train[1:which.min(wl_initial_cm)])
  # lm(quad_prime(WA, b1, b2)[1:which.min(WL)] ~ negWL[1:which.min(WL)] + I(negWL[1:which.min(WL)]^0.5))
  # quantreg::rq(quad_prime(WA, b1, b2)[1:which.min(WL)] ~ WL[1:which.min(WL)])
  # mle2(esy ~ dlst(a - (a - b) * exp (c * wl), sigma = sigma, df = nu), 
  #      start = list(a = 10, b = 0, c = -0.3, sigma = 1, nu = 3), 
  #      data = data.frame(esy = quad_prime(WA, b1, b2)[1:which.min(WL)], wl= negWL[1:which.min(WL)])[!is.na(data.frame(esy = quad_prime(WA, b1, b2)[1:which.min(WL)], wl= negWL[1:which.min(WL)])$wl), ],
  #      method = "L-BFGS-B",
  #      lower = list(a = -Inf, b = 0, c = -Inf, sigma = .Machine$double.eps, nu = 1),
  #      upper = list(a = Inf, b = Inf, c = 0, sigma = 1, nu = sum(!is.na(WL[1:which.min(WL)]))-1))
esy_fun <- 
  as.function(list(wl = NULL,
                   bquote(pmax(1, .(a) - (.(a) - .(b)) * exp(-.(c) * wl)),
                          where = as.list(coef(esy_mod)))))
  # as.function(list(wl = NULL,
  #                  offset = NULL,
  #                  bquote(ifelse(wl > offset, 
  #                                1,
  #                                pmax(1, .(b) + .(m1) * (-wl + offset) + .(m2) * (-wl + offset)^0.5)),
  #                         where = list(b = coef(esy_mod)[[1]],
  #                                      m1 = coef(esy_mod)[[2]],
  #                                      m2 = coef(esy_mod)[[3]]))))
  # as.function(list(wl = NULL, 
  #                  bquote(pmax(.(b), .(b) + .(m)*wl),
  #                         where = list(b = coef(esy_mod)[[1]],
  #                                      m = coef(quantreg::rq(quad_prime(WA, b1, b2)[1:which.min(WL)] ~ WL[1:which.min(WL)]))[[2]]))))

last_plot() +
  geom_function(fun = esy_fun)

param_train <- 
  dat[water_year == 2014]

(mod_opt <- 
    optimize_params(par = list(MP = 1.5, MM = 1.5, MQ = 0.5),
                    fixed = list(MPET = 1)))

# Test Fit Model ----------------------------------------------------------

test <- 
  dat[water_year %in% 2013]

{PET2 <- test$pet_cm
P2 <- test$rain_cm
M2 <- test$melt_cm
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
  numeric(n)

# Initialize model at full water level
wl_hat2[1] <- 
  max.wl

# Loop through weather data
for(t in 2:n){
  
  
  # Calculate gradient2 of drawdown 
  gradient2[t] <- 
    esy_fun(wl_hat2[t-1])
  
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
  wl_hat2[t] <-
    wl_hat2[t] + mod_opt$par[["MM"]] * M2[t] * gradient2[t]
  
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
lplot(wl_hat2 ~ D2); lines(WL2 ~ D2, col = 'gray40')
hydroGOF::md(wl_hat2[!is.na(WL2)], WL2[!is.na(WL2)])
}




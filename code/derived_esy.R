# Need to get IDW data for all external met sites

# Derived ESY Method

# 1. Fit asymptotic regression model Esy ~ a - (a - b) * exp (- c * wl) with 
#    dlst errors using optim. This supersedes having to fit a drawdown model and
#    then fit a second model to the gradient of that function to estimate ESy
# 2. Optimize MPET, MP, MM, and MQ for the training year
# 3. Validate optimized model against test datasets

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
  "157"

kenton <- 
  external_met[station_name == "kenton"]

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

# Optimze ESy -------------------------------------------------------------

esy_fun <- 
  as.function(list(wl = NULL,
                   min.esy = 1,
                   asym = FALSE,
                   bquote(if(asym){return(.(a))} else{pmax(min.esy, .(a) - (.(a) - .(b)) * exp (.(c) * wl))},
                          where = as.list(optim(par = list(a = 13, b = 2, c = 0.02,
                                                           sigma = 1, nu = 3),
                                                # method = "L-BFGS-B",
                                                # lower = c(MP = 0, MM = 0, MQ = 0, MPET = 1),
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
  esy_fun(asym = TRUE)

# Parametize Water Balance Coefs ------------------------------------------

ggplot(dat,
       aes(x = sample_date,
           y = wl_initial_cm)) +
  geom_line() +
  facet_wrap(~water_year,
             scales = "free_x")

param_train <- 
  dat[water_year == c(2012)]

(mod_opt <- 
    optimize_params(par = list(MP = 1.5, MQ = 0.5, ESY = 1, MPET = 1),
                    method = "L-BFGS-B",
                    lower = list(MP = 0, MQ = 0, ESY = 0, MPET = 0),
                    upper = list(MP = 10, MQ = 10, ESY = 20, MPET = 10)))

# Test Fit Model ----------------------------------------------------------

test <- 
  dat[water_year %in% 2016]

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
  wl_hat2[t] <-
    wl_hat2[t] + mod_opt$par[["MP"]] * M2[t] * gradient2[t]
  
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




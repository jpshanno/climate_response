# Derived ESY Method

# 1. Fit drawdown model of wl_initial_cm ~ quad(ytd_water_balance) for selected
#    year (largest drawdown year is best)
# 2. Use b1, b2 to calculate quad_prime() as an approximation of Esy
# 3. Fit asymptotic regression model Esy ~ a - (a - b) * exp (- c * wl)
# 4. Optimize MPET, MP, MM, and MQ for the training year
# 5. Validate optimized model agains test datasets

ggplot(data = data.frame(esy = quad_prime(WA, b1, b2)[1:which.min(WL)], 
                         wl= WL[1:which.min(WL)])) + 
  aes(x = wl, 
      y = esy) +
  geom_point() 

esy_mod <- 
  robustbase::nlrob(esy ~ (a - (a - b) * exp (- c * wl)),
                    start = c(a = 10, b = 0, c = -0.3), na.action = na.exclude,
                    data = data.frame(esy = quad_prime(WA, b1, b2)[1:which.min(WL)], wl= WL[1:which.min(WL)]))
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
                   bquote(.(a) - (.(a) - .(b)) * exp(-.(c) * wl),
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

mod_opt <- 
  optim(par = list(MP = 1.5, MPET = 1, MM = 1.5, MQ = 0.5),
        control = list(fnscale = -1),
        fn = 
          function(params){
            
            MPET <- 
              params[["MPET"]]
            
            MP <- 
              params[["MP"]]
            
            MM <- 
              params[["MM"]]
            
            MQ <- 
              params[["MQ"]]
            
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
                esy_fun(wl_hat[t-1])
              
              # Use net input to determine if water level increases or decreases
              if((P[t] + PET[t]) <= 0){
                
                # Water level drawdown = PET2, if P2 <= PET2 then it can be assumed to be
                # less than interception (not necessarily true, but works as a 
                # simplifying assumption)
                wl_hat[t] <-
                  wl_hat[t-1] + (MPET*PET[t]) * gradient[t]
                
              } else {
                
                # Water rise is P2 - PET2, which fits better, but need to rethink my
                # justification for this
                wl_hat[t] <-
                  wl_hat[t-1] + (MP*P[t]) * gradient[t]
                
              }
              
              # Directly add melt to water level. This should probably have some sort of
              # multiplier
              wl_hat[t] <-
                wl_hat[t] + MM * M[t] * gradient[t]
              
              # If WL is above spill point threshold then lose some to streamflow. 
              # This could probably be improved using the morphology models to determine
              # streamflow
              
              if(wl_hat[t-1] > max.wl){
                
                q_hat[t] <- 
                  MQ * (wl_hat[t-1] - max.wl)
                
                wl_hat[t] <- 
                  wl_hat[t] - q_hat[t]
              }
              
            }
            
            hydroGOF::md(wl_hat[!is.na(WL)], WL[!is.na(WL)])
          })


# Test Fit Model ----------------------------------------------------------

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
      wl_hat2[t-1] + (1.0405217*PET2[t]) * gradient2[t]
    
  } else {
    
    # Water rise is P2 - PET2, which fits better, but need to rethink my
    # justification for this
    wl_hat2[t] <-
      wl_hat2[t-1] + (1.1575809*P2[t]) * gradient2[t]
    
  }
  
  # Directly add melt to water level. This should probably have some sort of
  # multiplier
  wl_hat2[t] <-
    wl_hat2[t] + 0.8901357 * M2[t] * gradient2[t]
  
  # If WL is above spill point threshold then lose some to streamflow. 
  # This could probably be improved using the morphology models to determine
  # streamflow
  
  if(wl_hat2[t-1] > max.wl){
    
    q_hat2[t] <- 
      0.1219465 * (wl_hat2[t-1] - max.wl)
    
    wl_hat2[t] <- 
      wl_hat2[t] - q_hat2[t]
  }
  
}

lplot(wl_hat2 ~ D2); lines(WL2 ~ D2, col = 'gray40')
hydroGOF::md(wl_hat2[!is.na(WL2)], WL2[!is.na(WL2)])
# lplot(wa_hat ~ D2)
# lplot(q_hat2 ~ D2)
lplot(g_hat2 ~ D2)
lplot(s_hat ~ D2)




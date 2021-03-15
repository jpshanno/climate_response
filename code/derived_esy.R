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
        wl_hat[t] + m_hat[t]
      
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
                 control = list(fnscale = -1, 
                                maxit = 2000),
                 fn = 
                   function(params, fixed = NULL){
                     
                     params <- 
                       c(params, fixed)
                     
                     wl_hat <- 
                       wetland_model(data = data, params)$wl_hat
                     
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

kenton <- 
  kenton[water_year %in% 2012:2019 & format(sample_date, "%m%d") != "0229",
         .(sample_date, 
           water_year,
           pet_cm = -pet_cm,
           rain_cm,
           melt_cm,
           snow_cm,
           total_input_cm,
           tmax_c,
           tmin_c,
           tmean_c)
  ]

# Drop anomalous melt
kenton[between(melt_cm, 0, 0.1) & month(sample_date) %in% 6:9, 
      `:=`(total_input_cm = total_input_cm - melt_cm,
           melt_cm = 0)]

# Drop October snowfall to avoid melt contributions to wetland water level 
# rebound
# kenton[month(sample_date) %in% c(9, 10),
#        `:=`(total_input_cm = total_input_cm + snow_cm - melt_cm,
#             rain_cm = rain_cm + snow_cm,
#             melt_cm = 0,
#             snow_cm = 0)]

# Drop anomalous pet
kenton[tmax_c <= 0,
      pet_cm := 0]

# Calculate YTD after dropping above points
kenton[, 
      `:=`(ytd_pet_cm = cumsum(pet_cm),
           ytd_p_cm = cumsum(rain_cm),
           ytd_m_cm = cumsum(melt_cm),
           ytd_input_cm = cumsum(rain_cm + melt_cm)),
      by = .(water_year)]

kenton[, ytd_water_balance := ytd_pet_cm + ytd_p_cm + ytd_m_cm]
kenton[, dWB := c(0, diff(ytd_water_balance))]

train <-
  kenton[water_year == 2012][
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


# NEED MODELS FOR SITES WITHOUT 2012 DATA
# CHECK FOR TREATEMENT SITES LISTED AS CONTORL IN 2014

test <- 
  kenton[CJ(site = unique(train$site),
            sample_date = kenton$sample_date),
         on = "sample_date"]

test[water_budget,
     `:=`(site_status = i.site_status,
          wl_initial_cm = wl_initial_cm,
          best_precip_cm = best_precip_cm,
          idw_precip_cm = idw_precip_cm),
     on = c("site", "sample_date")]  


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
      save_model = "output/models/esy_model.stan", 
      data = train[hydro_period == "drawdown" & !is.na(wl_initial_cm)])

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
                                        bquote(pmax(min.esy,
                                                    .(a) - (.(a) - .(b)) * exp (.(c) * wl)),
                                               where = .SD[1, .(a, b, c)])))),
            max_wl = density(na.omit(wl_initial_cm))$x[which.max(density(na.omit(wl_initial_cm))$y)]), 
        keyby = .(site)]

# 119 isn't great (has an early peak in spring). Could also work to identify 
# better local min(wl_initial_)
ggplot(train[hydro_period == "drawdown"],
       aes(x = wl_initial_cm,
           y = esy_emp)) +
  geom_point(aes(color = site)) +
  geom_line(aes(y = esy_hat)) +
  geom_vline(data = esy_functions,
             aes(xintercept = max_wl)) +
  facet_wrap(~site)



# Optimize Control Period Parameters --------------------------------------

optim_res <- 
  train[, .(res = 
              list(optimize_params(
                data = .SD, 
                par = list(MPET = 1,
                           MP = 1.5,
                           MM = 1,
                           MQ = 0.5,
                           minESY = 1,
                           phiM = 0.9,
                           phiP = 0.5),
                fixed = list(maxWL = esy_functions[.BY[[1]], max_wl],
                             funESY = esy_functions[.BY[[1]], pred_fun[[1]]]),
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
                             minESY = Inf, 
                             phiM = 1-.Machine$double.neg.eps, 
                             MQ = 1-.Machine$double.neg.eps, 
                             phiP = 1-.Machine$double.neg.eps)))),
keyby = .(site)]


optim_res[, modified_index_of_agreement := map(res, pluck, "value")]
optim_res[, params := map(res, pluck, "par")]
optim_res[, names(optim_res$res[[1]]$par) := map_dfr(res, ~ as.data.table(.x[["par"]]))]

train <- 
  train[, cbind(.SD, wetland_model(.SD, optim_res[.BY[[1]], params[[1]]])),
      by = .(site)]

ggplot(train,
       aes(x = sample_date)) +
  geom_line(aes(y = wl_initial_cm,
                color = "Observed",
                linetype = "Observed")) +
  geom_line(aes(y = wl_hat,
                color = "Modeled",
                linetype = "Modeled")) +
  geom_text(data = optim_res,
             aes(x = as.Date("2012-04-01"),
                 y = maxWL,
                 label = sprintf("md: %1.2f", modified_index_of_aggreement)),
             vjust = 2,
             hjust = -0.1) +
  scale_color_manual(name = "legend",
                     values = c(Observed = "gray40",
                                Modeled = "black")) +
  scale_linetype_manual(name = "legend",
                        values = c(Observed = "dashed",
                                   Modeled = "solid")) +
  ylab("Water Level (cm)") +
  theme_classic() +
  theme(panel.grid.major.y = element_line(),
        legend.position = "top",
        axis.title.x = element_blank(),
        legend.title = element_blank()) +
  facet_wrap(~site,
             scales = "free")

ggplot(train,
       aes(x = sample_date)) +
  geom_col(aes(y = m_hat + p_hat,
               fill = "Precip")) +
  geom_col(aes(y = pet_hat - q_hat,
               fill = "Streamflow")) +
  geom_col(aes(y = pet_hat,
               fill = "ET")) +
  geom_col(aes(y = m_hat,
               fill = "Melt")) +
  scale_fill_manual(breaks = c("Precip", "Melt", "ET", "Streamflow"),
                    values = c(Precip = "#7aa49f", 
                               Melt = "#207972",
                               ET = "#bc956e",
                               Streamflow = "#a47138")) +
  ylab("Change in Water Level (cm)") +
  xlab(NULL) +
  guides(fill = guide_legend(title = "Contributing Source", ncol = 2)) +
  theme_classic() +
  theme(panel.grid.major.y = element_line(),
        legend.position = c(0.95, 0.05),
        legend.justification = c(1, 0)) +
  facet_wrap(~site,
             scales = "free")

# Test Models -------------------------------------------------------------

test_control <-
  test[site_status == "Control", 
       cbind(.SD, wetland_model(.SD, optim_res[.BY[[1]], params[[1]]])),
       by = .(site)]

ggplot(test_control[water_year == 2017],
       aes(x = sample_date)) +
  geom_line(aes(y = wl_initial_cm,
                color = "Observed",
                linetype = "Observed")) +
  geom_line(aes(y = wl_hat,
                color = "Modeled",
                linetype = "Modeled")) +
  scale_color_manual(name = 'legend',
                     values = c(Observed = "gray40",
                                Modeled = "black")) +
  scale_linetype_manual(name = 'legend',
                        values = c(Observed = "dashed",
                                   Modeled = "solid")) +
  ylab("Water Level (cm)") +
  theme_classic() +
  theme(panel.grid.major.y = element_line(),
        legend.position = "top",
        axis.title.x = element_blank(),
        legend.title = element_blank()) +
  facet_wrap(~site,
             scales = "free")

ggplot(test,
       aes(x = sample_date)) +
  geom_line(aes(y = wl_initial_cm,
                color = "Observed",
                linetype = "Observed")) +
  geom_line(aes(y = wl_hat,
                color = "Modeled",
                linetype = "Modeled")) +
  scale_color_manual(values = c(Observed = "gray40",
                                Modeled = "black")) +
  scale_linetype_manual(values = c(Observed = "dashed",
                                   Modeled = "solid")) +
  ylab("Water Level (cm)") +
  theme_classic() +
  theme(panel.grid.major.y = element_line(),
        legend.position = "top",
        axis.title.x = element_blank(),
        legend.title = element_blank()) +
  facet_wrap(~site,
             scales = "free")

metrics_control <- 
  test[!is.na(wl_initial_cm) & site_status == "Control",
       data.table(t(hydroGOF::gof(wl_hat, wl_initial_cm))),
       by = .(site, water_year)]

ggplot(test[water_year == 2013],
       aes(x = sample_date)) +
  geom_col(aes(y = m_hat + p_hat,
               fill = "Precip")) +
  geom_col(aes(y = pet_hat - q_hat,
               fill = "Streamflow")) +
  geom_col(aes(y = pet_hat,
               fill = "ET")) +
  geom_col(aes(y = m_hat,
               fill = "Melt")) +
  scale_fill_manual(breaks = c("Precip", "Melt", "ET", "Streamflow"),
                    values = c(Precip = "#7aa49f", 
                               Melt = "#207972",
                               ET = "#bc956e",
                               Streamflow = "#a47138")) +
  ylab("Change in Water Level (cm)") +
  xlab(NULL) +
  guides(fill = guide_legend(title = "Contributing Source", ncol = 2)) +
  theme_classic() +
  theme(panel.grid.major.y = element_line(),
        legend.position = c(0.95, 0.05),
        legend.justification = c(1, 0)) +
  facet_wrap(~site,
             scales = "free")

# Fit Treatment Conditions ------------------------------------------------

# This could be done in a different way that let's regrowth show up. Instead of
# fitting a new MPET, I could fit a % of MPET variable to multiple years and 
# see if it increases over time. Training the model on 2018 data, seems to 
# overestimate the drawdown of immediate post-treatment periods (2014 & 2015)
# indicating that I may be able to detect the increase in transpiration 
# following vegetation regrowth

test_treat <- 
  test[site_status == "Treated"]

ggplot(test_treat,
       aes(x = sample_date,
           y = wl_initial_cm,
           color = factor(water_year))) +
  geom_line() +
  facet_wrap(~site,
             scales = "free")

# Not worth using water year for this optimization because it's only the last 
# month of sample year 2017's data
train_treat <- 
  test_treat[water_year == 2014 & year(sample_date) == 2014]

optim_treat <- 
  train_treat[, .(res = 
              list(optimize_params(
                data = .SD, 
                par = list(MPET = optim_res[.BY[[1]], MPET]),
                fixed = list(MP = optim_res[.BY[[1]], MP],
                             MM = optim_res[.BY[[1]], MM],
                             MQ = optim_res[.BY[[1]], MQ],
                             minESY = optim_res[.BY[[1]], minESY],
                             phiM = optim_res[.BY[[1]], phiM],
                             phiP = optim_res[.BY[[1]], phiP],
                             maxWL = esy_functions[.BY[[1]], max_wl],
                             funESY = esy_functions[.BY[[1]], pred_fun[[1]]]),
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
                             minESY = Inf, 
                             phiM = 0.9, 
                             MQ = 0.9, 
                             phiP = 0.9)))),
        keyby = .(site)]


optim_treat[, modified_index_of_aggreement := map(res, pluck, "value")]
optim_treat[, params := map(res, pluck, "par")]
optim_treat[, names(optim_treat$res[[1]]$par) := map_dfr(res, ~ as.data.table(.x[["par"]]))]

train_treat <- 
  train_treat[, c("wl_hat", "q_hat", "m_hat", "p_hat", "pet_hat") := wetland_model(.SD, optim_treat[.BY[[1]], params[[1]]]),
        by = .(site)]

ggplot(train_treat,
       aes(x = sample_date)) +
  geom_line(aes(y = wl_initial_cm,
                color = "Observed",
                linetype = "Observed")) +
  geom_line(aes(y = wl_hat,
                color = "Modeled",
                linetype = "Modeled")) +
  geom_text(data = optim_treat,
            aes(x = as.Date("2014-04-01"),
                y = maxWL,
                label = sprintf("md: %1.2f", modified_index_of_aggreement)),
            vjust = 2,
            hjust = -0.1) +
  scale_color_manual(name = "legend",
                     values = c(Observed = "gray40",
                                Modeled = "black")) +
  scale_linetype_manual(name = "legend",
                        values = c(Observed = "dashed",
                                   Modeled = "solid")) +
  ylab("Water Level (cm)") +
  theme_classic() +
  theme(panel.grid.major.y = element_line(),
        legend.position = "top",
        axis.title.x = element_blank(),
        legend.title = element_blank()) +
  facet_wrap(~site,
             scales = "free")


# Test Treated Models -----------------------------------------------------

test_treat <- 
  test_treat[, c("wl_hat", "q_hat", "m_hat", "p_hat", "pet_hat") := wetland_model(.SD, optim_treat[.BY[[1]], params[[1]]]),
       by = .(site)]


ggplot(test_treat[water_year == 2015],
       aes(x = sample_date)) +
  geom_line(aes(y = wl_initial_cm,
                color = "Observed",
                linetype = "Observed")) +
  geom_line(aes(y = wl_hat,
                color = "Modeled",
                linetype = "Modeled")) +
  scale_color_manual(name = "legend",
                     values = c(Observed = "gray40",
                                Modeled = "black")) +
  scale_linetype_manual(name = "legend",
                        values = c(Observed = "dashed",
                                   Modeled = "solid")) +
  ylab("Water Level (cm)") +
  theme_classic() +
  theme(panel.grid.major.y = element_line(),
        legend.position = "top",
        axis.title.x = element_blank(),
        legend.title = element_blank()) +
  facet_wrap(~site,
             scales = "free")

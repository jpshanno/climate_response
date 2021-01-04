##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param data
##' @return A two component list for fixed and random effects coefficients and prediction functions
##' @author Joe Shannon
##' @export
model_climate_flow <- 
  function(data,
           flow.models) {
    
    # Save original column names to make dropping temporary columns easier
    orig_cols <- 
      copy(names(data))
    
    # Calculate expected morphology flow
    data[, morphology_flow_cm := flow.models[.BY[[1]], f_predict[[1]]](wl_initial_cm, outflow.only = FALSE),
         by = .(site)]
    
    # Determine residual flow explained by weather variation
    data[, climate_flow_cm := morphology_flow_cm - net_flow_cm]
    
    # Plots showing mixed-model structure
    # ggplot(data[best_precip_cm == 0 & between(climate_flow_cm, -7.5, 2.5)],
    #        aes(x = pet_cm,
    #            y = climate_flow_cm,
    #            color = site)) +
    #   geom_point(alpha = 0.25) +
    #   geom_smooth(method = lmrob,
    #               formula = y ~ x,
    #               method.args = list(setting = "KS2014")) +
    #   facet_wrap(~site_status,
    #              scales = "free")
    # 
    # ggplot(data[best_precip_cm == 0 & between(climate_flow_cm, -7.5, 2.5)],
    #        aes(x = water_availability_cm,
    #            y = climate_flow_cm,
    #            color = site)) +
    #   geom_point(alpha = 0.25) +
    #   geom_smooth(method = lmrob,
    #               formula = y ~ cos(2 * pi * x / diff(range(x))) + sin(2 * pi * x / diff(range(x))),
    #               method.args = list(setting = "KS2014")) +
    #   facet_wrap(~site_status,
    #              scales = "free")
    
    # Convert water availability to range spanning a total of 1 (not necessarily 0-1)
    data[, water_availability_decimal := water_availability_cm / diff(range(water_availability_cm)),
         by = .(site)]
    
    # Convert water availability to a radian expression
    data[, water_availability_rad := 2 * pi * water_availability_decimal]
    
    # Create site status models with site as a random effect
    # Use `method = "DASvar"` because "In .rlmerInit(lcall, pf, formula, data, method, rho.e, rho.b, rho.sigma.e,  :
    # Method 'DAStau' does not support blocks of size larger than 2. Falling back to method 'DASvar'."
    mods <- 
      data[best_precip_cm == 0,
           .(mod = list(rlmerRcpp(climate_flow_cm ~ pet_cm*(cos(water_availability_rad) + sin(water_availability_rad)) + (1 + cos(water_availability_rad) + sin(water_availability_rad) | site),
                                  data = .SD,
                                  method = "DASvar"))),
           keyby = .(site_status)]
    
    # # Model evaluation plots
    # # Predict using random and fixed effects
    # data[, `:=`(pred_random = predict(mods[CJ(.BY[[1]]), mod[[1]]],
    #                                   newdata = .SD),
    #             pred_fixed = predict(mods[CJ(.BY[[1]]), mod[[1]]],
    #                                  re.form = NA,
    #                                  newdata = .SD)),
    #      by = .(site_status)]
    #
    # ggplot(data[best_precip_cm == 0],
    #        aes(x = climate_flow_cm,
    #            y = pred_random)) +
    #   geom_point(shape = 20) +
    #   geom_abline() +
    #   facet_wrap(~site,
    #              scales = "free")
    # 
    # ggplot(data[best_precip_cm == 0],
    #        aes(x = pet_cm,
    #            y = pred_random - climate_flow_cm)) +
    #   geom_point(shape = 20) +
    #   facet_wrap(~site,
    #              scales = "free")
    # 
    # ggplot(data[best_precip_cm == 0],
    #        aes(x = water_availability_cm,
    #            y = pred_random - climate_flow_cm)) +
    #   geom_point(shape = 20) +
    #   facet_wrap(~site,
    #              scales = "free")
    
    # # Compare mixed models and site-level models
    # # Site-level model with no treatment status information
    # site_mods <-
    #   data[best_precip_cm == 0,
    #        .(mod_crossed = list(lmrob(climate_flow_cm ~ pet_cm*(water_availability_rad + cos(water_availability_rad) + sin(water_availability_rad)),
    #                                   setting = "KS2014",
    #                                   data = .SD))),
    #        keyby = .(site)]
    # 
    # data[, `:=`(pred_site = predict(site_mods[CJ(.BY[[1]]), mod_crossed[[1]]],
    #                                    newdata = .SD)),
    #      by = .(site)]
    # 
    # # Compare models graphically
    # ggplot(data[best_precip_cm == 0],
    #        aes(x = climate_flow_cm,
    #            y = pred_random)) +
    #   geom_point(aes(y = pred_site),
    #              color = "red",
    #              shape = 1) +
    #   geom_point(shape = 20) +
    #   geom_abline() +
    #   facet_wrap(~site,
    #              scales = "free")
    # 
    # # Compare models by RMSE
    # evals <-
    #   data[best_precip_cm == 0,
    #        lapply(.SD,
    #               hydroGOF::rmse,
    #               obs = climate_flow_cm,
    #               na.rm = TRUE),
    #        by = .(site, site_status),
    #        .SDcols = patterns("pred_")]
    # 
    # melt(evals,
    #      id.vars = c("site", "site_status"))[, .(RMSE = median(value)),
    #                                          by = .(variable)][order(RMSE)]
    # 
    # ggplot(melt(evals,
    #             id.vars = c("site", "site_status"),
    #             value.name = "RMSE"),
    #        aes(x = variable, y = RMSE)) +
    #   geom_boxplot() +
    #   geom_jitter(width = 0.15)
    
    
# Create prediction functions ---------------------------------------------
    
    # Fixed effects predictions:
    fixed_pred <-
      mods[, c(.SD, map_dfr(mod, lme4::fixef))]
    
    # Clean up coefficient names
    setnames(fixed_pred,
             c("site_status", "mod", "intercept", "m_pet_cm", "m_cos", 
               "m_sin", "m_pet_cos", "m_pet_sin"))
    
    # For fixed effects predictions water_availability_rad will be calculated
    # using the global min & max of observed water availability values
    # Calculate water_availability_range for converting to radians
    fixed_pred[, water_availability_range_cm := diff(range(data$water_availability_cm))]
    
    # Created prediction functions
    fixed_predict_functions <- 
      split(fixed_pred,
            by = "site_status") %>%
      map(~as.function(list(pet = NULL,
                            water.availability = NULL,
                            convert.to.radians = TRUE,
                            substitute({
                              
                              if(convert.to.radians){
                                water_rad <- 
                                  2 * pi * water.availability / water_availability_range_cm
                              } else {
                                water_rad <- 
                                  water.availability
                              }
                              
                              intercept + 
                                m_pet_cm * pet + 
                                m_cos * cos(water_rad) +
                                m_sin * sin(water_rad) +
                                m_pet_cos * pet * cos(water_rad) +
                                m_pet_sin * pet * sin(water_rad)},
                              env = .x))))
    
    # Add prediction functions to data.table
    fixed_pred[, f_predict := fixed_predict_functions]
    
    # Calculate predicted fixed-effect values with predict()
    true_predict <- 
      data[, predict(fixed_pred[.BY[[1]], mod[[1]]],
                     re.form = NA,
                     newdata = .SD),
           by = .(site_status)]$V1
      
    # Calculate predicted fixed-effect values with function
    test_predict <- 
      data[, fixed_pred[.BY[[1]], f_predict[[1]]](pet_cm, water_availability_rad, FALSE),
         by = .(site_status)]$V1
    
    # Check that manual and predict() predictions match
    if(isFALSE(all.equal(true_predict,test_predict))){
      stop("The fixed effects prediction function does not match the output from predict()")
    }
    
    # Random effects predictions:
    # Create data.table of site-site_status combinations
    random_pred <- 
      unique(data[, .(site, site_status)])
    
    # Carry over fixed effects and models from fixed_pred
    random_pred <- 
      random_pred[fixed_pred[, -c("f_predict", "water_availability_range_cm")], 
                  on = "site_status"]

    # Extract the random effects from the site_status models
    random_effects <- 
      mods[, map_dfr(mod, function(.x){
        df <- lme4::ranef(.x)
        data.table(site = row.names(df$site), 
                   df$site)}),
        by = .(site_status)]
    
    # Adjust the fixed-effects coefficients for random effects
    random_pred[random_effects,
                `:=`(intercept = intercept + `(Intercept)`,
                     m_cos = m_cos + `cos(water_availability_rad)`,
                     m_sin = m_sin + `sin(water_availability_rad)`),
                on = c("site", "site_status")]
    
    # Get min and max water level by site to calculate radians and check for out
    # of range predictions in climate simulations
    random_pred[data[, .(min_obs_water_availability_cm = min(water_availability_cm),
                         max_obs_water_availability_cm = max(water_availability_cm)),
                     by = .(site)],
                `:=`(min_obs_water_availability_cm = i.min_obs_water_availability_cm,
                     max_obs_water_availability_cm = i.max_obs_water_availability_cm),
                on = c("site")]
    
    # Set key to ensure ordering and look-up for prediction
    setkey(random_pred, "site", "site_status")
    
    # Create prediction functions
    random_predict_functions <- 
      split(random_pred,
            by = c("site", "site_status")) %>%
      map(~as.function(list(pet = NULL,
                            water.availability = NULL,
                            convert.to.radians = TRUE,
                            substitute({
                              
                              if(convert.to.radians){
                                water_rad <- 
                                  2 * pi * water.availability / (max_obs_water_availability_cm - min_obs_water_availability_cm)
                              } else {
                                water_rad <- 
                                  water.availability
                              }
                              
                              intercept + 
                                m_pet_cm * pet + 
                                m_cos * cos(water_rad) +
                                m_sin * sin(water_rad) +
                                m_pet_cos * pet * cos(water_rad) +
                                m_pet_sin * pet * sin(water_rad)},
                              env = .x))))
    
    # Add prediction functions back into data.table
    random_pred[, f_predict := random_predict_functions]
    
    # Calculate predicted random-effect values with predict()
    true_rand_predict <- 
      data[, predict(random_pred[CJ(.BY[[1]], .BY[[2]]), mod[[1]]],
                     newdata = cbind(.SD, site = .BY[[1]])),
           by = .(site, site_status)]$V1
    
    # Calculate predicted random-effect values with function
    test_rand_predict <- 
      data[, random_pred[CJ(.BY[[1]], .BY[[2]]), f_predict[[1]]](pet_cm, water_availability_rad, FALSE),
           by = .(site, site_status)]$V1
    
    # Check that manual and predict() predictions match for random predictions
    if(isFALSE(all.equal(true_rand_predict,test_rand_predict))){
      stop("The random effects prediction function does not match the output from predict()")
    }
    
    # Remove columns added to `data`
    data[,names(data)[!(names(data) %in% orig_cols)] := NULL]
    
    # Return a list with elements for fixed & full predictions
    list(fixed = fixed_pred,
         random = random_pred)
      }

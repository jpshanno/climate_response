##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param data
##' @return
##' @author Joe Shannon
##' @export
model_residual_flow <- 
  function(data,
           flow.models) {
    
    orig_cols <- 
      copy(names(data))
    
    data[, expected_flow_cm := flow.models[.BY[[1]], f_predict[[1]]](wl_initial_cm, outflow.only = FALSE),
         by = .(site)]
    
    data[, residual_flow_cm := expected_flow_cm - net_flow_cm]
    
    ggplot(data[best_precip_cm == 0 & between(residual_flow_cm, -7.5, 2.5)],
           aes(x = pet_cm,
               y = residual_flow_cm,
               color = site)) +
      geom_point(alpha = 0.25) +
      geom_smooth(method = lmrob,
                  formula = y ~ x,
                  method.args = list(setting = "KS2014")) +
      facet_wrap(~site_status,
                 scales = "free")
    
    ggplot(data[best_precip_cm == 0 & between(residual_flow_cm, -7.5, 2.5)],
           aes(x = water_availability_cm,
               y = residual_flow_cm,
               color = site)) +
      geom_point(alpha = 0.25) +
      geom_smooth(method = lmrob,
                  formula = y ~ cos(2 * pi * x / diff(range(x))) + sin(2 * pi * x / diff(range(x))),
                  method.args = list(setting = "KS2014")) +
      facet_wrap(~site_status,
                 scales = "free")
    
    data[, water_availability_m := water_availability_cm / 100]
    
    data[, water_availability_decimal := water_availability_cm / diff(range(water_availability_cm)),
         by = .(site_status)]
    
    data[, water_availability_rad := 2 * pi * water_availability_cm]
    
    mods <- 
      data[best_precip_cm == 0,
           .(mod_pet = list(lmrob(residual_flow_cm ~ pet_cm,
                                  setting = "KS2014",
                                  data = .SD)),
             mod_crossed = list(lmrob(residual_flow_cm ~ pet_cm*water_availability_m,
                                      setting = "KS2014",
                                      data = .SD)),
             mod_wa = list(lmrob(residual_flow_cm ~ cos(water_availability_rad) + sin(water_availability_rad),
                                 setting = "KS2014",
                                 data = .SD)),
             mod_add = list(lmrob(residual_flow_cm ~ pet_cm + cos(water_availability_rad) + sin(water_availability_rad),
                                  setting = "KS2014",
                                  data = .SD)),
             mod_int = list(lmrob(residual_flow_cm ~ pet_cm + pet_cm:water_availability_m,
                                  setting = "KS2014",
                                  data = .SD)),
             mod_no_int = list(lmrob(residual_flow_cm ~ pet_cm + cos(water_availability_rad) + sin(water_availability_rad),
                                     setting = "KS2014",
                                     data = .SD)),
             mod_full = list(lmrob(residual_flow_cm ~ pet_cm + pet_cm:water_availability_m + cos(water_availability_rad) + sin(water_availability_rad),
                                   setting = "KS2014",
                                   data = .SD))),
           keyby = .(site)]
    
    # pet*water_avilability_cm is a worse fit
    # Need to test this harmonic method
    mmods <- 
      data[best_precip_cm == 0,
           .(
             mmod = list(rlmerRcpp(residual_flow_cm ~ pet_cm + pet_cm:water_availability_m + (pet_cm + pet_cm:water_availability_m || site),
                              data = .SD)),
             mmod_harmonic = list(rlmerRcpp(residual_flow_cm ~ pet_cm + pet_cm:water_availability_m + cos(water_availability_rad) + sin(water_availability_rad) + (1 | site) + (0 + pet_cm | site) + (0 + cos(water_availability_rad) + sin(water_availability_rad) | site),
                                            data = .SD))),
           keyby = .(site_status)]
    
    data[, `:=`(pred_pet = predict(mods[CJ(.BY[[1]]), mod_pet[[1]]],
                                   newdata = .SD),
                pred_crossed = predict(mods[CJ(.BY[[1]]), mod_crossed[[1]]],
                                   newdata = .SD),
                pred_wa = predict(mods[CJ(.BY[[1]]), mod_wa[[1]]],
                                   newdata = .SD),
                pred_add = predict(mods[CJ(.BY[[1]]), mod_add[[1]]],
                                   newdata = .SD),
                pred_int = predict(mods[CJ(.BY[[1]]), mod_int[[1]]],
                                   newdata = .SD),
                pred_no_int = predict(mods[CJ(.BY[[1]]), mod_no_int[[1]]],
                                   newdata = .SD),
                pred_full = predict(mods[CJ(.BY[[1]]), mod_full[[1]]],
                                   newdata = .SD)),
         by = .(site)]
    
    data[, `:=`(pred_mmod_full = predict(mmods[CJ(.BY[[1]]), mmod[[1]]],
                                         newdata = .SD),
                pred_mmod_harmonic = predict(mmods[CJ(.BY[[1]]), mmod_harmonic[[1]]],
                                         newdata = .SD),
                pred_mmod_fixed = predict(mmods[CJ(.BY[[1]]), mmod[[1]]],
                                          re.form = NA,
                                          newdata = .SD)),
         by = .(site_status)]
    
    # ggplot(data[best_precip_cm == 0],
    #        aes(x = residual_flow_cm,
    #            y = pred_pet)) +
    #   geom_point() +
    #   geom_abline() +
    #   facet_wrap(~site,
    #              scales = "free")
    # 
    # ggplot(data[best_precip_cm == 0],
    #        aes(x = residual_flow_cm,
    #            y = pred_add)) +
    #   geom_point() +
    #   geom_abline() +
    #   facet_wrap(~site,
    #              scales = "free")
    # 
    # ggplot(data[best_precip_cm == 0],
    #        aes(x = residual_flow_cm,
    #            y = pred_wa)) +
    #   geom_point() +
    #   geom_abline() +
    #   facet_wrap(~site,
    #              scales = "free")
    # 
    # ggplot(data[best_precip_cm == 0],
    #        aes(x = residual_flow_cm,
    #            y = pred_int)) +
    #   geom_point() +
    #   geom_abline() +
    #   facet_wrap(~site,
    #              scales = "free")
    # 
    # ggplot(data[best_precip_cm == 0],
    #        aes(x = residual_flow_cm,
    #            y = pred_no_int)) +
    #   geom_point() +
    #   geom_abline() +
    #   facet_wrap(~site,
    #              scales = "free")
    # 
    # ggplot(data[best_precip_cm == 0],
    #        aes(x = residual_flow_cm,
    #            y = pred_full)) +
    #   geom_point() +
    #   geom_abline() +
    #   facet_wrap(~site,
    #              scales = "free")

# 
#     evals <-
#       data[best_precip_cm == 0,
#            lapply(.SD,
#                   hydroGOF::mae,
#                   obs = residual_flow_cm,
#                   na.rm = TRUE),
#            by = .(site),
#            .SDcols = patterns("pred_")]
# 
# 
#     ggplot(melt(evals, id.vars = "site"),
#            aes(x = variable, y = value)) +
#       geom_boxplot() +
#       geom_jitter(width = 0.15)
    # 
    # data[,
    #      `:=`(resid_pet = as.numeric(scale(pred_pet - residual_flow_cm)),
    #           resid_full = as.numeric(scale(pred_full - residual_flow_cm)),
    #           resid_int = as.numeric(scale(pred_int - residual_flow_cm))),
    #      by = .(site)]
    # 
    # 
    # ggplot(data[best_precip_cm == 0],
    #        aes(x = pet_cm,
    #            y = resid_full)) +
    #   geom_point() +
    #   geom_smooth(method = lmrob,
    #               formula = y ~ x,
    #               se = FALSE) +
    #   geom_hline(aes(yintercept = 0),
    #              color = "red") +
    #   facet_wrap(~site,
    #              scales = "free")
    # 
    # ggplot(data[best_precip_cm == 0],
    #        aes(x = water_availability_rad,
    #            y = resid_full)) +
    #   geom_point() +
    #   geom_smooth(method = lmrob,
    #               formula = y ~ x,
    #               se = FALSE) +
    #   geom_hline(aes(yintercept = 0),
    #              color = "red") +
    #   facet_wrap(~site,
    #              scales = "free")
    
    
    # Using pred_int is the simpler model and performance is really similar
    # between that and pred_full
    
    mods[, c("intercept", "slope", "interaction_slope") := map_dfr(mod_int, coef)]
    
    predict_functions <- 
      split(mods,
            by = "site") %>%
      map(~as.function(list(x = NULL,
                            substitute({
                              intercept + slope * x + interaction_slope * x},
                              env = .x))))
    
    mods[, f_predict := predict_functions]
    
    mods[, .(site, 
             mod = mod_int,
             intercept,
             slope,
             interaction_slope,
             f_predict)]
    
    # Remove columns added to `data`
    data[,names(data)[!(names(data) %in% orig_cols)] := NULL]
    
    mods
      }

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param data
##' @return
##' @author Joe Shannon
##' @export
model_net_flow <- 
  function(data) {
    
    # ggplot(data,
    #        aes(x = wl_initial_cm,
    #            y = Ds_cm + pet_cm - net_precip_cm - nafill(melt_cm, "const", 0))) +
    #   geom_point(alpha = 0.25) +
    #   geom_smooth(method = "loess",
    #               method.args = list(family = "gaussian"),
    #               span = 0.75,
    #               color = "red") +
    #   geom_smooth(method = "loess",
    #               method.args = list(family = "symmetric"),
    #               span = 0.75) +
    #   facet_wrap(~site,
    #              scales = "free")
    
    flow_mods <- 
      data[,
           .(mod = list(autoloess(net_flow_cm ~ wl_initial_cm, 
                                  data = .SD,
                                  model = TRUE,
                                  family = "gaussian",
                                  span = c(0.3, 1)))),
           keyby = .(site)]
    
    flow_mods[, span := map(mod, pluck, "pars", "span")]
    
    # ggplot(data,
    #        aes(x = wl_initial_cm,
    #            y = net_flow_cm)) +
    #   geom_point(alpha = 0.25) +
    #   geom_smooth(method = "loess",
    #               method.args = list(family = "gaussian"),
    #               span = 0.75,
    #               color = "red") +
    #   facet_wrap(~site,
    #              scales = "free")
    # 
    # data[, loess_pred := predict(flow_mods[.BY[[1]], mod[[1]]], newdata = .SD), 
    #      by = .(site)]
    # 
    # ggplot(data,
    #        aes(x = wl_initial_cm,
    #            y = net_flow_cm)) +
    #   geom_point(alpha = 0.25) +
    #   geom_line(aes(y = loess_pred),
    #             color = "red") +
    #   facet_wrap(~site,
    #              scales = "free")
    
    smoothed_dat <- 
      unique(data[!is.na(wl_initial_cm), .(site, wl_initial_cm = round(wl_initial_cm, 1))])
    
    smoothed_dat[,loess_pred := predict(flow_mods[.BY[[1]], mod[[1]]], newdata = .SD), 
                 by = .(site)]
    
    # Could use MARS from earth package as well
    
    seg_mods <- 
      smoothed_dat[!is.na(loess_pred), 
              .(lm_mod = list(lm(loess_pred ~ wl_initial_cm))),
              keyby = .(site)]
    
    seg_mods[, seg_mod := map(lm_mod, 
                              selgmented, 
                              seg.Z=~wl_initial_cm, 
                              type = "davies")]
    
    # data[, seg_pred := predict(seg_mods[.BY[[1]], seg_mod[[1]]], newdata = .SD),
    #      by = .(site)]
    # 
    # ggplot(data,
    #        aes(x = wl_initial_cm,
    #            y = net_flow_cm)) +
    #   geom_point(alpha = 0.25) +
    #   geom_line(aes(y = test_pred,
    #                 color = "loess")) +
    #   geom_line(aes(y = seg_pred,
    #                 color = "segmented")) +
    #   facet_wrap(~site,
    #              scales = "free")
    
    seg_coefs <- 
      seg_mods[, map_dfr(seg_mod, 
                         ~as.list(c(intercept(.x)$wl_initial_cm[, "Est."], 
                                    slope(.x)$wl_initial_cm[, "Est."], 
                                    .x$psi[, "Est."]))), 
               by = .(site)]
    
    setnames(seg_coefs,
             gsub("[\\(\\)\\.]+", "_", tolower(names(seg_coefs))))
    
    setnames(seg_coefs,
             gsub("(^_+|_+$)", "", names(seg_coefs)))

    # Could make this more flexible to account for models with more or less than
    # two break points by pasting together an expression from the names
    # grep("[0-9](?=_)", names(seg_coefs), perl = TRUE, value = TRUE)
    predict_functions <- 
      split(seg_coefs,
            by = "site") %>%
      map(~as.function(list(x = NULL,
                            outflow.only = NULL,
                            substitute({q <- fcase(x < psi1_wl_initial_cm, 
                                                   intercept1 + slope1 * x,
                                                   x >= psi1_wl_initial_cm & x < psi2_wl_initial_cm,
                                                   intercept2 + slope2 * x,
                                                   x >= psi2_wl_initial_cm,
                                                   intercept3 + slope3 * x)
                            if(!outflow.only){return(q)}
                            pmin(q, 0)},
                                       env = .x))))

    seg_mods[, f_predict := predict_functions]
    
    # Check of prediction function    
    # data[, 
    #      test_pred := seg_mods[.BY[[1]], f_predict[[1]]](wl_initial_cm, outflow.only = FALSE),
    #      by = .(site)]
    # 
    # all.equal(data$test_pred, data$seg_pred)
    # 
    # ggplot(data,
    #        aes(x = wl_initial_cm,
    #            y = net_flow_cm)) +
    #   geom_point(alpha = 0.25) +
    #   geom_line(aes(y = test_pred,
    #                 color = "manual")) +
    #   geom_line(aes(y = seg_pred,
    #                 color = "segmented"),
    #             linetype = "dashed") +
    #   facet_wrap(~site,
    #              scales = "free")
    
    seg_mods[, lm_mod := NULL]
    
    seg_mods
    
}

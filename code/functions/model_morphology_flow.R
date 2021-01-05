##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param data
##' @return
##' @author Joe Shannon
##' @export
model_morphology_flow <- 
  function(data) {
    
    # Create water level bins of equal width. n is set to be the same number as
    # rounding to 2 cm bins. Doing it this way provides a better distribution of
    # points to model net flow. Attempted doing it as rounding to nearest 
    # centimeter or 2 centimeter, or with doing equal number of observations 
    # within each bin. 
    
    data[, wl_bin := cut_interval(wl_initial_cm,
                                  n = uniqueN(round(wl_initial_cm, 0))/2), 
         by = .(site)]
    
    # Reduced data to median flow for each centimeter of water level to remove 
    # noise
    reduced_dat <- 
      data[net_flow_cm < 5, 
           .(net_flow_cm = median_na(net_flow_cm),
             wl_initial_cm = median_na(wl_initial_cm)), 
           by = .(site, wl_bin)]
    
    # Remove wl_bin from original data
    data[, wl_bin := NULL]
    
    # Remove NAs for modeling & flows > 5 cm as these all seem to occur at high
    # water levels and likely represent melt/rain on snow or missed precip 
    # events. Removing data from before May removes many of these too, but then 
    # the most extreme high water levels are missed
    
    reduced_dat <- 
      subset(x = reduced_dat,
             subset = 
               !is.na(net_flow_cm) & 
               !is.na(wl_initial_cm))
    
    mcp_mods <- 
      reduced_dat[, .(mod = list(mcp(list(net_flow_cm ~ wl_initial_cm,
                                          ~0 + wl_initial_cm + I(wl_initial_cm^2)),
                                     data = .SD))),
                  keyby = .(site)]
    
    mcp_mods <- 
      mcp_mods[, cbind(.SD,
                       map_dfr(mod,
                               ~dcast(setDT(fixef(.x)[, c("name", "mean")]), 
                                      . ~ name, 
                                      value.var = "mean")[, -c(".")]))]
    
    # Drop model sigma estimate
    mcp_mods[, sigma_1 := NULL]
    
    setnames(mcp_mods,
             c("cp_1", "int_1", "wl_initial_cm_1", "wl_initial_cm_2", "wl_initial_cm_2_E2"),
             c("changepoint", "intercept", "linear_slope_1", "linear_slope_2", 
               "quadratic_slope_2"))
    
    predict_functions <- 
      split(mcp_mods,
            by = "site") %>%
      map(~as.function(list(x = NULL,
                            outflow.only = NULL,
                            substitute({
                              q <- 
                                ifelse(x <= changepoint,
                                       intercept + linear_slope_1 * x,
                                       intercept + linear_slope_1 * (x - changepoint) + linear_slope_2 * (x - changepoint) + quadratic_slope_2 * (x - changepoint)^2)
                            if(!outflow.only){return(q)}
                            pmin(q, 0)},
                                       env = .x))))

    
    mcp_mods[, f_predict := predict_functions]
    
    # # Check of prediction function
    # # They are not exactly the same. I think the issue may be that the fixef()
    # # function gives mean coefficient estimates from the posterior distribution.
    # # The more skewed the posterior estimates, the more the manual prediction and
    # # the predict() prediction will differ.
    # 
    # reduced_dat[,
    #      func_pred := mcp_mods[.BY[[1]], f_predict[[1]]](wl_initial_cm, outflow.only = FALSE),
    #      by = .(site)]
    # 
    # reduced_dat[,
    #      predict_pred := predict(mcp_mods[.BY[[1]], mod[[1]]],
    #                           newdata = .SD)$predict,
    #      by = .(site)]
    # 
    # map2(mcp_mods$site,
    #      mcp_mods$mod,
    #      ~plot(.y, q_predict = TRUE) +
    #        ggtitle(.x) +
    #        geom_line(data = reduced_dat[site == .x],
    #                  aes(x = wl_initial_cm,
    #                      y = predict_pred),
    #                  linetype = "dashed",
    #                  color = "red") +
    #        geom_line(data = reduced_dat[site == .x],
    #                  aes(x = wl_initial_cm,
    #                      y = func_pred))) %>%
    #   reduce(`+`)
    # 
    # all.equal(reduced_dat$func_pred, reduced_dat$predict_pred)
    
    # Add Observed Water Level Range
    
    mcp_mods[data[, .(wl_range = list(range(wl_initial_cm, na.rm = TRUE))), 
                  by = .(site)],
             wl_range_cm := i.wl_range,
             on = "site"]
             
    mcp_mods
    
}

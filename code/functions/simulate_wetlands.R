##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param data
##' @param model.params
##' @return
##' @author Joe Shannon
##' @export
simulate_wetlands <- 
  function(data, model.params, out.path) {
    
    start_time <- Sys.time()
    
    if(fs::file_exists(out.path)){
      
      fs::file_delete(out.path)
      
      # writeLines doesn't automatically create connection to gz file, so using
      # fwrite (or any other function that does recognize and create gz connection)
      fwrite(list("simulation_id","site","site_status","gcm","scenario","simulation_date","wl_hat","q_hat","m_hat","p_hat","pet_hat","gradient"),
                 out.path)

    }
    
    # Memory is a problem with this. Limiting memory usage by going rowwise
    # through the site/treatment combinations. Going through the gcm/scenario 
    # combinations takes more memory (365 * 6 * 10000 * 34)
    
    simulation_intervals <-
      vector("list", nrow(model.params))
    
    simulation_proportions <- 
      vector("list", nrow(model.params))
    
    for(i in seq_len(nrow(model.params))){
      SITE <- model.params$site[i]
      STATUS <- model.params$site_status[i]
      PARAMS <- model.params$params[[i]]
      ASH_PERCENT <- model.params$future.forest.change[i]
      
      simulations <- 
        data[, 
             cbind(site = SITE,
                   site_status = STATUS,
                   gcm,
                   scenario,
                   simulation_date = .SD[["simulation_date"]],
                   wetland_model(.SD, params = PARAMS, future.forest.change = ASH_PERCENT)),
             by = .(simulation_id)]
      
      fwrite(simulations, 
             out.path,
             append = TRUE)
      
      cat(paste("Simulations for", SITE, STATUS, "complete\n"))
      
      cat(paste("Summarizing Simulations for", SITE, STATUS, "\n"))
      
      simulation_intervals[[i]] <- 
        simulations[,
                    summarize_doy(wl_hat),
                    keyby = .(site, site_status, gcm, scenario, simulation_date)]
      
      simulation_proportions[[i]] <- 
        simulations[,
                    .(prop_above_max_wl = sum(wl_hat >= PARAMS$maxWL) / .N,
                      prop_within_5cm_max_wl = sum(wl_hat >= (PARAMS$maxWL - 5)) / .N,
                      prop_above_0 = sum(wl_hat >= 0) / .N,
                      prop_above_neg_10 = sum(wl_hat >= -10) / .N,
                      prop_above_neg_25 = sum(wl_hat >= -25) / .N,
                      prop_above_neg_50 = sum(wl_hat >= -50) / .N,
                      prop_above_neg_100 = sum(wl_hat >= -100)  / .N),
                    keyby = .(site, site_status, gcm, scenario, simulation_date)]
      
      cat(paste("Simulation summary for", SITE, STATUS, "complete\n"))
      if("RPushbullet" %in% .packages(TRUE) && suppressWarnings(suppressMessages(RPushbullet::pbValidateConf()))){
        RPushbullet::pbPost("note",
                            paste(SITE, STATUS, "complete"))
      }
      
    }
    
    if("RPushbullet" %in% .packages(TRUE) && suppressWarnings(suppressMessages(RPushbullet::pbValidateConf()))){
      RPushbullet::pbPost("note",
                          paste("All Simulations complete", trunc(difftime(Sys.time(), start_time, units = 'mins')), "minutes."))
    }
    
    list(intervals = rbindlist(simulation_intervals),
         proportions = rbindlist(simulation_proportions))
    
    # ggplot(simulations[gcm == "ccsm4" & scenario == "rcp45"]) +
    #   aes(x = simulation_date,
    #       y = wl_hat) +
    #   stat_lineribbon()
  }


summarize_doy <- 
  function(x, interval.widths = c(0.5, 0.66, 0.95)){
    
    qis <- 
      data.table(ggdist::median_qi(x, 
                                   .width = interval.widths))
    
    hdcis <- 
      data.table(ggdist::median_hdci(x, 
                                     .width = interval.widths))
    
    intervals <- 
      data.table(median = median(x),
                 mean = mean(x),
                 mode = ggdist::Mode(x),
                 interval_name = c(sprintf("qi_%0.3f", c(0.5 * (1 - interval.widths), 0.5 * ( 1 +interval.widths))),
                                   sprintf("hdci_%0.3f", c(0.5 * (1 - interval.widths), 0.5 * ( 1 +interval.widths)))),
                 interval_value = c(qis$ymin, qis$ymax, hdcis$ymin, hdcis$ymax))
    
    dcast(intervals, 
          median + mean + mode ~ interval_name, 
          value.var = "interval_value")
  }

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param data
##' @param interception
##' @return
##' @author Joe Shannon
##' @export
model_ecosystem_sy <- 
  function(data, 
           interception,
           precip.col) {
    
    sy_dat <- 
      calculate_net_precip(data, 
                           interception, 
                           precip.col)[
                             net_precip_cm > 0, 
                             .(site, 
                               sample_date,
                               wl_min_cm, 
                               Dl_signed_cm, 
                               Dl_time_range_s, 
                               Ds_cm, 
                               site_status, 
                               obs_precip_cm,
                               idw_precip_cm,
                               best_precip_cm,
                               net_precip_cm,
                               i_cm)]
    
    sy_dat[, sy := best_precip_cm / Dl_signed_cm]
    
    sy_dat[, precip_intensity_cm_hr := best_precip_cm / (Dl_time_range_s / 3600)]
    
    # Low intensity storms have extremely low Sy, even at high water levels
    # sy_dat <-
    #   sy_dat[24*precip_intensity_cm_hr > i_cm]
    
    # Remove April (dominated by snowmelt not rain)
    sy_dat <- 
      sy_dat[month(sample_date) > 4]
    
    # Artificially high sy values (may be from well bottoming out)] 
    sy_dat <- 
      sy_dat[sy < 2]
    
    # Artificially high sy values (may be from well bottoming out)] 
    sy_dat <- 
      sy_dat[!(wl_min_cm < -60 & sy > 0.3)]
    
    # High water levels at two sites with no outflow
    sy_dat <- 
      sy_dat[!(site %in% c("119", "156") & wl_min_cm > 10)]
    
    # Drop remaining negative Esy values
    sy_dat <- 
      sy_dat[sy > 0]
    
    sy_mods <- 
      sy_dat[, 
             .(mod = list(nlrob(sy ~ b + m ** (wl_min_cm - c),
                                data = .SD,
                                start = list(b = min(.SD$sy, na.rm = TRUE), m = 1.25, c = -50),
                                algorithm = "port",
                                maxit = 100,
                                lower = list(b = min(.SD$sy, na.rm = TRUE), m = 0, c =-100000),
                                control = nls.control(warnOnly = TRUE, maxiter = 100)))),
             keyby = .(site)]
    
    sy_mods[, c("b", "m", "c") := map_dfr(mod, coef)]
    
    prediction_functions <- 
      split(sy_mods, 
            by = "site") %>% 
      map(~as.function(list(water.level = NULL, 
                            substitute(b + m ** (water.level - c), 
                                       env = .x))))
    
    sy_mods[, f_predict := prediction_functions]
    
    sy_mods
    
    }

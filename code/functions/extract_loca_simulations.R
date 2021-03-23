##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param out.path
##' @return
##' @author Joe Shannon
##' @export
##' 
##' Originally used tidync to get the data from the dodsC server, but it's much
##' faster to build the url and then download from the ncss server
##' https://cida.usgs.gov/thredds/catalog.html?dataset=cida.usgs.gov/loca_future
##' https://cida.usgs.gov/thredds/catalog.html?dataset=cida.usgs.gov/loca_historical
extract_loca_simulations <- 
  function(nc.files,
           coords) {
    
    nc <- 
      nc.files
    
    dat <- 
      data.table(file_name = nc)
    
    dat[, c("period", "variable") := tstrsplit(basename(file_name), "_")[c(1, 3)]]
    
    dat[, variable := str_remove(variable, ".nc")]
    
    dat[, nc := map(file_name, tidync)]
    dat[, start_date := ifelse(period == "future", 
                               as_loca_date("2070-01-01"),
                               as_loca_date("1980-01-01"))]
    dat[, end_date := ifelse(period == "future", 
                             as_loca_date("2099-12-31"),
                             as_loca_date("2009-12-31"))]
    dat[, projections := ifelse(period == "future",
                                list(c(1:4)),
                                list(c(1,3)))]
    dat[, projection_names := map(file_name,
                                  ~nc_att(.x, "NC_GLOBAL", "Projections") %>% 
                                    setDT() %>% 
                                    .[, value[[1]]] %>% 
                                    strsplit(", ") %>% 
                                    .[[1]])]
    dat[period == "present",
        projection_names := map(projection_names,
                                ~str_replace_all(.x, "rcp(45|85)", "historical"))]
    
    dat[, extracted := list(list(pmap_dfr(coords[, .(ID = station_name, LON = lon, LAT = lat)],
                                          function(ID, LON, LAT){
                                            nc[[1]] %>%
                                              hyper_tibble(select_var = variable,
                                                           time = between(time, start_date, end_date),
                                                           lon = index == which.min(abs(lon - LON)),
                                                           lat = index == which.min(abs(lat - LAT))) %>%
                                              setDT() %>% 
                                              set(x = ., j = "station_name", value = ID) %>% 
                                              setnames(variable, "value") %>% 
                                              set(x = ., j = "variable", value = variable) %>% 
                                              subset(projection %in% projections[[1]]) %>%
                                              set(x = ., j = "projection", value = projection_names[[1]][.$projection]) %>% 
                                              setcolorder(c("station_name", "lon", "lat", "projection", "time", "variable", "value"))
                                          }))), 
        by = .(file_name)]
    
    simulations <- 
      rbindlist(dat$extracted)
    
    simulations[, c("gcm", "model_run", "scenario") := tstrsplit(projection, "\\.")]
    
    simulations[variable == "pr", value := 0.1*value]
    simulations[, variable := fcase(variable == "pr", "precip_cm",
                                    variable == "tasmin", "tmin_c",
                                    variable == "tasmax", "tmax_c")]
    
    cast_names <- 
      names(simulations)[which(!(names(simulations) %in% c("value", "variable")))]
    
    cast_form <- 
      paste(paste(cast_names, collapse = "+"), "~ variable")
    
    simulations <- 
      dcast(simulations,
            cast_form, 
            value.var = "value")
    
    setnames(simulations, "time", "sample_date")
    simulations[, sample_date := from_loca_date(sample_date)]
    
    simulations
    
  }



as_loca_date <- 
  function(x){
    stopifnot(is.character(x))
    as.numeric(as.Date(x) - as.Date("1900-01-01"))
  }

# This needs to convert to character and then to date to avoid date values that
# are not integers. Could probably do it by truncating x or doing as.integer x
# but I haven't had a chance to experiment with that
from_loca_date <- 
  function(x){
    stopifnot(is.numeric(x))
    as.Date(as.character(as.Date("1900-01-01") + x))
  }

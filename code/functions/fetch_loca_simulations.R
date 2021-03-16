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
fetch_loca_simulations <- function(out.path) {

  require(tidync)
  require(ncmeta)
  require(furrr)
  require(purrr)
  require(dplyr)
  require(tidyr)
  require(stringr)
  require(readr)
  require(httr)
  
  # Get the variable names from the thredds server
  
  future_vars <- 
    tidync('http://cida.usgs.gov/thredds/dodsC/loca_future') %>% 
    hyper_vars() %>% 
    filter(str_detect(name, "_(CCSM4|GFDL-CM3)")) %>% 
    pull(name) %>% 
    as.list() %>% 
    set_names(rep("var", times = length(.)))
  
  historical_vars <- 
    tidync('http://cida.usgs.gov/thredds/dodsC/loca_historical') %>% 
    hyper_vars() %>% 
    filter(str_detect(name, "_(CCSM4|GFDL-CM3)")) %>% 
    pull(name) %>% 
    as.list() %>% 
    set_names(rep("var", times = length(.)))
  
  stations <- 
    targets::tar_read(external_met) %>% 
    distinct(station_name, lat, lon) %>% 
    mutate(lon = 360 + lon) %>% 
    as_tibble()
  
  extract_params <- 
    bind_rows(stations %>% mutate(target = "historical", vars = list(historical_vars), start_date = "1980-01-01", end_date = "2009-12-31"),
              stations %>% mutate(target = "future", vars = list(future_vars), start_date = "1980-01-01", end_date = "2009-12-31"),
              stations %>% mutate(target = "future", vars = list(future_vars), start_date = "2070-01-01", end_date = "2099-12-31")) %>% 
    mutate(ncss_url = pmap_chr(list(lat, lon, start_date, end_date, vars, target),
                               ~build_ncss_url(c(latitude = ..1,
                                                 longitude = ..2,
                                                 time_start = ..3,
                                                 time_end = ..4,
                                                 ..5),
                                               target = ..6)))
  
  plan(multisession,
       workers = 4)
  
  ncss_responses <- 
    future_imap(set_names(extract_params$ncss_url,
                          extract_params$station_name)[1],
                ~GET(.x, 
                     write_disk(tempfile(tmpdir = "tmp", fileext = ".html"), overwrite = TRUE)),
                .progress = TRUE)
  
  plan(sequential)
  
  problems <- 
    map_dfr(ncss_responses, http_status) %>% 
    filter(reason != "OK")
  
  if(nrow(problems) != 0){
    stop("There we download problems")
  }
  
  ncss_content <- 
    map(ncss_responses,
        content,
        type = "text")
  
  ncss_data <- 
    map(ncss_content,
        read_csv)
  
  
  
  %>% 
                  content(type = "text") %>% 
                  read_csv() %>% 
                  mutate(station_name = .y))
  
  
  
# tidync & dodSc approach
time_start <-
  nc_atts(loca_url, "time") %>%
  mutate(value = map_chr(value, ~.x)) %>%
  filter(name == "units") %>%
  pull(value) %>%
  str_extract("[0-9]{4}-[0-9]{2}-[0-9]{2}") %>%
  as.Date()

present_start <-
  as.numeric(as.Date("1980-01-01") - time_start)

present_end <-
  as.numeric(as.Date("2009-12-31") - time_start)

future_start <-
  as.numeric(as.Date("2070-01-01") - time_start)

future_end <-
  as.numeric(as.Date("2099-12-31") - time_start)

  # future_simulations <- 
  #   pmap_dfr(extract_params, 
  #            ~nc %>% 
  #              hyper_tibble(select_var = my_vars$name,
  #                           time = between(time, ..4, ..5),
  #                           lon = index %in% map_dbl(..3, ~which.min(abs(lon - .x))),
  #                           lat = index %in% map_dbl(..2, ~which.min(abs(lat - .x)))) %>% 
  #              mutate(station_name = ..1,
  #                     sample_date = time_start + time))
  # 
  # 
  historical_simulations <-
    pmap_dfr(stations,
             ~nc %>%
               hyper_tibble(select_var = my_vars$name,
                            time = between(time, present_start, present_end),
                            lon = index %in% map_dbl(..3, ~which.min(abs(lon - .x))),
                            lat = index %in% map_dbl(..2, ~which.min(abs(lat - .x)))) %>%
               mutate(station_name = ..1,
                      sample_date = time_start + time))
  
  
  simulations <- 
    simulations %>% 
    pivot_longer(cols = -c(station_name, time, starts_with("lon"), starts_with('lat')),
                 names_to = "variable")
  
  
  simulations %>% 
    extract(variable,
            into = c("variable", "unit"),
            regex = '^([A-z0-9_-]+)\\[unit="([A-z0-9]+)"\\]$') %>% 
    separate(variable,
             into = c("variable", "gcm", "run", "scenario"),
             sep = "_") %>% 
    mutate(variable = case_when(variable == "pr" ~ "precip_cm",
                                variable == "tasmax" ~ "tmax_c",
                                variable == "tasmin" ~ "tmin_c")) %>% 
    select(-unit) %>% 
    pivot_wider(names_from = c(variable, unit),
                )
  
  saveRDS(simulations, out.path)

}


build_ncss_url <- 
  function(query.vars, target){

    match.arg(target, c("future", "historical"))
    stopifnot(inherits(query.vars, "list"))
    
    base_url <- 
      ifelse(target == "future",
             "https://cida.usgs.gov/thredds/ncss/loca_future",
             "https://cida.usgs.gov/thredds/ncss/loca_historical")
      
    query_vars <- 
      c(query.vars,
        list(req="station",
             timeStride=1,
             addLatLon="true",
             accept="csv"))
    
    base_url <- 
      parse_url(base_url)
    
    base_url$query <- 
      query_vars
    
    build_url(base_url)
   
    # query_spatial <- 
    #   list(latitude=46.5345,
    #        longitude=-89.1839)
    # 
    # query_time <- 
    #   list(time_start="2070-01-01",
    #        time_end="2070-01-02")
    # 
    # query_vars <- 
    #   list(var = "pr_CCSM4_r6i1p1_rcp45",
    #        var = "pr_CCSM4_r6i1p1_rcp85",
    #        var = "pr_GFDL-CM3_r1i1p1_rcp45",
    #        var = "pr_GFDL-CM3_r1i1p1_rcp85",
    #        var = "tasmax_CCSM4_r6i1p1_rcp45",
    #        var = "tasmax_CCSM4_r6i1p1_rcp85",
    #        var = "tasmax_GFDL-CM3_r1i1p1_rcp45",
    #        var = "tasmax_GFDL-CM3_r1i1p1_rcp85",
    #        var = "tasmin_CCSM4_r6i1p1_rcp45",
    #        var = "tasmin_CCSM4_r6i1p1_rcp85",
    #        var = "tasmin_GFDL-CM3_r1i1p1_rcp45",
    #        var = "tasmin_GFDL-CM3_r1i1p1_rcp85")
  }

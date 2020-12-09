##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author Joe Shannon
##' @export
fetch_ncei_data <- 
  function(out.path, # Character vector for HPDN and GHCND files
           coords, # named numeric vector c(lat = , lon = )
           radius, # in km
           start.date,
           end.date,
           force.download = FALSE) {

    datasets <- 
      c("hpdn", "ghcnd")
    
    dat <- 
      data.table(dataset = datasets,
                 out_path = file.path(out.path, 
                                      datasets),
                 readme_url = c("https://www.ncei.noaa.gov/data/coop-hourly-precipitation/v2/doc/readme.csv.txt",
                                "https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/readme.txt"),
                 inventory_url = c("https://www.ncei.noaa.gov/data/coop-hourly-precipitation/v2/station-inventory/HPD_v02r02_stationinv_c20201129.csv",
                                   "https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/ghcnd-stations.txt"),
                 citation_text = c(paste0("Hourly Precipitation Data (HPD) Network, Version 2.r2, NOAA National Centers for Environmental Information. ",
                                          Sys.Date(), "."),
                                   paste0("Menne, M.J., I. Durre, B. Korzeniewski, S. McNeal, K. Thomas, X. Yin, S. Anthony, R. Ray, R.S. Vose, B.E.Gleason, and T.G. Houston, 2012: Global Historical Climatology Network - Daily (GHCN-Daily), Version ",
                                          grep("^[0-9]{1,}\\.[A-z0-9-]{1,}$", 
                                               scan("https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/ghcnd-version.txt", what = "character", quiet = TRUE),
                                               value = TRUE),
                                          ". NOAA National Climatic Data Center. http://doi.org/10.7289/V5D21VHZ ",
                                          Sys.Date(),
                                          ".")),
                 key = "dataset")
    
    # Create download directory if it doesn't exist
    walk(dat$out_path, ~if(!dir.exists(.x) | force.download){dir.create(.x)})
    
    # Download readme
    walk2(dat$out_path, 
          dat$readme_url, 
          ~if(!file.exists(file.path(.x, "readme.txt")) | force.download){
            download.file(.y, 
                          file.path(.x, "readme.txt"))})
    
    # Download station inventories
    walk2(dat$out_path, 
          dat$inventory_url, 
          ~if(!file.exists(file.path(.x, "station_inventory.txt")) | force.download){
            download.file(.y, 
                          file.path(.x, "station_inventory.txt"))})
    
    # Write Citation
    walk2(dat$out_path, 
          dat$citation_text, 
          ~if(!file.exists(file.path(.x, "citation.txt")) | force.download){
            writeLines(.y, 
                       file.path(.x, "citation.txt")) })
    
    # Read in HPDN stations
    hpdn <- 
      fread(file.path(dat["hpdn", out_path], "station_inventory.txt"))

    # Convert Period of Record to start and end dates    
    hpdn[, c("record_start", "record_end") := tstrsplit(POR_Date_Range, "-")]
    
    hpdn[, c("record_start", "record_end") := lapply(list(record_start, record_end),
                                                     ymd)]
    
    # Calculate distance from coordinates
    hpdn[, distance_km := geodist(hpdn, coords, measure = "geodesic")[, 1] / 1000]
    
    # Select bounded sites
    selected_stations <- 
      hpdn[distance_km < radius & record_start <= start.date & record_end >= end.date]
    
    # Download HPD Data
    
    # Create url for the CSV
    selected_stations[, hpdn_url := paste0("https://www.ncei.noaa.gov/data/coop-hourly-precipitation/v2/access/",
                                           StnID,
                                           ".csv")]
    
    # Clean up station names
    selected_stations[, station_name := gsub("[^A-z09]{1,}", "_", tolower(Name))]
    
    # Create the output file name
    selected_stations[, hpdn_file := file.path(dat["hpdn", out_path],
                                               paste0(station_name, ".csv"))]
    
    # Check for missing files to download
    missing_hpdn_files <- 
      !file.exists(selected_stations$hpdn_file) | force.download
    
    # Download any missing files
    if(sum(missing_hpdn_files) > 0){
      download.file(selected_stations$hpdn_url[missing_hpdn_files],
                    selected_stations$hpdn_file[missing_hpdn_files])
    }
    
    # Download GHCND Data -----------------------------------------------------
    # Not all HPD sites have GHCND data
    
    ghcnd_stations <- 
      read.fortran(file.path(dat["ghcnd", out_path], "station_inventory.txt"),
                   format = c( "A11", "F9", "F10", "F7", "X1","A2",
                               "X1","A30", "X1", "A3", "X1", "A3", "X1", "A5"), 
                   col.names = c("ID", "LAT", "LON", "ELEV", "ST", "NAME","GSN", "HCN", "WMOID"),
                   comment.char="")
    
    # Find sites withouth GHCND data
    selected_stations[, has_ghcnd := ifelse(StnID %in% ghcnd_stations$ID,
                                            TRUE,
                                            FALSE)]
    
    # Create the ghcnd url
    # These are available as CSVs too, which are daily rows (much tidier than .dly)
    # files, but the readme and version files are in the same folder as the .dly 
    # files, not the .csv files
    # https- https://www.ncei.noaa.gov/data/global-historical-climatology-network-daily/access/
    
    selected_stations[(has_ghcnd), 
                      ghcnd_url := paste0("https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/all/",
                                          StnID,
                                          ".dly")]
    
    # Create the output file name
    selected_stations[(has_ghcnd), 
                      ghcnd_file := file.path(dat["ghcnd", out_path],
                                              paste0(station_name, ".dly"))]
    
    # Check for missing files to download
    missing_ghcnd_files <- 
      (!file.exists(selected_stations$ghcnd_file) & selected_stations$has_ghcnd) | force.download
    
    # Download any missing files
    if(sum(missing_ghcnd_files) > 0){
      # Not working as a vector
      # 7: In download.file(selected_stations$ghcnd_url[missing_ghcnd_files],  :
      # cannot open URL 'ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily/all/USC00476939.dly': FTP status was '530 Not logged in'
      # for(i in selected_stations[(has_ghcnd) & missing_ghcnd_files, station_name]){
      download.file(selected_stations$ghcnd_url[missing_ghcnd_files],
                    selected_stations$ghcnd_file[missing_ghcnd_files])
      # }
    }

    list.files(out.path,
               pattern = "(csv|dly)$",
               recursive = TRUE,
               full.names = TRUE)
}


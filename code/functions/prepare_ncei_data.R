##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author Joe Shannon
##' @export
prepare_ncei_data <- 
  function(path, start.date, end.date, ghcnd.units) {
    
    hpdn_data <- 
      lapply(path[grepl("csv$", path)],
             prepare_hpdn_data,
             start.date = start.date,
             end.date = end.date) %>% 
      rbindlist()
             
    ghcnd_data <- 
      lapply(path[grepl("dly$", path)],
             prepare_ghcnd_data,
             start.date = start.date,
             end.date = end.date,
             ghcnd.units = ghcnd.units) %>% 
      rbindlist(fill = TRUE)
    
    hpdn_data[ghcnd_data[, .(station_name, sample_date, tmin_c, tmax_c)],
              `:=`(tmin_c = i.tmin_c,
                   tmax_c = i.tmax_c),
              on = c("station_name", "sample_date")]
    
    hpdn_data[, sample_year := year(sample_date)]
    
    # Remove days where tmax <= tmin as they are likely errors, but keep missing
    # values
    hpdn_data[tmax_c > tmin_c | is.na(tmin_c) | is.na(tmax_c)]  
  }

prepare_hpdn_data <- 
  function(x, start.date, end.date) {
    
    # Read each file
    dat <- 
      fread(x)
    
    # Expand the dataset to include all possible dates in range
    dat <- 
      dat[data.table(DATE = seq(start.date, end.date, by = 1)),
          on = c("DATE")]
    
    # Extract the daily summary columns
    daily_precip <- 
      dat[, c(list(sample_date = DATE),
              .SD),
          .SDcols = patterns("Dly")]
    
    # Convert precipitation from hundredths of an inch to cm
    daily_precip[, precip_cm := 2.54 * DlySum / 100]
    
    # Set -9999 values to NA
    daily_precip[DlySum == -9999,
                 precip_cm := NA_real_]
    
    # Set water year based on USGS approach (begins October 1)
    daily_precip[, water_year := as.water_year(sample_date, wy.month = 11)]
    
    # Set day of water year (dowy)
    daily_precip[, dowy := as.dowy(sample_date, wy.month = 11)]
    
    # Final output of daily data
    daily_precip <- 
      daily_precip[, .(sample_date,
                       water_year, 
                       dowy, 
                       precip_cm)]
    
    # Output data
    # na.exclude removes NAs induced by expanding data to include all 
    # timesteps
    data.table(station_name = sub(".csv", "", basename(x)), 
               station_id = unique(na.exclude(dat$STATION)), 
               lon = unique(na.exclude(dat$LONGITUDE)), 
               lat = unique(na.exclude(dat$LATITUDE)), 
               elevation_m = ifelse(unique(na.exclude(dat$ELEVATION)) == -999.9, 
                                    NA_real_, 
                                    unique(na.exclude(dat$ELEVATION))),
               daily_precip)
    
  }

prepare_ghcnd_data <- 
  function(x, start.date, end.date, ghcnd.units) {
    
    # Specify column classes and widths
    ghcnd_specs <- 
      c( "A11", "I4", "I2", "A4",
         rep( c( "I5", "A1", "A1", "A1"), 31))
    
    # Specify column names
    ghcnd_names <- 
      c("station_id",
        "year",
        "month",
        "element",
        paste0(c("value", "mflag", "qflag", "sflag"), rep(1:31, each = 4)))

    # Processing flags right now, but not saving them
    
    # Read in data see note above for source for formatting fortran-style
    # columns
    dat <- 
      setDT(read.fortran(file = x,
                         format = ghcnd_specs,
                         col.names = ghcnd_names,
                         na.strings = -9999))
    
    # Get columsn with flag information
    flags <- 
      melt(dat, 
           id.vars = c("element", "year", "month"),
           measure.vars = patterns("flag"),
           variable.name = "flag_type")
    
    # Get day of month from column name
    flags[, dom := gsub("[^0-9]", "", flag_type)]
    
    # Get flag type from column name
    flags[, flag_type := gsub("flag[0-9]{1,2}$", "_flag", flag_type)]
    
    # Create column for each flag_type
    flags <- 
      dcast(flags, 
            element + year + month + dom ~ flag_type,
            value.var = "value")
    
    # Get value columns for each day of the month
    values <- 
      melt(dat, 
           id.vars = c("element", "year", "month"),
           measure.vars = patterns("value"),
           variable.name = "dom")
    
    # Get day of month from column name
    values[, dom := gsub("[^0-9]", "", dom)]
    
    # Join flags and values
    ghcnd <- 
      values[flags, 
             on = c("element", "year", "month", "dom")]
    
    # Create sample date
    ghcnd[, sample_date := as.Date(paste(year, month, dom, sep = "-"),
                                   quiet = TRUE)]
    
    # Remove days that won't parse correction (eg Feb. 30). These are 
    # created because of the melting
    ghcnd <- 
      ghcnd[!is.na(sample_date)]
    
    # Convert value from an integer to a real to do unit conversions
    ghcnd[, value := as.numeric(value)]
    
    # Join units
    ghcnd[ghcnd.units,
          unit := i.units, 
          on = c("element")]
    
    # Convert any unit that is stored in tenths
    ghcnd[grepl("tenths", unit), 
          `:=`(value = 0.1 * value,
               unit = gsub(".* ([A-z]+$)", "\\1", unit))]
    
    # Convert mm to cm
    ghcnd[unit == "mm", 
          `:=`(value = 0.1 * value,
               unit = "cm")]
    
    # Convert 'C' and any others to lowercase
    ghcnd[, unit := tolower(unit)]
    
    ghcnd[, element := paste(tolower(element), unit, sep = "_")]
    
    ghcnd[, station_name := sub(".dly", "", basename(x))]
    
    ghcnd <- 
      dcast(ghcnd, 
            station_name + sample_date ~ element, 
            value.var = "value")
  }

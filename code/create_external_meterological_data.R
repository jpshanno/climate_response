# Download and format data from the NOAA NCEI Hourly Precipitation Data Network
# Limit sites by radius from a pair of coordinates and by bounding dates
# See the README for the dataset for more information
# https://www.ncei.noaa.gov/data/coop-hourly-precipitation/v2/doc/readme.csv.txt

# The script will download the CSVs limited spatially or temporally and will return
# three datatables containing the hourly dataset, daily dataset, and the station 
# coordinates and elevation.
# Downloaded data are saved to the provided directory using the station name,
# not station ID
# This script does not retain any data measurement flags in the output

# Formatting for reading in GHCND fixed with data from 
# http://spatialreasoning.com/wp/20170307_1244_r-reading-filtering-weather-data-from-the-global-historical-climatology-network-ghcn
# Would also work with read.fwf setting colum widths and column classes 
# separately

library(data.table)
library(lubridate)
library(geodist)

# Set geographic and temporal bounds --------------------------------------

center_coords <- 
  data.table(lon = -89.61332, 
             lat = 46.43133)

radius_km <- 
  100

start_date <- 
  ymd("2011-10-01")

end_date <- 
  ymd("2020-09-30")

ncei_dir <- 
  "data/ncei_hourly_precipitation_data/"

ghcnd_dir <- 
  "data/ghcnd_weather_observations/"


# Download Metadata -------------------------------------------------------

# Create the hpd download directory if it doesn't exist
if(!dir.exists(ncei_dir)){
  dir.create(ncei_dir)
}

if(!file.exists(paste0(ncei_dir, "readme.txt"))) {
  download.file("https://www.ncei.noaa.gov/data/coop-hourly-precipitation/v2/doc/readme.csv.txt",
                paste0(ncei_dir, "readme.txt"))
}

if(!file.exists(paste0(ncei_dir, "hpd_station_inventory.txt"))) {
  download.file("https://www.ncei.noaa.gov/data/coop-hourly-precipitation/v2/station-inventory/HPD_v02r02_stationinv_c20201129.csv",
                paste0(ncei_dir, "hpd_station_inventory.txt"))
}

if(!file.exists(paste0(ncei_dir, "citation.txt"))){
  hpd_citation <- 
    paste0("Hourly Precipitation Data (HPD) Network, Version 2.r2, NOAA National Centers for Environmental Information. ",
           Sys.Date(), ".")
  
  writeLines(hpd_citation, 
             paste0(ncei_dir, "citation.txt")) 
}


# Create the hpd download directory if it doesn't exist
if(!dir.exists(ghcnd_dir)){
  dir.create(ghcnd_dir)
}


if(!file.exists(paste0(ghcnd_dir, "ghcnd_stations.txt"))){
  download.file("https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/ghcnd-stations.txt",
                paste0(ghcnd_dir, "ghcnd_stations.txt"))
}

if(!file.exists(paste0(ghcnd_dir, "readme.txt"))) {
  download.file("https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/readme.txt",
                paste0(ghcnd_dir, "readme.txt"))
}

if(!file.exists(paste0(ghcnd_dir, "version.txt"))) {
  download.file("https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/ghcnd-version.txt",
                paste0(ghcnd_dir, "version.txt"))
}

if(!file.exists(paste0(ghcnd_dir, "citation.txt"))){
  ghcnd_citation <- 
    paste0("Menne, M.J., I. Durre, B. Korzeniewski, S. McNeal, K. Thomas, X. Yin, S. Anthony, R. Ray, R.S. Vose, B.E.Gleason, and T.G. Houston, 2012: Global Historical Climatology Network - Daily (GHCN-Daily), Version ",
           grep("^[0-9]{1,}\\.[A-z0-9-]{1,}$", 
                scan(paste0(ghcnd_dir, "version.txt"), what = "character", quiet = TRUE),
                value = TRUE),
           ". NOAA National Climatic Data Center. http://doi.org/10.7289/V5D21VHZ ",
           Sys.Date(),
           ".")
  
  writeLines(ghcnd_citation, 
             paste0(ghcnd_dir, "citation.txt")) 
}

# Select HPDN Stations ----------------------------------------------------

hpdn <- 
  fread("data/ncei_hourly_precipitation_data/hpd_station_inventory.txt")

hpdn[, c("record_start", "record_end") := tstrsplit(POR_Date_Range, "-")]

hpdn[, c("record_start", "record_end") := lapply(list(record_start, record_end),
                                                 ymd)]

hpdn[, distance_km := geodist(hpdn, center_coords, measure = "geodesic")[, 1] / 1000]


# # Visualze the Matching Sites
# # Requires leaflet
library(leaflet)
# leaflet(hpdn[distance_km < radius_km & record_start <= start_date & record_end >= end_date]) %>%
#   addCircles(lat = ~Lat,
#              lng = ~Lon,
#              popup = ~StnID) %>%
#   addMarkers(data = center_coords,
#              lat = ~lat,
#              lng = ~lon) %>%
#   addTiles()

# Select bounded sites
selected_stations <- 
  hpdn[distance_km < radius_km & record_start <= start_date & record_end >= end_date]

# Download HPD Data -------------------------------------------------------

# Create url for the CSV
selected_stations[, ncei_url := paste0("https://www.ncei.noaa.gov/data/coop-hourly-precipitation/v2/access/",
                             StnID,
                             ".csv")]

# Clean up station names
selected_stations[, station_name := gsub("[^A-z09]{1,}", "_", tolower(Name))]

# Create the output file name
selected_stations[, ncei_file := paste0(ncei_dir,
                                   station_name,
                                   ".csv")]

# Check for missing files to download
missing_ncei_files <- 
  !file.exists(selected_stations$ncei_file)

# Download any missing files
if(sum(missing_ncei_files) > 0){
  download.file(selected_stations$ncei_url[missing_ncei_files],
                selected_stations$ncei_file[missing_ncei_files])
}

# Download GHCND Data -----------------------------------------------------
# Not all HPD sites have GHCND data

ghcnd_stations <- 
  read.fortran("data/ghcnd_weather_observations/ghcnd_stations.txt",
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
                  ghcnd_file := paste0(ghcnd_dir,
                                         station_name,
                                         ".dly")]

# Check for missing files to download
missing_ghcnd_files <- 
  !file.exists(selected_stations$ghcnd_file) & selected_stations$has_ghcnd

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

# Load & Format HPD Data --------------------------------------------------

# KENTON HAS TWO ROWS WITH AN EXTRA COLUMN NEAR THE BEGINNING 1984-06-01 & 1984-09-01
# Manually adjusted

hpd_dat <- 
  lapply(selected_stations$ncei_file,
         function(.x){
          
           # Read each file
           dat <- 
             fread(.x)
           
           # Expand the dataset to include all possible dates in range
           dat <- 
             dat[data.table(DATE = seq(start_date, end_date, by = 1)),
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
           daily_precip[, water_year := ifelse(month(sample_date) > 9,
                                               year(sample_date) + 1,
                                               year(sample_date))]
           
           # Set day of water year (dowy)
           daily_precip[, dowy := as.numeric(1 + sample_date - min(sample_date)),
                        by = .(water_year)]
           
           # Final output of daily data
           daily_precip <- 
             daily_precip[, .(sample_date,
                              water_year, 
                              dowy, 
                              precip_cm)]
           
           # Take columns of hourly data and melt it into a long tidy format
           hourly_precip <- 
             melt(dat, 
                  id.vars = c("DATE"), 
                  measure.vars = patterns("HR[0-9]{2}Val"))

           # Create sample_time from date and hour column header
           hourly_precip[, sample_time := ymd_h(paste(DATE, gsub("[^0-9]", "", variable)))]
           
           # Expand to include all hourly timesteps
           hourly_precip <- 
             hourly_precip[data.table(sample_time = seq(ymd_hms(paste(start_date, "00:00:00")), ymd_hms(paste(end_date, "23:00:00")), by = 3600)),
                           on = c("sample_time")]
           
           # Convert precip from hundredths of an inch to cm
           hourly_precip[, precip_cm := 2.54 * value / 100]
            
           # Set water year based on USGS approach (begins October 1)
           hourly_precip[, water_year := ifelse(month(sample_time) > 9,
                                                year(sample_time) + 1,
                                                year(sample_time))]
           
           # Set -9999 values to NA
           hourly_precip[value == -9999,
                        precip_cm := NA_real_]

           # Set date-time of water year (dtowy)
           hourly_precip[, dtowy := as.numeric(difftime(sample_time, min(sample_time), units = "days")),
                         by = .(water_year)]
           
           # Final hourly output
           hourly_precip <- 
             hourly_precip[, .(sample_time, water_year, dtowy, precip_cm)]
           
           # Metadata
           # na.exclude removes NAs induced by expanding data to include all 
           # timesteps
           meta <- 
             data.table(station = unique(na.exclude(dat$STATION)), 
                        lon = unique(na.exclude(dat$LONGITUDE)), 
                        lat = unique(na.exclude(dat$LATITUDE)), 
                        elevation_m = ifelse(unique(na.exclude(dat$ELEVATION)) == -999.9, 
                                             NA_real_, 
                                             unique(na.exclude(dat$ELEVATION))))
           
           # Output hourly and daily data as a nested datatable
           data.table(station_name = sub(".csv", "", basename(.x)), 
                      metadata = list(meta), 
                      hourly_precip = list(hourly_precip), 
                      daily_precip = list(daily_precip))
           
         }) 
  
hpd_dat <- 
  rbindlist(hpd_dat)

# Unnest daily data
daily_hpd <- 
  hpd_dat[, daily_precip[[1]], by = .(station_name)]

# Unnest hourly data
hourly_hpd <- 
  hpd_dat[, hourly_precip[[1]], by = .(station_name)]

# Unnest station metadata separately
station_info <- 
  hpd_dat[, metadata[[1]], by = .(station_name)]

# Load & Format GHCND -----------------------------------------------------

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

# Elements units
# Taken from GHCND readme
ghcnd_units <- 
  fread("resources/ghcnd_element_codes.csv")


# Processing flags right now, but not saving them
ghcnd_dat <- 
  lapply(selected_stations[!is.na(ghcnd_file), ghcnd_file],
         function(.x){
           
           # Read in data see note above for source for formatting fortran-style
           # columns
           dat <- 
             setDT(read.fortran(file = .x,
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
                   ... ~ flag_type,
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
           ghcnd[, sample_date := ymd(paste(year, month, dom),
                                      quiet = TRUE)]
           
           # Remove days that won't parse correction (eg Feb. 30). These are 
           # created because of the melting
           ghcnd <- 
             ghcnd[!is.na(sample_date)]
           
           # Convert value from an integer to a real to do unit conversions
           ghcnd[, value := as.numeric(value)]
           
           # Join units
           ghcnd[ghcnd_units,
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
           
           ghcnd[, station_name := sub(".dly", "", basename(.x))]
           
           dcast(ghcnd, 
                 station_name + sample_date ~ element, 
                 value.var = "value")
           
         })
  
ghcnd_dat <- 
  rbindlist(ghcnd_dat, 
            fill = TRUE)

ghcnd_dat <- 
  ghcnd_dat[tmin_c < tmax_c]

# Create Weights for IDW --------------------------------------------------

mesowest_stations <- 
  fread("output/tabular/bristow-campbell_solar_radiation_coefficients.csv")

# Calculate distances between sites
station_distances <- 
  geodist(rbind(station_info, mesowest_stations, fill = TRUE),
          measure = "geodesic") / 1000

# Set row and column names to station names for easy lookup
dimnames(station_distances) <- 
  list(c(station_info$station_name, mesowest_stations$station_name),
       c(station_info$station_name, mesowest_stations$station_name))

# Cacluate IDW weights
station_weights <- 
  ifelse(station_distances == 0,
         0,
         1 / station_distances^3)


# Fill Missing Precip Data ------------------------------------------------
# Tested with GLM -> predicted precip -> weighted mean of predicted precip and
# saw slight decline in model metrics
# Create nested datatable to performed IDW for each station using all other 
# stations. This requires nesting all data except for a given station
idw <- 
  station_info[, 
               .(daily_dat = list(set(x = copy(daily_hpd[station_name != .BY[[1]]]),
                                      j = "site",
                                      value = .BY[[1]])),
                 hourly_dat = list(set(x = copy(hourly_hpd[station_name != .BY[[1]]]),
                                       j = "site",
                                       value = .BY[[1]]))), 
               by = .(station_name)]

idw[, daily_dat := 
      lapply(daily_dat,
             function(x){
               x$weights <- 
                 station_weights[, first(x$site)][x$station_name]
               x
             })]

idw[, hourly_dat := 
      lapply(hourly_dat,
             function(x){
               x$weights <- 
                 station_weights[, first(x$site)][x$station_name]
               x
             })]


idw[, daily_dat := 
      lapply(daily_dat, 
             function(x){
               x[, .(precip_idw_cm = weighted.mean(precip_cm, weights, na.rm = TRUE)), 
                 by = .(sample_date)]})]

idw[, hourly_dat := 
      lapply(hourly_dat, 
             function(x){
               x[, .(precip_idw_cm = weighted.mean(precip_cm, weights, na.rm = TRUE)), 
                 by = .(sample_time)]})]

daily_hpd[idw[, daily_dat[[1]], by = .(station_name)],
          precip_idw_cm := i.precip_idw_cm,
          on = c("station_name", "sample_date")]

hourly_hpd[idw[, hourly_dat[[1]], by = .(station_name)],
          precip_idw_cm := i.precip_idw_cm,
          on = c("station_name", "sample_time")]

daily_hpd[, lapply(.SD, modeval, measured = precip_cm, stat = "RMSE"),
          by = .(station_name),
          .SDcols = patterns("idw_cm")]

if(interactive()){
  library(ggplot2)
  ggplot(hourly_hpd,
         aes(x = precip_cm,
             y = precip_idw_cm)) +
    geom_point() +
    facet_wrap(~station_name)
}

if(interactive()){  
  ggplot(daily_hpd,
         aes(x = precip_cm,
             y = precip_idw_cm)) +
    geom_point() +
    facet_wrap(~station_name)
}

if(interactive()){
  ggplot(hourly_hpd[, .(dtowy, precip_cm, cum_precip_cm = cumsum(fcoalesce(precip_cm, precip_idw_cm))),
                    by = .(station_name, water_year)]) +
    geom_line(aes(x = dtowy,
                  y = cum_precip_cm,
                  color = station_name),
              alpha = 0.4) +
    geom_line(data = daily_hpd[, .(dowy, precip_cm, cum_precip_cm = cumsum(fcoalesce(precip_cm, precip_idw_cm))),
                               by = .(station_name, water_year)],
              aes(x = dowy,
                  y = cum_precip_cm,
                  color = station_name),
              linetype = "dashed") +
    facet_wrap(~water_year)
}


if(interactive()){
  ggplot(hourly_hpd[, .(dtowy, precip_cm, cum_precip_cm = cumsum(nafill(precip_cm, "const", 0))),
                    by = .(station_name, water_year)][is.na(precip_cm), cum_precip_cm := NA_real_]) +
    geom_line(aes(x = dtowy,
                  y = cum_precip_cm,
                  color = station_name),
              alpha = 0.4) +
    geom_line(data = daily_hpd[, .(dowy, precip_cm, cum_precip_cm = cumsum(nafill(precip_cm, "const", 0))),
                               by = .(station_name, water_year)][is.na(precip_cm), cum_precip_cm := NA_real_],
              aes(x = dowy,
                  y = cum_precip_cm,
                  color = station_name),
              linetype = "dashed") +
    facet_wrap(~water_year)
}


# Fill Missing Temperatures -----------------------------------------------

meso_dat <- 
  fread("output/tabular/mesowest_daily_temperature_and_radiation.csv")

meso_dat[, sample_date := as.Date(sample_date)]

temperature_dat <- 
  rbind(ghcnd_dat[station_info, 
                  .(station_name, sample_date, lat, lon, tmin_c, tmax_c), 
                  on = "station_name"],
        meso_dat[mesowest_stations,
                 .(station_name, sample_date, lat, lon, tmin_c, tmax_c), 
                 on = "station_name"])

if(interactive()){
  ggplot(data = temperature_dat,
         aes(x = tmax_c,
             y = ..scaled..)) +
    geom_density() +
    geom_density(color = NA,
                 alpha = 0.25,
                 aes(fill = fcase(month(sample_date) %in% 3:5, "mam",
                                  month(sample_date) %in% 6:8, "jja",
                                  month(sample_date) %in% 9:11, "son", 
                                  default = "djf")))  
}


idw_temp <- 
  station_info[, .(daily_dat = list(set(x = copy(temperature_dat[station_name != .BY[[1]]]),
                                        j = "site",
                                        value = .BY[[1]]))),
               by = .(station_name)]

idw_temp[, daily_dat := lapply(daily_dat, 
                               function(x){
                                 x$weights <- 
                                     station_weights[, first(x$site)][x$station_name]
                                 x})]

idw_temp[, daily_dat := lapply(daily_dat, 
                          function(x){x[, lapply(.SD, 
                                                 weighted.mean, 
                                                 w = weights,
                                                 na.rm = TRUE), 
                                        by = .(sample_date),
                                        .SDcols = patterns("t(min|max)_c")]})]

daily_hpd[ghcnd_dat,
          `:=`(obs_tmin_c = i.tmin_c,
               obs_tmax_c = i.tmax_c),
          on = c("station_name", "sample_date")]

daily_hpd[idw_temp[, daily_dat[[1]], by = .(station_name)],
          `:=`(idw_tmin_c = i.tmin_c,
               idw_tmax_c = i.tmax_c),
          on = c("station_name", "sample_date")]

if(interactive()){
  ggplot(daily_hpd,
         aes(x = obs_tmin_c,
             y = idw_tmin_c)) +
    geom_point() +
    geom_abline(color = "red") +
    facet_wrap(~station_name)
}

if(interactive()){
  ggplot(daily_hpd,
         aes(x = obs_tmax_c,
             y = idw_tmax_c)) +
    geom_point() +
    geom_abline(color = "red") +
    facet_wrap(~station_name)
}

if(interactive()){
  ggplot(daily_hpd[!is.na(obs_tmin_c)],
         aes(x = dowy,
             y = obs_tmin_c)) +
    geom_line() +
    geom_line(aes(y = idw_tmin_c),
              color = "red",
              size = rel(0.5)) +
    facet_grid(water_year~station_name)
}

# Create Final Product ----------------------------------------------------

daily_met <- 
  daily_hpd[, .(station_name,
                sample_date,
                water_year,
                dowy,
                precip_cm = fcoalesce(precip_cm, precip_idw_cm),
                tmin_c = fcoalesce(obs_tmin_c, idw_tmin_c),
                tmax_c = fcoalesce(obs_tmax_c, idw_tmax_c),
                idw_precip = is.na(precip_cm),
                idw_tmin = is.na(obs_tmin_c),
                idw_tmax = is.na(obs_tmax_c))]

fwrite(daily_met, "output/tabular/external_daily_meterology.csv")
fwrite(daily_hpd, "output/tabular/ncei_hpd_daily_precipitation.csv")
fwrite(hourly_hpd, "output/tabular/ncei_hpd_hourly_precipitation.csv")
fwrite(station_info, "output/tabular/ncei_hpd_site_info.csv")

# Precipitation Comparison ------------------------------------------------

ggplot(daily_met[ghcnd_dat, on = c("station_name", "sample_date"), nomatch = NULL][, .(dowy, 
                     cum_ghcnd = cumsum(nafill(prcp_cm, "const", 0)),
                     cum_precip_cm = cumsum(nafill(precip_cm, "const", 0))),
                 by = .(station_name, water_year),
                 ][station_name == "ironwood"],
       aes(x = dowy)) +
  geom_line(aes(y = cum_precip_cm,
                linetype = "hpd")) +
  geom_line(aes(y = cum_ghcnd,
                linetype = "ghcnd")) +
  scale_linetype_manual(values = c(hpd = "solid", 
                                   ghcnd = "dashed")) +
  facet_wrap(~water_year)

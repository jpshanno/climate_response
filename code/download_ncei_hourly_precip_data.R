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

library(data.table)
library(lubridate)
library(leaflet)
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

output_dir <- 
  "data/hourly_precipitation_data/"

hpdn <- 
  fread("data/HPD_station_inventory.csv")
# https://www.ncei.noaa.gov/data/coop-hourly-precipitation/v2/station-inventory/HPD_v02r02_stationinv_c20201129.csv

hpdn[, c("record_start", "record_end") := tstrsplit(POR_Date_Range, "-")]

hpdn[, c("record_start", "record_end") := lapply(list(record_start, record_end),
                                                 ymd)]

hpdn[, distance_km := geodist(hpdn, center_coords, measure = "geodesic")[, 1] / 1000]


# Visualze the Matching Sites ---------------------------------------------
# Requires leaflet
# leaflet(hpdn[distance_km < radius_km & record_start <= start_date & record_end >= end_date]) %>%
#   addCircles(lat = ~Lat,
#              lng = ~Lon,
#              popup = ~StnID) %>%
#   addMarkers(data = center_coords,
#              lat = ~lat,
#              lng = ~lon) %>%
#   addTiles()


# Download Raw Data -------------------------------------------------------

# Select bounded sites
selected_hpdn <- 
  hpdn[distance_km < radius_km & record_start <= start_date & record_end >= end_date]

# Create url for the CSV
selected_hpdn[, access_url := paste0("https://www.ncei.noaa.gov/data/coop-hourly-precipitation/v2/access/",
                             StnID,
                             ".csv")]

# Clean up station names
selected_hpdn[, station_name := gsub("[^A-z09]{1,}", "_", tolower(Name))]

# Create the output file name
selected_hpdn[, filename := paste0(output_dir,
                                   station_name,
                                   ".csv")]

# Create the download directory if it doesn't exist
if(!dir.exists(output_dir)){
  dir.create(output_dir)
}

# Check for existing files
existing_files <- 
  !file.exists(selected_hpdn$filename)

# Download any missing files
if(sum(existing_files) > 0){
  download.file(selected_hpdn$access_url[existing_files],
                selected_hpdn$filename[existing_files])
}


# Load & Format Data ------------------------------------------------------

# KENTON HAS TWO ROWS WITH AN EXTRA COLUMN NEAR THE BEGINNING 1984-06-01 & 1984-09-01
# Manually adjusted

hpd_dat <- 
  lapply(selected_hpdn$filename,
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
           
         }) %>% 
  rbindlist()

# Final Products ----------------------------------------------------------

# Unnest daily data
daily_hpd <- 
  hpd_dat[, daily_precip[[1]], by = .(station_name)]
  
# Unnest hourly data
hourly_hpd <- 
  hpd_dat[, hourly_precip[[1]], by = .(station_name)]

# Unnest station metadata separately
station_info <- 
  hpd_dat[, metadata[[1]], by = .(station_name)]

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


station_distances <- 
  geodist(station_info, 
          measure = "geodesic") / 1000

dimnames(station_distances) <- 
  list(station_info$station_name,
       station_info$station_name)

station_weights <- 
  ifelse(station_distances == 0,
         0,
         1 / station_distances^3)



# IDW Fill missing data ---------------------------------------------------
library(purrr)
idw <- 
  station_info[, 
               .(daily_dat = list(copy(daily_hpd[station_name != .BY[[1]]])[, site := .BY[[1]]]),
                 hourly_dat = list(copy(hourly_hpd[station_name != .BY[[1]]])[, site := .BY[[1]]])), 
               by = .(station_name)]

idw[, daily_dat := lapply(daily_dat, function(x){x[, weights := station_weights[, first(site)][x$station_name]]})]
idw[, hourly_dat := lapply(hourly_dat, function(x){x[, weights := station_weights[, first(site)][x$station_name]]})]

idw[, daily_dat := lapply(daily_dat, function(x){x[, .(precip_idw_cm = weighted.mean(precip_cm, weights, na.rm = TRUE)), by = .(sample_date)]})]
idw[, hourly_dat := lapply(hourly_dat, function(x){x[, .(precip_idw_cm = weighted.mean(precip_cm, weights, na.rm = TRUE)), by = .(sample_time)]})]

daily_hpd[idw[, daily_dat[[1]], by = .(station_name)],
          precip_idw_cm := i.precip_idw_cm,
          on = c("station_name", "sample_date")]

hourly_hpd[idw[, hourly_dat[[1]], by = .(station_name)],
          precip_idw_cm := i.precip_idw_cm,
          on = c("station_name", "sample_time")]

if(interactive()){
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

fwrite(daily_hpd, "output/tabular/ncei_hpd_daily_precipitation.csv")
fwrite(hourly_hpd, "output/tabular/ncei_hpd_hourly_precipitation.csv")
fwrite(station_info, "output/tabular/ncei_hpd_site_info.csv")

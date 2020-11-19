# Retreived From
# https://gdo-dcp.ucllnl.org/downscaled_cmip_projections/dcpInterface.html

# Want to switch to using geoknife
# Need to check metadata for missing value id and valid value range

library(stars)
library(purrr)
library(data.table)
library(ggplot2)

# Future Projections ------------------------------------------------------


loca_precip_2070 <- 
  read_ncdf("data/loca5_gui_GFDL-CM3/projections_2070-2099/Extraction_pr.nc",
            var = "pr")

loca_tmax_2070 <- 
  read_ncdf("data/loca5_gui_GFDL-CM3/projections_2070-2099/Extraction_tasmax.nc",
            var = "tasmax")

loca_tmin_2070 <- 
  read_ncdf("data/loca5_gui_GFDL-CM3/projections_2070-2099/Extraction_tasmin.nc",
            var = "tasmin")

# There are some NAs in the projections over Lake Superior, just remove the top
# 6 rows of the grid
# min(which(is.na(st_apply(loca_precip, c("lat", "lon"), sum)[[1]]), arr.ind = TRUE)[,"row"])

init_date <- 
  as.Date("2070-01-01")

set.seed(1234)

# Over Sample (200 each) so that I can drop duplicates and still have 100 pairs
sampled_rows <- 
  sample(1:min(which(is.na(st_apply(loca_precip_2070, c("lat", "lon"), sum)[[1]]), arr.ind = TRUE)[,"row"]),
         200,
         replace = TRUE)

sampled_columns <- 
  sample(1:dim(loca_precip_2070)[["lon"]],
         200,
         replace = TRUE)

sampled_points <- 
  data.table(lat = sampled_rows,
             lon = sampled_columns)

# Remove Duplicates
sampled_points <- 
  unique(sampled_points)

# Keep first 100
sampled_points <- 
  sampled_points[1:100]

future_sequences <- 
  map2_dfr(.x = sampled_points$lon,
           .y = sampled_points$lat,
           ~data.frame(sample_date = init_date + 0:(dim(loca_precip_2070)[["projection"]]-1),
                       lon = st_bbox(loca_precip_2070)[["xmin"]] + (.x-1)*0.0625 - 360,
                       lat = st_bbox(loca_precip_2070)[["ymin"]] + (.y-1)*0.0625,
                       precip = as.numeric(loca_precip_2070[1, .x, .y, ][[1]]),
                       tmin = as.numeric(loca_tmin_2070[1, .x, .y, ][[1]]),
                       tmax = as.numeric(loca_tmax_2070[1, .x, .y, ][[1]]))) %>% 
  setDT()

# Get DOY
future_sequences[, doy := as.numeric(format(sample_date, "%j"))]

# Assign Sequence ID
future_sequences[, seq_ID := rleid(year(sample_date), lon, lat)]

fwrite(future_sequences,
       "output/tabular/weather_sequences_gfdl-cm3_2070-2099.csv")



# Present Projections -----------------------------------------------------


# Future Projections ------------------------------------------------------


loca_precip_1980 <- 
  read_ncdf("data/loca5_gui_GFDL-CM3/projections_1980-2009/Extraction_pr.nc",
            var = "pr")

loca_tmax_1980 <- 
  read_ncdf("data/loca5_gui_GFDL-CM3/projections_1980-2009/Extraction_tasmax.nc",
            var = "tasmax")

loca_tmin_1980 <- 
  read_ncdf("data/loca5_gui_GFDL-CM3/projections_1980-2009/Extraction_tasmin.nc",
            var = "tasmin")

# There are some NAs in the projections over Lake Superior, just remove the top
# 6 rows of the grid
# min(which(is.na(st_apply(loca_precip, c("lat", "lon"), sum)[[1]]), arr.ind = TRUE)[,"row"])

present_init_date <- 
  as.Date("1980-01-01")

# Over Sample (200 each) so that I can drop duplicates and still have 100 pairs
sampled_rows <- 
  sample(1:min(which(is.na(st_apply(loca_precip_1980, c("lat", "lon"), sum)[[1]]), arr.ind = TRUE)[,"row"]),
         200,
         replace = TRUE)

sampled_columns <- 
  sample(1:dim(loca_precip_1980)[["lon"]],
         200,
         replace = TRUE)

sampled_points <- 
  data.table(lat = sampled_rows,
             lon = sampled_columns)

# Remove Duplicates
sampled_points <- 
  unique(sampled_points)

# Keep first 100
sampled_points <- 
  sampled_points[1:100]

present_sequences <- 
  map2_dfr(.x = sampled_points$lon,
           .y = sampled_points$lat,
           ~data.frame(sample_date = present_init_date + 0:(dim(loca_precip_1980)[["projection"]]-1),
                       lon = st_bbox(loca_precip_1980)[["xmin"]] + (.x-1)*0.0625 - 360,
                       lat = st_bbox(loca_precip_1980)[["ymin"]] + (.y-1)*0.0625,
                       precip = as.numeric(loca_precip_1980[1, .x, .y, ][[1]]),
                       tmin = as.numeric(loca_tmin_1980[1, .x, .y, ][[1]]),
                       tmax = as.numeric(loca_tmax_1980[1, .x, .y, ][[1]]))) %>% 
  setDT()

# Get DOY
present_sequences[, doy := yday(sample_date)]

# Assign Sequence ID
present_sequences[, seq_ID := rleid(year(sample_date), lon, lat)]

# Check for sanity
ggplot(data = present_sequences,
       aes(x = doy,
           y = tmin,
           color = factor(year(sample_date)))) +
  geom_line()

ggplot(data = present_sequences,
       aes(x = doy,
           y = tmax,
           color = factor(year(sample_date)))) +
  geom_line()

fwrite(present_sequences,
       "output/tabular/weather_sequences_gfdl-cm3_1980-2009.csv")

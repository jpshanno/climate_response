library(targets)
# This is an example target script.
# Read the tar_script() help file for details.

# Load functions
for(i in list.files("code/functions", full.names = TRUE)){
  source(i, verbose = FALSE)
  rm(i)
} 

# Set target-specific options such as packages.
tar_option_set(packages = c("lutz", 
                            "data.table", 
                            "purrr", 
                            "robustbase", 
                            "sirad", 
                            "lubridate",
                            "geodist"),
               format = "fst_dt",
               resources = list(compress = 100))

# Define targets
targets <- list(
  
  tar_target(
    mesowest_dir,
   "data/mesowest_met/",
    format = "file"
  )
  
  , tar_target(
    mesowest_met,
    prepare_mesowest_data(dir = mesowest_dir, 
                          start_date = as.Date("2011-10-01"),
                          end_date = as.Date("2020-09-30"))
  )
  
  , tar_target(
    solrad_coefs,
    calculate_solrad_coefs(data = mesowest_met)
  )
  
  , tar_target(
    # Elements units
    # Taken from GHCND readme
    ghcnd_units,
    fread("resources/ghcnd_element_codes.csv")
  )
  
  , tar_target(
    ncei_files,
    fetch_ncei_data(out.path = "data/ncei_data",
                    coords = c(lon = -89.61332,
                               lat = 46.43133),
                    radius = 100,
                    start.date = as.Date("2011-10-01"),
                    end.date = as.Date("2020-09-30"),
                    force.download = FALSE),
    format = "file"
  )
  
  
  # IDW does not currently include mesowest precip, could add in for daily data
  , tar_target(
    external_met,
    prepare_ncei_data(path = ncei_files,
                      start.date = as.Date("2011-10-01"),
                      end.date = as.Date("2020-09-30"),
                      ghcnd.units = ghcnd_units) %>%
      fill_missing_data(additional.data = mesowest_met,
                        source = "ncei") %>%
      calculate_solar_radiation(coefs = solrad_coefs)
  )
  
)

# End with a call to tar_pipeline() to wrangle the targets together.
# This target script must return a pipeline object.
tar_pipeline(targets)

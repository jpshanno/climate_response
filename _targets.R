library(targets)
# This is an example target script.
# Read the tar_script() help file for details.

# Load functions
for(i in list.files("code/functions", full.names = TRUE)){
  source(i, verbose = FALSE)
  rm(i)
} 

# Set target-specific options such as packages.
tar_option_set(packages = readLines("code/climate_packages.R"),
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
      calculate_mean_temp() %>% 
      calculate_solar_radiation(coefs = solrad_coefs) %>% 
      calculate_hargreaves_pet(lambda.MJ.kg = 2.45)
  )

  , tar_target(
    study_data_paths,
    c("data/well_levels.csv",
      "data/precipitation_daily.csv",
      "data/daily_snowmelt.csv"),
    format = "file"
  )
  
  , tar_target(
    daily_water_levels,
    prepare_water_levels(path = study_path("well_levels", study_data_paths)) %>% 
      append_study_precip(precip.path = study_path("precip", study_data_paths),
                          snow.path = study_path("snowmelt", study_data_paths))
  )
  
  , tar_target(
    mod_interception,
    model_interception_loss(data = daily_water_levels,
                            y = "Ds_cm",
                            x = c("best_precip_cm", "best_precip_cm:cos(doy_decimal*2*pi)"),
                            g = "site_status"),
    format = "rds"
  )

  , tar_target(
    mod_esy,
    model_ecosystem_sy(data = daily_water_levels,
                       interception = mod_interception),
    format = "rds"
  )
  # 
  # , tar_target(
  #   mod_flow,
  #   NULL
  # )
  # 
  # , tar_target(
  #   mod_rise,
  #   NULL
  # )
  # 
  # , tar_target(
  #   mod_drawdown,
  #   NULL
  # )
  # 
  # , tar_target(
  #   water_balance,
  #   NULL
  # )
)

# End with a call to tar_pipeline() to wrangle the targets together.
# This target script must return a pipeline object.
tar_pipeline(targets)

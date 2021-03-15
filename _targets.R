library(targets)
# This is an example target script.
# Read the tar_script() help file for details.

# Load functions
for(i in list.files("code/functions", full.names = TRUE)){
  source(i, verbose = FALSE)
  rm(i)
} 

# Set target-specific options such as packages.
tar_option_set(packages = readLines("code/climate_packages.txt"),
               format = "fst_dt",
               resources = list(compress = 100))

# Define targets
targets <- list(
  
  
  # Set Data Paths ----------------------------------------------------------
  
  tar_target(
    study_data_paths,
    c("data/well_levels.csv",
      "data/precipitation_daily.csv"),
    format = "file"
  )
  
  , tar_target(
    mesowest_dir,
    "data/mesowest_met/",
    format = "file"
  )
  
  # The station inventory for the HPD data is not a stable address, the date
  # changes. Need to update the function to account for that.
  , tar_target(
    ncei_files,
    fetch_ncei_data(out.path = "data/ncei_data",
                    coords = c(lon = -89.61332,
                               lat = 46.43133),
                    radius = 100,
                    start.date = as.Date("2005-11-01"),
                    end.date = as.Date("2020-10-30"),
                    force.download = FALSE),
    format = "file"
  )
  
  # Prepare Met Data ------------------------------------------------------
  
  , tar_target(
    mesowest_met,
    prepare_mesowest_data(dir = mesowest_dir, 
                          start.date = as.Date("2011-11-01"),
                          end.date = as.Date("2020-10-30"))
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
    external_met,
    prepare_ncei_data(path = ncei_files,
                      start.date = as.Date("2005-11-01"),
                      end.date = as.Date("2020-10-31"),
                      ghcnd.units = ghcnd_units) %>%
      fill_missing_data(additional.data = mesowest_met,
                        source = "ncei") %>%
      calculate_mean_temp() %>% 
      calculate_solar_radiation(coefs = solrad_coefs[station_name == "PIEM4"],
                                return.vector = FALSE) %>% 
      calculate_hargreaves_pet(lambda.MJ.kg = 2.45) %>% 
      subset(!(station_name %in% c("ironwood", "alberta_ford_for_center", "stambaugh_sse"))) %>% 
      calculate_snowmelt() %>% 
      calculate_water_availability()
  )
  
  # Prepare Water Budget Data ---------------------------------------------
  , tar_target(
    water_budget,
    prepare_water_levels(path = study_path("well_levels", study_data_paths)) %>% 
      append_study_precip(precip.path = study_path("precip", study_data_paths)) %>% 
      merge_external_met(external.met = external_met)
  )
  
  # Split Data --------------------------------------------------------------

  # , tar_target(
  #   validation_data,
  #   select_validation_data(data = water_budget,
  #                           n.in.group = 1,
  #                           groups = c("site", "site_status"))
  # )
  # 
  # , tar_target(
  #   training_data,
  #   water_budget[!validation_data]
  # )  

  # Build Water Level Models ------------------------------------------------

  # Make priors
    # Besy priors from min of empirical Sy
    # Mesy prior from papers?
    # 
  # Run Models
  # Evaluate Models
  
  # Simulate Wetland Runs ---------------------------------------------------
  
  # simulate_water_levels()
  # , tar_target(
  #   simulations,
  #   run_simulations(n, ...) %>% 
  #   simulate_water_levels(site_status, 
  #                       hydrology_model = NULL, 
  #                       initial.wl = NA,
  #                       weather, 
  #                       return.components = FALSE)

)

# End with a call to tar_pipeline() to wrangle the targets together.
# This target script must return a pipeline object.
tar_pipeline(targets)

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
      "data/precipitation_daily.csv",
      "data/daily_snowmelt.csv"),
    format = "file"
  )
  
  , tar_target(
    mesowest_dir,
    "data/mesowest_met/",
    format = "file"
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
  
  # Prepare Data ------------------------------------------------------------
  
  , tar_target(
    mesowest_met,
    prepare_mesowest_data(dir = mesowest_dir, 
                          start.date = as.Date("2011-10-01"),
                          end.date = as.Date("2020-09-30"))
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
                      start.date = as.Date("2011-10-01"),
                      end.date = as.Date("2020-09-30"),
                      ghcnd.units = ghcnd_units) %>%
      fill_missing_data(additional.data = mesowest_met,
                        source = "ncei") %>%
      calculate_mean_temp() %>% 
      calculate_solar_radiation(coefs = solrad_coefs) %>% 
      calculate_hargreaves_pet(lambda.MJ.kg = 2.45) %>% 
      calculate_water_availability(precip.col = "precip_cm",
                                   pet.col = "pet_cm",
                                   group.cols = c("station_name", "water_year"))
  )
  
  , tar_target(
    daily_water_levels,
    prepare_water_levels(path = study_path("well_levels", study_data_paths)) %>% 
      append_study_precip(precip.path = study_path("precip", study_data_paths),
                          snow.path = study_path("snowmelt", study_data_paths))
  )
  
  , tar_target(
    water_budget,
    prepare_water_budget(data = daily_water_levels,
                         external.met = external_met,
                         interception.mods = mod_interception,
                         sy.mods = mod_esy)
  )
  
  # Model Physical Components -----------------------------------------------
  
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
                       interception = mod_interception,
                       precip.col = "best_precip_cm"),
    format = "rds"
  )
  
  # Need formal model validations (save residual checks, etc)
  , tar_target(
    mod_morphology_flow,
    model_morphology_flow(data = water_budget),
    format = "rds"
  )
  
  # Model Response Components -----------------------------------------------
  # Split data
  
  , tar_target(
    validation_data,
    select_validation_years(data = water_budget,
                            date.col = "sample_date",
                            n.in.group = 1,
                            groups = c("site", "site_status"))
  )
  
  , tar_target(
    training_data,
    water_budget[!validation_data]
  )
  
  # Model Water level Rise

  , tar_target(
    mod_rise,
    model_precip_rise(training_data),
    format = "rds"
  )
    
  # Model Residual Flow
  # Flow not explained by water level fluctuations
  # Does not allow random slope, only random intercept. Should probably involve
  # both, but run times took forever
  , tar_target(
    mod_residual,
    model_residual_flow(training_data,
                        mod_morphology_flow),
    format = "rds"
  )
    
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

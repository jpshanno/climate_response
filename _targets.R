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
                    end.date = as.Date("2020-10-31"),
                    force.download = TRUE),
    format = "file"
  )
  
  # Prepare Met Data ------------------------------------------------------
  
  , tar_target(
    mesowest_met,
    prepare_mesowest_data(dir = mesowest_dir, 
                          start.date = as.Date("2011-11-01"),
                          end.date = as.Date("2020-10-31"))
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
      subset(site != "006") %>% 
      append_study_precip(precip.path = study_path("precip", study_data_paths)) %>% 
      merge_external_met(external.met = external_met)
  )
  
  # Split Data --------------------------------------------------------------
  # Training data for the control period for each site comes from the first 
  # availabile year of data. For most sites this is 2012 (which is ideal because
  # it has the largest drawdown period, which means it provides the best data
  # for ESy models). For treatment data just choosing the first post-treatment
  # year
  
  , tar_target(
    training_data,
    list(control = select_training_data(water_budget),
         treated = water_budget[water_budget[!is.na(wl_initial_cm) & site_status == "Treated", 
                                .(training_year = min(water_year)),
                                by = .(site)],
                                on = c("site", water_year = "training_year")]),
    format = "rds"
  )
  

  # Train Wetland Models ----------------------------------------------------

  , tar_target(
    esy_functions,
    build_esy_functions(training_data[["control"]]),
    format = "rds"
  )
  
  , tar_target(
    control_optimization,
    fit_models(training_data[["control"]],
               esy_functions,
               par = list(MPET = 1,
                          MP = 1.5,
                          MM = 1,
                          MQ = 0.5,
                          minESY = 1,
                          phiM = 0.9,
                          phiP = 0.5)),
    format = "rds"
  )
  
  , tar_target(
    treated_optimization,
    refit_model(training_data[["treated"]],
                control_optimization,
                refit = list(MPET = 1)),
    format = "rds"
  )
  
  , tar_target(
    model_params,
    setkey(rbind(control_optimization[, .(site, site_status = "Control", params)],
                 treated_optimization[, .(site, site_status = "Treated", params)]),
           site, site_status)[],
    format = "rds"
  )
 
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

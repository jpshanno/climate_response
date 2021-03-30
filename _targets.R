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
  
  # Paths to Study Field Data
  tar_target(
    study_data_paths,
    c("data/well_levels.csv",
      "data/precipitation_daily.csv"),
    format = "file"
  )
  
  # Paths to Mesowest Meterological Stations
  , tar_target(
    mesowest_dir,
    "data/mesowest_met/",
    format = "file"
  )
  
  # Paths to HPN & GHCND Stations
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

  # Path to Climate Normals
  , tar_target(
    climate_norms_file,
    "data/bergland_climate_normals.csv",
    format = "file")
    
  # Paths to LOCA Daily Projections
  , tar_target(
    loca_files,
    list.files("data/loca_ccsm4_gfdlcm3/", 
               pattern = "\\.nc$", 
               full.names = TRUE),
    format = "file"
  )
  
  # Prepare Data ----------------------------------------------------------
  
  # Create Climate Normals Dataset
  , tar_target(
    climate_norms,
    prepare_climate_norms(norms.file = climate_norms_file,
                          ncei.path = ncei_files,
                          ghcnd.units = ghcnd_units,
                          solrad.coefs = solrad_coefs[station_name == "PIEM4"]),
  )
  
  # Units for GHCND columns to automate conversion & standardization
  # taken from GHCND README
  , tar_target(
    ghcnd_units,
    fread("resources/ghcnd_element_codes.csv")
  )
  
  # Prepare GHCND data for use in SWG creation
  , tar_target(
    swg_data,
    prepare_ghcnd_swg_data(ncei.paths = ncei_files[grepl("\\.dly$", ncei_files)], 
                           ghcnd.units = ghcnd_units,
                           additional.stations = c("USC00200718", "USW00014858", "USC00208706", 
                                                   "USR0000MWAT", "USC00206220"))
  )
  
  # Combine and format Mesowest stations
  , tar_target(
    mesowest_met,
    prepare_mesowest_data(dir = mesowest_dir, 
                          start.date = as.Date("2011-11-01"),
                          end.date = as.Date("2020-10-31"))
  )
  
  # Calculate Bristow-Campbell Solar Radiation Coefficients
  , tar_target(
    solrad_coefs,
    calculate_solrad_coefs(data = mesowest_met)
  )

  # Create Combined External Met Data for Wetland Model Creation
  # The removed stations may be okay looking at the raw GHCND data, but they 
  # showed very strange trends in water availability, likely the result of an 
  # instrumentation change or a the backfilling process
  # For some reason ncei_data is not saving sample_date as Date, but as IDate
  # 250 is taken from 2010 climate normals @ Bergland Dam. Running 
  # CreateRunOptions() using Bergland Dam station 1980-2009 gives 341, but CMX
  # gives 245 and Ironwood gives 233
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
      calculate_snowmelt(mean.snowfall = 274,
                         cn.params = c(0.4, 2.9)) %>% 
      calculate_water_availability()
  )
  
  # Extract LOCA daily projections at various meteorological stations
  , tar_target(
    loca_simulations,
    extract_loca_simulations(nc.files = loca_files,
                             coords = unique(rbind(swg_data[, .SD[1, .(lon = 360 + lon, lat)], by = .(station_name)],
                                                   external_met[, .SD[1, .(lon = 360 + lon, lat)],  by = .(station_name)]),
                                             by = "station_name")) %>% 
    .[, swg_format_data(.SD), by = .(gcm, scenario)]
  )
  
  # Combine Study Data and External Met Data
  , tar_target(
    water_budget,
    prepare_water_levels(path = study_path("well_levels", study_data_paths)) %>% 
      subset(site != "006") %>% 
      set(i = which(.$water_year %in% c(2016, 2020)), j = "wl_initial_cm", value = NA_real_) %>% 
      append_study_precip(precip.path = study_path("precip", study_data_paths)) %>% 
      merge_external_met(external.met = external_met)
  )
  
  # Split Data --------------------------------------------------------------
  # Training data for the control period for each site comes from the first 
  # availabile year of data. For most sites this is 2012 (which is ideal because
  # it has the largest drawdown period, which means it provides the best data
  # for ESy models). The same process is used for the treatment period
  
  , tar_target(
    training_data,
    list(control = select_training_data(water_budget),
         treated = selected_treatment_training_data(water_budget)),
    format = "rds"
  )
  
  # Train Wetland Models ----------------------------------------------------

  # Create Esy models for each wetland
  , tar_target(
    esy_functions,
    build_esy_functions(training_data[["control"]]),
    format = "rds"
  )
  
  # Optimize Control Models
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
  
  # Reoptimize models for the control period
  , tar_target(
    treated_optimization,
    refit_model(training_data[["treated"]],
                control_optimization,
                refit = list(MPET = 1, 
                             MP = 1,
                             maxWL = 10)),
    format = "rds"
  )
  
  # Combine Control & Treated Model Parameters
  , tar_target(
    model_params,
    setkey(rbind(control_optimization[, .(site, site_status = "Control", params)],
                 treated_optimization[, .(site, site_status = "Treated", params)]),
           site, site_status)[],
    format = "rds"
  )

  # Evaluate Models ---------------------------------------------------------
  # This should probably be done breaking the water years up rather than doing 
  # them sequentially. Need to also consider adjusting max.wl for different 
  # water years within a given site.
  

  # Simulate Weather Series -------------------------------------------------

  # Create & Run SWG for Observed Data  
  , tar_target(
    swg_simulations_observed,
    swg_single_site(con.data = swg_data[station_name == "bergland_dam"],
                    n.simulations = 10000,
                    n.workers = 4,
                    solar.coefs = solrad_coefs[station_name == "PIEM4"],
                    simulation.dates = seq(as.Date("2008-11-01"), as.Date("2009-10-31"), by = "days")),
  )

  # Create & Run SWG for LOCA Data
  , tar_target(
    swg_simulations_loca,
    loca_simulations[station_name == "bergland_dam",
                     swg_single_site(.SD, 
                                     simulation.dates = seq(as.Date("2008-11-01"),
                                                            as.Date("2009-10-31"),
                                                            by = "days"),
                                     n.simulations = 10000,
                                     n.workers = 4,
                                     solar.coefs = solrad_coefs[station_name == "PIEM4"]),
                     by = .(gcm, scenario)]
  )
  
)

# End with a call to tar_pipeline() to wrangle the targets together.
# This target script must return a pipeline object.
tar_pipeline(targets)

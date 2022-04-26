library(targets)
# After the original creation of this project and writing the dissertation
# chapter a large model overhaul was conducted to improve the wetland models,
# including partial pooling using a Bayesian VI model. The data munging and
# other processing have not been adapated to the updated model, but rather the
# model inputs and outputs are maniuplated to match the existing pipeline. This
# may lead to some strange design choices below.

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
  
  # Path to raw wetland simulations
# , tar_target(
#   wetland_simulations_path,
#   "/run/media/jpshanno/Joe_Shannon/wetland_simulations.csv.gz",
#   format = "file"
# )
  
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
      merge_external_met(external.met = external_met) %>%
      # Adjusting site 053 water levels because of the well being moved
      .[site == "053" & sample_year < 2016, wl_initial_cm := wl_initial_cm + 16]
  )
  
  , tar_target(
    treatment_sites,
    water_budget[site_status == "Treated", unique(site)],
    format = "rds"
  )

  # Split Data --------------------------------------------------------------
  # Training data for the control period for each site comes from the first 
  # availabile year of data. For most sites this is 2012 (which is ideal because
  # it has the largest drawdown period, which means it provides the best data
  # for ESy models). The same process is used for the treatment period

  , tar_target(
    training_data,
    list(control = select_training_data(water_budget[site %in% treatment_sites]),
         treated = selected_treatment_training_data(water_budget[site %in% treatment_sites])),
    format = "rds"
  )

  , tar_target(
    testing_data,
    water_budget[!training_data[["control"]]][!training_data[["treated"]]][site %in% treatment_sites]
  )


  # Train Wetland Models ----------------------------------------------------

  # Create Esy models for each wetland
  , tar_target(
    esy_model,
    fit_esy_model(training_data[["control"]]),
    format = "rds"
  )

  , tar_target(
    esy_coefs,
    get_esy_coefs(esy_model),
  )

  , tar_target(
    esy_functions,
    build_esy_functions(esy_coefs),
    format = "rds"
  )

  , tar_target(
    wetland_model_code,
    "code/wetland_model.stan",
    format = "file"
  )

  , tar_target(
    wetland_model,
    fit_wetland_model(wetland_model_code, training_data, esy_coefs, "output/models/wetland_model.rds"),
    format = "file"
  )

  # Optimize Control Models
  # TODO: It it possible to set the upper limit MPET to 1. May likely require 
  # evaluating Esy within the optimization model?
  , tar_target(
    control_optimization,
    format_model_parameters("Control", wetland_model, esy_functions),
    format = "rds"
  )

  # Reoptimize models for the treatment period
  , tar_target(
    treated_optimization,
    format_model_parameters("Treated", wetland_model, esy_functions),
    format = "rds"
  )

  # Combine Control & Treated Model Parameters
  , tar_target(
    model_params,
    concatenate_model_params(Control = control_optimization, 
                             Treated = treated_optimization),
    format = "rds"
  )

  # Evaluate Models ---------------------------------------------------------
  # This should probably be done breaking the water years up rather than doing 
  # them sequentially. Need to also consider adjusting max.wl for different 
  # water years within a given site.

  , tar_target(
    train_data_fits,
    predict_train_period(data = training_data,
                         model.params = model_params),
    format = "rds"
  )

  , tar_target(
    test_data_fits,
    predict_test_period(data = testing_data[site %in% treatment_sites], 
                        model.params = model_params)
  )

  , tar_target(
    test_data_future_forest_fits,
    predict_test_period_future_forest(data = testing_data[site %in% treatment_sites], 
                        model.params = model_params)
  )

  , tar_target(
    wetland_model_metrics,
    calculate_wetland_model_metrics(data = test_data_fits, max_wl_data = esy_functions))

  , tar_target(
    wetland_model_metrics_plot,
    create_wetland_model_metrics_plot(
      data = wetland_model_metrics,
      metrics = c("r2", "med_err", "rmedse", "rmedse_range"),
      output.file = "output/figures/boxplot_wetland_model_metrics.tiff",
      type = "cairo",
      compression = "lzw",
      dpi = 600,
      width = 4,
      height = 6,
      bg = "white"
      ),
    format = "file")

  , tar_target(
    wetland_model_metrics_table,
    create_wetland_model_metrics_table(wetland_model_metrics),
    format = "rds"
  )

  , tar_target(
    wetland_model_predicted_probabilities,
    calculate_predicted_probabilities(test_data_fits, max_wl_data = esy_functions)
  )

  , tar_target(
    wetland_model_predicted_probs_tests,
    test_predicted_probabilities(wetland_model_predicted_probabilities)
  )

  , tar_target(
    wetland_model_predicted_probs_table,
    create_wetland_model_probability_table(wetland_model_predicted_probs_tests),
    format = "rds"
  )
  
  , tar_target(
    wetland_model_predicted_probs_plot,
    create_wetland_model_probability_plot(
      data = wetland_model_predicted_probabilities,
      tests = wetland_model_predicted_probs_tests,
      output.file = "output/figures/jitterpoint_and_pointrange_wetland_model_predicted_probability.tiff",
      type = "cairo",
      compression = "lzw",
      dpi = 600,
      width = 7.5,
      height = 3,
      bg = "white"),
    format = "file"
  )

  # Simulate Weather Series -------------------------------------------------

  # Create & Run SWG for Observed Data  
  , tar_target(
    swg_simulations_observed,
    swg_single_site(con.data = swg_data[station_name == "bergland_dam"],
                    n.simulations = 10000,
                    n.workers = 4,
                    solar.coefs = solrad_coefs[station_name == "PIEM4"],
                    simulation.dates = seq(as.Date("2009-01-01"), as.Date("2009-12-31"), by = "days")),
  )

  # Create & Run SWG for LOCA Data
  # Stopped using water year start and end because that really works better with
  # continuous data than with the simulation data. Had lots of days
  # (specifically) historical gfdl-cm3 in early water year (Nov-Dec) that had
  # low water levels. It may have to do with the seasonal conditioning of the
  # swg with Sept & Oct bringing November values up
  , tar_target(
    swg_simulations_loca,
    loca_simulations[station_name == "bergland_dam",
                     swg_single_site(.SD, 
                                     simulation.dates = seq(as.Date("2009-01-01"),
                                                            as.Date("2009-12-31"),
                                                            by = "days"),
                                     n.simulations = 10000,
                                     n.workers = 4,
                                     solar.coefs = solrad_coefs[station_name == "PIEM4"]),
                     by = .(gcm, scenario)]
  )

  # Evalute GCM & SWG ---------------------------------------------------------

    # TODO: Run and compare all 4 future climate scenarios. I can then 
    # empirically decide which ones are the less and more sensitive scenario.
    , tar_target(
    gcm_check_plot,
    create_gcm_check_plot(
      observed.data = tar_read(swg_data)[station_name == "bergland_dam"],
      gcm.data = tar_read(loca_simulations)[station_name == "bergland_dam" & scenario == "historical"],
      output.file = "output/figures/density_ridges_gcm_and_observed_climatology.tiff",
      type = "cairo",
      compression = "lzw",
      dpi = 600,
      width = 10,
      height = 7.5),
    format = "file"
  )
  
  , tar_target(
    swg_check_table,
    create_swg_table(
      swg.data = swg_simulations_loca[scenario == "historical"],
      loca.data = loca_simulations[station_name == "bergland_dam" & scenario == "historical"]
    ),
    format = "rds"
  )
  
  , tar_target(
    swg_check_plot,
    create_swg_plot(
      swg.data = swg_simulations_loca[scenario == "historical"],
      loca.data = loca_simulations[station_name == "bergland_dam" & scenario == "historical"],
      output.file = "output/figures/density_lines_loca_and_swg_climatology.tiff",
      type = "cairo",
      compression = "lzw",
      dpi = 600,
      width = 10,
      height = 6,
      bg = "white"
    ),
    format = "file"
  )

  # Run Wetland Models on Synthetic Weather ---------------------------------
  
  # TODO: Save output as individual disk.frames by treatment sharded by site and
  # site status
  , tar_target(
    wetland_simulation_summaries,
    simulate_wetlands(data = simplify_scenarios(swg_simulations_loca),
                      model.params = model_params[site %in% treatment_sites],
                      out.dir = "output/tabular/"),
    format = "rds"
  )


  , tar_target(
    pet_effects,
    summarize_pet(
      files = fs::dir_ls("output/tabular/", regexp = "wetland_simulations_(Control|Treated|Future_Forested)"),
      dummy = wetland_simulation_summaries
    )
  )

  , tar_target(
    analysis_simulations,
    map(wetland_simulation_summaries, ~ simplify_scenarios(.x[site %in% treatment_sites])),
    format = "rds"
  )

  # Evaluate EAB and CC Impacts ---------------------------------------------

  , tar_target(
    total_impact_plot,
    create_total_impact_plot(
      proportions = analysis_simulations[["proportions"]],
      output.file = "output/figures/pointrange_and_slabinterval_ecohydrological_level_total_impact.tiff",
      type = "cairo",
      compression = "lzw",
      dpi = 600,
      width = 7,
      height = 8,
      bg = "white"
    ),
    format = "file"
  )

  , tar_target(
    eab_impact_plot,
    create_eab_impact_plot(
      proportions = analysis_simulations[["proportions"]],
      output.file = "output/figures/pointrange_and_slabinterval_ecohydrological_level_eab_impact.tiff",
      type = "cairo",
      compression = "lzw",
      dpi = 600,
      width = 10,
      height = 5,
      bg = "white"
    ),
    format = "file"
  )

  , tar_target(
    climate_impact_plot,
    create_climate_impact_plot(
      proportions = analysis_simulations[["proportions"]],
      output.file = "output/figures/pointrange_and_slabinterval_ecohydrological_level_climate_impact.tiff",
      type = "cairo",
      compression = "lzw",
      dpi = 600,
      width = 10,
      height = 5,
      bg = "white"
    ),
    format = "file"
  )

  , tar_target(
    impact_attribution_plot,
    create_impact_attribution_plot(
      proportions = analysis_simulations[["proportions"]],
      output.file = "output/figures/pointrange_impact_attribution.tiff",
      type = "cairo",
      compression = "lzw",
      dpi = 600,
      width = 6,
      height = 8,
      bg = "white"
    ),
    format = "file"
  )

  , tar_target(
    impact_attribution_table,
    create_impact_attribution_table(analysis_simulations[["proportions"]]),
    format = "rds"
  )

  # Discussion Figures ------------------------------------------------------

  , tar_target(
    water_balance_plot,
    create_water_balance_plot(
      data = external_met,
      output.file = "output/figures/line_plot_study_period_water_availability.tiff",
      type = "cairo",
      compression = "lzw",
      dpi = 600,
      width = 6,
      height = 3,
      bg = "white"
    ),
    format = "file"
  )

  , tar_target(
    esy_impact_plot,
    create_esy_impact_plot(
      data = training_data[["control"]],
      esy_data = esy_functions,
      parameters = control_optimization,
      output.file = "output/figures/scatterplot_and_zoom_delta_water_level_esy_predictions.tiff",
      type = "cairo",
      compression = "lzw",
      dpi = 600,
      width = 6,
      height = 8,
      bg = "white"
    ),
    format = "file"
  )
    
  , tar_target(
    future_climate_plot,
    create_future_climate_plot(
      observed.data = swg_data[station_name == "bergland_dam"],
      gcm.data = simplify_scenarios(loca_simulations[station_name == "bergland_dam"]),
      solrad.coefs = solrad_coefs[station_name == "PIEM4"],
      output.file = "output/figures/density_ridges_future_and_observed_climatology.tiff",
      type = "cairo",
      compression = "lzw",
      dpi = 600,
      width = 10,
      height = 7.5,
      bg = "white"),
    format = "file"
  )

  , tar_target(
    pet_impact_plot,
    create_pet_impact_plot(
      data = pet_effects,
      output.file = "output/figures/lines_pet_demand_by_cover_and_climate.tiff",
      type = "cairo",
      compression = "lzw",
      dpi = 600,
      width = 8,
      height = 4,
      bg = "white"
    ),
    format = "file"
  )

  , tar_target(
    future_climate_table,
    create_future_climate_table(
      observed.data = swg_data[station_name == "bergland_dam"],
      gcm.data = simplify_scenarios(loca_simulations[station_name == "bergland_dam"]),
      solrad.coefs = solrad_coefs[station_name == "PIEM4"]
    ),
    format = "rds"
  )

)

# End with a call to tar_pipeline() to wrangle the targets together.
# This target script must return a pipeline object.
# tar_pipeline(targets)
targets

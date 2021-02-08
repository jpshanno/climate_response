source('code/load_project.R')

norms <- 
  fread("data/bergland_climate_normals.csv", na.strings = c("", "-7777"))

norms <- 
  setnafill(norms[, .(placeholder_date =as.Date(paste0("2020-", DATE)),
                      doy = yday(as.Date(paste0("2020-", DATE))),
                      lat = LATITUDE,
                      grdd_10 = `DLY-GRDD-BASE50`,
                      cldd_10 = `DLY-CLDD-BASE50`,
                      htdd_10 = `DLY-HTDD-BASE50`,
                      tmax_c = (0.1 * `DLY-TMAX-NORMAL` - 32) / 1.8,
                      tmean_c = (0.1 * `DLY-TAVG-NORMAL` - 32) / 1.8,
                      tmin_c = (0.1 * `DLY-TMIN-NORMAL` - 32) / 1.8,
                      ytd_precip_cm = `YTD-PRCP-NORMAL` / 39.37)], 
            "locf", 
            cols = "ytd_precip_cm")

bc_coefs <- 
  tar_read(solrad_coefs)[station_name == "PIEM4"]

norms[, et_solrad := extrat(doy, first(radians(lat)))$ExtraTerrestrialSolarRadiationDaily]
norms[, solrad_MJ_m2 := sirad::bc(days = placeholder_date, 
                                  lat = first(lat),
                                  BCb = bc_coefs$BC_B,
                                  extraT = et_solrad,
                                  Tmax = tmax_c,
                                  Tmin = tmin_c,
                                  BCc = bc_coefs$BC_C,
                                  tal = bc_coefs$BC_A)]
calculate_hargreaves_pet(norms, 2.45)
norms[order(doy), ytd_pet_cm := cumsum(pet_cm)]
norms[, dowy := ifelse((1:366) - 274 >= 1, (1:366) - 274, 366 + (1:366) - 274)]

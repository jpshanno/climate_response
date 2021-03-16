##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param data
##' @return
##' @author Joe Shannon
##' @export
calculate_snowmelt <- 
  function(data) {
    
    # as.POSIXct gets tripped up with converting date to POSIX without specifying
    # a time, it uses EDT even with tz="UTC" what so that 2020-10-31 becomes 2020-10-30.
    # as.POSIXlt avoids that
    cn_mods <- 
      data[, 
           .(n_records = .N,
             input_mod = list(CreateInputsModel(FUN_MOD = RunModel_CemaNeigeGR4J, 
                                                DatesR = as.POSIXlt(sample_date), 
                                                Precip = precip_cm * 10, 
                                                PotEvap = pmax(0, pet_cm * 10),
                                                TempMean = tmean_c,
                                                TempMin = tmin_c,
                                                TempMax = tmax_c,
                                                ZInputs = 380,
                                                HypsoData = rep(380, 101),
                                                NLayers = 1))),
           by = .(station_name)]
    
    cn_mods[,
            run_opts := map2(input_mod, n_records,
                             ~CreateRunOptions(FUN_MOD = RunModel_CemaNeigeGR4J, 
                                               InputsModel = .x,
                                               IndPeriod_WarmUp = 1:365,
                                               IndPeriod_Run = (365 + 1):.y,
                                               IsHyst = FALSE,
                                               MeanAnSolidPrecip = 250)
            )]
    
    cn_mods[, mod := map2(input_mod,
                          run_opts,
                          ~RunModel_CemaNeige(InputsModel = .x,
                                              RunOptions = .y,
                                              Param = c(0.25, 3.74)))]
    
    daily_melt <- 
      cn_mods[, map_dfr(mod, 
                        ~data.table(sample_date = as.Date(.x[["DatesR"]]),
                                    rain_cm = 0.1 * .x$CemaNeigeLayers$Layer01$Pliq,
                                    snow_cm = 0.1 * .x$CemaNeigeLayers$Layer01$Psol,
                                    swe_cm = 0.1 * .x$CemaNeigeLayers$Layer01$SnowPack,
                                    melt_cm = 0.1 * .x$CemaNeigeLayers$Layer01$Melt,
                                    total_input_cm = 0.1 * .x$CemaNeigeLayers$Layer01$PliqAndMelt)),
              by = .(station_name)]
    
    data[daily_melt,
         `:=`(rain_cm = i.rain_cm,
              snow_cm = i.snow_cm,
              swe_cm = i.swe_cm,
              melt_cm = i.melt_cm,
              total_input_cm = i.total_input_cm),
         on = c("station_name", "sample_date")]
    
    # The uncalibrated CN model shows anamolous melt in the early/mid summer 
    # (really long tail on melt)
    data[between(melt_cm, 0, 0.1) & month(sample_date) %in% 6:9, 
         `:=`(total_input_cm = total_input_cm - melt_cm,
              melt_cm = 0)]

    # Remove WY 2006 data as it has no CN model output    
    data[water_year != 2006]
  }

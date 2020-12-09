get_gmt_tzone <- 
  function(lat, lon){
    
    stopifnot(is.numeric(lat))
    stopifnot(is.numeric(lon))

    tz <- 
      tz_lookup_coords(lat = lat,
                       lon = lon,
                       method = "accurate")
    
    tzones <- 
      setDT(tz_list(), 
            key = "tz_name")[!(is_dst)]
    
    std_time_offset <- 
      tzones[tz, utc_offset_h]
    
    offset_sign <- 
      ifelse(std_time_offset > 0,
             "-",
             "+")
    
    paste0("Etc/GMT", offset_sign, abs(std_time_offset))
  }

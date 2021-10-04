##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param nameme1
##' @return
##' @author Joe Shannon
##' @export
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

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param nameme1
##' @return
##' @author Joe Shannon
##' @export
sum_na <- 
  function(x, class.out = class(x)){
    if(sum(is.na(x)) == length(x)){
      res <- NA
      class(res) <- class.out
      return(res)
    }
    
    sum(x, na.rm = TRUE)
  }

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param nameme1
##' @return
##' @author Joe Shannon
##' @export
mean_na <- 
  function(x, class.out = class(x)){
    if(sum(is.na(x)) == length(x)){
      res <- NA
      class(res) <- class.out
      return(res)
    }
    
    mean(x, na.rm = TRUE)
  }

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param nameme1
##' @return
##' @author Joe Shannon
##' @export
median_na <- 
  function(x, class.out = class(x)){
    if(sum(is.na(x)) == length(x)){
      res <- NA
      class(res) <- class.out
      return(res)
    }
    
    median(x, na.rm = TRUE)
  }

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param nameme1
##' @return
##' @author Joe Shannon
##' @export
max_na <- 
  function(x, class.out = class(x)){
    if(sum(is.na(x)) == length(x)){
      res <- NA
      class(res) <- class.out
      return(res)
    }
    
    max(x, na.rm = TRUE)
  }

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param nameme1
##' @return
##' @author Joe Shannon
##' @export
min_na <- 
  function(x, class.out = class(x)){
    if(sum(is.na(x)) == length(x)){
      res <- NA
      class(res) <- class.out
      return(res)
    }
    
    min(x, na.rm = TRUE)
  }

span_na <- 
  function(x, class.out = class(x)){
    if(sum(is.na(x)) == length(x)){
      res <- NA
      class(res) <- class.out
      return(res)
    }
    
    diff(range(x, na.rm = TRUE))
  }


##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param nameme1
##' @return
##' @author Joe Shannon
##' @export
study_path <- 
  function(dataset,
           paths = study_data_paths) {
  
    match.arg(dataset,
              choices = c("well_levels", 
                          "precip",
                          "snowmelt"))
    
    grep(pattern = dataset,
         x = paths,
         fixed = TRUE,
         value = TRUE)
  
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param y
##' @param ... can be a character vector [c("x1", "x2")] or comma separated x values "x1", "x2"
##' @return
##' @author Joe Shannon
##' @export
make_formula <- 
  function(y, ...) {
    
    form <- 
      paste(y, "~", Reduce(function(x, y){paste(x, y, sep = " + ")}, ...))
    
    as.formula(form, env = parent.env(environment()))
    
  }

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param x a vector of class Date/Datetime or a numeric vector representing month
##' @return
##' @author Joe Shannon
##' @export
as.climate_season <- 
  function(x, return.factor){
    stopifnot(is.numeric(x) | is.Date(x))
    
    if(is.numeric(x) & !all(data.table::between(x, 1, 12))){
      stop("If x is numeric then all values must be between 1 and 12")
    }
    
    if(is.Date(x)){
      x <- 
        month(x)
    } 
    
    seas <- 
      data.table::fcase(x %in% c(12, 1, 2), "djf", 
                        x %in% 3:5, "mam", 
                        x %in% 6:8, "jja", 
                        x %in% 9:11, "son")
  
    if(!return.factor){
      return(seas)
    }
      
    factor(seas, 
           levels = c("djf", "mam", "jja", "son"),
           ordered = TRUE)
  }

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param x a vector of class Date/Datetime
##' @return
##' @author Joe Shannon
##' @export
as.water_year <- 
  function(x, wy.month){
    stopifnot(any(grepl("(Date|POSIXt)", class(x))))
    
    ifelse(month(x) > (wy.month - 1),
           year(x) + 1,
           year(x))
  }

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param x a vector of class Date/Datetime
##' @return
##' @author Joe Shannon
##' @export
as.dowy <- 
  function(x, wy.month){
    stopifnot(any(grepl("(Date|POSIXt)", class(x))))
    
    1 + as.numeric(x - as.Date(paste(as.water_year(x, wy.month)-1, wy.month, "01", sep = "-")))
    
  }

lyne_hollick <- 
  function(q,
           alpha = 0.925,
           reflection = 0.5,
           passes = 3,
           na.rm = TRUE){
    stopifnot(is.numeric(q),
              is.numeric(alpha),
              is.numeric(reflection),
              is.numeric(passes))
    
    if(passes %% 2 != 1){
      stop("passes must be an odd integer.")
    }
    
    leading_nas <- trailing_nas <- min_q <- 0
    
    n_q <- 
      length(q)
    
    if(any(is.na(q))){
      
      if(!na.rm){
        stop("There are NA's in 'q', please use na.rm = TRUE to proceed.")
      }
      
      if(sum(!is.na(q)) == 0){return(q)}
      
      q <- 
        zoo::na.approx(q, na.rm = FALSE)
      
      if(is.na(q[1])){
        leading_nas <- 
          min(which(!is.na(q))) - 1
      }
      
      if(is.na(q[n_q])){
        trailing_nas <- 
          n_q - max(which(!is.na(q)))
      }
      
      q <- 
        q[(1 + leading_nas):(length(q) - trailing_nas)]
      
    }
    
    if(any(q < 0)){
      
      min_q <- 
        min(q, na.rm = TRUE)
      
      q <- 
        q - min_q
    }
    
    if(reflection != 0){
      reflection_n <- 
        floor(reflection*length(q))
      
      q_pad <- 
        c(rev(q[1:reflection_n]), q, rev(q)[1:reflection_n])
    } else {
      q_pad <- 
        q
    }
    
    q_f <- 
      rep(NA_real_, length(q_pad))
    
    q_f[1] <-
      q_pad[1+reflection_n]
    
    q_b <- 
      q_pad - q_f
    
    pass <- 0
    
    q_prev <- 
      q_pad
    
    while(pass < passes){
      for (i in 2:length(q_pad)){
        
        q_f[i] <-
          # lyne_hollick_eq(alpha, q_f, q_prev, i)
          alpha * q_f[i -1] + ((1 + alpha) / 2) * (q_prev[i] - q_prev[i - 1])
        
        if(q_f[i] < 0) q_f[i] <- 0
        if(q_f[i] > q_prev[i]) q_f[i] <- q_prev[i]
        
        q_b[i] <-
          q_prev[i] - q_f[i]
      }
      
      pass <- 
        pass + 1
      
      if(pass < passes){
        
        q_prev <- 
          rev(q_b)
        
        q_f[1] <-
          rev(q_b)[1]
        
        q_b[1] <- 
          q_prev[1] - q_f[1]
        
      }
      
    }
    c(rep(NA_real_, leading_nas), 
      q_b[(reflection_n + 1):(length(q_b) - reflection_n)] + min_q,
      rep(NA_real_, trailing_nas))
  }


ytd_sum <- 
  function(x){
    
    stopifnot(is.numeric(x))
    
    if(all(!is.na(x))){return(cumsum(x))}
    
    indx <- 
      which(!is.na(x))
    
    # x <- 
    #   zoo::na.approx(x, rule = 2)
    
    ytd_x <- 
      cumsum(na.omit(x))
    
    out <- 
      rep(NA, times = length(x))
    
    out[indx] <- 
      ytd_x
    
    out
  }

diff_na <- 
  function(x){
    
    stopifnot(is.numeric(x))
    
    if(all(!is.na(x))){return(c(NA, diff(x)))}
    
    indx <- 
      which(!is.na(x))
    
    dx <- 
      c(NA, diff(na.omit(x)))
    
    out <- 
      rep(NA, length(x))
    
    out[indx] <- 
      dx
    
    out
  }


pretty_round <- function(data, n_digits) {
  
  format(round(data, n_digits), nsmall = n_digits)
  
}

filled_rolling_sum <- function(x, na.rm = TRUE) {
  
  # Mirror first and last 15 rows to get rid of NA values from rolling mean
  x <- c(
    head(x, 15),
    x,
    tail(x, 15)
  )
  
  # Get 30 day rolling mean of precip variables
  x <- data.table::frollsum(x = x, n = 30, align = "center", na.rm = na.rm)
  
  # Drop mirrored values
  x <- head(x, -15)
  x <- tail(x, -15)
  
}

filled_rolling_mean <- function(x, na.rm = TRUE) {
  
  # Mirror first and last 15 rows to get rid of NA values from rolling mean
  x <- c(
    head(x, 15),
    x,
    tail(x, 15)
  )
  
  # Get 30 day rolling mean of precip variables
  x <- data.table::frollmean(x = x, n = 30, align = "center", na.rm = na.rm)
  
  # Drop mirrored values
  x <- head(x, -15)
  x <- tail(x, -15)
  
}
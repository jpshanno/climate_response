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
calculate_climate_season <- 
  function(x){
    stopifnot(is.numeric(x) | is.Date(x))
    
    if(is.numeric(x) & any(data.table::between(x, 1, 12))){
      stop("If x is numeric then all values must be between 1 and 12")
    }
    
    if(is.Date(x)){
      mon <- 
        month(x)
    }
    
    data.table::fcase(mon %in% c(12, 1, 2), "djf", 
                      mon %in% 3:5, "mam", 
                      mon %in% 6:8, "jja", 
                      mon %in% 9:11, "son")
    
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
  function(x){
    stopifnot(any(grepl("(Date|POSIXt)", class(x))))
    
    ifelse(month(x) > 9,
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
  function(x){
    stopifnot(any(grepl("(Date|POSIXt)", class(x))))
    
    1 + as.numeric(x - as.Date(paste(as.water_year(x)-1, "10-01", sep = "-")))
    
  }

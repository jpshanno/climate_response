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

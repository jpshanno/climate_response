# autoloess.R: compute loess metaparameters automatically
# Kyle Gorman <gormanky@ohsu.edu>
# https://gist.githubusercontent.com/kylebgorman/6444612/raw/04ee28214e1ee2f959c97f2c87c7d4733eae77f0/autoloess.R

aicc.loess <- function(fit) {
  # compute AIC_C for a LOESS fit, from:
  # 
  # Hurvich, C.M., Simonoff, J.S., and Tsai, C. L. 1998. Smoothing 
  # parameter selection in nonparametric regression using an improved 
  # Akaike Information Criterion. Journal of the Royal Statistical 
  # Society B 60: 271–293.
  # 
  # @param fit        loess fit
  # @return           'aicc' value
  stopifnot(inherits(fit, 'loess'))
  # parameters
  n <- fit$n
  trace <- fit$trace.hat
  sigma2 <- sum(resid(fit) ^ 2) / (n - 1)
  return(log(sigma2) + 1 + (2 * (trace + 1)) / (n - trace - 2))
}

# Approach from gist did not work in nested data frame using data = .SD. Use the
# autoloess below
# autoloess <- function(fit, span=c(.1, .9)) {
#   # compute loess fit which has span minimizes AIC_C
#   # 
#   # @param fit        loess fit; span parameter value doesn't matter
#   # @param span       a two-value vector representing the minimum and 
#   #                   maximum span values
#   # @return           loess fit with span minimizing the AIC_C function
#   stopifnot(inherits(fit, 'loess'), length(span) == 2)
#   # loss function in form to be used by optimize
#   f <- function(span) aicc.loess(update(fit, span=span))
#   # find best loess according to loss function
#   return(update(fit, span=optimize(f, span)$minimum))
# }


autoloess <- 
  function(...) {
    
    # @param ... Arguments to loess
    
    args <- 
      list(...)
    
    if(length(args[["span"]]) == 1){
      return(loess(...))
    }
    
    span_range <- 
      args[["span"]]
    
    args[["span"]] <- 
      args[["span"]][1]
    
    fit <- 
      do.call("loess", args)
    
    # loss function in form to be used by optimize
    f <- function(span) aicc.loess(update(fit, span=span))
    # find best loess according to loss function
    return(update(fit, span=optimize(f, span_range)$minimum))
  }

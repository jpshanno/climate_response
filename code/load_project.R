# Load packages
library(targets)
pkgs <- purrr::map(readLines("code/climate_packages.txt"), strsplit, " ")
pkgs <- purrr::map(pkgs, unlist)
purrr::map_if(pkgs, 
              ~length(.x) == 1, 
              ~library(.x[[1]], 
                       character.only = TRUE))
purrr::map_if(pkgs, 
              ~length(.x) > 1, 
              ~library(.x[[1]], 
                       character.only = TRUE, 
                       include.only = .x[-c(1)]))
rm(pkgs)
# Load functions
for(i in list.files("code/functions", full.names = TRUE)){
  source(i, verbose = FALSE)
  rm(i)
} 

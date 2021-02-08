# Load packages
library(targets)
purrr::walk(readLines("code/climate_packages.txt"), library, character.only = TRUE)

# Load functions
for(i in list.files("code/functions", full.names = TRUE)){
  source(i, verbose = FALSE)
  rm(i)
} 

# Title: Import CA base map 
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Sat Oct  3 01:33:59 2020

# preparation --------
rm(list = ls())
cat("\014")
setwd("~/Lab/caledna_transect/")

library(ggmap)
today = format(Sys.Date(), "%Y%m%d")

# def functions --------

# def variables --------
# output dir 
dir.create("./derive_data/step0_prepare_data/sample_map/", recursive = T)
# start plotting
ggmap::register_google(key = "YOUR API KEY")

# load data --------

# main --------
cabase <- get_map(location = "california", maptype = "roadmap", zoom = 6) # base map
save(cabase, file = 
         paste0("./derive_data/step0_prepare_data/sample_map/CA_ggplot_basemap_", today, ".RData"))

# cleanup --------

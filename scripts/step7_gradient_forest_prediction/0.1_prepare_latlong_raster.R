# Title : prepare raster for generate prediction --------
# stack the variables 
# Author: Meixi Lin
# Date: Sun Jun  9 19:44:04 2019

# preparation --------
rm(list = ls())
cat("\014")
# here set the directory on UCLA hoffman2 layers
# setwd("/u/home/m/meixilin/project-rwayne/abiotic_transect/")
# # if you need to set things locally 
setwd("~/UCLA/Lab/abiotic_transect/")
library(rgeos)
library(dplyr)
library(raster)
library(gstat)
library(rgdal)
library(rasterVis)

# resolution 
myres = 100

# CRS TO use
mycrs = "+proj=longlat +datum=WGS84"

# some functions 
source("./r_codes/step0_prepare_data/metadata/generate/function_raster.R")
functions.file <- list.files("./r_codes/step4_gradient_forest/functions/")
lapply(functions.file, function(xx) {source(file = paste0("./r_codes/step4_gradient_forest/functions/", xx))})

# load a file with metadata name and directory --------
rasterdir <- read.csv(file = "./final_data/metadata/new_raster_dir.csv", stringsAsFactors = F)

# load raster layers --------
myvar.files <- rasterdir$final_dir[rasterdir$raster %in% myvar]

# load ca boundary -------
ca0 = readOGR("./maps_vectors_rasters/vectors/CA_boundary_TIGER/CA_boundary_TIGER.shp")
ca = spTransform(ca0, CRS(mycrs)); rm(ca0)

# stack, crop and mask the layers --------
# get the values that was not Longitude or Latitude 
myvar.ras <- stack(myvar.files)
# get the Latitude and Longitude layers
Longitude = init(myvar.ras$EE_NDVI32_CA_100m, 'x')
Latitude = init(myvar.ras$EE_NDVI32_CA_100m, 'y')

# write out the latitude and longitude variables --------
writeRaster(x = Longitude, filename = "./maps_vectors_rasters/projected/x100m/longitude_ca_100m_wgs84.tif", format = "GTiff")
writeRaster(x = Latitude, filename = "./maps_vectors_rasters/projected/x100m/latitude_ca_100m_wgs84.tif", format = "GTiff")

# clean up --------
closeAllConnections()

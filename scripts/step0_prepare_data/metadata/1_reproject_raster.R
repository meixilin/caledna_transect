# Title: redo extraction of metadata layers step one: clip all layers to same resolution and projection --------
# Author: Meixi Lin
# Date: Thu Mar  7 22:52:33 2019
# Author: Meixi Lin 
# Date: Tue Mar 26 12:46:16 2019
# Modification: change all resolution to 100m and see 
# Date: Mon Apr  1 11:29:01 2019
# Modifications: match local and hoffman2 directory format 

# This script is a major overhaul to match up gradient forest analysis 
# Most of these data will be extracted by first clipping to the same resolution 
# calculate the distance to the points, and extract values from the closest points

# preparation --------------------------------------------------------------------------------
rm(list = ls())
cat("\014")
# here set the directory on UCLA hoffman2 layers
setwd("/u/home/m/meixilin/project-rwayne/abiotic_transect/")
library(dplyr)
library(raster)
library(gstat)
library(rgdal)
library(rasterVis)

# CRS TO use
mycrs = "+proj=longlat +datum=WGS84"

# target resolution in meter 
myres = 100

# some functions 
source("./scripts/step0_prepare_data/metadata/function_raster.R")

# get index to perform the analysis ------------------------------------------------------------
ii = as.integer(commandArgs(trailingOnly=TRUE)[1])

# load a file with metadata name and directory ------------------------------------------------------------
rasterdir <- read.csv(file = "./final_data/raster_dir.csv", stringsAsFactors = F)

# load the earth engine layer as reference -------------------------------------------------------------
rasterref <- raster(x = "./maps_vectors_rasters/rasters/earth_engine/EE_B1_CA_100m.tif")
rastertemp <- projectExtent(rasterref, crs = mycrs)

# load raster layers by input value --------
# 1. load individual rasters ######################
print(paste("processing raster:", rasterdir[ii,'raster'], date()))
rawrr = raster(x = rasterdir[ii,'raster_dir'])
print(rawrr)

# 2. check if need aggregate ######################
rawres = rasterdir[ii, 'org_res']
if (rawres < myres) {
    print("Raster had lower resolution, need aggregate")
    fa = floor(myres/rawres)
    aggrrfile = paste0("./maps_vectors_rasters/projected/aggregate/", rasterdir[ii,'raster'], "_ca_fact_", fa, ".tif")
    aggrr = raster::aggregate(x = rawrr, fact = fa, fun = mean, expand = T, na.rm = T,
                              filename = aggrrfile)
}else {
    print("Raster had higher resolution, no need to aggregate")
    aggrr = rawrr
}
rm(rawrr)

# 3. check if need reprojection ###################
print(compareCRS(aggrr, rastertemp))
print(res(aggrr) == res(rastertemp)) 

if (rasterdir[ii, 'raster'] == 'cecsol') {
    # the range of CELSOL was too large and somehow project raster couldn't handle it
    aggrr = raster::crop(x = aggrr, y = rastertemp, filename = "./maps_vectors_rasters/rasters/soilgrid_x100m/cecsol_ca_250m_wgs84.tif")
}

if (rasterdir[ii, 'transform'] == T) {
    print("Raster need transformation to a slightly different resolution")
    prorrfile = paste0("./maps_vectors_rasters/projected/x100m/", rasterdir[ii,'raster'], "_ca_100m_wgs84.tif")
    if (rasterdir[ii, 'data_type'] == 'factor') {
        prmethod = 'ngb'
    } else {
        prmethod = 'bilinear'
    }
    prorr = projectRaster(from = aggrr, to = rastertemp, 
                          filename = prorrfile, 
                          method = prmethod,
                          overwrite = T)
    finaldir = prorr@file@name
} else {
    print("Raster is the same as raster template, no need to project.")
    finaldir = aggrr@file@name
}

# 4. write out the final directory for every file ###
print(paste("final directory for", rasterdir[ii, 'raster'], "is", finaldir))
print(paste("finished raster:", rasterdir[ii,'raster'], date()))
newdir = cbind(rasterdir[ii,], finaldir)
write.table(newdir, file = "./derive_data/metadata/20190403/new_raster_dir.tsv", append = T, sep = "\t", col.names = F)

# close all device --------
closeAllConnections()
# Title: redo extraction of metadata layers step one: clip all layers to same resolution and projection --------
# Author: Meixi Lin
# Date: Thu Mar  7 22:52:33 2019
# Author: Meixi Lin 
# Date: Tue Mar 26 12:46:16 2019
# Modification: change all resolution to 100m and see 
# Date: Mon Apr  1 11:29:01 2019
# Modifications: match local and hoffman2 directory format 
# Date: Tue Apr 30 07:58:14 2019
# Modifications: match local and hoffman2 directory format 

# This script is a major overhaul to match up gradient forest analysis 
# Most of these data will be extracted by first clipping to the same resolution 
# calculate the distance to the points, and extract values from the closest points

# preparation --------
rm(list = ls())
cat("\014")
# here set the directory on UCLA hoffman2 layers
setwd("/u/home/m/meixilin/project-rwayne/abiotic_transect/")
# # if you need to set things locally 
# setwd("~/UCLA/Lab/abiotic_transect/")
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
source("./scripts/step0_prepare_data/metadata/function_raster.R")

# load and prepare georeferenced points --------
# metadata and latlong
metadata <- read.csv("./final_data/Final_metadata_v0.csv", stringsAsFactors = F) 
# ALWAYS REMEMBER THAT LONGITUDE COMES FIRST 
latlong <- metadata %>% dplyr::select(Longitude, Latitude)
rownames(latlong) <- metadata$MatchName
# define points
pts <- SpatialPoints(latlong, proj4string = CRS(mycrs)) # code for WGS1984
# plot(pts)

# load a file with metadata name and directory --------
rasterdir <- read.csv(file = "./final_data/new_raster_dir.csv", stringsAsFactors = F)

# start a dataframe to store the data --------
newmeta <- as.data.frame(metadata[,c('MatchName', 'Longitude', 'Latitude')])

# load raster layers and perform batch extraction in a for loop --------
for (ii in  1:nrow(rasterdir)) {
# load individual rasters 
print(paste("processing raster:", rasterdir[ii,'raster'], date()))
newrr = raster(x = rasterdir[ii,'final_dir'])
print(newrr)
# now extract values from newrr
ptsextract <- raster::extract(x = newrr, y = pts, method = "simple", sp = T)
ptsvalue = cbind(row.names(ptsextract@coords), ptsextract@data,
                 stringsAsFactors = F)
colnames(ptsvalue) = c("MatchName", "raster_value")
# now fix values that was na 
if (length(which(is.na(ptsvalue))) > 0) {
  print(paste("current raster:", rasterdir[ii,'raster'], 
              "have NA values, finding the closest neighbor."))
  if (rasterdir[ii, 'category'] == "topology") {
    nafixed <- fix.napts(newrr, pts, myex = ptsvalue, mywidth = 0.5, maxdist = myres)
  } else {
    nafixed <- fix.napts(newrr, pts, myex = ptsvalue, mywidth = 0.5)
  }
  write.csv(nafixed, 
            file = paste0("./derive_data/metadata/20190430/x100m/", 
                          rasterdir[ii,'raster'], "_detail_na_move.csv"))
  # now join the nafixed and ptsvalue
  naid = which(is.na(ptsvalue[,'raster_value']))
  # check if nafixed worked 
  ptsvalue[naid,'raster_value'] = nafixed[,'raster_value']
} else {
  print(paste("current raster:", rasterdir[ii,1], "does not have NA values, moving on."))
}
# get data out 
colnames(ptsvalue) <- c("MatchName", rasterdir[ii,'raster'])
write.csv(ptsvalue, 
          file = paste0("./derive_data/metadata/20190430/x100m/", 
                        rasterdir[ii,'raster'], "_metadata_value.csv"))
# add the value to new metadata 
newmeta <- left_join(x = newmeta, y = ptsvalue, by = 'MatchName')
print(paste("finished raster:", rasterdir[ii,1], date()))
} 

# write out the data to csv --------
print("all done for")
print(colnames(newmeta))
write.csv(newmeta, file = "./derive_data/metadata/20190430/x100m/Final_newmeta.csv")
date()
closeAllConnections()

# Title: scale the variables as comparison for gradient forest predictions  
# Store the scaled CA_grid into raster formats (step7_gradient_forest_prediction/all_100m/4.1_step7_ref_scale_grid_20201014.R) --------
# Author: Meixi Lin
# Date: Thu Oct 31 14:04:57 2019
# Author: Meixi Lin
# Date: Mon Oct 19 01:21:07 2020
# Modification: check before publications

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE)
setwd("/u/project/rwayne/meixilin/caledna_transect/")
library(data.table)
library(sp)
library(raster)
source("./scripts/function_transect.R")
print(sessionInfo())

# setting variables --------
# CRS TO use
mycrs = "+proj=longlat +datum=WGS84"
# mark the transformation as "nobio_CAgrid"
nobio = "nobio_CAgrid"
# extent for plotting 
myextent <- raster::extent(CAlimit)

# the gradient forest setting
gfdir = "./derive_data/step6_gradient_forest/"
gfname = "gf_deco_all_Family_Presence_2000_2_0.05_FALSE_17"

# output directory 
outdir = paste0("./derive_data/step7_gradient_forest_prediction/", gfname, "/", nobio, "/")
# the Scld_CA filename (step7_gradient_forest_prediction/all_100m/4.1_step7_ref_scale_grid_20201014.R)
Scld_CA_filepath = paste0(outdir, "predict_nogf_Scld_CAgrid.RData")

# load data --------
# load the CAgrid
print(paste(date(), "Loading Scld_CA ...", Scld_CA_filepath))
load(file = Scld_CA_filepath) # should be named as "Scld_CA"
print(str(Scld_CA))

# write data --------
print(paste(date(), "Convert Scld_CA to raster rsScld_CA ..."))
# Remove data.table class
Scld_CA = data.table::setDF(Scld_CA)
data.table::is.data.table(Scld_CA)
sp::coordinates(Scld_CA) <- ~ lon + lat
sp::proj4string(Scld_CA) <- CRS(mycrs)
# coerce to SpatialPixelsDataFrame
sp::gridded(Scld_CA) <- TRUE
# coerce to raster
rsScld_CA <- raster::stack(Scld_CA)
print(rsScld_CA)
rm(Scld_CA) # release space 

print(paste(date(), "Exporting rsScld_CA ...", paste0(outdir, "predict_nogf_rsScld_CAgrid.tif")))
raster::writeRaster(rsScld_CA, filename = paste0(outdir, "predict_nogf_rsScld_CAgrid.tif"), format = "GTiff")
gc()
print(paste(date(), "Exporting rsScld_CA ...", paste0(outdir, "predict_nogf_rsScld_CAgrid.grd")))
raster::writeRaster(rsScld_CA, filename = paste0(outdir, "predict_nogf_rsScld_CAgrid.grd"), format = "raster")

# cleanup --------
print(paste(date(), "All done. Job finished successfully."))
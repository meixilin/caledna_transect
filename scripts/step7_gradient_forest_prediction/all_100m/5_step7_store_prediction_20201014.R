# Title: store the predictions of gradient forest output --------
# Author: Meixi Lin
# Date: Sun Nov  3 16:38:59 2019
# Author: Meixi Lin
# Date: Sun Oct 18 00:08:01 2020
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

# get arguments ---------
args <-  commandArgs(trailingOnly=TRUE)

extrap <- as.logical(args[1]) # extrapolation 
Trns_CA_filepath <- as.character(args[2]) # Path to file "predict_gf_Trns_CAgrid_extrap.RData/predict_gf_Trns_CAgrid_noextrap.RData"
outdir <- as.character(args[3])

dir.create(outdir, recursive = FALSE, showWarnings = FALSE)

# setting variables --------
# CRS TO use
mycrs = "+proj=longlat +datum=WGS84"

# def functions --------
getrsTrns_CA_filename <- function(extrap, outdir) {
    if (extrap == TRUE) {
        rsTrns_CA_filename = paste0(outdir, "predict_gf_rsTrns_CAgrid_extrap")
    } else { 
        if (extrap == FALSE) {
            rsTrns_CA_filename = paste0(outdir, "predict_gf_rsTrns_CAgrid_noextrap")
        } else {
            stop("Wrong extrap value!")
        }
    } 
    return(rsTrns_CA_filename)
}

# load data --------
print(paste(date(), "Loading ...", Trns_CA_filepath))
load(Trns_CA_filepath) # should be named as Trns_CA
str(Trns_CA)

# write data --------
print(paste(date(), "Convert Trns_CA to raster rsTrns_CA ..."))
Trns_CA <- data.table::setDF(Trns_CA)
str(Trns_CA)
sp::coordinates(Trns_CA) <- ~ lon + lat
sp::proj4string(Trns_CA) <- CRS(mycrs)
# coerce to SpatialPixelsDataFrame
sp::gridded(Trns_CA) <- TRUE
# coerce to raster
rsTrns_CA <- raster::stack(Trns_CA)
print(rsTrns_CA)
rm(Trns_CA)

rsTrns_CA_filename=getrsTrns_CA_filename(extrap, outdir)
print(paste(date(), "Exporting ...", paste0(rsTrns_CA_filename,".tif")))
raster::writeRaster(rsTrns_CA, filename = paste0(rsTrns_CA_filename,".tif"), format = "GTiff")
gc()
print(paste(date(), "Exporting ...", paste0(rsTrns_CA_filename,".grd")))
raster::writeRaster(rsTrns_CA, filename = paste0(rsTrns_CA_filename,".grd"), format = "raster")

# cleanup --------
print(paste(date(), "All done. Job finished successfully."))

# Title: Get soilgrid layers
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Fri Sep  4 17:24:18 2020

# preparation --------
rm(list = ls())
cat("\014")
setwd("~/Lab/caledna_transect/")

library(rgdal)
library(raster)

source("./scripts/function_transect.R")
today = format(Sys.Date(), "%Y%m%d")

# def functions --------
plot_raster <- function(rr, rrname) {
    plot(rr, main = rrname, xlab = "Longitude", ylab = "Latitude", cex = 1.5)
    return(0)
}

# def variables --------
mycrs = "+proj=longlat +datum=WGS84"
uncertainty_layers = list.files(path = "./maps_vectors_rasters/rasters/uncertainty_rasters/", pattern = ".tif")
uncertainty_headers1 = stringr::str_extract(uncertainty_layers, "(?<=EE_)(.+)(?=\\_CA)")
uncertainty_headers = c("Sentinel B1 stdDev", "Sentinel B10 stdDev", "Sentinel B11 stdDev", "Sentinel B12 stdDev", "Sentinel B2 stdDev",
                       "Sentinel B3 stdDev", "Sentinel B4 stdDev", "Sentinel B5 stdDev", "Sentinel B6 stdDev", "Sentinel B7 stdDev", 
                       "Sentinel B8 stdDev", "Sentinel B8A stdDev", "Sentinel B9 stdDev", "NLCD ptrcv Uncertainty %")
# reorder the layers and the headers
uncertainty_layers = uncertainty_layers[c(1,5:13,2:4,14)]
uncertainty_headers = uncertainty_headers[c(1,5:13,2:4,14)]
uncertainty_headers1 = uncertainty_headers1[c(1,5:13,2:4,14)]
# load data --------
# get ca boundary 
ca0 = readOGR("./maps_vectors_rasters/vectors/CA_boundary_TIGER/CA_boundary_TIGER.shp")
ca = spTransform(ca0, CRS(mycrs)); rm(ca0)

uncertainty_rasters = lapply(uncertainty_layers, function(xx) {
    rr = raster(paste0("./maps_vectors_rasters/rasters/uncertainty_rasters/", xx))
    return(rr)
})

# main --------
for (ii in 1:length(uncertainty_headers)) {
    rm(rr, rrname, rrhead)
    rrname = uncertainty_headers[ii]
    rr = uncertainty_rasters[[ii]]
    rrhead = uncertainty_headers1[ii]
    print(rr)
    print(rrname)
    print(rrhead)
    pdf(file = paste0("./plots/step0_prepare_data/", rrhead,"_uncertainty_",today,".pdf"), width = 6, height = 6)
    plot_raster(rr, rrname)
    dev.off()
}

# cleanup --------


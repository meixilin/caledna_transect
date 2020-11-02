# Title : prepare splitted raster for use of predictions --------
# for all potential gradient forest models 
# Author: Meixi Lin
# Date: Sun Jun  9 23:43:48 2019

# preparation --------
rm(list = ls())
cat("\014")
# here set the directory on UCLA hoffman2 layers
options(echo = T)
options(java.parameters = "-Xmx90G")
setwd("/u/project/rwayne/meixilin/caledna_transect/")
library(raster)
library(dplyr)

source("./scripts/function_transect.R")

# set variables -------
# CRS TO use
mycrs = "+proj=longlat +datum=WGS84"
mysize = 1000000
indir = "./derive_data/step7_gradient_forest_prediction/"
outdir = "./derive_data/step7_gradient_forest_prediction/split_all_grid/"
dir.create(outdir, recursive = T)

# load raster layers --------
myvar_ras = raster::stack(paste0(indir, "myvar_ras.grd"))
print(myvar_ras)
print(ncell(myvar_ras)) 
# [1] 121613742

# sample all sites by spliting them into tiles  --------
grids = raster::as.data.frame(myvar_ras, xy = TRUE, na.rm = TRUE) # removing any sites that might contain missing data 
colnames(grids) = c(contlist, "lon", "lat") 
endgrid = nrow(grids)

# write the entire csv 
print(paste("Write grid data frame ...", date()))
save(grids, file = paste0(indir, "myvar_ras_dataframe.RData"))
write.csv(grids , file=paste0(indir, "myvar_ras_dataframe.csv"), row.names = FALSE)
print(paste("Done grid data frame ...", date()))

# for each mysize lines write to a rdata for downstream processing 
start = 1
end = (start + mysize - 1)
ii = 0

while (start < endgrid) {
    # subset grid
    end = min(end, endgrid)
    id = sprintf("%02d", ii)
    print(paste("Subsetting...", start, end))
    grid = grids[start:end, ]
    outname = paste0(outdir, "grid_ca_", id, ".RData")
    print(paste(date(), "Output...", outname))
    save(grid, file = outname)
    # reset values 
    start = end + 1
    end = (start + mysize - 1)
    ii = ii + 1
}

# end --------
closeAllConnections()

# Title: make predictions for all variables without foreach --------
# Author: Meixi Lin
# Date: Thu Oct 31 19:30:17 2019
# Author:
# Date: Wed Oct 14 13:25:35 2020
# Modification: check before publications
# Author:
# Date: Wed Oct 14 14:40:05 2020
# Modification: No foreach package

# preparation --------
rm(list = ls())
cat("\014")
options(echo = T)
options(java.parameters = "-Xmx10G")
setwd("/u/project/rwayne/meixilin/caledna_transect/")
library(gradientForest)
source("./scripts/function_transect.R")

# get arguments ---------
args <-  commandArgs(trailingOnly=TRUE)

extrap <- as.logical(args[1])
gfdir <- as.character(args[2])
gfname <- as.character(args[3]) 
outdir <- as.character(args[4]) 
ii <- as.character(args[5])

# setting variables --------
# CRS TO use
mycrs = "+proj=longlat +datum=WGS84"
griddir = "./derive_data/step7_gradient_forest_prediction/split_all_grid/grid_ca_"

print(ii)
dir.create(outdir,recursive = FALSE, showWarnings = FALSE)

# load data --------
gf = loadRData(file = paste0(gfdir, gfname, ".RData"))
imp.vars = names(importance(gf))
 
# predict grid --------
# set file names 
if (extrap == TRUE) {
    Trns_grid_filename = paste0(outdir, "predict_gf_cagrid_", ii, "_extrap.RData")
} else { 
    if (extrap == FALSE) {
        Trns_grid_filename = paste0(outdir, "predict_gf_cagrid_", ii, "_noextrap.RData")
    } else {
        stop("Wrong extrap value!")
    }
} 

mygrid = paste0(griddir, ii, ".RData")
print(paste(date(), "Loading ...", mygrid))
load(mygrid)
# this generates a normalized importance based on the value of your grids
# how the grid is ordered does not affect the output, here i reordered columns from the most important to the least variables 
print(paste(date(), "Transforming ...", Trns_grid_filename))
Trns_grid = cbind(grid[,c("lon", "lat")], 
                  gradientForest::predict.gradientForest(gf,grid[,imp.vars], extrap = extrap))
save(Trns_grid, file = Trns_grid_filename)
print(paste(date(), "Transforming done ...", Trns_grid_filename))


# Title: scale the variables as comparison for gradient forest predictions --------
# Author: Meixi Lin
# Date: Thu Oct 31 14:04:57 2019
# Author:
# Date:
# Modification:

# preparation --------
rm(list = ls())
cat("\014")
# here set the directory on UCLA hoffman2 layers
options(echo = T)
options(java.parameters = "-Xmx90G")

library(sp)
library(raster)
library(gradientForest)
setwd("/u/project/rwayne/meixilin/caledna_transect/")

source("./scripts/function_transect.R")

# setting variables --------
# CRS TO use
mycrs = "+proj=longlat +datum=WGS84"
# mark the transformation as "nobio_CAgrid"
# seed to use
myseed = 17

myextent <- raster::extent(CAlimit)

# the gradient forest setting
gfdir <- "./derive_data/step6_gradient_forest/"
gfname <- "gf_deco_all_Family_Presence_2000_2_0.05_FALSE_17"
rand_grids_file = "./derive_data/step7_gradient_forest_prediction/gf_rand_grids.csv"
# output directories
# make a directory for each gf because the PC vectors will be different 
outdir <- paste0("./derive_data/step7_gradient_forest_prediction/", gfname, "/nobio_CAgrid/")
dir.create(outdir, recursive = T)

# load data --------
gf = loadRData(file = paste0(gfdir, gfname, ".RData"))
imp.vars = names(importance(gf))
print(gf)
print(imp.vars)

# load the entire CA grids # should be a variable named "grids"
load(file = "./derive_data/step7_gradient_forest_prediction/myvar_ras_dataframe.RData")

# load the rand_grid for plotting 
grid = read.csv(rand_grids_file)

# predict site and random grid --------
# transform the site and random grid environmental predictors 
print(paste("Scaling site and random grids ...", date()))
centers <- apply(grids[,imp.vars], 2, mean)
print(centers)
scales <- apply(grids[,imp.vars], 2, sd)
print(scales)
save(centers, scales, file = paste0(outdir, "centers_scales_CAgrid.RData"))

# scale the sites and the random points as well 
Scld_site = scale(gf$X[,imp.vars], center = centers, scale = scales)
Scld_randgrid = cbind(grid[,c("lon", "lat")], 
                       scale(grid[,imp.vars], center = centers, scale = scales))

save(Scld_site, file = paste0(outdir, "predict_nogf_Scld_site.RData"))
save(Scld_randgrid, file = paste0(outdir, "predict_nogf_Scld_randgrid.RData"))

closeAllConnections()
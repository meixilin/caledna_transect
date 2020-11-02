# Title : make predictions for the random sites in CA for plotting --------
# gf_deco_10_all_FALSE_0.05_Presence_100_17
# Author: Meixi Lin
# Date: Tue Oct 29 15:11:53 2019
# Author: Meixi Lin
# Date: Tue Apr 21 17:19:08 2020

# preparation --------
rm(list = ls())
cat("\014")
options(echo = T)
setwd("~/Lab/caledna_transect/")

library(gradientForest)

# some functions 
source("./scripts/function_transect.R")

# setting variables --------
# CRS TO use
mycrs = "+proj=longlat +datum=WGS84"
# the rasternames
predictors = contlist
rand_grids_file = "./derive_data/step7_gradient_forest_prediction/gf_rand_grids.csv"
gfdir <- "./derive_data/step6_gradient_forest/"
gfname <- "gf_deco_all_Family_Presence_2000_2_0.05_FALSE_17"
outdir <- paste0("./derive_data/step7_gradient_forest_prediction/", gfname, "/")
outdir1 <- paste0(outdir, "extrap_TRUE/")
outdir2 <- paste0(outdir, "extrap_FALSE/")
dir.create(outdir1, recursive = T)
dir.create(outdir2, recursive = T)

# load data --------
gf = loadRData(file = paste0(gfdir, gfname, ".RData"))
imp.vars = names(importance(gf))

# read the random grids (1% of all data points in CA, 505012 sites)
grid = read.csv(rand_grids_file)

# predict grid --------
# this generate a normalized importance based on the value of your random grids
Trns_randgrid_extrap = cbind(grid[,c("lon", "lat")], 
                  predict.gradientForest(gf,grid[,imp.vars], extrap = T))
Trns_randgrid_noextrap = cbind(grid[,c("lon", "lat")], 
                         predict.gradientForest(gf,grid[,imp.vars], extrap = F))

save(Trns_randgrid_extrap, file = paste0(outdir1, "predict_gf_Trns_randgrid_extrap.RData"))
save(Trns_randgrid_noextrap, file = paste0(outdir2, "predict_gf_Trns_randgrid_noextrap.RData"))

# ending --------
closeAllConnections()
# Title : make predictions for the sites we sampled and plot the PCs --------
# gf_deco_10_all_FALSE_0.05_Presence_100_17
# Author: Meixi Lin
# Date: Tue Jun 26 15:12:01 2018
# Author: Meixi Lin
# Date: Sat Apr 18 12:16:38 2020

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

# predict site --------
# extrapolate or not extrapolate were slightly different 
Trns_site_extrap = gradientForest::predict.gradientForest(gf, extrap = T)
Trns_site_noextrap = gradientForest::predict.gradientForest(gf, extrap = F)

save(Trns_site_extrap, file = paste0(outdir1, "predict_gf_Trns_site_extrap.RData"))
save(Trns_site_noextrap, file = paste0(outdir2, "predict_gf_Trns_site_noextrap.RData"))

# ending --------
closeAllConnections()

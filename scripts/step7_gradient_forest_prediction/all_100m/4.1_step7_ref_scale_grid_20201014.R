# Title: scale the variables as comparison for gradient forest predictions --------
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
library(gradientForest)
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

# the CAgrid (step7_gradient_forest_prediction/0.3_prepare_raster_splits.R)
indir = "./derive_data/step7_gradient_forest_prediction/"
gridsfile = "myvar_ras_dataframe.RData"

# output directory 
outdir = paste0(indir, gfname, "/", nobio, "/")

# load data --------
# load the CAgrid
print(paste(date(), "Loading CAgrid ...", paste0(indir, gridsfile)))
load(file = paste0(indir, gridsfile)) # should be named as "grids"
print(str(grids))

print(paste(date(), "Loading gradient forest object ...", paste0(gfdir, gfname, ".RData")))
gf = loadRData(file = paste0(gfdir, gfname, ".RData"))
print(summary(gf))
imp.vars = names(importance(gf))
print(imp.vars)
rm(gf)

# center and scale grids --------
print(paste(date(), "Scaling CAgrid ..."))
Scld_CA = cbind(grids[,c("lon", "lat")], scale(grids[,imp.vars]))
print(str(Scld_CA))
rm(grids)
# Add data.table class to match Trns_CA
Scld_CA = data.table::setDT(Scld_CA)
data.table::is.data.table(Scld_CA)
print(paste(date(), "Exporting Scld_CA ...", paste0(outdir, "predict_nogf_Scld_CAgrid.RData")))
save(Scld_CA, file = paste0(outdir, "predict_nogf_Scld_CAgrid.RData"))

# cleanup --------
print(paste(date(), "All done. Job finished successfully."))
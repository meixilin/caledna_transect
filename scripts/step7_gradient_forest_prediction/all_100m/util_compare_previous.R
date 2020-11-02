# Title: Compare the two versions of rasters 
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Tue Oct 20 15:23:37 2020

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

# def functions --------
check_diff <- function (rr1, rr2) {
    output = (rr2-rr1)/rr2
    return(output)
}

# def variables --------
today = format(Sys.Date(), "%Y%m%d")
rrbeforefile = "/u/home/m/meixilin/project-rwayne/abiotic_transect/derive_data/step6_gf_prediction/gf_deco_all_Family_Presence_2000_2_0.05_FALSE_17/cagrid/FALSE/rsgrid/rsTrns_CA_grid.grd"
rrafterfile = "/u/home/m/meixilin/project-rwayne/caledna_transect/derive_data/step7_gradient_forest_prediction/gf_deco_all_Family_Presence_2000_2_0.05_FALSE_17/extrap_FALSE/predict_gf_rsTrns_CAgrid_noextrap.grd"
outdir = "/u/project/rwayne/meixilin/caledna_transect/derive_data/step7_gradient_forest_prediction/gf_deco_all_Family_Presence_2000_2_0.05_FALSE_17/extrap_FALSE/compare_previous/"

dir.create(outdir, recursive = FALSE, showWarnings = FALSE)

# load data --------
rrbefore = raster::stack(rrbeforefile)
rrafter = raster::stack(rrafterfile)

all(names(rrbefore) == names(rrafter))
# main --------
varlist = names(rrafter)
for (ii in 1:length(varlist)) {
    print(paste(date(), "Checking", varlist[ii], "..."))
    rrbeforetemp = raster::subset(rrbefore, subset = ii)
    print(rrbeforetemp)
    rraftertemp = raster::subset(rrafter, subset = ii)
    print(rraftertemp)
    rrdiffname = paste0(outdir, "diff_predict_gf_", varlist[ii], ".tif")
    rrdifftemp = raster::overlay(rrbeforetemp, rraftertemp, fun = check_diff, 
                                 filename = rrdiffname)
    print(paste(date(), "Percent difference (rrafter - rrbefore)/rrafter raster stored in ...", rrdiffname))
    print(rrdifftemp)
    print(raster::hist(rrdifftemp, breaks = c(-1, -5e-1, -1e-1, -1e-2, -1e-3, -1e-5, 0, 1e-5, 1e-3, 1e-2, 1e-1, 5e-1, 1), plot = FALSE))
    pdf(file = paste0(outdir, "plotHIST_diff_predict_gf_", varlist[ii], "_", today, ".pdf"), width = 8, height = 4) 
        pp=par(mfrow=c(1,2))
        raster::hist(rrdifftemp, main = paste0("plotHIST_diff_predict_gf_", varlist[ii]))
        raster::plot(rrdifftemp)
        par(pp)
    dev.off()
    print(paste(date(), "Checking", varlist[ii], "SUCCESS ..."))
} 

# cleanup --------
print(paste(date(), "All done. Job finished successfully."))


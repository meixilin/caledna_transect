# Title: combine the prediction --------
# Author: Meixi Lin
# Date: Thu Oct 31 14:50:12 2019
# Author: Meixi Lin
# Date: Wed Oct 14 16:38:39 2020
# Modification: check before publications

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE)
setwd("/u/project/rwayne/meixilin/caledna_transect/")
library(data.table)
source("./scripts/function_transect.R")

# get arguments ---------
args <-  commandArgs(trailingOnly=TRUE)

extrap <- as.logical(args[1])
rdatadir <- as.character(args[2]) # before files "predict_gf_cagrid_<idx>_noextrap.RData"
outdir <- as.character(args[3]) # should be a full directory

# def functions --------
getTrns_grid_filename <- function(extrap, ii, rdatadir) {
    if (extrap == TRUE) {
        Trns_grid_filename = paste0(rdatadir, "predict_gf_cagrid_", ii, "_extrap.RData")
    } else { 
        if (extrap == FALSE) {
            Trns_grid_filename = paste0(rdatadir, "predict_gf_cagrid_", ii, "_noextrap.RData")
        } else {
            stop("Wrong extrap value!")
        }
    } 
    return(Trns_grid_filename)
}

getTrns_CA_filename <- function(extrap, outdir) {
    if (extrap == TRUE) {
        Trns_CA_filename = paste0(outdir, "predict_gf_Trns_CAgrid_extrap.RData")
    } else { 
        if (extrap == FALSE) {
            Trns_CA_filename = paste0(outdir, "predict_gf_Trns_CAgrid_noextrap.RData")
        } else {
            stop("Wrong extrap value!")
        }
    } 
    return(Trns_CA_filename)
}

# main --------
dir.create(outdir, recursive = FALSE, showWarnings = FALSE)

# load the split transform files
allid = sprintf("%02d", 0:50)  
mylist = lapply(allid, function(ii) {
    rdata  = getTrns_grid_filename(extrap, ii, rdatadir)
    print(paste(date(), "Loading ...", rdata))
    Trns_grid = loadRData(rdata)
    return(Trns_grid)
})

# combine the files 
print(paste(date(), "Combining ..."))
Trns_CA = data.table::rbindlist(mylist)
print(dim(Trns_CA)) # check dimension 
# SHOULD BE: 50449570 35

# output --------
Trns_CA_filename <- getTrns_CA_filename(extrap, outdir)
print(paste(date(), "Exporting ...", Trns_CA_filename))
save(Trns_CA, file = Trns_CA_filename)
print(paste(date(), "All done. Job finished successfully."))




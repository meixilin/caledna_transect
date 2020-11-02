# Title: envfit calculation --------
# Author: Meixi Lin
# Date: Mon Aug  5 19:37:31 2019
# Author:
# Date:
# Modification:

# preparation --------
rm(list = ls())
cat("\014")

setwd("~/Lab/caledna_transect/")
library(dplyr)
library(vegan)
library(phyloseq)
library(ggplot2)
date() # the execution date

source("./scripts/function_transect.R")

# # source function to correct for multiple testing
# source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/p.adjust.envfit.r')

# set some variable --------
outdir = "./derive_data/step4_beta_diver/"
sink(paste0(outdir, "logs/pcoa_envfit_calc", Sys.Date(), ".log"))

mysig <- 0.05
# here I will only use major habitat as the category for plotting
groupv = "majorhab"

# load data -------
load(file = "./derive_data/phy_rare/phyrare.RData")
load(file = "./derive_data/step4_beta_diver/jadiss.RData")

# start a list to store envfit result -------
eflist <- vector(length = length(primers_commeco), mode = "list")

# for each variables --------
for (ii in 1:length(primers_commeco)) {
    # get phyloseq objects 
    primer <- primers_commeco[ii]
    physeq <- phyrare[[primer]]
    mydiss <- jadiss[[primer]]
    ord.res <- stats::cmdscale(d = mydiss, eig = T)
    print(primer)
    print(physeq)
    
    # get sample dataframe 
    sampledf <- data.frame(sample_data(physeq)) %>%
        dplyr::select(contlist) 
    
    # now start the envfit 
    ef <- envfit(ord.res, sampledf, permu = 1999, na.rm = T)
    print(ef)
    eflist[[ii]] <- ef
    names(eflist)[ii] <- primer
}

# store the ef result -------
# Please note that the p.value will change slightly since no seed was set
efdir = paste0(outdir, "envfit_ordisurf/")
dir.create(efdir, recursive = T)
save(eflist, file = paste0(efdir, "envfit_pcoa_",Sys.Date(),".RData"))

# end ========
sink()
closeAllConnections()
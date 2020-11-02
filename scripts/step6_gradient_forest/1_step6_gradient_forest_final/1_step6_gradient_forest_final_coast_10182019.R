# Title:  Perform final gradient forest analysis on all primer sets with tested parameters  --------
# Delete entire coastal sites 
# Author: Meixi Lin
# Date: Fri Oct 18 23:17:57 2019
# Modification: Clean up 

# preparation --------
rm(list = ls())
cat("\014")

# setwd("~/Lab/caledna_transect/")
setwd("/u/project/rwayne/meixilin/abiotic_transect/")
# load packages
library(dplyr)
library(phyloseq)
library(gradientForest)
library(ggplot2)
source("./scripts/function_transect.R")
functions.file <- list.files("./r_codes/step4_gradient_forest/functions/")
lapply(functions.file, function(xx) {source(file = paste0("./r_codes/step4_gradient_forest/functions/", xx))})

# load and preprocess dataset --------
####################
# all
primer <- primers[6]
bothokflag <- FALSE
txlevel <- "Family"
gfmethod <- "Presence"
abun <- 2 # minimum sequence abundance per site per family 
cutoff <- 0.05 # minimum occurrence of sites per family as a percentage of the total number of sites (this is a two sided cutoff, as in, the sample should not occur too often, but also not to )
droptopo <- F # should it drop aspect, CTI and DAH
ntrees <- 2000
my.seed <- 17
###################

# for ONLY "primer == all" parameters --------
physeq <- loadRData(paste0("./derive_data/phy_deco_", primer, ".RData"))
physeq <- subset_samples(physeq, transect == "Coastal")
gfinput <- gfprep(physeq = physeq, glomtax = txlevel, metavar = myvar,
                  subsettax = NULL, subsettarget = NULL,
                  bothok = bothokflag, droptopo = droptopo, abun = abun, cutoff = cutoff)
aa <- Sys.time()
# perform gradient forest 
gf <- gfdefault(Phys_site = gfinput[[1]], Sp_mat = gfinput[[2]], ntrees = ntrees, lev = gfinput[[3]], seed = my.seed)
# write out the data
filename <- paste("coast", primer, txlevel, gfmethod, ntrees, abun, cutoff, droptopo, my.seed, sep = "_")
directory <- "./derive_data/step4_gradient_forest/2_final/gf_"
save(gf, file = paste0(directory, filename, ".RData"))
print(Sys.time() - aa)

# here you evaluate what is the most important predicting factor --------
# get basic parameters for gf
print_gf_summary(gf)


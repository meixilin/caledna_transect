# Title: Correlation test and something else  --------
# Author: Meixi Lin
# Date: Oct 26 2018
# Author: Meixi Lin
# Date: Tue Mar 17 13:56:46 2020
# Modification: finish up on generating the code used 

# preparation --------
rm(list = ls())
cat("\014")

setwd("~/Lab/caledna_transect")
library(ggplot2)
library(ggcorrplot)
date() # the execution date

# read data ---------
load("./final_data/Final_metadata.RData")
myvar <- c("Longitude", "Latitude", "hfp", 
           "bio1", "bio2", "bio3", "bio4", "bio5", "bio6",
           "bio7", "bio8", "bio9", "bio10", "bio11", "bio12",
           "bio13", "bio14", "bio15", "bio16", "bio17","bio18", "bio19",
           "phihox", "orcdrc", "cecsol", "sndppt", "clyppt", "bldfie", "ntot",
           "elev", "Slope", "aspect","CTI", "Roughness", "Ruggedness", "DAH",
           "B1", "B2", "B3", "B4", "B5", "B6", "B7",
           "B8", "B8A", "B9", "B10", "B11", "B12",
           "NDVIS2", "NDVI32", "EVI", "NBRT", "greenness",
           "imprv", "ptrcv")

# calculate correlation --------
corr <- matrix(nrow = length(myvar), ncol = length(myvar), dimnames = list(myvar, myvar))
for (ii in 1:length(myvar)) {
    for (jj in 1:length(myvar)) {
        corr[ii,jj] <- cor(biom[,myvar[ii]], biom[,myvar[jj]], use = "complete.obs")
    }
}

save(corr, file = "./derive_data/step0_prepare_data/metadata_correlation.RData")
p1 <- ggcorrplot(corr)
corr.a <- abs(corr)
p2 <- ggcorrplot(corr.a)
# This is figure S1
ggsave(filename = "./plots/step0_prepare_data/metadata_correlation.pdf", plot = p1, width = 10, height = 10, device = "pdf")
ggsave(filename = "./plots/step0_prepare_data/metadata_correlation_abs.pdf", plot = p2, width = 10, height = 10, device = "pdf")


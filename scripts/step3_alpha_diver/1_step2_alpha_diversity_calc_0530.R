# Title: basic alpha diversity calculation --------
# Author: Meixi Lin
# Date: Tue Jun 26 15:12:01 2018
# Author: Meixi Lin
# Date: Fri Jul 27 20:36:37 2018
# Modification:

# preparation --------
rm(list = ls())
cat("\014")

setwd("~/Lab/caledna_transect")
library(vegan)
library(tibble)
library(dplyr)
library(ggplot2)
library(phyloseq)
library(microbiomeSeq)

source("./scripts/function_transect.R")

date() # the execution date

# import phyloseq object --------
load("./derive_data/phyrare.RData")

# make a list for data storage --------
adivrare <- vector("list", length = length(phyrare))
names(adivrare) <- names(phyrare)

mymeasures <- c("Observed", "Shannon")

# calculate alpha diversity for rarefied --------
for (ii in 1:length(phyrare)) {
    psdata <- phyrare[[ii]]
    psdata
    primer <- names(phyrare)[ii]
    adiv <- phyloseq::estimate_richness(physeq = psdata, measures = mymeasures)
    adivrare[[ii]] <- adiv
    write.csv(x = adiv, file = paste0("./derive_data/step2_alpha_diver/adivrare/", primer, "_adiv_rare_all_meas_06012019.csv"))
}

save(adivrare, file = "./derive_data/step2_alpha_diver/adivrare/adivrare_06012019.RData")


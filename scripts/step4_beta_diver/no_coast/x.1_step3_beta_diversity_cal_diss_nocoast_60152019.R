# Title: calculate beta diversity measures --------
# Author: Meixi Lin
# Date: Wed Aug  7 20:38:57 2019
# Author:
# Date:
# Modification:

# preparation --------
rm(list = ls())
cat("\014")

setwd("~/UCLA/Lab/abiotic_transect/")
library(dplyr)
library(vegan)
library(phyloseq)
library(ggplot2)
date() # the execution date

source("./scripts/function_transect.R")

# load the data -------
load(file = "./derive_data/phyrare.RData")
phyrare
phyrare <- lapply(phyrare, function(xx) {
    table(sample_sums(xx))
    xx <- subset_samples(physeq = xx, transect != "Coastal")
    return(xx)
    })
phyrare

# rank indices --------
# looks like bray and jaccard are equally good 
lapply(phyrare, function(phy) {
    spmat <- t(as.data.frame(otu_table(phy)))
    physsite <-  as.data.frame(sample_data(phy))[,myvar]
    
    print(rankindex(grad = scale(physsite), veg = spmat, indices = c("euc","man","bray","jac","kul")))
})

# calculate jaccard dissimilarity for all phyloseq objects and store ----------
jadiss_nocoast <- lapply(phyrare, function(phy) {
    yy <- phyloseq::distance(physeq = phy, method = "jaccard", binary = T)
})
save(jadiss_nocoast, file = "./derive_data/step3_beta_diver/jadiss_nocoast.RData")

# perform adonis -------

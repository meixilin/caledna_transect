# Title: calculate beta diversity measures --------
# Author: Meixi Lin
# Date: Sat Jun 15 01:00:09 2019
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

# load the data -------
load(file = "./derive_data/phy_rare/phyrare.RData")
phyrare
lapply(phyrare, function(xx) {table(sample_sums(xx))})

# rank indices --------
# looks like bray and jaccard are equally good 
lapply(phyrare, function(phy) {
    spmat <- t(as.data.frame(otu_table(phy)))
    physsite <-  as.data.frame(sample_data(phy))[,contlist]
    # not binary
    print(rankindex(grad = scale(physsite), veg = spmat, indices = c("euc","man","bray","jac","kul")))
    # binary 
    spmat2 = apply(spmat, 2, function(xx) {
        xx = as.integer(as.logical(xx))
        return(xx)
    })
    print(rankindex(grad = scale(physsite), veg = spmat2, indices = c("euc","man","bray","jac","kul")))
})

# calculate jaccard dissimilarity for all phyloseq objects and store ----------
jadiss <- lapply(phyrare, function(phy) {
    yy <- phyloseq::distance(physeq = phy, method = "jaccard", binary = T)
})
save(jadiss, file = "./derive_data/step3_beta_diver/jadiss.RData")
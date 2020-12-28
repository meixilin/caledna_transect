# Title: generate summary statistics --------
# Author: Meixi Lin
# Date: Mon Jan  6 00:18:24 2020

# preparation --------
rm(list = ls())
cat("\014")

setwd("~/Lab/caledna_transect/")
# load packages
library(dplyr)
library(ggplot2)
library(reshape2)
library(phyloseq)
library(ranacapa)
source("./scripts/function_transect.R")

# load data --------
load(file = "./derive_data/phy_deco/phy_deco_all.RData")
load(file = "./raw_data/transect_ASV/taxonomy_detail_all.RData")

# get total reads --------
sum(sample_sums(phy_deco_all))
# 16,157,425

# get mean taxonomic entry --------
test <- phyloseq::transform_sample_counts(phy_deco_all, function(abun){1*(abun > 0)}) 
sample_sums(test)
mean(sample_sums(test)) 

# In discussion: 
# For instance, the CALeDNA program samples each contained an average of 778 identified taxonomic lineages.

# get the total amount of phyla --------
phylum <- unique(tax_table(phy_deco_all)@.Data[,'Phylum'])
# 88 - "" - "NA" = 86

# get the raw read depth --------
get_raw_depth <- function(taxdetail) {
    asv = taxdetail$asv %>%
        dplyr::select(starts_with("K0")) %>%
        colSums()
    return(summary(asv))
}
lapply(asvdb, get_raw_depth)

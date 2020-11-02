# Title: compare ucnrs sites --------
# Author: Meixi Lin
# Date: Fri May 17 15:02:56 2019
# Author:
# Date:
# Modification:

# preparation --------
rm(list = ls())
cat("\014")

setwd("~/Lab/caledna_transect/")
library(dplyr)
library(phyloseq)
library(ggplot2)
date() # the execution date
source("./scripts/function_transect.R")

# define variables --------
myprimer <- "all"
myphylum <- c("Chordata", "Arthropoda", "Streptophyta")
mylev <- c("Order", "Family", "Genus")
mycutoff <- 0
outdir <- "./derive_data/step2_data_description/ucnrs/bysite/"
dir.create(outdir, recursive = T)

# load data -------
fauna <- read.csv(file = "./final_data/ucnrs/UCNRS_fauna_list.csv", stringsAsFactors = F)
flora <- read.csv(file = "./final_data/ucnrs/UCNRS_flora_list.csv", stringsAsFactors = F)
load(file = "./derive_data/phy_deco/phydeco_uc.RData")

psdata <- phydeco_uc[[myprimer]]
myloc <- sample_data(psdata)$loc %>% unique() %>% as.character()
# 
# sort(table(sample_data(psdata)$loc))
# # fort ord had most samples (26)

# loop through all sites --------
for (ii in myloc) {
    # this is preliminary in trying to keep the distinct row per reserve, but each Species, Family entries could still overlap, not keeping the replicate count since it is mostly database error
    dtfauna <- fauna[fauna$loc == ii,] %>% 
        dplyr::distinct()
    dtflora <- flora[flora$loc == ii,] %>% 
        dplyr::distinct()
    
    dtedna <- subset_samples(psdata, loc == ii)
    leftotu <- names(taxa_sums(dtedna))[(taxa_sums(dtedna) > mycutoff)]
    dtedna <- prune_taxa(leftotu, x = dtedna)
    dtedna <- subset_taxa(dtedna, Phylum %in% myphylum)
    
    dtedna_otu <- otu_table(dtedna)@.Data %>% as.data.frame()
    dtedna_tax <- tax_table(dtedna)@.Data %>% as.data.frame()
    
    write.csv(dtfauna, file = paste0(outdir, ii, "_fauna.csv"), row.names = F)
    write.csv(dtflora, file = paste0(outdir, ii, "_flora.csv"), row.names = F)
    write.csv(dtedna_otu, file = paste0(outdir, ii, "_edna_otu.csv"), row.names = F)
    write.csv(dtedna_tax, file = paste0(outdir, ii, "_edna_tax.csv"), row.names = F)
}

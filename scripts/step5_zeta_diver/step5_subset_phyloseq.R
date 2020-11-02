# Title: subset phyloseq by taxonomy abundance at class level --------
# Author: Meixi Lin
# Date: Tue Jun 26 15:12:01 2018
# Author: Meixi Lin
# Date: Fri Jul 27 20:36:37 2018
# Modification:

# preparation --------
rm(list = ls())
cat("\014")

setwd("/Users/linmeixi/UCLA/Lab/abiotic_transect")
library(tibble)
library(dplyr)
library(ggplot2)
library(phyloseq)
library(microbiomeSeq)

source("./r_codes/function_beta_0725.R")

date() # the execution date

# delete the primers 12S
primers <- primers[-1]

# subset by each class level --------
for (primer in primers) {
    # load different primers
    physeq <- get(load(paste0("./derive_data/phy_deco_", primer, ".RData")))
    # get taxonomic info for each asv at certain level, here, class
    taxinfo <- as.character(tax_table(physeq)[,'Class'])
    # get taxonomy label 
    taxlist <- unique(taxinfo)
    taxlist <- taxlist[!is.na(taxlist)]
    # get approximate summary table
    classsum <- sort(table(taxinfo)); print(classsum)
    # get how many NA is in the database 
    classna <- sum(is.na(taxinfo))
    # TO LEVI: you can change which class you want to by the asv abundance, here I subsetted all 
    for (ii in taxlist) {
        new_phy <- subset_taxa(physeq, Class == ii)
        otu_phy <- as.data.frame(otu_table(new_phy))
        write.csv(otu_phy, file = paste0("./derive_data/levi/asv_", primer, "_",ii,".csv"))
    }
}

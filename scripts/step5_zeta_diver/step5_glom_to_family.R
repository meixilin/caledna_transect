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

# define taxonomic level 
glomtax <- "Family"

# subset by each family level (all families) --------
for (primer in primers) {
    # load different primers
    physeq <- get(load(paste0("./derive_data/phy_deco_", primer, ".RData")))
    # get taxonomic info for each asv at certain level, here, class
    taxinfo <- as.character(tax_table(physeq)[, glomtax])
    # get taxonomy label 
    taxlist <- unique(taxinfo)
    taxlist <- taxlist[!is.na(taxlist)]
    # get approximate summary table
    taxsum <- sort(table(taxinfo)); print(taxsum)
    # get how many NA is in the database 
    taxna <- sum(is.na(taxinfo)); print(taxna)
    # TO LEVI: you can change which class you want to by the asv abundance, here I subsetted all 
    for (ii in taxlist) {
        new_phy <- subset_taxa(physeq, Family == ii)
        otu_phy <- as.data.frame(otu_table(new_phy))
        write.csv(otu_phy, file = paste0("./derive_data/levi/asv_fam_", primer, "_",ii,".csv"))
    }
}

# glom to family level and not subsetting each family --------
for (primer in primers) {
    # load different primers
    physeq <- get(load(paste0("./derive_data/phy_deco_", primer, ".RData")))
    # get them to family level
    physeq <- glom_tax(physeq, taxlevel = glomtax)
    otu_phy <- as.data.frame(otu_table(physeq))
    write.csv(otu_phy, file = paste0("./derive_data/levi/asv_fam_", primer, "_all", ".csv"))
}

# exactly what I used for gradient forest analysis --------
# here is the link to data 
# https://github.com/meixilin/abiotic_transect/blob/master/derive_data/new_gf/gf_Y_deco_10_all_fam.csv


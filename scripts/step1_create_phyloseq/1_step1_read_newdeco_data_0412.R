# Title: Import data in phyloseq and specify factor files --------
# Author: Meixi Lin
# Date: Sun Oct 14 17:43:00 2018
# Author:
# Date:
# Modification:

# preparation --------
rm(list = ls())
cat("\014")

setwd("~/Lab/caledna_transect")
library(dplyr)
library(ranacapa)
library(phyloseq)
library(reshape2)

source("./scripts/function_transect.R")

date() # the execution date


# import data and biom --------
asv.16S <- read.csv("./final_data/deco_3/asv_deco_dedup_16S.csv")
asv.18S <- read.csv("./final_data/deco_3/asv_deco_dedup_18S.csv")
asv.CO1 <- read.csv("./final_data/deco_3/asv_deco_dedup_CO1.csv")
asv.FITS <- read.csv("./final_data/deco_3/asv_deco_dedup_FITS.csv")
asv.PITS <- read.csv("./final_data/deco_3/asv_deco_dedup_PITS.csv")

load(file = "./final_data/Final_metadata.RData")

# convert data to phyloseq --------
# phy_deco_12S <- convert_anacapa_to_phyloseq(asv.12S, biom)
phy_deco_16S <- convert_anacapa_to_phyloseq(asv.16S, biom)
phy_deco_18S <- convert_anacapa_to_phyloseq(asv.18S, biom)
phy_deco_CO1 <- convert_anacapa_to_phyloseq(asv.CO1, biom)
phy_deco_FITS <- convert_anacapa_to_phyloseq(asv.FITS, biom)
phy_deco_PITS <- convert_anacapa_to_phyloseq(asv.PITS, biom)

# save(phy_deco_12S, file = "./derive_data/phy_deco_12S.RData")
save(phy_deco_16S, file = "./derive_data/phy_deco/phy_deco_16S.RData")
save(phy_deco_18S, file = "./derive_data/phy_deco/phy_deco_18S.RData")
save(phy_deco_CO1, file = "./derive_data/phy_deco/phy_deco_CO1.RData")
save(phy_deco_FITS, file = "./derive_data/phy_deco/phy_deco_FITS.RData")
save(phy_deco_PITS, file = "./derive_data/phy_deco/phy_deco_PITS.RData")

# merge phyloseq objects ---------
phy_deco_all <- phyloseq::merge_phyloseq(phy_deco_16S, phy_deco_18S, phy_deco_CO1, phy_deco_FITS, phy_deco_PITS)
save(phy_deco_all, file = "./derive_data/phy_deco/phy_deco_all.RData")

# make a list --------
phydeco <- list(phy_deco_16S, phy_deco_18S, phy_deco_CO1, phy_deco_FITS, phy_deco_PITS, phy_deco_all)
names(phydeco) <- primers
save(phydeco, file = "./derive_data/phy_deco/phydeco.RData")

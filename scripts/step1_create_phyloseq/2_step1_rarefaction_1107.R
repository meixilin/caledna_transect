# Title: rarefy decontaminated phyloseq --------
# Author: Meixi Lin
# Date: Tue Jun 26 15:12:01 2018
# Author: Meixi Lin
# Date: Wed Nov 7 14:10:00 2018
# Modification: use new decontaminated datasets
# Author: Meixi Lin
# Date: Mon Jan 28 10:19:05 2019
# Modification: rarefy by different depth 
# NOTE: the phyrare also included a rarefied "all" primers combined dataset that was not used. 

# preparation --------
rm(list = ls())
cat("\014")

setwd("~/Lab/caledna_transect/")
library(tibble)
library(dplyr)
library(ggplot2)
library(phyloseq)
library(microbiomeSeq)

source("./scripts/function_transect.R")

date() # the execution date

# import phyloseq object --------
load("./derive_data/phy_deco/phy_deco_all.RData")
load("./derive_data/phy_deco/phy_deco_16S.RData")
load("./derive_data/phy_deco/phy_deco_18S.RData")
load("./derive_data/phy_deco/phy_deco_CO1.RData")
load("./derive_data/phy_deco/phy_deco_FITS.RData")
load("./derive_data/phy_deco/phy_deco_PITS.RData")

# rarefy samples --------
rarefaction_reps  <- 10
# myseeds 
# seeds <- as.integer(runif(rarefaction_reps, min = 0, max = 1000)) 
seeds <- c(825,598,744,276,705,872,895,901,759,544)

phy_10k_all <- custom_rarefaction_2(phy_deco_all, sample_size = 10000, replicates = rarefaction_reps, seeds = seeds)
phy_2000_16S <- custom_rarefaction_2(phy_deco_16S, sample_size = 2000, replicates = rarefaction_reps, seeds = seeds)
phy_4000_18S <- custom_rarefaction_2(phy_deco_18S, sample_size = 4000, replicates = rarefaction_reps, seeds = seeds)
phy_1000_CO1 <- custom_rarefaction_2(phy_deco_CO1, sample_size = 1000, replicates = rarefaction_reps, seeds = seeds)
phy_4000_FITS <- custom_rarefaction_2(phy_deco_FITS, sample_size = 4000, replicates = rarefaction_reps, seeds = seeds)
phy_1000_PITS <- custom_rarefaction_2(phy_deco_PITS, sample_size = 1000, replicates = rarefaction_reps, seeds = seeds)

# write out the rarefied object --------
save(phy_10k_all, file = "./derive_data/phy_rare/phy_10k_all.RData")
save(phy_2000_16S, file = "./derive_data/phy_rare/phy_2000_16S.RData")
save(phy_4000_18S, file = "./derive_data/phy_rare/phy_4000_18S.RData")
save(phy_1000_CO1, file = "./derive_data/phy_rare/phy_1000_CO1.RData")
save(phy_4000_FITS, file = "./derive_data/phy_rare/phy_4000_FITS.RData")
save(phy_1000_PITS, file = "./derive_data/phy_rare/phy_1000_PITS.RData")

# write out the rarefied data frame -------
write.csv(x = phy_10k_all@otu_table@.Data, file = "./derive_data/step1_create_phyloseq/physeq_rare_otu/phy_10k_all.csv")
write.csv(x = phy_2000_16S@otu_table@.Data, file = "derive_data/step1_create_phyloseq/physeq_rare_otu/phy_2000_16S.csv")
write.csv(x = phy_4000_18S@otu_table@.Data, file = "derive_data/step1_create_phyloseq/physeq_rare_otu/phy_4000_18S.csv")
write.csv(x = phy_1000_CO1@otu_table@.Data, file = "derive_data/step1_create_phyloseq/physeq_rare_otu/phy_1000_CO1.csv")
write.csv(x = phy_4000_FITS@otu_table@.Data, file = "derive_data/step1_create_phyloseq/physeq_rare_otu/phy_4000_FITS.csv")
write.csv(x = phy_1000_PITS@otu_table@.Data, file = "derive_data/step1_create_phyloseq/physeq_rare_otu/phy_1000_PITS.csv")

# make a list of it -------
phyrare <- list(phy_2000_16S, phy_4000_18S, phy_1000_CO1, phy_4000_FITS, phy_1000_PITS, phy_10k_all)
names(phyrare) <- primers
save(phyrare, file = "./derive_data/phy_rare/phyrare.RData")


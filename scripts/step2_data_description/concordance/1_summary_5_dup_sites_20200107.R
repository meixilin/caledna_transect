# Title: generate R object for duplicated sites --------
# Author: Meixi Lin
# Date: Tue Jan  8 12:49:38 2019
# Author:
# Date: Mon Jan 28 10:19:05 2019
# Modification: add 10 random seeds 

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

seeds <- c(825,598,744,276,705,872,895,901,759,544)

outdir <- "./derive_data/step2_data_description/concordance/"
dir.create(outdir, recursive = T)

# read in data --------
asv.16S <- read.csv("./other_version_data/Oct_25_2018_CALeDNA_Transect_Results/16S_table_3_Oct25_2018.csv")
asv.18S <- read.csv("./other_version_data/Oct_25_2018_CALeDNA_Transect_Results/18S_table_3_Oct25_2018.csv")
asv.CO1 <- read.csv("./other_version_data/Oct_25_2018_CALeDNA_Transect_Results/CO1_table_3_Oct25_2018.csv")
asv.FITS <- read.csv("./other_version_data/Oct_25_2018_CALeDNA_Transect_Results/FITS_table_3_Oct25_2018.csv")
asv.PITS <- read.csv("./other_version_data/Oct_25_2018_CALeDNA_Transect_Results/PITS_table_3_Oct25_2018.csv")

# rename scrub data --------
matchname <- as.matrix(colnames(asv.16S)[-1])
newnames <- apply(matchname, 1, function(x) {
    if (!is.na(colsplit(x, "[.]",1:2))[,2]) {
        x <- paste0(colsplit(x, "[.]",1:2)[,1], colsplit(x, "[.]",1:2)[,2])
    } else {
        x <- x
    }
    return(x)
})
newnames <- c("sum.taxonomy", unlist(newnames))
colnames(asv.16S) <- newnames
colnames(asv.18S) <- newnames
colnames(asv.CO1) <- newnames
colnames(asv.FITS) <- newnames
colnames(asv.PITS) <- newnames

droptran <- colsplit(newnames[-1], "[_]", 1:2)[,1] # should be 283


# find the samples that were sequenced twice --------
dup <- droptran[duplicated(droptran)]
primers <- c("16S", "18S", "CO1", "FITS", "PITS")

# make the phyloseq objects -------
asvlist <- list(asv.16S, asv.18S, asv.CO1, asv.FITS, asv.PITS)
names(asvlist) <- primers
phylist <- lapply(asvlist, convert_asv_to_phyloseq)
# add sample data for describe 5 duplicate sites 
phylist <- lapply(phylist, function(xx) {
    sampledata = sample_data(data.frame(
        MatchName = substr(sample_names(xx), 1, 7),
        indup = substr(sample_names(xx), 1, 7),
        Transect = substr(sample_names(xx), 8, 9), 
        row.names=sample_names(xx),
        stringsAsFactors=FALSE
    ))
    xx <- merge_phyloseq(xx, sampledata)
})

# merge the phylist 
phyall <- phyloseq::merge_phyloseq(phylist[[1]], phylist[[2]],  phylist[[3]], phylist[[4]], phylist[[5]])
phylist[[6]] <- phyall
primers <- c(primers, "all")
names(phylist) <- primers

# perform rarefaction on the entire phyloseq list and select the duplicated sites
# for rarefaction depth 
phylist.r <- vector("list", length = 6)
for (ii in 1:length(rare_depth)) {
    phylist.r[[ii]] <- custom_rarefaction_2(phylist[[ii]], sample_size = rare_depth[ii], replicates = 10, seeds = seeds)
}
names(phylist.r) <- primers

# get sample names and sample data 
dupsample <- c("K0058C2_S", "K0117C2_S", "K0142A2_S", "K0142C2_S", "K0177A1_S", "K0058C2_F", "K0117C2_F", "K0142A2_F", "K0142C2_F", "K0177A1_C")

subphy <- lapply(phylist, function(xx) {
    hh <- prune_samples(dupsample, x = xx)
})
names(subphy) <- primers

subphy.r <- lapply(phylist.r, function(xx) {
    hh <- prune_samples(dupsample, x = xx)
})
names(subphy.r) <- primers

# save the output 
save(phylist, file = paste0(outdir, "phydeco_wdup.RData"))
save(phylist.r, file = paste0(outdir, "phyrare_wdup.RData"))
save(subphy, file =paste0(outdir, "phy_dup.RData"))
save(subphy.r, file =paste0(outdir, "phy_dup_rare.RData"))


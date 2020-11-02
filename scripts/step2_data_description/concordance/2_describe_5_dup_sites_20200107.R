# Title: generate summary for duplicated sites --------
# Author: Meixi Lin
# Date: Mon Jan  6 00:50:15 2020

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
library(stringr)
source("./scripts/function_transect.R")

indir <- "./derive_data/step2_data_description/concordance/"
outdir <- "./derive_data/step2_data_description/concordance/raw_tables/"

dir.create(outdir, recursive = T)
# get sample names 
dupsample <- c("K0058C2_S", "K0117C2_S", "K0142A2_S", "K0142C2_S", "K0177A1_S", "K0058C2_F", "K0117C2_F", "K0142A2_F", "K0142C2_F", "K0177A1_C")
dup <- c("K0058C2", "K0117C2", "K0142A2", "K0142C2", "K0177A1")

# get arguments --------
args <-  commandArgs(trailingOnly=TRUE)
rare <- args[1] # rare/deco
cutoff <- as.integer(args[2])
txlevel <- as.character(args[3])

# define function --------
# NOTE: here, we deleted "absence" to reduce the denominator to taxonomy entries that occurred in either replicates
count_overlap <- function(physeq1, ii, prune = T) {
    if (!(ii %in% sample_data(physeq1)$MatchName)) {
        output <- "neither"
        names(output) <- "percent"
    } else {
        phy = subset_samples(physeq = physeq1, MatchName == ii)
        if (all(otu_table(object = phy)@.Data == 0)) {
            output <- "neither"
            names(output) <- "percent"
        } else {
            if (prune == T) {
                phy <- prune_taxa(taxa_sums(phy) > 0, phy)
            }
            mydata <- otu_table(object = phy)@.Data
            if (dim(mydata)[2] == 1) {
                output <- sample_names(phy)
                names(output) <- "percent"
            } else {
                true <- sum(mydata[,1] == mydata[,2])
                false <- sum(mydata[,1] != mydata[,2])
                percent <- true/(true + false); names(percent) = "percent"
                output <- c(true, false, percent)
            }
        }
    }
    return(output)
}

# for plotting
filter_chr <- function(xx) {
    if (str_detect(xx, "neither") | str_detect(xx, "K")) {
        xx <- NA
    } else {
        xx <- as.numeric(xx)
    }
    return(xx)
}

# read in data --------
if (rare == "rare") {
    rawdata <- loadRData(file = paste0(indir, "phy_dup_rare.RData"))
} else {
    if (rare == "deco") {
        rawdata <- loadRData(file = paste0(indir, "phy_dup.RData"))
    }
}

# add taxonomy summary --------
# Here species is essentially the full path 
if (txlevel != "Species") {
    subphy1 <- lapply(rawdata, function(xx, var1 = txlevel) {
        xx <- glom_tax(xx, taxlevel = var1)
    })
} else {
    subphy1 <- rawdata
}

# add cutoff -------
subphy <-  lapply(subphy1, function(xx, ct = cutoff) {
    xx <- transform_sample_counts(xx, function(abund) {1*(abund > ct)})
    # delete the taxonomy entries that didn't occur 
    # NOTE: here, we deleted "absence" to reduce the denominator to taxonomy entries that occurred in entire dataset 
    xx <- prune_taxa(taxa_sums(xx) > 0, xx)
    return(xx)
})

# evaluate overlapping entries --------
summary.percent <- data.frame(matrix(NA, nrow = length(dup), ncol = length(primers)))
colnames(summary.percent) <- primers
row.names(summary.percent) <- dup

for (ii in dup) {
    for (jj in primers) {
        summary.percent[ii, jj] <- count_overlap(subphy[[jj]], ii, prune = T)['percent']
    }
}

write.csv(summary.percent, file = paste0(outdir, txlevel, "_", rare, "_", cutoff, ".csv"))
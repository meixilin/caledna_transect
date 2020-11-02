# Title: load adonis and beta dispersion result --------
# Author: Meixi Lin
# Date: Mon Aug  5 15:37:40 2019
# Author:
# Date:
# Modification:


# preparation --------
rm(list = ls())
cat("\014")

setwd("~/Lab/caledna_transect")
date() # the execution date
library(dplyr)
library(phyloseq)
library(ggplot2)
library(vegan)

source("./scripts/function_transect.R")

# set some variables -------
mycutoff <- 4
indir = "./derive_data/step4_beta_diver/nocoast/"
outdir = indir
nnperm = 2999

# load the adonis and beta dispersion result --------
adres <- loadRData(paste0(indir, "adonis_result_nocoast_cutoff_", mycutoff, ".RData"))
bdspres <- loadRData(paste0(indir, "betadisper_result_nocoast_cutoff_", mycutoff, ".RData"))

# make a df ---------
addf <- data.frame(matrix(nrow = length(primers) * length(catlist), ncol = 9))
colnames(addf) <- c("name", "var.df", "res.df", "var.ssq", "res.ssq", "f.stat", "var.r2", "nperm","p.value")
bdspdf <- data.frame(matrix(nrow = length(primers) * length(catlist), ncol = 8))
colnames(bdspdf) <- c("name", "var.df", "res.df", "var.ssq", "res.ssq", "f.stat", "nperm", "p.value")

# write out the adonis --------
for (ii in 1:length(adres)) {
    xx <- adres[[ii]]
    addf[ii, "name"] <- names(adres)[ii]
    addf[ii, "var.df"] <- xx$aov.tab$Df[1]
    addf[ii, "res.df"] <- xx$aov.tab$Df[2]
    addf[ii, "var.ssq"] <- xx$aov.tab$SumsOfSqs[1]
    addf[ii, "res.ssq"] <- xx$aov.tab$SumsOfSqs[2]
    addf[ii, "f.stat"] <- xx$aov.tab$F.Model[1]
    addf[ii, "var.r2"] <- xx$aov.tab$R2[1]
    addf[ii, "nperm"] <- nnperm
    addf[ii, "p.value"] <- xx$aov.tab$`Pr(>F)`[1]
}

write.csv(addf, file = paste0(outdir, "adonis_nocoast_beta_div.csv"), row.names = F)

# write out the betadispersion ---------
for (ii in 1:length(bdspres)) {
    colnames(bdspdf) <- c("name", "var.df", "res.df", "var.ssq", "res.ssq", "f.stat", "nperm", "p.value")
    xx <- permutest(bdspres[[ii]], permutations = how(nperm = nnperm))
    bdspdf[ii, "name"] <- names(bdspres)[ii]
    bdspdf[ii, "var.df"] <- xx$tab$Df[1]
    bdspdf[ii, "res.df"] <- xx$tab$Df[2]
    bdspdf[ii, "var.ssq"] <- xx$tab$`Sum Sq`[1]  
    bdspdf[ii, "res.ssq"] <- xx$tab$`Sum Sq`[2]
    bdspdf[ii, "f.stat"] <- xx$tab$`F`[1]
    bdspdf[ii, "nperm"] <- xx$tab$N.Perm[1]
    bdspdf[ii, "p.value"] <- xx$tab$`Pr(>F)`[1]
}

write.csv(bdspdf, file = paste0(outdir, "beta_dispersion_nocoast_beta_div.csv"), row.names = F)


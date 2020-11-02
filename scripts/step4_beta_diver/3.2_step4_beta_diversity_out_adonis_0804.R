# Title: load adonis and beta dispersion result --------
# Author: Meixi Lin
# Date: Mon Aug  5 15:37:40 2019
# Author:
# Date:
# Modification:


# preparation --------
rm(list = ls())
cat("\014")

setwd("~/Lab/caledna_transect/")
date() # the execution date
library(dplyr)
library(phyloseq)
library(ggplot2)
library(vegan)

source("./scripts/function_transect.R")

# set some variables -------
mycutoff <- 4
mysig <- 0.05
# sysdate <- Sys.Date() # 2019-10-11
sysdate <- "2020-03-30"
nnperm <- 2999
nnobs <- 9 # significance correction
outdir = "./derive_data/step4_beta_diver/adonis_bdsp/"

# load the adonis and beta dispersion result --------
load(paste0("./derive_data/step4_beta_diver/adonis_bdsp/adonis_result_cutoff_", mycutoff, "_",sysdate, ".RData"))
load(paste0("./derive_data/step4_beta_diver/adonis_bdsp/betadisper_result_cutoff_", mycutoff, "_",sysdate, ".RData"))

# make a df ---------
addfcolnames <- c("Metabarcode","Variable", "var.df", "res.df", "var.ssq", "res.ssq", "f.stat", "var.r2", "nperm", "P.value", "Adj.P.value", "Significant")
addf <- data.frame(matrix(nrow = length(primers_commeco) * length(catlist), ncol = length(addfcolnames)))
colnames(addf) <- addfcolnames

bdspdfcolnames <-  c("Metabarcode","Variable", "var.df", "res.df", "var.ssq", "res.ssq", "f.stat", "nperm", "P.value", "Adj.P.value", "Significant")
bdspdf <- data.frame(matrix(nrow = length(primers_commeco) * length(catlist), ncol = length(bdspdfcolnames)))
colnames(bdspdf) <- bdspdfcolnames

# write out the adonis --------
for (ii in 1:length(adres)) {
    xx <- adres[[ii]]
    addf[ii, "Metabarcode"] <- base::strsplit(names(adres)[ii],"[|]")[[1]][1]
    addf[ii, "Variable"] <- base::strsplit(names(adres)[ii],"[|]")[[1]][2]
    addf[ii, "var.df"] <- xx$aov.tab$Df[1]
    addf[ii, "res.df"] <- xx$aov.tab$Df[2]
    addf[ii, "var.ssq"] <- xx$aov.tab$SumsOfSqs[1]
    addf[ii, "res.ssq"] <- xx$aov.tab$SumsOfSqs[2]
    addf[ii, "f.stat"] <- xx$aov.tab$F.Model[1]
    addf[ii, "var.r2"] <- xx$aov.tab$R2[1]
    addf[ii, "nperm"] <- nnperm
    addf[ii, "P.value"] <- xx$aov.tab$`Pr(>F)`[1]
    addf[ii, "Adj.P.value"] <- min(1, addf[ii, "P.value"]*nnobs)
    addf[ii, "Significant"] <-  (addf[ii, "Adj.P.value"] < mysig)
}

write.csv(addf, file = paste0(outdir, "adonis_beta_div_", sysdate, ".csv"), row.names = F, quote = F)

# write out the betadispersion ---------
for (ii in 1:length(bdspres)) {
    xx <- permutest(bdspres[[ii]], permutations = how(nperm = nnperm))
    bdspdf[ii, "Metabarcode"] <- base::strsplit(names(bdspres)[ii],"[|]")[[1]][1]
    bdspdf[ii, "Variable"] <- base::strsplit(names(bdspres)[ii],"[|]")[[1]][2]
    bdspdf[ii, "var.df"] <- xx$tab$Df[1]
    bdspdf[ii, "res.df"] <- xx$tab$Df[2]
    bdspdf[ii, "var.ssq"] <- xx$tab$`Sum Sq`[1]
    bdspdf[ii, "res.ssq"] <- xx$tab$`Sum Sq`[2]
    bdspdf[ii, "f.stat"] <- xx$tab$`F`[1]
    bdspdf[ii, "nperm"] <- xx$tab$N.Perm[1]
    bdspdf[ii, "P.value"] <- xx$tab$`Pr(>F)`[1]
    bdspdf[ii, "Adj.P.value"] <- min(1, bdspdf[ii, "P.value"]*nnobs)
    bdspdf[ii, "Significant"] <- (bdspdf[ii, "Adj.P.value"] < mysig)
}

write.csv(bdspdf, file = paste0(outdir, "beta_dispersion_beta_div_", sysdate, ".csv"), row.names = F, quote = F)

# note that the R2 is unadjusted and here is an interesting thread on interpretation:
# https://github.com/vegandevs/vegan/issues/316
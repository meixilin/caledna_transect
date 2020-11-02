# Title: look at the category of overlap --------
# Author: Meixi Lin
# Date: Mon Oct  7 14:34:44 2019
# Author:
# Date:
# Modification:

# preparation --------
rm(list = ls())
cat("\014")

setwd("~/UCLA/Lab/abiotic_transect/")
library(dplyr)
library(vegan)
library(phyloseq)
library(ggplot2)
date() # the execution date

source("./scripts/function_transect.R")

# load data -------
load("./final_data/Final_metadata.RData")

# write a function -------
table_factor <- function(var1, var2, biom, outdir) {
    hh = table(biom[,c(var1,var2)])
    write.csv(hh, file = paste0(outdir, var1, "_", var2, ".csv"))
    return(hh)
}
chisq_biom <- function(var1, var2, biom) {
    hh = chisq.test(biom[,var1], biom[,var2], simulate.p.value = T)
    return(hh)
}

fisher_biom <- function(var1, var2, biom) {
    hh = fisher.test(biom[,var1], biom[,var2], simulate.p.value = T)
    return(hh)
}

# set paramters ------
mycutoff = 0.01

# export --------
outdir <- "./derive_data/step0_prepare_data/corr_factor/"
dir.create(outdir, recursive = T)
table_factor("majorhab", "minorhab", biom, outdir)
table_factor("NLCD", "SoS", biom, outdir)

# start a correlation test matrix --------
chisq = data.frame(matrix(nrow = length(catlist), ncol = length(catlist) - 1))
row.names(chisq) = catlist
colnames(chisq) = catlist[-1]
for (ii in catlist) {
    for (jj in catlist[-1]) {
        hh = chisq_biom(ii, jj, biom)
        hh$data.name = paste(ii, jj)
        # print(hh)
        # chisq[ii, jj] = hh$p.value
        chisq[ii, jj] = hh$p.value > mycutoff
    }
}


# result: 
# only majorhab -- taxousda had a non-significant pair of 

# test of independence using fisher's exact test --------
# start a correlation test matrix 
fisher = data.frame(matrix(nrow = length(catlist), ncol = length(catlist) - 1))
row.names(fisher) = catlist
colnames(fisher) = catlist[-1]
for (ii in catlist) {
    for (jj in catlist[-1]) {
        hh = fisher_biom(ii, jj, biom)
        hh$data.name = paste(ii, jj)
        # print(hh)
        fisher[ii, jj] = hh$p.value
        # fisher[ii, jj] = hh$p.value > mycutoff
    }
}

write.csv(fisher, file = paste0(outdir, "fisher_exact.csv"))


# only test within habitat --------
catlist = catlist[c(2,3,4,5,9)]
fisher.1 = data.frame(matrix(nrow = length(catlist), ncol = length(catlist) - 1))
row.names(fisher.1) = catlist
colnames(fisher.1) = catlist[-1]
for (ii in catlist) {
    for (jj in catlist[-1]) {
        hh = fisher_biom(ii, jj, biom)
        hh$data.name = paste(ii, jj)
        # print(hh)
        fisher.1[ii, jj] = hh$p.value
        # fisher[ii, jj] = hh$p.value > mycutoff
    }
}

 vwrite.csv(fisher, file = paste0(outdir, "fisher_exact_habitat.csv"))


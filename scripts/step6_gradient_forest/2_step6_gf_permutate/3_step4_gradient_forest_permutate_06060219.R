#!/usr/bin/env Rscript

# Title:  shuffle gradient forest in hoffman2 --------
# Author: Meixi Lin
# Credit: Ryan Harrigan
# Date: Sun Sep 16 16:11:50 2018

# This script was used for reading a gradient forest object, shuffle the explaining variables and make predictions again.
# It was used as a "null model".

# Parent script:
# gradient_forest_tune_hoff_qsub.sh

# preparation --------
rm(list = ls())
cat("\014")

setwd("/u/project/rwayne/meixilin/abiotic_transect")
# setwd("/Users/linmeixi/UCLA/Lab/abiotic_transect") # this was for local mode
outdir <- "/u/flashscratch/m/meixilin/3_permutate/"
options(echo=TRUE)

# load packages
library(dplyr)
library(phyloseq)
library(gradientForest)

source("./scripts/function_transect.R")
functions.file <- list.files("./r_codes/step4_gradient_forest/functions/")
lapply(functions.file, function(xx) {source(file = paste0("./r_codes/step4_gradient_forest/functions/", xx))})

# get arguments
args <-  commandArgs(trailingOnly=TRUE)

# this script allow 3 inputs:
# args[1]: the absolute path to the original gradient forest object
# args[2]: output file pathes
# args[3]: number of replication
# args[4]: SGE_TASK_ID

load(as.character(args[1])) # all the gradient forests should be named as "gf"
outname <- as.character(args[2]) # this should be the name of the "gf" object tested
nrep <- as.integer(args[3])
sgeid <- as.integer(args[4])

# load and preprocess dataset --------
# the predictors
Phys_site <- gf$X
dim(Phys_site)

# the response
Sp_mat <- gf$Y
dim(Sp_mat)

# the model parameters
ntrees <- gf$ntree
# ntrees = 10

lev <- gf_lev(Sp_mat)
lev

# do the random runs --------
for (ii in 1:nrep) {
    # shuffle the predictors
    predsR <- Phys_site[sample(nrow(Phys_site)),]
    # repeat the gradient forest predictions
    speciesforestR <- gfdefault(predsR, Sp_mat, ntrees, lev)

    randtotal <- speciesforestR$species.pos.rsq
    randaverage <- sum(speciesforestR$result)/randtotal

    write.table(randtotal,file=paste0("./derive_data/step4_gradient_forest/3_permutate/", outname, "_randtotal.tsv"),row.names=FALSE,col.names=FALSE,append=T)
    write.table(randaverage,file=paste0("./derive_data/step4_gradient_forest/3_permutate/", outname, "_randaverage.tsv"),row.names=FALSE,col.names=FALSE,append=T)
    
    dir.create(path = paste0(outdir, outname,"/"))
    save(speciesforestR, file = paste0(outdir, outname,"/", sgeid, "_", ii,"_permu.RData"))
}

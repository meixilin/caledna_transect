# Title: Plot the result from "step4_gradient_forest_1018" --------
# Author: Meixi Lin
# Date: Mon Sep  3 11:55:04 2018
# Author: Meixi Lin
# Date:
# Modification: Pair with the step4_gradient_forest_1018.R

# preparation --------
options(echo = TRUE)
rm(list = ls())
cat("\014")

setwd("~/Lab/caledna_transect/")
# setwd("/u/home/m/meixilin/project-rwayne/abiotic_transect/")
# load packages
library(dplyr)
library(phyloseq)
library(gradientForest)
source("./scripts/function_transect.R")
functions.file <- list.files("./r_codes/step4_gradient_forest/functions/")
lapply(functions.file, function(xx) {source(file = paste0("./r_codes/step4_gradient_forest/functions/", xx))})

# get arguments --------
args <-  commandArgs(trailingOnly=TRUE)

# this script allow 3 inputs:
# args[1]: the absolute path to the original gradient forest object
# args[2]: filename
# args[3]: plotting directory

load(as.character(args[1])) # RENAME WHAT EVER GF to "gf"
filename <- as.character(args[2]) # this should be the name of the "gf" object tested
plotdir <- as.character(args[3])

# # get basic parameters --------
# print(gf)
# print(summary(gf))
# # print(gf$X)
# most_imp <- importance(gf)[1:5]
# most_imp_write <- t(c(most_imp, filename))
# write.table(most_imp_write, file = "./derive_data/new_gf/most_imp_1107.tsv", sep = "\t", append = T, row.names = F, col.names = T)
# rsq <- gf_res(gf)
# rsq_write <- t(c(rsq, filename))
# write.table(rsq_write, file = "./derive_data/new_gf/rsq_1107.tsv", sep = "\t", append = T, row.names = F, col.names = F)

# plot the standard plots --------
# (accuracy importance for each parameters) =================================
dir.create(plotdir)
pdf(file = paste0(plotdir, filename, "_O", ".pdf"), height = 6, width = 8)
plot(gf, plot.type = "O")
dev.off()

# (importance density for 8 most important parameters) =================================
most_important <- names(importance(gf))[1:8]
pdf(file = paste0(plotdir, filename, "_S", ".pdf"), height = 6, width = 8)
par(mgp = c(2, 0.75, 0))
plot(gf, plot.type = "S", imp.vars = most_important,
     leg.posn = "topright", cex.legend = 1, cex.axis = 0.6,
     cex.lab = 0.7, line.ylab = 0.9,
     par.args = list(mgp = c(1.5,0.5, 0), mar = c(3.1, 1.5, 0.1, 1)))
dev.off()

# (species cumulative plot) =================================
# These show cumulative change in abundance of individual species, where changes occur on the gradient, and the species changing most on each gradient.
pdf(file = paste0(plotdir, filename, "_C", ".pdf"), height = 8, width = 6)
plot(gf, plot.type = "C", imp.vars = most_important,
     show.overall = F, legend = T, leg.posn = "topleft",
     leg.nspecies = 5, cex.lab = 0.7, cex.legend = 1,
     cex.axis = 0.6, line.ylab = 0.9,
     par.args = list(mgp = c(1.5,0.5, 0), mar = c(2.5, 1, 0.1, 0.5), omi = c(0, 0.3, 0, 0)))
dev.off()

# (predictor cumulative plot, not show speciess) =================================
# predictor cumulative plot (plot.type="C", show.species=F),
# which for each predictor shows cumulative importance distributions of splits improvement scaled by R2 weighted importance, and standardised by density of observations, averaged over all species.
pdf(file = paste0(plotdir, filename, "_C2", ".pdf"), height = 8, width = 6)
plot(gf, plot.type = "C", imp.vars = most_important,
     show.species = F, common.scale = T, cex.axis = 0.6,
     cex.lab = 0.7, line.ylab = 0.9,
     par.args = list(mgp = c(1.5,0.5, 0), mar = c(2.5, 1, 0.1, 0.5), omi = c(0,0.3, 0, 0)))
dev.off()

# (ordered R2 for each species) =================================
# the R2 measure of the fit of the random forest model for each species, ordered in various ways.
pdf(file = paste0(plotdir, filename, "_P", ".pdf"), height = 6, width = 10)
plot(gf, plot.type = "P", show.names = T, horizontal = F,
     cex.axis = 1, cex.labels = 0.6, line = 2.5, ncutoff = 50)
dev.off()

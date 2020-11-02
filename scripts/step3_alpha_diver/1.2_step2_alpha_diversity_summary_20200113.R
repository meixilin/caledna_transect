# Title: alpha diversity descriptive statistics --------
# Author: Meixi Lin
# Date: Tue May  7 12:34:27 2019
# Author:
# Date:
# Modification:

# preparation --------
rm(list = ls())
cat("\014")
setwd("~/Lab/caledna_transect")

# create output directory
outdir <- "./derive_data/step3_alpha_diver/description/"
dir.create(outdir, recursive = T)
sink(file = paste0(outdir, "summary_stats.log"))

library(vegan)
library(tibble)
library(dplyr)
library(ggplot2)
library(phyloseq)
require(ggpubr)

source("./scripts/function_transect.R")

date() # the execution date

# import phyloseq object --------
load("./derive_data/phy_rare/phyrare.RData")
load("./derive_data/step3_alpha_diver/adivrare_06012019.RData")

# get summary --------
summary <- lapply(adivrare, summary)
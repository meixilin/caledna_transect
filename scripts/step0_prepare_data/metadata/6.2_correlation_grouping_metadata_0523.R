# Title: grouping correlation test result and generate vargroup variable --------
# Author: Meixi Lin
# Date: Thu May 23 10:25:50 2019
# Author:
# Date:
# Modification:

# preparation --------
rm(list = ls())
cat("\014")

setwd("~/Lab/caledna_transect")
date() # the execution date
library(dplyr)
library(ggplot2)
library(ggcorrplot)

# load data --------
myvar <- c("Longitude", "Latitude", "hfp", 
           "bio1", "bio2", "bio3", "bio4", "bio5", "bio6",
           "bio7", "bio8", "bio9", "bio10", "bio11", "bio12",
           "bio13", "bio14", "bio15", "bio16", "bio17","bio18", "bio19",
           "phihox", "orcdrc", "cecsol", "sndppt", "clyppt", "bldfie", "ntot",
           "elev", "Slope", "aspect","CTI", "Roughness", "Ruggedness", "DAH",
           "B1", "B2", "B3", "B4", "B5", "B6", "B7",
           "B8", "B8A", "B9", "B10", "B11", "B12",
           "NDVIS2", "NDVI32", "EVI", "NBRT", "greenness",
           "imprv", "ptrcv")
mycorrcutoff <- 0.80
load(file = "./derive_data/step0_prepare_data/metadata_correlation.RData")

# use neighbor joining --------
# generated "vargroup" variable in Table S1
mydist <- as.dist(1 - abs(corr))
varclust2 <- hclust(d = mydist)
plot(varclust2)
dd <- 1- mycorrcutoff
clust <- as.data.frame(cutree(varclust2, h = dd)); colnames(clust) <- "vargroup"
table(clust)
# from each group, pick one variable
reducevar <- clust %>% 
    tibble::rownames_to_column(var = "Column.name") %>%
    dplyr::arrange(vargroup) 
write.csv(reducevar, file = paste0("./derive_data/step0_prepare_data/reduce_variable_corr_", mycorrcutoff, ".csv"))

# make some minor changes to pick a value more interpretable
# vargroup 29: NDVIS2 -- greenness: pick greenness
# vargroup 30: NDVI32 -- EVI: NDVI32

# generate reduced 33 environmental variables
myvar2 <- c("Longitude", "hfp", "bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio8", "bio14", "bio15", "phihox", "orcdrc", "cecsol", "sndppt", "bldfie", "ntot", "elev", "Slope", "aspect","CTI", "DAH", "B1", "B4", "B6", "B9", "B10", "B11", "NDVI32", "NBRT", "greenness", "imprv", "ptrcv")
reducevar[reducevar$Column.name %in% myvar2, 'vargroup']

# source the function_transect.R 
source("./scripts/function_transect.R")
myvar2 == contlist

# 
# # read and combine reduce var 
# metaexpl <- read.csv(file = "./final_data/metadata/Metadata_explanation_05012019.csv")
# metaexpl <- dplyr::left_join(metaexpl, reducevar, by = "Column.name")
# 
# write.csv(metaexpl, file = "./final_data/metadata/Metadata_explanation_05232019.csv", row.names = F)

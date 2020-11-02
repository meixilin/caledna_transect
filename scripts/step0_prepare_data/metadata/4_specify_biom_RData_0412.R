# Title: Import data in phyloseq and specify factor files --------
# Author: Meixi Lin
# Date: Sun Oct 14 17:43:00 2018
# Author:
# Date:
# Modification:

# preparation --------
rm(list = ls())
cat("\014")

setwd("/Users/linmeixi/UCLA/Lab/abiotic_transect")
library(dplyr)
library(reshape2)

source("./scripts/function_transect.R")

date() # the execution date


# import data and biom --------
biom <- read.csv("./final_data/metadata/Final_metadata_05012019.csv")

# read label definition for taxousda 
taxousda <- read.csv(file = "./final_data/metadata/USDA_levs.csv")
labs_taxousda <- taxousda %>% 
    dplyr::filter(Number %in% biom$taxousda)
write.csv(labs_taxousda, file = "./derive_data/step1_mk_phyloseq/taxousda_labs.csv")

# read label definition for NLCD 
nlcd <- read.csv(file = "./final_data/metadata/NLCD_levs.csv")
labs_nlcd <- nlcd %>% 
    dplyr::filter(Value %in% biom$NLCD)
write.csv(labs_nlcd, file = "./derive_data/step1_mk_phyloseq/NLCD_labs.csv")

# specify the biom structure 
biom <- biom %>% 
    dplyr::mutate(clust = as.factor(clust), 
                  Zeta_4 = as.factor(Zeta_4), 
                  taxousda = factor(x = taxousda, levels = labs_taxousda$Number, labels = labs_taxousda$Group), 
                  NLCD = factor(x = NLCD, levels = labs_nlcd$Value, labels = labs_nlcd$Definition)) 

save(biom, file = "./final_data/Final_metadata.RData")
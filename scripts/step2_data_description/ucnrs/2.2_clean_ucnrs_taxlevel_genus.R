# Title: generate matching taxonomy hierachy for UCNRS and CALeDNA --------
# Author: Meixi Lin
# Date: Thu Jan 16 11:31:13 2020
# Author:
# Date:
# Modification:

# match the taxonomy hierachy --------
rm(list = ls())
cat("\014")

setwd("~/Lab/caledna_transect/")
library(dplyr)
library(taxonomizr)
library(taxize)
library(taxizedb)
library(RSQLite)

date() # the execution date

# # make a sql database for query --------
# nodes.sql <- read.nodes.sql(nodeFile = "./raw_data/TAXO/nodes.dmp", sqlFile = paste0(outdir, "nodes.sqlite"))
# names.sql <- read.names.sql(nameFile = "./final_data/TAXO/names.dmp", sqlFile = paste0(outdir, "nodes.sqlite"))

# define pathes --------
indir <- "./derive_data/step2_data_description/ucnrs/"
outdir <- "./derive_data/step2_data_description/ucnrs/"
dir.create(outdir, recursive = T)
mysql <- paste0(outdir, "nodes.sqlite")

# load data -------
fauna <- read.csv(file = paste0(indir, "UC-NRS-Species-List_cleaned.csv"), stringsAsFactors = F)
flora <- read.csv(file = paste0(indir, "reserve_plant_list_cleaned.csv"), stringsAsFactors = F)

# test the fauna --------
fauna <- cbind(fauna, reshape2::colsplit(fauna$Species, " ", names = c("Genus", "suffix"))[, "Genus"]) %>%
    dplyr::mutate_if(is.factor, as.character)
colnames(fauna)[ncol(fauna)] <- "Genus"

# get the family name from NCBI and see if it matches 
family_id <- getId(fauna$Family, sqlFile = mysql)
genus_id <- getId(fauna$Genus, sqlFile = mysql)
taxa_id <- getId(fauna$Species, sqlFile = mysql)
fauna <- cbind(fauna, family_id, genus_id, taxa_id)

# test the flora -------
family_idp <- getId(flora$Family, sqlFile = mysql)
genus_idp <- getId(flora$Genus, sqlFile = mysql)
taxa_idp <- getId(flora$Species, sqlFile = mysql)
flora <- cbind(flora, family_idp, genus_idp, taxa_idp)

# write out the fauna with species id and try to fix it --------
write.csv(x = fauna, file = paste0(outdir, "fauna_ncbi_bind_v0.csv"))
write.csv(x = flora, file = paste0(outdir, "flora_ncbi_bind_v0.csv"))

#######################################################
# then fixed them manually, see README_ucnrs.md for details 
#######################################################

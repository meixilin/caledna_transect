# Title: compare ucnrs sites --------
# Author: Meixi Lin
# Date: Fri May 17 15:02:56 2019
# Author:
# Date:
# Modification:

# preparation --------
rm(list = ls())
cat("\014")

setwd("~/Lab/caledna_transect/")
library(dplyr)
library(phyloseq)
library(ggplot2)
date() # the execution date
source("./scripts/function_transect.R")

# define variables --------
myprimer <- "all"
myphylum <- c("Chordata", "Arthropoda", "Streptophyta")
mylev <- c("Order", "Family", "Genus")
mycutoff <- 0
outdir <- "./derive_data/step2_data_description/ucnrs/"
dir.create(outdir, recursive = T)

# define functions --------
strip_na_dup <- function(data) {
    data = as.character(data)
    id = which(is.na(data) | data == "")
    if (length(id) > 0) {
        data = data[-id]
    }
    data = sort(unique(data))
    return(data)
}

# get rid of human sequences in eDNA 
strip_human <- function(data) {
    id = which(data$Order == "Primates")
    if (length(id) > 0) {
        data = data[-id,]
    }
    return(data)
}

# load data -------
fauna <- read.csv(file = "./final_data/ucnrs/UCNRS_fauna_list.csv", stringsAsFactors = F)
flora <- read.csv(file = "./final_data/ucnrs/UCNRS_flora_list.csv", stringsAsFactors = F)
load(file = "./derive_data/phy_deco/phydeco_uc.RData")

psdata <- phydeco_uc[[myprimer]]
myloc <- sample_data(psdata)$loc %>% unique() %>% as.character()
# 
# sort(table(sample_data(psdata)$loc))
# # fort ord had most samples (26)

# add overall comparisons --------
leftotu <- names(taxa_sums(psdata))[(taxa_sums(psdata) > mycutoff)]
dtedna <- prune_taxa(leftotu, x = psdata)
dtedna_otu <- otu_table(dtedna)@.Data %>% as.data.frame() 
dtedna_tax <- tax_table(dtedna)@.Data %>% as.data.frame()


internames <- c("Phylum", "taxlevel", "intersect", "unique_edna", "unique_trad", "count_intersect", "count_unique_edna", "count_unique_trad")
inter <- data.frame(matrix(nrow = (length(mylev) * length(myphylum)), 
                           ncol = length(internames)))
colnames(inter) <- internames
inter$Phylum <- rep(myphylum, each = length(mylev))
inter$taxlevel <- rep(mylev, times = length(myphylum))

# start iteration --------
for (jj in myphylum) {
    print(jj)
    edna <- dtedna_tax %>% dplyr::filter(Phylum == jj) 
    edna <- strip_human(edna)
    
    if (jj != "Streptophyta") {
        trad <- fauna[fauna$Phylum == jj,] 
    } else {
        trad <- flora
    }
    
    for (kk in mylev) {
        # subset a finer scale
        ftrad <- strip_na_dup(trad[, kk])
        fedna <- strip_na_dup(edna[, kk])
        # strip the unknown and NA
        
        index <- (inter$Phylum == jj) & (inter$taxlevel == kk)
        if (length(intersect(ftrad, fedna)) == 0) {
            inter[index, 'intersect'] <- NA
            inter[index, 'count_intersect'] <- 0
        } else{
            chars <- paste(intersect(ftrad, fedna), collapse = ";")
            inter[index, 'intersect'] <- chars
            inter[index, 'count_intersect'] <- length(intersect(ftrad, fedna))
        }
        
        if (length(setdiff(ftrad, fedna)) == 0) {
            inter[index, 'unique_trad'] <- NA
            inter[index, 'count_unique_trad'] <- 0
        } else{
            chars <- paste(setdiff(ftrad, fedna), collapse = ";")
            inter[index, 'unique_trad'] <- chars
            inter[index, 'count_unique_trad'] <- length(setdiff(ftrad, fedna))
        }
        
        if (length(setdiff(fedna,ftrad)) == 0) {
            inter[index, 'unique_edna'] <- NA
            inter[index, 'count_unique_edna'] <- 0
        } else{
            chars <- paste(setdiff(fedna, ftrad), collapse = ";")
            inter[index, 'unique_edna'] <- chars
            inter[index, 'count_unique_edna'] <- length(setdiff(fedna, ftrad))
        }
    }
}

# percent intersect 
inter = inter %>%
    dplyr::mutate(percent_inter_trad = count_intersect/(count_intersect + count_unique_trad)) %>%
    dplyr::mutate(percent_inter_edna = count_intersect/(count_intersect + count_unique_edna))

# write the output --------
write.csv(inter, file = paste0(outdir, "overall_species_intersect.csv"))



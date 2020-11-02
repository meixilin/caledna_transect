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
mycutoff <- 0
mylev <- c("Order", "Family", "Genus")
names(mylev) <- myphylum

indir <- "./derive_data/step2_data_description/ucnrs/bysite/"
outdir <- "./derive_data/step2_data_description/ucnrs/"

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
rsnames <- read.csv(file = "./raw_data/ucnrs/reserve_name_binder.csv", stringsAsFactors = F)
myloc <- rsnames$loc

# create a dataframe to store the intersection info --------
internames <- c("loc", "Phylum", "taxlevel", "intersect", "unique_edna", "unique_trad", "count_intersect", "count_unique_edna", "count_unique_trad")
inter <- data.frame(matrix(nrow = (length(myloc) * length(mylev) * length(myphylum)), 
                           ncol = length(internames)))
colnames(inter) <- internames
inter$loc <- rep(myloc, each = length(mylev) * length(myphylum))
inter$Phylum <- rep(myphylum, each = length(mylev), times = length(myloc))
inter$taxlevel <- rep(mylev, times = length(myloc) * length(myphylum))

# loop through all sites --------
for (ii in myloc) {
    print(ii)
    # load data ########
    # here each taxonomy entry is unique 
    dtfauna <- read.csv(file = paste0(indir, ii, "_fauna.csv"), stringsAsFactors = F)
    dtflora <- read.csv(file = paste0(indir, ii, "_flora.csv"), stringsAsFactors = F)
    dtedna_tax <- read.csv(file = paste0(indir, ii, "_edna_tax.csv"), stringsAsFactors = F)
    for (jj in myphylum) {
        print(jj)
        edna <- dtedna_tax %>% dplyr::filter(Phylum == jj)
        edna <- strip_human(edna)
        if (jj != "Streptophyta") {
            trad <- dtfauna[dtfauna$Phylum == jj,] 
        } else {
            trad <- dtflora
        }
        
        for (kk in mylev) {
            # subset a finer scale
            ftrad <- strip_na_dup(trad[, kk])
            fedna <- strip_na_dup(edna[, kk])
            # strip the unknown and NA
            
            index <- (inter$loc == ii) & (inter$Phylum == jj) & (inter$taxlevel == kk)
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
}

# percent intersect 
inter = inter %>%
    dplyr::mutate(percent_inter_trad = count_intersect/(count_intersect + count_unique_trad)) %>%
    dplyr::mutate(percent_inter_edna = count_intersect/(count_intersect + count_unique_edna))

# write the output --------
write.csv(inter, file = paste0(outdir, "sites_species_intersect.csv"))

# end analysis --------
closeAllConnections()


# Title: generate matching taxonomy hierachy for UCNRS and CALeDNA --------
# Author: Meixi Lin
# Date: Thu Jan 16 11:31:13 2020
# Author:
# Date:
# Modification:

#######################################################
# genus had been fixed  
#######################################################

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

# define variables --------
# define pathes 
indir <- "./derive_data/step2_data_description/ucnrs/"
outdir <- "./final_data/ucnrs/"
dir.create(outdir, recursive = T)
mysql <- paste0(indir, "nodes.sqlite")

# define functions --------
# load the binded data
load_bind_v1 <- function(indir, filename){
    data <- read.csv(file = paste0(indir, filename), row.names = 1, stringsAsFactors = F) %>% 
        dplyr::select(loc, CRUX_Genus) %>% 
        dplyr::distinct()
    # fill in the genus and genus id
    genus <- data$CRUX_Genus
    genusid <- getId(genus, sqlFile = mysql)
    data <- cbind(data, genusid) %>% 
        dplyr::filter(!is.na(genusid)) %>% 
        dplyr::mutate_if(is.factor, as.character) %>% 
        dplyr::arrange(loc, CRUX_Genus)
    # remove the data that still had NA in genus level id 
    colnames(data)[2:3] = c("Genus", "Genus_ID")
    return(data)
}

# define a function to fix multiple taxon id problem 
fix_dup_genus <- function(data, ref = 2, refcan) {
    # get duplicate id indix
    dupid = grep(",", data[, "Genus_ID"])
    # define desired taxonomy level
    mytaxa = c("superkingdom", "phylum", "class", "order", "family")
    # define a new data frame 
    newdata = data.frame(matrix(nrow = nrow(data), ncol = length(mytaxa)))
    for (ii in 1:nrow(data)) {
        if (ii %in% dupid) {
            # check the right genus entry
            ids = reshape2::colsplit(data[ii,"Genus_ID"], pattern = ",", names = c("Genus_id1", "Genus_id2"))
            test1 <- getTaxonomy(ids = ids$Genus_id1, sqlFile = mysql, desiredTaxa = mytaxa)
            test2 <- getTaxonomy(ids = ids$Genus_id2, sqlFile = mysql, desiredTaxa = mytaxa)
            testlist = list(test1, test2)
            if (all(test1 == test2, na.rm = T)) {
                test = test1 
                rightid = 1 # take the first id 
            } else {
                list = c(test1[ref], test2[ref])
                rightid = which(list %in% refcan)
                test = testlist[[rightid]]
            }
            # fix the right genus id 
            data[ii, "Genus_ID"] = as.character(ids[rightid])
            newdata[ii, ] = test
        } else {
            newdata[ii, ] = getTaxonomy(ids = data[ii,"Genus_ID"], sqlFile = mysql, desiredTaxa = mytaxa)
        }
    }
    
    colnames(newdata) = c("Superkingdom", "Phylum", "Class", "Order", "Family")
    # bind the data 
    data = cbind(data, newdata)

    # change the heading and the order 
    neworder = c("loc", "Genus_ID", "Superkingdom", "Phylum", "Class", "Order", "Family", "Genus")
    data = data[,neworder] 
    return(data)
}


# load cleaned data --------
fauna <- load_bind_v1(indir, "fauna_ncbi_bind_v1.csv")
flora <- load_bind_v1(indir, "flora_ncbi_bind_v1.csv")

# main --------
date()
fauna_fix <- fix_dup_genus(data = fauna, refcan = c("Chordata", "Arthropoda"))
flora_fix <- fix_dup_genus(data = flora, refcan = c("Streptophyta", "Chlorophyta"))
date()

write.csv(fauna_fix, file = paste0(outdir, "UCNRS_fauna_list.csv"))
write.csv(flora_fix, file = paste0(outdir, "UCNRS_flora_list.csv"))

#######################################################
# fixed one line manually since it had three ids 
# 3063	Sweeney_Granite_Mountains_Desert_Research_Center	56373	Eukaryota	Arthropoda	Insecta	Lepidoptera	Nolidae	Baileya
# 3063	Sweeney_Granite_Mountains_Desert_Research_Center	128736	Eukaryota	Streptophyta	NA	Asterales	Asteraceae	Baileya
#######################################################



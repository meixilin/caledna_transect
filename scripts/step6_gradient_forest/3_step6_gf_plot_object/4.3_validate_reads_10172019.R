# Title: dig into reads tested positive with gf --------
# Author: Meixi Lin
# Date: Thu Oct 17 10:36:44 2019
# Author:
# Date:
# Modification:

# preparation --------
rm(list = ls())
cat("\014")

setwd("~/UCLA/Lab/abiotic_transect/")
# setwd("/u/home/m/meixilin/project-rwayne/abiotic_transect/")
# load packages
library(dplyr)
library(phyloseq)
library(stringr)
library(gradientForest)

# define function --------
# take in asv as a confidence level
get_conf <- function(asv, mylevel = "family") {
    # match for minimum confidence level for family 
    taxlevels <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
    myconf <- asv$taxonomy_confidence %>% 
        reshape2::colsplit(., ";", taxlevels) 
    myconf <- myconf[,mylevel] %>% 
        reshape2::colsplit(., ":", c("txlevel", "conf")) %>% 
        dplyr::tbl_df() %>%
        dplyr::mutate(conf = as.numeric(conf)) 
    return(myconf)
}

# get read count, asv is an asv that had freq table after "accessions" 
count_reads <- function (asv) {
    hh <- asv %>%
        dplyr::select(-(1:accessions))
    mycount <- rowSums(hh)
    return(mycount)
}

query_name <- function(xx, query, conf, mylevel = "family") {
    # first match for the string 
    matchid <- str_which(xx$asv$taxonomy, pattern = query)
    if (length(matchid) == 0) {
        print(paste("no match for", query, "in db", xx$primer, xx$transect))
        myhit <- NULL
    } else {
        # output the sequences 
        myhit <- xx$asv[matchid, ] %>% 
            dplyr::tbl_df(.) 
        # get confidence level
        myconf <- get_conf(myhit)
        confid <- myconf$conf >= conf
        if (all(confid == FALSE)) {
            print(paste(length(matchid), "matches for", query, "in db", xx$primer, xx$transect,".", "None passed confidence cutoff at", conf))
            myhit <- NULL
        } else {
            myhit <- myhit[confid,] %>% 
                dplyr::mutate(primer = xx$primer,
                              transect = xx$transect) %>% 
                dplyr::select(primer, transect, everything())
            # add myhit summary reads 
            nreads <- count_reads(myhit)
            myhit <- cbind(myhit, nreads) %>% 
                dplyr::select(primer, transect, nreads, everything()) %>% 
                dplyr::arrange(., desc(nreads))
        }
    }
    return(myhit)
}

# main --------
confcut <- 60
# load gf, get families ========
load("./derive_data/step4_gradient_forest/2_final/gf_deco_all_Family_Presence_2000_2_0.05_FALSE_17.RData")
topfamilies <- sort(gf$result, decreasing = T) %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var = "family")
colnames(topfamilies)[2] <- "imp"
write.csv(topfamilies, file ="./derive_data/step4_gradient_forest/2_final/gf_deco_all_Family_Presence_2000_2_0.05_FALSE_17_result.csv")
rm(gf)

interest = topfamilies[1:10, 'family']

# load gf Y ==========
gf_Y <- read.csv(file = "./derive_data/step4_gradient_forest/2_final/gf_deco_all_Family_Presence_2000_2_0.05_FALSE_17_Y.csv")
View(gf_Y)
gf.in <- gf_Y[,interest]


# load the taxonomy detailed reads ========
load("./final_data/transect_ASV/taxonomy_detail_all.RData")
# load the gf of interest

# define output directory 
outdir <- "./derive_data/step4_gradient_forest/4_validate_reads/"
dir.create(outdir, recursive = T)

## starts check the taxonomy assignments ========
for (ii in interest) {
    query = ii
    myhit <- lapply(asvdb, function(xx) {
        xxhit <- query_name(xx, query, conf = confcut)
        return(xxhit)
    })
    
    # remove null elements 
    myhit <- myhit[-which(sapply(myhit, is.null))]
    
    # write to a csv file 
    lapply(myhit, function(xx) {
        write.table(x = xx, sep = "," , file = paste0(outdir, query, "_taxonomy_detailed_conf_", confcut, ".csv"), append = T)
        return(0)
    })
}



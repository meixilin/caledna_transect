# Title: basic statistics, e.g. read depth for each species, taxonomic distribution, etc. --------
# For rarefied dataset 
# Author: Meixi Lin
# Date: Thu Aug 30 17:51:56 2018
# Author: 
# Date: Wed Mar 18 13:25:04 2020
# Modification: Clean up

# preparation --------
rm(list = ls())
cat("\014")

setwd("~/Lab/caledna_transect")
library(dplyr)
library(ranacapa)
library(phyloseq)
library(reshape2)

source("./scripts/function_transect.R")

date() # the execution date

# define variables --------
txlevel <- "Phylum"
outdir <- "./derive_data/step2_data_description/summary_stats/"
indir = outdir 
# define functions --------
# This function returns the number of LCA entries within the selected category 
# i.e. how many different taxonomy entries (species) were accounted as a child of "Proteobacteria" 
get_taxa_entry_count <- function(phyloseq1, taxlevel = c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
    txcount <- tax_table(phyloseq1)@.Data %>%
        as.data.frame %>%
        # select by taxlevel
        group_by(eval(as.name(taxlevel))) %>% 
        # get count for each value
        summarise(count = n())
    colnames(txcount) <- c(taxlevel, "count")
    return(txcount)
}

# This function returns the number of seqeuencing reads in each phyloseq object  
# i.e. how many reads were accounted as ";Chrysophyceae;Chromulinales;Chromulinaceae;NA;" 
get_taxa_read_count <- function(phyloseq1) {
    # get taxonomy reads
    txcount <- taxa_sums(phyloseq1) %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "taxa")
    colnames(txcount)[ncol(txcount)] = "count"
    return(txcount)
}

join_dtlist <- function(dtlist, length, by) {
    mydt = dtlist[[1]]
    for (ii in 2:length) {
        mydt = dplyr::left_join(mydt, dtlist[[ii]], by = by, suffix = as.character(c(ii-1, ii)))
    }
    return(mydt)
}

# import data --------
load(file = "./derive_data/phy_rare/phyrare.RData")
load(file = "")
# use the same order as the decontaminated datasets
alltaxa1 = read.csv(file = paste0(indir, "phy_deco_taxa_entry_count_Phylum.csv"), row.names = 1, stringsAsFactors = F) %>%
    select(Phylum)
alltaxa2 = read.csv(file = paste0(indir, "phy_deco_taxa_read_sums.csv"), row.names = 1, stringsAsFactors = F) %>%
    select(-starts_with("primer"))

# 1. taxonomy entries by phylum by primers ---------
# get taxonomy counts in the sense of the samples
txcount <- lapply(phyrare, get_taxa_entry_count, taxlevel = txlevel)
names(txcount) <- primers

# merge the different taxonomy by "alltaxa1", use the original names 
taxcount <- lapply(txcount, function(x) {
    hh <- merge(alltaxa1, x, by = txlevel, all.x = T, sort = F)
    return(hh)
})

names(taxcount) <- primers

# get a more beautifully shaped dataframe
taxcounts <- cbind(taxcount[[1]],taxcount[[2]][,2], taxcount[[3]][,2], taxcount[[4]][,2], taxcount[[5]][,2], taxcount[[6]][,2])
colnames(taxcounts)[-1] <- primers
# replace na
taxcounts <- taxcounts %>% 
    mutate(Phylum = ifelse(is.na(Phylum), "NA", Phylum)) 
taxcounts[is.na(taxcounts)] <- 0

write.csv(taxcounts, file = paste0(outdir, "phy_rare_taxa_entry_count_", txlevel, ".csv"))

# 2. read depth each taxonomy each primer --------
# get taxonomy read count in each metabarcode 
taxasums <- lapply(phyrare, get_taxa_read_count)

# combine the taxonomy table with the read count in each metabarcode 
tsdf <- join_dtlist(dtlist = c(list(alltaxa2), taxasums), length = 7, by = "taxa") 
colnames(tsdf)[8:13] <- paste0("primer_", primers)
tsdf <- tsdf %>%
    mutate_at(vars(starts_with("primer")), funs(ifelse(is.na(.),0,.)))

write.csv(tsdf, file = paste0(outdir, "phy_rare_taxa_read_sums.csv"))

# 3. after rarefaction, read depth each sample each primer ---------
samplenames <- as.data.frame(sample_names(phyrare[['all']])) 
samplenames <- as.data.frame()
colnames(samplenames) <- "sample"
samplesums <- lapply(phyrare, function(xx){
    hh <- as.data.frame(sample_sums(xx)) %>% 
        tibble::rownames_to_column(var = "sample") %>%
        mutate(sample = as.character(sample))
    }) 
ssdf <- lapply(samplesums, function(xx) {
    hh <- left_join(samplenames, xx, by = "sample", all.x = T)
}) 
ssdf1 <- join_dtlist(ssdf, length = 6, by = "sample")
colnames(ssdf1)[-1] <- primers
ssdf1[is.na(ssdf1)] <- 0

write.csv(ssdf1, file = paste0(outdir, "phy_rare_sample_read_sums.csv"))

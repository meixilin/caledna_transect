# Title: generate summary PCoA plots for duplicated sites --------
# Author: Meixi Lin
# Date: Mon Jan  6 00:18:24 2020

# preparation --------
rm(list = ls())
cat("\014")

setwd("~/Lab/caledna_transect/")
# load packages
library(dplyr)
library(ggplot2)
library(reshape2)
library(phyloseq)
source("./scripts/function_transect.R")
source("./scripts/step6_gradient_forest/functions/function_gfprep.R")

# def functions --------
add_Replicate <- function(myphyseq, dup) {
    Replicate = sample_data(myphyseq)$indup
    for (ii in 1:length(Replicate)) {
        if (Replicate[ii] %in% dup) {
            # keep it 
        } else {
            Replicate[ii] <- NA
        }
    }
    sample_data(myphyseq)$Replicate = Replicate
    return(myphyseq)
}
# def variables --------
plotdir <- "./plots/step2_data_description/concordance/"
dir.create(plotdir, recursive = T)
# abundance cutoff 
cutoff = 2
txlevel = "Family"
# duplicate sites 
dup <- c("K0058C2", "K0117C2", "K0142A2", "K0142C2", "K0177A1")

# load data --------
load(file = "./derive_data/step2_data_description/concordance/phydeco_wdup.RData")

myphyseq = phylist[['all']]
rm(phylist)

# add duplication marks --------
myphyseq = add_Replicate(myphyseq = myphyseq, dup = dup)

# process data in the same way as the gradient forest process --------
# add taxonomy summary 
myphyseq = glom_tax(myphyseq, taxlevel = txlevel)
# add cutoff 
myphyseq <- transform_sample_counts(myphyseq, function(abund) {1*(abund > cutoff)})
# delete the taxonomy entries that didn't occur 
# NOTE: here, we deleted "absence" to reduce the denominator to taxonomy entries that occurred in entire dataset 
myphyseq <- prune_taxa(taxa_sums(myphyseq) > 0, myphyseq)

# plot the ordination from the processed phyloseq -------
myordination = ordinate(myphyseq, method = "MDS", distance = "jaccard", binary = T)
pp1 <- plot_ordination(myphyseq, myordination, color = "Replicate") +
    geom_point(size = 1.5) +
    theme_bw() +
    theme(text = element_text(size = 16)) 
ggsave(filename = paste0(plotdir, "ord_dup_jaccard_MDS_supp_fig6.pdf"), plot = pp1, height = 6, width = 7)

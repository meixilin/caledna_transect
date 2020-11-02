# Title: iNext evaluation of the rarefied items
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Thu Mar  5 13:57:38 2020

# preparation --------
rm(list = ls())
cat("\014")

setwd("~/Lab/caledna_transect/")

outdir = "./derive_data/step1_create_phyloseq/eval_rarefaction/"
sink(file = paste0(outdir, "iNEXT_report.log"))

library(iNEXT)
library(ggplot2)
library(ranacapa)
library(phyloseq)
date()


# def functions --------
get_iNEXT <- function(physeq, cutoff = 0, hill = 0) {
    # convert to presence
    physeq = phyloseq::transform_sample_counts(physeq, function(x) {1*(x > cutoff)})
    physeq = phyloseq::prune_taxa(taxa_sums(physeq) > 0, physeq)
    otu = t(ranacapa::vegan_otu(physeq))
    otu = as.incfreq(otu)
    # perform iNEXT 
    out.inc <- iNEXT(otu, q = hill, datatype="incidence_freq")
    return(out.inc)
}

get_SC <- function(iNEXT) {
    xx = iNEXT$iNextEst
    xx = xx[xx$method == "observed", 'SC']
    return(xx)
}

# def variables --------

# load data --------
load(file = "./derive_data/phy_deco/phydeco.RData")
load(file = "./derive_data/phy_rare/phyrare.RData")

# main --------
iNEXT_deco <- lapply(phydeco, get_iNEXT)
iNEXT_rare <- lapply(phyrare, get_iNEXT)


# output --------
save(iNEXT_deco, file = paste0(outdir, "iNEXT_deco.RData"))
save(iNEXT_rare, file = paste0(outdir, "iNEXT_rare.RData"))

# get sample coverage --------
lapply(iNEXT_deco, get_SC)
lapply(iNEXT_rare, get_SC)

# get all info -------
print(iNEXT_deco)
print(iNEXT_rare)


# ending --------
sink()
closeAllConnections()



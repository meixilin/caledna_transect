# Title: Basic beta diversity: LCBD --------
# Author: Meixi Lin
# Date: Tue Jun 26 15:12:01 2018
# Author: Meixi Lin
# Date: Fri Jul 27 20:36:37 2018
# Modification:
# Date: Fri Aug 31 23:31:26 2018

# preparation --------
rm(list = ls())
cat("\014")

setwd("/Users/linmeixi/UCLA/Lab/abiotic_transect")
library(tibble)
library(dplyr)
library(ggplot2)
library(phyloseq)
library(microbiomeSeq)
library(reshape2)
library(vegan) # microbiomeSeq calls a function "ordiellipse"
library(adespatial) # microbiomeSeq calls a function "beta.div"
library(grid) # microbiomeSeq calls a function "textGrob"
library(gridExtra) # microbiomeSeq calls a function "ttheme_minimal"
library(fso)

source("./r_codes/function_beta_0725.R")
txlevel <- "Phylum"
bothokflag <- F
mycutoff <- 4

date() # the execution date

# load --------
load(file = "./derive_data/phyrare.RData")

# import phyloseq object --------
for (ii in 1:length(primers)) {
    
primer <- primers[ii]
physeq <- phyrare[[ii]]

if (primer %in% c("FITS", "PITS")) {
    txlevel <- "Class"
} else {
    txlevel <- "Phylum"
}

# beta diversity: local contribution to beta diversity --------
if (bothokflag == T) {
    physeq1 <- subset_samples(physeq, bothok == bothokflag)
} else {
    physeq1 <- physeq
}
sampledf <- data.frame(sample_data(physeq1))
phyglom <- glom_tax(phyloseq1 = physeq1, taxlevel = txlevel)
mpg0 <- normalise_data(t(phyglom), norm.method = "relative") # calculate relative density 

# (by major habitat) =================================
for (jj in 1:length(catlist)) {
groupv = catlist[jj]
category <- get_variable(physeq1, varName = groupv)
sumby <- table(category) %>% as.data.frame() %>% 
    filter(Freq > mycutoff)
id <- as.character(sampledf[(sampledf[,groupv] %in% sumby$category), 
                            'MatchName'])
mpg <- subset_samples(mpg0, MatchName %in% id)
pp <- plot_taxa(mpg, grouping_column = groupv, method = "hellinger", number.taxa = 10) +
    labs(title = paste("grouping =", groupv, "primer =", primer, "bothok =", bothokflag))

plotdir <- paste0("./plots_important/step3_beta_diver/", primer)
dir.create(plotdir, recursive = T)
ggsave(paste0("beta_diver_lcbd_", groupv, "_hellinger.pdf"), plot = pp, path = plotdir, height = 6, width = 20)
}
}

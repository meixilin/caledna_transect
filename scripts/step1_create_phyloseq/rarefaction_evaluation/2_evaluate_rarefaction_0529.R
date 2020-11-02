# Title: evaluate rarefaction depth --------
# Author: Meixi Lin
# Date: Wed May 29 10:58:17 2019
# Author:
# Date:
# Modification:

# preparation --------
rm(list = ls())
cat("\014")

setwd("~/Lab/caledna_transect/")
library(tibble)
library(dplyr)
library(ggplot2)
library(ranacapa)
library(phyloseq)
library(microbiomeSeq)
source("./scripts/function_transect.R")

date() # the execution date

# define the variables --------
indir = "./derive_data/step1_create_phyloseq/eval_rarefaction/"
outdir = indir
plotdir = "./plots/step1_create_phyloseq/eval_rarefaction/"

# load the object --------
raredata <- lapply(primers_commeco, function(primer) {
    xx <- get(load(file = paste0(indir, primer, "_rarefaction_summary_verbose.RData"))) %>% 
        dplyr::select(Sample, Depth, Measure, Alpha_diversity_mean, Alpha_diversity_sd) %>% 
        dplyr::filter(Measure == "Observed") %>% 
        dplyr::count(Depth) %>% 
        dplyr::mutate(Metabarcode = primer) %>%
        dplyr::mutate(Chosen = (Depth == rare_depth[primer])) %>%
        dplyr::mutate(Depth = as.factor(Depth))
    return(xx)
})

forplot <- rbind(raredata[[1]], raredata[[2]], raredata[[3]], raredata[[4]], raredata[[5]])

pp <- ggplot(forplot, aes(x = Depth, y = n, fill = Metabarcode, color = Chosen)) + 
    scale_color_manual(values = c("white", "black")) + 
    labs(x = "Rarefaction Depth",
         y = "Number of samples kept") +
    geom_bar(stat = "identity", size = 1) + 
    facet_wrap(. ~ Metabarcode, nrow = 3, ncol = 2) + 
    geom_hline(yintercept = 200) +
    geom_hline(yintercept = 250) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90))

# output --------
ggsave(filename = paste0(plotdir, "rarefaction_sample_count.pdf"), plot = pp, device = "pdf", width = 8, height = 8)
write.csv(forplot, file = paste0(outdir, "rarefaction_sample_count.csv"), row.names = F)
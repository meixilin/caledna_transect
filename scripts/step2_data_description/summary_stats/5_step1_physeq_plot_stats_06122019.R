# Title: summary stats plotting --------
# Author: Meixi Lin
# Date: Wed Jun 12 15:06:31 2019
# Author:
# Date:
# Modification:

# preparation --------
rm(list = ls())
cat("\014")

setwd("~/Lab/caledna_transect")
library(dplyr)
library(ggplot2)
library(ranacapa)
library(phyloseq)
library(reshape2)
library(ggsci)

source("./scripts/function_transect.R")

date() # the execution date

# define variables --------
indir = "./derive_data/step2_data_description/summary_stats/"
plotdir = "./plots/step2_data_description/summary_stats/"

# load the data ------
mydata <- read.csv(file = paste0(indir, "phy_deco_taxa_read_sums.csv"), row.names = 1)

# reformat for phylum --------
byphylum <- mydata %>%
    dplyr::select(Phylum, starts_with("primer_")) %>%
    dplyr::mutate(Phylum = as.character(Phylum)) %>%
    dplyr::mutate(Phylum = ifelse(is.na(Phylum) | Phylum == "", "unknown", Phylum)) %>%
    dplyr::group_by(Phylum) %>%
    dplyr::summarise(count_16S = sum(primer_16S), 
                     count_18S = sum(primer_18S),
                     count_CO1 = sum(primer_CO1),
                     count_FITS = sum(primer_FITS),
                     count_PITS = sum(primer_PITS),
                     count_all = sum(primer_all))

# forplot --------
txlevel <- "Phylum"
pplist <- vector("list", length = length(primers))
for (ii in primers) {
    forplot <- byphylum %>%
        dplyr::select(Phylum, contains(ii))
    colnames(forplot) <- c("variable", "value")
    maxid <- min(sum(forplot$value > 0), 9)
    forplot1 <- dplyr::top_n(forplot, n = maxid, wt = value) %>%
        dplyr::mutate(variable = as.character(variable)) %>%
        dplyr::arrange(desc(value)) 
    forplot2 <- dplyr::top_n(forplot, n = (-nrow(forplot) + maxid), wt = value)
    forplot2 <- c("Other", sum(forplot2[,'value']))
    forplot <- rbind(forplot1, forplot2) %>%
        dplyr::mutate(variable = factor(variable, levels = unique(variable)), 
                      primer = ii,
                      value = as.numeric(value))
    pplist[[which(primers %in% ii)]] <- forplot

}

# for the entire plot 
all_top10_sum <- rbind(pplist[[1]],pplist[[2]], pplist[[3]], pplist[[4]], pplist[[5]], pplist[[6]])

# plot entire thing 
mycategories <- names(sort(table(all_top10_sum$variable), decreasing = T))
all_top10_sum$variable <- factor(all_top10_sum$variable, levels = mycategories)
mycolours <- unname(pals::polychrome(n = 24))

forplot2 <- all_top10_sum %>%
    dplyr::filter(primer != "all")

pp2 <- ggplot(forplot2, aes(x = primer, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = "stack") + 
    scale_fill_manual(values = mycolours, 
                      name = txlevel) +
    labs(x = "Metabarcode", y = "Read count") +
    theme_bw() +
    theme(text = element_text(size = 16), 
          axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave(filename = paste0(plotdir, "taxa_sum_top10_combine_phylum.pdf"), width = 8, height = 6, plot = pp2)

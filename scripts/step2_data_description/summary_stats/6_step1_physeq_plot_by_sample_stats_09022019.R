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

# load the data ------
mydata <- read.csv(file = "./derive_data/step2_data_description/summary_stats/phy_deco_sample_read_sums.csv", row.names = 1) %>%
    dplyr::select(-all) 
# get forplot
forplot <- mydata %>%
    reshape2::melt(id.vars = "sample", value.name = "Reads")
# get labels 
xlabs <- mydata$sample[c((0:floor(nrow(mydata)/5)) *5 + 1, nrow(mydata))]

# get median 
mymedian <- median(forplot$Reads) # 7717

# plot the result -------
pp <- ggplot(forplot, aes(x = sample, y = Reads, fill = variable)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_x_discrete(breaks = xlabs) +
    scale_fill_discrete(name="Metabarcode") +
    labs(x = "Sample Name", y = "Read counts") +
    theme_bw() +
    theme(text = element_text(size = 20), 
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.title = element_text(size = 22)) 
    

ggsave(filename = paste0("./plots/step2_data_description/summary_stats/read_sum_sample_by_primer.pdf"), width = 20, height = 8, plot = pp)

# 


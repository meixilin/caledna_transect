# Title: PCoA plots for each category and primers --------
# combine the PCoA plots into one plot 
# Author: Meixi Lin
# Date: Mon Aug  5 19:37:31 2019
# Author:
# Date:
# Modification:

# preparation --------
rm(list = ls())
cat("\014")

setwd("~/UCLA/Lab/abiotic_transect/")
library(dplyr)
library(vegan)
library(grid)
library(gridExtra)
library(phyloseq)
library(ggplot2)
date() # the execution date

source("./scripts/function_transect.R")

# set some variable --------
mycutoff <- 4 # min items in one category 
# define plotting theme
mytheme <- theme_bw() +
    theme(text = element_text(size = 16))

# load data -------
load(file = "./derive_data/step3_beta_diver/pcoa_nocoast_category_plot.RData")

# for each variables --------
for (jj in 1:length(catlist)) {
    catid <- jj + length(catlist)*(0:(length(primers) - 1))
    groupv <- catlist[jj]
    names(pplist)[catid]
    mydata1 <- pplist[catid]
    
    mydata1 <- lapply(mydata1, function(xx) {
        xx <- xx + 
            theme(legend.position = "none")
    })
    
    allpp <- grid.arrange(mydata1[[1]], mydata1[[2]], mydata1[[3]],
                          mydata1[[4]], mydata1[[5]], mydata1[[6]],
                          nrow = 2, ncol = 3)
    plotdir <- "./plots_important/step3_beta_diver/combined/pcoa_nocoast/"
    dir.create(plotdir, recursive = T)
    ggsave(paste0("beta_diver_pcoa_", catlist[jj], ".pdf"), plot = allpp, path = plotdir, height = 12, width = 16)
}


# Title: plot permutation result --------
# Author: Meixi Lin
# Date: Mon Sep 16 11:05:05 2019
# Author:
# Date:
# Modification:

# preparation --------
options(echo = TRUE)
rm(list = ls())
cat("\014")
sink(file = paste0("./derive/3.2_plot_permutation_result_", Sys.Date(), ".log"))
setwd("~/UCLA/Lab/Coastal/coastal_gf/")
# setwd("/u/home/m/meixilin/project-rwayne/abiotic_transect/")
# load packages
library(ggplot2)
library(dplyr)
library(phyloseq)
library(gradientForest)
functions.file <- list.files("./scripts/functions/")
lapply(functions.file, function(xx) {source(file = paste0("./scripts/functions/", xx))})
sessionInfo()
date() # the execution date

# list the target gf objects --------
gfnames <- c("all_nf", "uni_nf", "multi_nf", "all_fl", "uni_fl", "multi_fl")

# load the dataset -------
for (ii in gfnames) {
    randaverage <- read.delim(file = paste0("./derive/3_permutate_gf/gfout_", ii, "_randaverage.tsv"), header = F, col.names = "randaverage")
    randtotal <- read.delim(file = paste0("./derive/3_permutate_gf/gfout_", ii, "_randtotal.tsv"), header = F, col.names = "randtotal")
    gf <- loadRData(fileName = paste0("./derive/gfout_", ii, ".RData"))
    gftotal <- gf_res(gf)[1]
    gfaverage <- gf_res(gf)[2]
    # get pseudo-P value 
    ptotal <- sum(randtotal > gftotal)/nrow(randtotal)
    paverage <- sum(randaverage > gfaverage)/nrow(randaverage)
    # print output and permutation --------
    print(ii)
    print(gf_res(gf))
    print(gf)
    print(dim(gf$Y))
    print(importance(gf)[1:5], digits = 4)
    print(paste("P-total:", ptotal))
    print(paste("P-Mean R2:", paverage))
    
    # start plotting -------
    pp1 <- ggplot(data = randaverage, aes(x = randaverage)) + 
        geom_histogram(color = 'black', fill = 'grey') +
        xlim(min(randaverage) * 0.95, gfaverage * 1.05) +
        geom_vline(xintercept = gfaverage, color = 'red') + 
        labs(x = "Mean R^2", y = "Frequency") + 
        theme_bw() +
        theme(text = element_text(size = 16)
              ) 
    
    pp2 <- ggplot(data = randtotal, aes(x = randtotal)) + 
        geom_histogram(color = 'black', fill = 'grey') +
        xlim(min(randtotal) * 0.95, gftotal * 1.05) +
        geom_vline(xintercept = gftotal, color = 'red') + 
        labs(x = "Number of response with positive R2", y = "Frequency") + 
        theme_bw() +
        theme(text = element_text(size = 16)
        ) 
    
    pp3 <- gridExtra::grid.arrange(pp1, pp2, ncol = 2)
    ggsave(filename = paste0("./plots/3_permutation_gfout_", ii, ".pdf"), plot = pp3, width = 8, height = 4)
}
    

# ending --------
date()
sink(NULL)
closeAllConnections()
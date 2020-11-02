# Title: Plot permutation output 
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Sat Apr 18 23:04:34 2020

# preparation --------
options(echo = TRUE)
rm(list = ls())
cat("\014")

setwd("~/Lab/caledna_transect/")
# load packages
library(ggplot2)
library(gradientForest)

source(file = "./scripts/step6_gradient_forest/functions/function_gfprep.R")

# load the dataset -------
randaverage <- read.csv(file = "./derive_data/step6_gradient_forest/gf_deco_all_Family_Presence_2000_2_0.05_FALSE_17_randaverage.tsv", header = F, col.names = "randaverage")
randtotal <- read.csv(file = "./derive_data/step6_gradient_forest/gf_deco_all_Family_Presence_2000_2_0.05_FALSE_17_randtotal.tsv", header = F, col.names = "randtotal")

load(file = "./derive_data/step6_gradient_forest/gf_deco_all_Family_Presence_2000_2_0.05_FALSE_17.RData")
gftotal <- gf_res(gf)[1]
gfaverage <- gf_res(gf)[2]

# start plotting -------
xlab.text1 = expression(paste("Mean R"^"2"))
pp1 <- ggplot(data = randaverage, aes(x = randaverage)) + 
    geom_histogram(color = 'black', fill = 'grey') +
    xlim(min(randaverage) * 0.95, gfaverage * 1.05) +
    geom_vline(xintercept = gfaverage, color = 'red') + 
    labs(x = xlab.text1, y = "Frequency") + 
    theme_bw() +
    theme(text = element_text(size = 16)) 

xlab.text2 = expression(paste("Number of families with positive R"^"2"))
pp2 <- ggplot(data = randtotal, aes(x = randtotal)) + 
    geom_histogram(color = 'black', fill = 'grey') +
    xlim(min(randtotal) * 0.95, gftotal * 1.05) +
    geom_vline(xintercept = gftotal, color = 'red') + 
    labs(x = xlab.text2, y = "Frequency") + 
    theme_bw() +
    theme(text = element_text(size = 16)) 

pp3 <- gridExtra::grid.arrange(pp1, pp2, ncol = 2)
ggsave(filename = "./plots/step6_gradient_forest/for_publication/gf_deco_all_permutation.pdf", plot = pp3, width = 8, height = 4, device = "pdf")
    





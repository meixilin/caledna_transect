# Title: evaluate species response --------
# Author: Meixi Lin
# Date: Wed Oct 16 14:36:25 2019
# Author:
# Date:
# Modification:

# preparation --------
options(echo = TRUE)
rm(list = ls())
cat("\014")

setwd("~/Lab/caledna_transect/")
# setwd("/u/home/m/meixilin/project-rwayne/abiotic_transect/")
# load packages
library(dplyr)
library(phyloseq)
library(gradientForest)
library(ggplot2) # 
library(ggmap)
library(rgdal)
library(ggsn)
library(gridExtra)
source("./scripts/function_transect.R")
functions.file <- list.files("./r_codes/step4_gradient_forest/functions/")
lapply(functions.file, function(xx) {
    source(file = paste0("./r_codes/step4_gradient_forest/functions/", xx))
    return(0)
    })

# load data --------
# gf <- loadRData(as.character(args[1])) # RENAME WHAT EVER GF to "gf"
load(file = "./derive_data/step4_gradient_forest/2_final/gf_deco_all_Family_Presence_2000_2_0.05_FALSE_17.RData")

# start counting -------
otu <- gf$Y %>%
    as.data.frame() %>%
    dplyr::mutate_if(is.factor, as.numeric)
otu <- otu - 1
countotu <- colSums(otu) %>%
    named_num2df(., c("family", "count"))

taxlist <- sort(gf$result, decreasing = T) %>% 
    named_num2df(., c("family", "r2")) %>%
    dplyr::left_join(. ,countotu, by = 'family')

pp <- ggplot(taxlist, aes(x = count, y = r2, color = count)) +
    geom_point() +
    geom_smooth() +
    geom_vline(xintercept = floor(nrow(otu) * 0.025), colour = 'red') + 
    theme(legend.position = "none") +
    theme_bw() +
    labs(title = "families' predictability change with positive count")

ggsave(filename = "./plots_important/step4_gradient_forest/imp_species/gf_deco_all_Family_Presence_2000_2_0.05_FALSE_17_r2_count.pdf", width = 8, height = 6)

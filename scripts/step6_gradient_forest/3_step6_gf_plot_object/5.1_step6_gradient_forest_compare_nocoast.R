# Title: Compare the coastal sites
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Wed Dec 23 12:31:01 2020

# preparation --------
rm(list = ls())
cat("\014")
setwd("~/Lab/caledna_transect/")

library(dplyr)
library(gradientForest)
source("./scripts/function_transect.R")
# def functions --------

# def variables --------

# load data --------
gf_deco = loadRData("./derive_data/step6_gradient_forest/gf_deco_all_Family_Presence_2000_2_0.05_FALSE_17.RData")
gf_nocoast = loadRData("./derive_data/step6_gradient_forest/gf_nocoast_all_Family_Presence_2000_2_0.05_FALSE_17.RData")

# main --------
imp_deco = named_num2df(importance(gf_deco), c("var", "imp_deco")) %>% dplyr::mutate(rank_deco = 1:33)
imp_nocoast = named_num2df(importance(gf_nocoast), c("var", "imp_nocoast")) %>% dplyr::mutate(rank_nocoast = 1:33)

imp_comp = dplyr::left_join(x = imp_deco, y = imp_nocoast, by = 'var') %>%
    dplyr::mutate(imp_diff = imp_nocoast - imp_deco,
                  rank_diff = rank_nocoast - rank_deco) 
# difference denoted in ranks 
# cleanup --------

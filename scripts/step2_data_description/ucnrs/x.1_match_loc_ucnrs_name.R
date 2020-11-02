# Title: UCNRS site names match with metadata --------
# Author: Meixi Lin
# Date: Fri May 17 15:02:56 2019
# Author:
# Date:
# Modification:

# preparation --------
rm(list = ls())
cat("\014")

setwd("~/Lab/caledna_transect/")
library(dplyr)
library(phyloseq)
date() # the execution date

# load data -------
bycoord2 <- read.csv("./final_data/ucnrs/sites_ucnrs_final.csv", stringsAsFactors = F, row.names = 1) 
load("./derive_data/phy_deco/phydeco_uc.RData")
fauna <- read.csv(file = "./raw_data/ucnrs/UC-NRS-Species-List_v0.csv", stringsAsFactors = F)
flora <- read.csv(file = "./raw_data/ucnrs/reserve_plant_list_v0.csv", stringsAsFactors = F)

# first check reserve name --------
locedna <- sort(unique(bycoord2$loc))
rsedna <- sort(unique(bycoord2$Name))
rsflora <- sort(unique(flora$Reserve))
rsfauna <- sort(unique(fauna$X))
rsnames <- cbind(locedna, rsedna, rsflora, rsfauna)
write.csv(rsnames, file = "./derive_data/ucnrs/reserve_name_binder_v0.csv", row.names = F)
# changed by excel, now the name is matched

# Title: UCNRS site phyloseq build --------
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

# load data and remove duplicates -------
bycoord2 <- read.csv("./maps_vectors_rasters_local/ucnrs/Final_metadata_ucnrs_bycoord_0.00833.csv", stringsAsFactors = F) %>% 
    base::unique()
# write the final data out 
write.csv(bycoord2, file = "./final_data/ucnrs/sites_ucnrs_final.csv")

# subset by sample name and build phyloseq object --------
bycoord2 <- read.csv("./final_data/ucnrs/sites_ucnrs_final.csv", stringsAsFactors = F, row.names = 1) 
load("./derive_data/phy_deco/phydeco.RData")
ucnrs <- bycoord2$MatchName

phydeco_uc <- lapply(phydeco, function(xx) {
    xx <- subset_samples(xx, MatchName %in% ucnrs)
})

save(phydeco_uc, file = "./derive_data/phy_deco/phydeco_uc.RData")


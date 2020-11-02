# Title: compare the coastal final data and update the current version --------
# Author: Meixi Lin
# Date: Thu May 23 18:17:38 2019
# Author:
# Date:
# Modification:

# preparation --------
rm(list = ls())
cat("\014")
setwd("/Users/linmeixi/UCLA/Lab/abiotic_transect")

library(dplyr)

# load data --------
load(file = "./final_data/Final_metadata.RData")
coastmeta <- read.csv(file = "./other_version_data/dannise_coastal/Final_metadata_coast_May2019.csv", stringsAsFactors = F)

biom0 <- biom
biom <- biom %>%
    mutate(minorhab = as.character(minorhab))

# define minor habitat variable to change
tochange <- c("K0154A1","K0154A2","K0155A1","K0155A2","K0156A1","K0157A2","K0160B1","K0212B1","K0212C2","K0213C1","K0214B1","K0217A1","K0217B1","K0217C2","K0218A1")

# now perform the join -------
for (ii in tochange) {
    biomid <- which(biom$MatchName %in% ii)
    coastid <- which(coastmeta$MatchName %in% ii)
    print(biom[biomid,'MatchName'] == coastmeta[coastid, 'MatchName'])
    print(paste(biom[biomid,'MatchName'], 
          biom[biomid, 'minorhab'],
          coastmeta[coastid, 'new.minor.habitat'], sep = "||"))
    biom[biomid, 'minorhab'] <- coastmeta[coastid, 'new.minor.habitat']
} 

# look at the minorhabitat names 
table(biom$minorhab)
# convert back to factor 
biom$minorhab <- as.factor(biom$minorhab)

# save the change --------
save(biom, file = "./final_data/Final_metadata.RData")
write.csv(biom, file = "./final_data/metadata/Final_metadata_05312019.csv")




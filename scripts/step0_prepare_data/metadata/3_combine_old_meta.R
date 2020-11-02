# Title: Combine old and new metadata layers --------
# Author: Meixi Lin
# Date: Tue Apr  9 10:10:23 2019
# Author:
# Date:
# Modification:

# preparation --------
rm(list = ls())
cat("\014")

setwd("~/Lab/caledna_transect/")
library(dplyr)
library(ggplot2)
library(devtools)
source_gist("524eade46135f6348140") # this for generating the qqplot
date() # the execution date

# load data -------
oldmeta <- read.csv(file = "./final_data/Final_metadata_v0.csv", stringsAsFactors = F)
newmeta <- read.csv(file = "./derive_data/metadata/20190430/x100m/Final_newmeta.csv", stringsAsFactors = F, row.names = 1)
colnames(oldmeta)
colnames(newmeta)

# newmeta -----
toappend <- oldmeta %>% 
    dplyr::select(MatchName, date, loc, ecoregion, majorhab, minorhab, transect, SoS, clust, Zeta_4)

finalmeta <- dplyr::left_join(x = toappend, y = newmeta, by = "MatchName")

# add up the gps accuracy layer --------
meta061318 <- read.csv(file = "./raw_data/metadata/Metadata_versions/Final_metadata_06132018_all_samples.csv", stringsAsFactors = F) %>% 
    dplyr::select(MatchName, GPS_resolution)
# new meta 
finalmeta <- dplyr::left_join(x = finalmeta, y = meta061318, by = "MatchName")
# move up after longitude 
finalmeta <- finalmeta[,c(1:12, 69, 13:68)]
colnames(finalmeta)[13] <- "gpsres"
write.csv(x = finalmeta, file = "./final_data/Final_metadata_05012019.csv", row.names = F)

# total not na points 
hh <- is.na(finalmeta[,c(41:46)])
table(rowSums(hh))

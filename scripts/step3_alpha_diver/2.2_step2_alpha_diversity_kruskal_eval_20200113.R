# Title: evaluate the kruskal testing result and output table --------
# Author: Meixi Lin
# Date: Wed May 22 10:45:02 2019
# Author:
# Date: Mon Jan 13 10:47:02 2020
# Modification:

# preparation --------
rm(list = ls())
cat("\014")
setwd("~/Lab/caledna_transect")

# create output directory
indir <- "./derive_data/step3_alpha_diver/kruskalres/"
outdir <- "./derive_data/step3_alpha_diver/kruskalres/tables/"
dir.create(outdir, recursive = T)

library(vegan)
library(tibble)
library(dplyr)
library(ggplot2)
library(phyloseq)
require(FSA)
require(ggpubr)

source("./scripts/function_transect.R")

date() # the execution date

# parameter def ---------
mymeasures <- c("Observed", "Shannon")
mysig <- 0.05 
mycutoff <- 4
nobs <- length(catlist) * length(mymeasures)

# define function --------
adjust_p_bonf <- function(pvalue, nobs) {
    palter = pvalue*nobs 
    if (palter > 1) {palter = 1}
    return(palter)
}

# start a dataframe --------
kruskaldf <- data.frame(matrix(nrow = nobs * length(primers_commeco), ncol = 7))
colnames(kruskaldf) <- c("Metabarcode", "Category", "Measure", "Df", "P.value", "Adj.P.value", "Significant")
kruskaldf$Metabarcode <- rep(primers_commeco, each = length(catlist) * length(mymeasures))
kruskaldf$Category <- rep(catlist, each = length(mymeasures), times = length(primers_commeco))
kruskaldf$Measure <- rep(mymeasures, times = length(primers_commeco) * length(catlist))

# load and write the kruskal out --------
for (ii in  1:length(primers_commeco)) {
    print(primers_commeco[ii])
    for (jj in 1:length(catlist)) {
        print(catlist[jj])
        for (kk in 1:length(mymeasures)) {
            print(mymeasures[kk])
            load(file = paste0(indir, primers_commeco[ii], "_", catlist[jj], "_", mymeasures[kk], "_kruskal_20200330.RData"))
            index <- (kruskaldf$Metabarcode == primers_commeco[ii]) & (kruskaldf$Category == catlist[jj]) & (kruskaldf$Measure == mymeasures[kk])
            kruskaldf[index, 4:5] <- c(unname(kruskal$parameter), kruskal$p.value)
            kruskaldf[index, 6] <- adjust_p_bonf(kruskaldf[index,5], nobs = nobs)
            kruskaldf[index, 7] <- kruskaldf[index, 6] < mysig
            # if the result is significant, output the dunn test result 
            if (kruskaldf[index, 7] == TRUE) {
                load(file = paste0(indir, primers_commeco[ii], "_", catlist[jj], "_", mymeasures[kk], "_dunn_20200330.RData"))
                notation <- paste0(primers_commeco[ii], "_", catlist[jj], "_", mymeasures[kk], ",,,\n")
                cat(notation,file=paste0(outdir, primers_commeco[ii], "_", catlist[jj], "_", mymeasures[kk], "_dunn_res.csv"), append = F)
                write.table(dunn$res, file = paste0(outdir, primers_commeco[ii], "_", catlist[jj], "_", mymeasures[kk], "_dunn_res.csv"), sep = ",", quote = T, row.names = F, append = T)
            }
        }
    }
}

# write out value --------
write.csv(kruskaldf, file = paste0(outdir, "kruskal_alpha_div.csv"), row.names = F)

# # to generate the tables in S5.2 ~ later, used terminal command to concatenate results by variable --------
# # Do the commented lines in terminal 
# outdir="./derive_data/step3_alpha_diver/kruskalres/tables/"
# cd $outdir

# paste -d"\t" *_dunn_res.csv > summary_dunn_res.csv

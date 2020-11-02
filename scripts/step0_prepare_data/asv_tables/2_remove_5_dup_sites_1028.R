# Title: remove 5 duplicated sites, prepare final data --------
# Author: Meixi Lin
# Date: Sun Oct 28 23:49:50 2018
# Author:
# Date:
# Modification:

# preparation --------
rm(list = ls())
cat("\014")

setwd("~/Lab/caledna_transect/")
# load packages
library(dplyr)
library(ggplot2)
library(reshape2)
library(phyloseq)
library(ranacapa)
source("./scripts/function_transect.R")

# read in data --------
asv.16S <- read.csv("./other_version_data/Oct_25_2018_CALeDNA_Transect_Results/16S_table_3_Oct25_2018.csv")
asv.18S <- read.csv("./other_version_data/Oct_25_2018_CALeDNA_Transect_Results/18S_table_3_Oct25_2018.csv")
asv.CO1 <- read.csv("./other_version_data/Oct_25_2018_CALeDNA_Transect_Results/CO1_table_3_Oct25_2018.csv")
asv.FITS <- read.csv("./other_version_data/Oct_25_2018_CALeDNA_Transect_Results/FITS_table_3_Oct25_2018.csv")
asv.PITS <- read.csv("./other_version_data/Oct_25_2018_CALeDNA_Transect_Results/PITS_table_3_Oct25_2018.csv")


# rename scrub data --------
matchname <- as.matrix(colnames(asv.16S)[-1])
newnames <- apply(matchname, 1, function(x) {
    if (!is.na(colsplit(x, "[.]",1:2))[,2]) {
        x <- paste0(colsplit(x, "[.]",1:2)[,1], colsplit(x, "[.]",1:2)[,2])
    } else {
        x <- x
    }
    return(x)
})
newnames <- c("sum.taxonomy", unlist(newnames))
colnames(asv.16S) <- newnames
colnames(asv.18S) <- newnames
colnames(asv.CO1) <- newnames
colnames(asv.FITS) <- newnames
colnames(asv.PITS) <- newnames

droptran <- colsplit(newnames[-1], "[_]", 1:2)[,1] # should be 283

# compare the two sequencing run results from same sites using biplot --------
# find the samples that were sequenced twice
dup <- droptran[duplicated(droptran)]
primer <- c("16S", "18S", "CO1", "FITS", "PITS")

# plot the sample distribution for the same sample --------
png(filename = "./plots/1019/dup_sample_reads.png", width = 20, height = 20, unit ="in", res = 150)
pp <- par(mfrow = c(length(dup),length(primer)), cex = 0.7)
for (ii in 1:length(dup)) {
    for (jj in 1:length(primer)) {
        # take in the asv
        asv <- eval(as.name(paste0("asv.", primer[jj])))
        forplot <- asv %>%
            select(sum.taxonomy, starts_with(dup[ii]))
        # only plot both non-zero values 
        forplot <- forplot[rowSums(forplot[,2:3]) > 0,]
        xylim <- range(forplot[,2], forplot[,3])
        plot(forplot[,2], forplot[,3], xlab = colnames(forplot)[2], ylab = colnames(forplot)[3],
             main = primer[jj], xlim = xylim, ylim = xylim)
        abline(a = 0, b = 1, col = "red")
    }
}
par(pp)
dev.off()

# glom taxa by "family" level and plot the distribution again --------
txlevel <- "Family"
png(filename = "./plots/1019/dup_sample_reads_family.png", width = 20, height = 20, unit ="in", res = 150)
pp <- par(mfrow = c(length(dup),length(primer)), cex = 0.7)
for (ii in 1:length(dup)) {
    for (jj in 1:length(primer)) {
        # take in the asv
        asv <- eval(as.name(paste0("asv.", primer[jj])))
        forplot <- asv %>%
            glom_tax_df(taxlevel = txlevel) %>%
            select(starts_with(dup[ii])) 
        forplot <- forplot[-which(rownames(forplot) == "unknown"),]
        # only plot both non-zero values 
        # forplot <- forplot[rowSums(forplot) > 0,]
        xylim <- range(forplot[,1], forplot[,2])
        plot(forplot[,1], forplot[,2], xlab = colnames(forplot)[1], ylab = colnames(forplot)[2],
             main = primer[jj], xlim = xylim, ylim = xylim)
        abline(a = 0, b = 1, col = "red")
    }
}
par(pp)
dev.off()

# drop all the sites from scrub --------
todrop <- paste0(dup, "_S")

for (jj in 1:length(primer)) {
    asv <- eval(as.name(paste0("asv.", primer[jj])))
    new.asv <- asv %>% select(-one_of(todrop)) %>%
        select(sum.taxonomy, order(colnames(.)[-1]) + 1) # this may be a little counterintuitive, it reorder columns by its name 
    colnames(new.asv)[-1] <- colsplit(colnames(new.asv)[-1], "[_]", 1:2)[,1]
    write.csv(new.asv, 
              file = paste0("./final_data/deco_3/asv_deco_dedup_", primer[jj], ".csv"), 
              row.names = F)
} 
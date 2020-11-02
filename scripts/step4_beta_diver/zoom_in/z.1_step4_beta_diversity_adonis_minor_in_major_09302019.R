# Title: PCoA plots for minor habitat within major habitat --------
# Author: Meixi Lin
# Date: Wed Sep 18 11:05:55 2019

# preparation --------
options(echo = TRUE)
rm(list = ls())
cat("\014")

setwd("~/Lab/caledna_transect/")

library(dplyr)
library(vegan)
library(phyloseq)
date() # the execution date

source("./scripts/function_transect.R")

# define variables --------
nnperm = 2999
mycutoff = 0
majorlist <- c("Aquatic", "Herbaceous-Dominated", "Shrub-Dominated", "Tree-Dominated")
adres <- vector(length = length(primers_commeco) * length(majorlist), mode = "list")
bdspres <- vector(length = length(primers_commeco) * length(majorlist), mode = "list")
indir = "./derive_data/step4_beta_diver/"
outdir = "./derive_data/step4_beta_diver/zoom_in/"
dir.create(outdir, recursive = T)

# load data -------
load(file = "./derive_data/phy_rare/phyrare.RData")
load(file = paste0(indir, "jadiss.RData"))
load(file = "./final_data/Final_metadata.RData")

# for each variables --------
for (ii in 1:length(primers_commeco)) {
    # get phyloseq objects 
    primer <- primers_commeco[ii]
    physeq1 <- phyrare[[primer]]
    mydiss <- jadiss[[primer]]
    
    # subset phyloseq by major habitat 
    aqua <- subset_samples(physeq1, majorhab == "Aquatic")
    herb <- subset_samples(physeq1, majorhab == "Herbaceous-Dominated")
    shru <- subset_samples(physeq1, majorhab == "Shrub-Dominated")
    tree <- subset_samples(physeq1, majorhab == "Tree-Dominated")
    
    # perform PCoA on all of them
    phy.majorhab = list(aqua, herb, shru, tree)
    
    # perform stats testing 
    # perform adonis and betadispersion test
    for (jj in 1:length(majorlist)) {
        groupv = majorlist[jj]
        id <- (ii - 1) * length(majorlist) + jj
        myheader = paste(primer, majorlist[jj], sep = "|")
        print(myheader)
        myphyseq = phy.majorhab[[jj]]
        # calculate dissimilarity 
        diss <- phyloseq::distance(physeq = myphyseq, method = "jaccard", binary = T)
        sampledf <- data.frame(sample_data(myphyseq))
        # discard the NA values
        sampledf1 <- sampledf %>%
            tidyr::drop_na(., minorhab)
        diss1 <- usedist::dist_subset(diss, sampledf1$MatchName)
        # adonis 
        adres[[id]] <- adonis(diss1 ~ minorhab, data = sampledf1, permutations = nnperm)
        names(adres)[id] <- paste(primer, groupv, sep = "|")
        print(adres[[id]])
        # beta dispersion 
        bdspres[[id]] <- with(sampledf1, betadisper(diss1, minorhab))
        names(bdspres)[id] <- paste(primer, groupv, sep = "|")
        print(bdspres[[id]])
    }
}

# write out the result ---------
save(adres, file = paste0(outdir, "adonis_result_zoom_in_cutoff_", mycutoff, "_",Sys.Date(),".RData"))
save(bdspres, file = paste0(outdir, "betadisper_result_zoom_in_cutoff_", mycutoff, "_",Sys.Date(), ".RData"))

# transform to a csv --------
# make a df ---------
addf <- data.frame(matrix(nrow = length(primers_commeco) * length(majorlist), ncol = 8))
colnames(addf) <- c("name", "var.df", "res.df", "var.ssq", "res.ssq", "f.stat", "var.r2", "p.value")
bdspdf <- data.frame(matrix(nrow = length(primers_commeco) * length(majorlist), ncol = 8))
colnames(bdspdf) <- c("name", "var.df", "res.df", "var.ssq", "res.ssq", "f.stat", "nperm", "p.value")
# write out the adonis --------
for (ii in 1:length(adres)) {
    xx <- adres[[ii]]
    addf[ii, "name"] <- names(adres)[ii]
    addf[ii, "var.df"] <- xx$aov.tab$Df[1]
    addf[ii, "res.df"] <- xx$aov.tab$Df[2]
    addf[ii, "var.ssq"] <- xx$aov.tab$SumsOfSqs[1]
    addf[ii, "res.ssq"] <- xx$aov.tab$SumsOfSqs[2]
    addf[ii, "f.stat"] <- xx$aov.tab$F.Model[1]
    addf[ii, "var.r2"] <- xx$aov.tab$R2[1]
    addf[ii, "p.value"] <- xx$aov.tab$`Pr(>F)`[1]
}

write.csv(addf, file = paste0(outdir, "adonis_beta_div_zoom_in_", Sys.Date(), ".csv"), row.names = F)

# write out the betadispersion ---------
for (ii in 1:length(bdspres)) {
    colnames(bdspdf) <- c("name", "var.df", "res.df", "var.ssq", "res.ssq", "f.stat", "nperm", "p.value")
    xx <- permutest(bdspres[[ii]], permutations = how(nperm = nnperm))
    bdspdf[ii, "name"] <- names(bdspres)[ii]
    bdspdf[ii, "var.df"] <- xx$tab$Df[1]
    bdspdf[ii, "res.df"] <- xx$tab$Df[2]
    bdspdf[ii, "var.ssq"] <- xx$tab$`Sum Sq`[1]  
    bdspdf[ii, "res.ssq"] <- xx$tab$`Sum Sq`[2]
    bdspdf[ii, "f.stat"] <- xx$tab$`F`[1]
    bdspdf[ii, "nperm"] <- xx$tab$N.Perm[1]
    bdspdf[ii, "p.value"] <- xx$tab$`Pr(>F)`[1]
}

write.csv(bdspdf, file = paste0(outdir, "beta_dispersion_beta_div_zoom_in_", Sys.Date(), ".csv"), row.names = F)

# ending --------
date()
closeAllConnections()


# Title: Basic beta diversity: adonis and beta disper --------
# Author: Meixi Lin
# Date: Sun Aug  4 23:55:41 2019

# preparation --------
rm(list = ls())
cat("\014")

setwd("~/Lab/caledna_transect")

options(echo = TRUE)
outdir = "./derive_data/step4_beta_diver/logs/" 
sink(paste0("./derive_data/step4_beta_diver/logs/adonis_c4_", Sys.Date(),".log"), append = F)

library(tibble)
library(dplyr)
library(ggplot2)
library(phyloseq)
library(microbiomeSeq)
library(reshape2)
library(vegan) # microbiomeSeq calls a function "ordiellipse"
library(adespatial) # microbiomeSeq calls a function "beta.div"
library(grid) # microbiomeSeq calls a function "textGrob"
library(gridExtra) # microbiomeSeq calls a function "ttheme_minimal"
library(fso)
library(usedist)

source("./scripts/function_transect.R")
myseed <- 17
mycutoff <- 4
nnperm <- 2999

date() # the execution date

# load data --------
load("./derive_data/phy_rare/phyrare.RData")
load("./derive_data/step4_beta_diver/jadiss.RData")

# start a list for store the testing results ---------
adres <- vector(length = length(primers_commeco) * length(catlist), mode = "list")
bdspres <- vector(length = length(primers_commeco) * length(catlist), mode = "list")
    
# here it starts ---------
for (ii in 1:length(primers_commeco)) {
    # get value from other 
    primer <- primers_commeco[ii]
    physeq1 <- phyrare[[primer]]
    diss <- jadiss[[primer]]
    plotdir <- paste0("./plots/step4_beta_diver/", primer, "/bdsp/")
    dir.create(plotdir, recursive = T)
    
    print(dim(otu_table(physeq1)))
    sampledf <- data.frame(sample_data(physeq1))
    
    for (jj in 1:length(catlist)) {
        id <- (ii - 1) * length(catlist) + jj
        groupv <- catlist[jj]
        print(paste(primer, groupv, sep = "|"))
        # delete the sites that does not have at least 5 elements in that category
        category <- get_variable(physeq1, varName = groupv)
        sumby <- table(category) %>% as.data.frame() %>% 
            filter(Freq > mycutoff)
        samid <- as.character(sampledf[(sampledf[,groupv] %in% sumby$category), 'MatchName'])
        diss1 <- usedist::dist_subset(diss, samid)
        sampledf1 <- sampledf[(sampledf[,groupv] %in% sumby$category), ]

        # adonis 
        adres[[id]] <- adonis(diss1 ~ eval(as.name(groupv)), data = sampledf1, permutations = nnperm)
        names(adres)[id] <- paste(primer, groupv, sep = "|")
        print(adres[[id]])
        # beta dispersion 
        bdspres[[id]] <- with(sampledf1, betadisper(diss1, eval(as.name(groupv))))
        names(bdspres)[id] <- paste(primer, groupv, sep = "|")
        print(bdspres[[id]])
        # print(anova(bdspres[[id]]))
        print(permutest(bdspres[[id]], permutations = how(nperm = nnperm)))
    }
}

# save the data --------
dir.create(path = "./derive_data/step4_beta_diver/adonis_bdsp/", recursive = T)
save(adres, file = paste0("./derive_data/step4_beta_diver/adonis_bdsp/adonis_result_cutoff_", mycutoff, "_",Sys.Date(),".RData"))
save(bdspres, file = paste0("./derive_data/step4_beta_diver/adonis_bdsp/betadisper_result_cutoff_", mycutoff, "_",Sys.Date(), ".RData"))

sink()
closeAllConnections()

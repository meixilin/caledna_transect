# Title: alpha diversity and category variables kruskal wallis testing --------
# Author: Meixi Lin
# Date: Tue May  7 12:34:27 2019
# Author:
# Date:
# Modification:

# preparation --------
rm(list = ls())
cat("\014")
setwd("~/Lab/caledna_transect")

# create output directory
outdir <- "./derive_data/step3_alpha_diver/kruskalres/"
dir.create(outdir, recursive = T)
sink(paste0(outdir, "kruskal_rare_logs_20200330.log"))

library(vegan)
library(tibble)
library(dplyr)
library(ggplot2)
library(phyloseq)
require(FSA)
require(ggpubr)

source("./scripts/function_transect.R")

date() # the execution date
# import phyloseq object --------
load("./derive_data/phy_rare/phyrare.RData")
load("./derive_data/step3_alpha_diver/adivrare_06012019.RData")

mymeasures <- c("Observed", "Shannon")
mysig <- 0.05 
mycutoff <- 4

# calculate kruskal wallis based on the variable is categorical and rarefied dataset --------
for (ii in  1:length(primers_commeco)) {
    print(primers_commeco[ii])
    psdata <- phyrare[[ii]]
    print(psdata)
    adiv <- adivrare[[ii]]
    
    # change the jj value --------
    for (jj in 1:length(catlist)) {
        print(catlist[jj])
        category <- get_variable(psdata, varName = catlist[jj])
        sumby <- table(category) %>% as.data.frame() %>% 
            filter(Freq > mycutoff)
        # if all data are scatterred, 0 observation for sumby  
        if (nrow(sumby) < 3) {
            print(nrow(sumby))
            print(paste("for", catlist[jj], "no data for stats test."))
            next
        }
        print(sort(table(category)))
        for (kk in 1:length(mymeasures)) {
            print(mymeasures[kk])
            adivvar <- adiv[,kk]
            # filter out samples with less than 5 data in that category 
            fortest <- data.frame(adivvar, category) %>%
                filter(category %in% sumby$category)
            fortest$category <- droplevels(fortest$category)
            # sort the data by mean alpha diversity level
            aorder <- fortest %>% 
                dplyr::group_by(category) %>% 
                dplyr::summarise(aorder = median(adivvar)) %>% 
                dplyr::arrange(aorder)
            # reorder category level for plotting 
            fortest$category <- factor(fortest$category, 
                                       levels = levels(fortest$category)[as.integer(aorder$category)])
            # perform stats test 
            kruskal <- kruskal.test(adivvar ~ category, data = fortest)
            print(kruskal)
            # add post-hoc dunn test 
            dunn <- FSA::dunnTest(adivvar ~ category, data = fortest, method = "bonferroni")
            
            save(kruskal, file = paste0(outdir, primers_commeco[ii], "_", catlist[jj], "_", mymeasures[kk], "_kruskal_20200330.RData"))
            save(dunn, file = paste0(outdir, primers_commeco[ii], "_", catlist[jj], "_", mymeasures[kk], "_dunn_20200330.RData"))
            rm(kruskal)
            rm(dunn)
        } 
    }
}


# cleaning up --------
sink()
closeAllConnections()

# Title: Basic beta diversity: adonis and beta disper --------
# Author: Meixi Lin
# Date: Sun Aug  4 23:55:41 2019

# preparation --------
rm(list = ls())
cat("\014")

setwd("~/Lab/caledna_transect")
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

outdir = "./derive_data/step4_beta_diver/nocoast/"
sink(paste0(outdir, "adonis_nocoast_c4_08052019.log"), append = T)
source("./scripts/function_transect.R")
bothokflag <- F
myseed <- 17
mycutoff <- 4
nnperm = 2999

date() # the execution date

# load data --------
load("./derive_data/phy_rare/phyrare.RData")
jadiss <- loadRData(paste0(outdir, "jadiss_nocoast.RData"))

phyrare <- lapply(phyrare, function(xx) {
    # table(sample_sums(xx))
    xx <- subset_samples(physeq = xx, transect != "Coastal")
    return(xx)
})
phyrare

# start a list for store the testing results ---------
adres <- vector(length = length(primers) * length(catlist), mode = "list")
bdspres <- vector(length = length(primers) * length(catlist), mode = "list")
    
# here it starts ---------
for (ii in 1:length(primers)) {
    # get value from other 
    primer <- primers[ii]
    physeq <- phyrare[[primer]]
    diss <- jadiss[[primer]]
    plotdir <- paste0("./plots/step4_beta_diver/", primer, "/bdsp_nocoast/")
    dir.create(plotdir, recursive = T)
    
    # delete the not okay sites 
    if (bothokflag == T) {
        physeq1 <- subset_samples(physeq, bothok == bothokflag)
        physeq1 <- prune_samples(sample_sums(physeq1) > 0, physeq1)
    } else {
        physeq1 <- physeq
    }
    
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
        # adres[[id]] <- adonis(diss ~ eval(as.name(groupv)), data = sampledf)
        adres[[id]] <- adonis(diss1 ~ eval(as.name(groupv)), data = sampledf1, permutations = nnperm)
        names(adres)[id] <- paste(primer, groupv, sep = "|")
        print(adres[[id]])
        # beta dispersion 
        bdspres[[id]] <- with(sampledf1, betadisper(diss1, eval(as.name(groupv))))
        names(bdspres)[id] <- paste(primer, groupv, sep = "|")
        print(bdspres[[id]])
        # print(anova(bdspres[[id]]))
        print(permutest(bdspres[[id]], permutations = how(nperm = nnperm)))
        # print(TukeyHSD(bdspres[[id]]))
        
        # start plotting 
        pdf(file = paste0(plotdir, "beta_diver_bdsp_", groupv, "_cutoff_", mycutoff ,"_bija.pdf"), height = 9, width = 4)
        pp <- par(mfrow = c(3,1))
        plot(bdspres[[id]], label = F)
        plot(bdspres[[id]], label = T)
        boxplot(bdspres[[id]])
        par(pp)
        dev.off()
    }
}

# save the data 
save(adres, file = paste0(outdir, "adonis_result_nocoast_cutoff_", mycutoff, ".RData"))
save(bdspres, file = paste0(outdir, "betadisper_result_nocoast_cutoff_", mycutoff, ".RData"))

sink()
dev.off()
closeAllConnections()

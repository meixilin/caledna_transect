# Title: output for capscale --------
# Author: Meixi Lin
# Date: Fri Aug  9 12:39:06 2019
# Author: Meixi Lin 
# Date: 
# Modification: Change directories 

# preparation --------
rm(list = ls())
cat("\014")

setwd("~/Lab/caledna_transect/")

options(echo = T)
sink(paste0("./derive_data/step4_beta_diver/logs/cap_latlong_output_c4_", Sys.Date(), ".log"), append = F)
library(dplyr)
library(ggplot2)
date() # the execution date

source("./scripts/function_transect.R")

# set variables --------
mycutoff <- 4
sysdate <- "2020-03-31" # date when the data was generated 

# load data -------
load(file = paste0("./derive_data/step4_beta_diver/cap_latlong/cap_result_", sysdate, ".RData"))
load(file = paste0("./derive_data/step4_beta_diver/cap_latlong/varpar_result_", sysdate, ".RData"))

# for this analyses category list of interest is short, only habitat and soil property 
# catlist = catlist[c(2,3,4,5,6,8,9)]
catlist = c("ecoregion","majorhab","minorhab","transect","SoS","taxousda","NLCD")

capdfcolnames = c("Metabarcode", "Variable", "Condition", "N.obs","R^2", "Adj.R^2", "anova.var.df", "anova.res.df", "anova.var.ssq", "anova.res.ssq", "anova.f.stats", "anova.p.value", "CAP1", "CAP2", "MDS1", "MDS2")
capdf <- data.frame(matrix(nrow = length(primers_commeco) * length(catlist), ncol = length(capdfcolnames)))
colnames(capdf) <- capdfcolnames
if (length(caplist) != nrow(capdf)) {
    stop("data not of the same size, check! ")
}

varpardfcolnames = c("Metabarcode", "Variable_of_interest (x1)", "Spatial_Variable (x2)", "N.obs", "x1.df", "x2.df", "all.df", "x1.r2.adj", "x2.r2.adj", "all.r2.adj", "x1|x2.r2.adj", "x2|x1.r2.adj", "x1Ax2.r2.adj")
varpardf <- data.frame(matrix(nrow = length(primers_commeco) * length(catlist), ncol = length(varpardfcolnames)))
colnames(varpardf) <- varpardfcolnames
if (length(varparlist) != nrow(varpardf)) {
    stop("data not of the same size, check! ")
}

# initiate a list to store the anova.cca output
ancaplist <- vector(length = length(primers_commeco) * length(catlist), mode = "list")

# output result for capdf -------
for (ii in 1:length(primers_commeco)) {
    primer <- primers_commeco[[ii]]
    for (jj in 1:length(catlist)) {
        groupv <- catlist[jj]
        mycapname <- paste(primer, groupv, sep = "|")
        print(mycapname)
        # get the variables 
        capid <- (ii - 1) * length(catlist) + jj
        mycap <- caplist[[capid]]
        # check if it's the right cap 
        if (names(caplist)[capid] != mycapname) {
            stop("capname not match, check!")
        }
        # print(anova.cca(mycap, by = "terms", strata = sampledf1$clust))
        ancap <- anova.cca(mycap, by = "terms") # permutation test, could have slightly different results everytime
        ancaplist[[capid]]=ancap
        names(ancaplist)[capid]=mycapname
        capdf[capid, "Metabarcode"] <- primer
        capdf[capid, "Variable"] <- groupv
        capdf[capid, "Condition"] <- paste0(groupv, " + Condition(Longitude + Latitude)")
        capdf[capid, "N.obs"] <- nobs(mycap)
        capdf[capid, "R^2"] <- RsquareAdj(mycap)$r.squared
        capdf[capid, "Adj.R^2"] <- RsquareAdj(mycap)$adj.r.squared
        capdf[capid, "anova.var.df"] <- ancap$Df[1]
        capdf[capid, "anova.res.df"] <- ancap$Df[2]
        capdf[capid, "anova.var.ssq"] <- ancap$SumOfSqs[1]
        capdf[capid, "anova.res.ssq"] <- ancap$SumOfSqs[2]
        capdf[capid, "anova.f.stats"] <- ancap$`F`[1]
        capdf[capid, "anova.p.value"] <- ancap$`Pr(>F)`[1]
        capdf[capid, "CAP1"] <- vegan::eigenvals(mycap)["CAP1"]
        capdf[capid, "CAP2"] <- vegan::eigenvals(mycap)["CAP2"]
        capdf[capid, "MDS1"] <- vegan::eigenvals(mycap)["MDS1"]
        capdf[capid, "MDS2"] <- vegan::eigenvals(mycap)["MDS2"]
    }
}

# output results for varpart --------
varpardfcolnames = c("Metabarcode", "Variable_of_interest (x1)", "Spatial_Variable (x2)", "N.obs", "x1.df", "x2.df", "all.df", "x1.r2.adj", "x2.r2.adj", "all.r2.adj", "x1|x2.r2.adj", "x2|x1.r2.adj", "x1Ax2.r2.adj")
for (ii in 1:length(primers_commeco)) {
    primer <- primers_commeco[[ii]]
    for (jj in 1:length(catlist)) {
        groupv <- catlist[jj]
        myvarparname <- paste(primer, groupv, sep = "|")
        print(myvarparname)
        # get the variables 
        varparid <- (ii - 1) * length(catlist) + jj
        myvarpar <- varparlist[[varparid]]
        # check if it's the right cap 
        if (names(varparlist)[varparid] != myvarparname) {
            stop("varpar name not match, check!")
        }
        varpardf[varparid, "Metabarcode"] <- primer
        varpardf[varparid, "Variable_of_interest (x1)"] <- groupv
        varpardf[varparid, "Spatial_Variable (x2)"] <- myvarpar$tables[2]
        varpardf[varparid, "N.obs"] <- myvarpar$part$n
        varpardf[varparid, "x1.df"] <- myvarpar$part$fract[1, 'Df']
        varpardf[varparid, "x2.df"] <- myvarpar$part$fract[2, 'Df']
        varpardf[varparid, "all.df"] <- myvarpar$part$fract[3, 'Df']
        varpardf[varparid, "x1.r2.adj"] <- myvarpar$part$fract[1, 'Adj.R.squared']
        varpardf[varparid, "x2.r2.adj"] <- myvarpar$part$fract[2, 'Adj.R.squared']
        varpardf[varparid, "all.r2.adj"] <- myvarpar$part$fract[3, 'Adj.R.squared']
        varpardf[varparid, "x1|x2.r2.adj"] <- myvarpar$part$indfract[1, 'Adj.R.squared']
        varpardf[varparid, "x1Ax2.r2.adj"] <- myvarpar$part$indfract[2, 'Adj.R.squared']
        varpardf[varparid, "x2|x1.r2.adj"] <- myvarpar$part$indfract[3, 'Adj.R.squared']
    }
}

# report out the result --------
save(x=ancaplist, file = paste0("./derive_data/step4_beta_diver/cap_latlong/cap_result_anova_", sysdate, ".RData"))
write.csv(x = capdf, file = paste0("./derive_data/step4_beta_diver/cap_latlong/capscale_beta_div_", sysdate, ".csv"), row.names = F)
write.csv(x = varpardf, file = paste0("./derive_data/step4_beta_diver/cap_latlong/varpar_beta_div_", sysdate, ".csv"), row.names = F)
sink()
closeAllConnections()

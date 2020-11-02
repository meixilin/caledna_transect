# Title: Summary 20 replicated runs in gradient forest 
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Wed Mar 11 10:39:26 2020

# preparation --------
rm(list = ls())
cat("\014")

setwd("~/Lab/caledna_transect")
library(gradientForest)
library(dplyr)
source("./scripts/function_transect.R")

# def functions --------
join_dtlist <- function(dtlist, length, by) {
    mydt = dtlist[[1]]
    for (ii in 2:length) {
        mydt = dplyr::left_join(mydt, dtlist[[ii]], by = by, suffix = as.character(c(ii-1, ii)))
    }
    return(mydt)
}

# get the number of occurrences for each sample
get_family_n_gfY <- function(gfname, sdfam) {
    load(file = gfname)
    gfY=gf$Y
    rm(gf)
    familyname=colnames(gfY)
    for (jj in 1:ncol(gfY)) {
        gfY[,jj] <- as.integer(gfY[,jj])-1
    }
    countfam=colSums(gfY)
    countfam=named_num2df(countfam,yourcols = c("family", "Num_occurrence"))
    sdcountfam=dplyr::full_join(sdfam, countfam, by ='family')
    return(sdcountfam)
}
# def variables --------
indir = "./derive_data/step6_gradient_forest/"

# load data --------
load(paste0(indir, "stab_res_gf_deco_all_Family_Presence_2000_2_0.05_FALSE.RData")) # import an object called "gf_sumlist" 

# main --------
# get variable importance
impdt = lapply(gf_sumlist, function(xx) {
    named_num2df(xx$most_imp, c("var", "importance"))
})

impdt = join_dtlist(impdt, length = 20, by = "var") %>%
    tibble::column_to_rownames(var = "var") 
sd = apply(impdt, 1, sd)
write.csv(sd, file = paste0(indir, "stab_importance_sd.csv"))

# get total Rsq 
rsq = unlist(lapply(gf_sumlist, function(xx) {xx$rsq['ave.rsq']}))
sd(rsq)

# get family rsq
familydt = lapply(gf_sumlist, function(xx) {
    named_num2df(xx$imp_families, c("family", "rsq"))
})

familydt = join_dtlist(familydt, length = 20, by = "family") %>%
    tibble::column_to_rownames(var = "family") 
sdfam = named_num2df(apply(familydt, 1, sd), c("family", "sd"))
orderedfam = named_num2df(sort(gf_sumlist[[1]]$imp_families, decreasing = T), c("family", "rsq"))
sdfam = left_join(orderedfam, sdfam, by = "family")
write.csv(sdfam, file = paste0(indir, "stab_family_sd.csv"))

# additionally add the number of occurrences --------
sdfam=read.csv(file=paste0(indir, "stab_family_sd.csv"), row.names = 1, stringsAsFactors = FALSE)
countsdfam=get_family_n_gfY(gfname = paste0(indir, "gf_deco_all_Family_Presence_2000_2_0.05_FALSE_17.RData"),sdfam=sdfam)
write.csv(countsdfam, file = paste0(indir, "stab_family_sd_count_", Sys.Date(),".csv"))

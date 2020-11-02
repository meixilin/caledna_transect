# Title: Metadata summary (for Table S1)
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Sun Mar  1 18:23:03 2020

# preparation --------
rm(list = ls())
cat("\014")
setwd("~/Lab/caledna_transect/")

library(dplyr)
source("./scripts/function_transect.R")

# def functions --------
get_summary <- function(var) {
    sum = as.data.frame(table(biom[,var]))
    colnames(sum) = c(var, "Freq")
    write.csv(sum, file = paste0(outdir, var, "_distribution.csv"), row.names = F)
    return(sum)
}

get_summary_num <- function(var) {
    min = min(biom[,var], na.rm= T)
    max = max(biom[,var], na.rm= T)
    median = median(biom[,var], na.rm= T)
    mean = mean(biom[,var], na.rm= T)
    sd = sd(biom[,var], na.rm = T)
    missing = sum(is.na(biom[,var]))
    output = c(var, min, max, median, mean, sd, missing)
    return(output)
}
# def variables --------
outdir = "./derive_data/step2_data_description/summary_stats/"
dir.create(outdir, recursive = T)

# load data --------
load("./final_data/Final_metadata.RData")

# main --------
# get the summary tables for factor variables 
catlist2 = c("date", catlist)
for (ii in catlist2) {
    get_summary(ii)
}

# get the summary for numerical values 
sumnames = c("Variable", "min", "max", "median","mean", "sd", "missing")
biom_names = colnames(biom)[c(11:33, 35:66, 68,69)]
sumdt = data.frame(matrix(NA, nrow = length(biom_names), ncol = length(sumnames)))
colnames(sumdt) = sumnames

for (ii in 1:length(biom_names)) {
    sumdt[ii,] = get_summary_num(biom_names[ii])
}

write.csv(sumdt, file = paste0(outdir,"continuous_variable_distribution.csv"), row.names = F)

# cleanup --------
closeAllConnections()

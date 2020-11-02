# Title: generate summary plots for duplicated sites --------
# This is for generating summarized plots of 
# Author: Meixi Lin
# Date: Mon Jan  6 00:50:15 2020

# preparation --------
rm(list = ls())
cat("\014")

setwd("~/Lab/caledna_transect/")
# load packages
library(dplyr)
library(stringr)
library(ggplot2)
library(reshape2)
library(phyloseq)
library(ranacapa)
source("./scripts/function_transect.R")

outdir <- "./derive_data/step2_data_description/concordance/"
plotdir <- "./plots/step2_data_description/concordance/"

dir.create(outdir, recursive = T)
dir.create(plotdir, recursive = T)

# get arguments --------
rare = "deco"
txlevel = c("Family", "Genus", "Species")
cutoff = c(0,2,4,9,49,99)

# define functions --------
# for plotting
filter_chr <- function(xx) {
    if (str_detect(xx, "neither") | str_detect(xx, "K")) {
        xx <- NA
    } else {
        xx <- as.numeric(xx)
    }
    return(xx)
}

full_full_union <- function(listxx, len) {
    data = listxx[[1]]
    for (ii in 2:len) {
        data = dplyr::union(data, listxx[[ii]])
    }
    return(data)
}

# read in data --------
# this is for only presence 
txlevel1 <- paste(txlevel, rare, sep = "_")
dict <- expand.grid(txlevel1, cutoff) %>%
    as.data.frame() %>%
    mutate(name = paste0(Var1, "_", Var2,".csv"),
           fullname = paste0(Var1, "_", Var2,"_full.csv"))
files <- dict$name
# this is for both presence and total absence 
fullfiles <- dict$fullname

rawdata <- lapply(files, function(xx) {
    data <- read.csv(paste0(outdir, xx), row.names = 1) 
    data <- apply(data, c(1,2), filter_chr) # get rid of the annotation, convert to NA 
    data <- data %>%
        reshape2::melt(varnames = c("MatchName", "Metabarcode"), value.name = "Overlap") %>%
        mutate(Taxonomy_level = strsplit(xx, "_")[[1]][1],
               cutoff = as.integer(str_sub(strsplit(xx, "_")[[1]][3], start = 1, end = 1)),
               presence_only = T) 
})

fullrawdata <- lapply(fullfiles, function(xx) {
    data <- read.csv(paste0(outdir, xx), row.names = 1) 
    data <- apply(data, c(1,2), filter_chr) # get rid of the annotation, convert to NA 
    data <- data %>%
        reshape2::melt(varnames = c("MatchName", "Metabarcode"), value.name = "Overlap") %>%
        mutate(Taxonomy_level = strsplit(xx, "_")[[1]][1],
               cutoff = as.integer(str_sub(strsplit(xx, "_")[[1]][3], start = 1, end = 1)),
               presence_only = F)
})

rawdata <- full_full_union(rawdata, length(rawdata))
rawdata$Taxonomy_level[rawdata$Taxonomy_level == "Species"] = "Original_LCA"
fullrawdata <- full_full_union(fullrawdata, length(rawdata))
fullrawdata$Taxonomy_level[fullrawdata$Taxonomy_level == "Species"] = "Original_LCA"

# plot the overlap by taxonomy level --------
pp1 <- ggplot(data = rawdata[rawdata$Taxonomy_level == "Family",], aes(x = cutoff, y = Overlap, 
                                  color = MatchName)) + 
    geom_line(aes(group = MatchName)) +
    facet_wrap(. ~  Metabarcode) +
    geom_point() +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90))

pp2 <- ggplot(data = fullrawdata, aes(x = Taxonomy_level, y = Overlap, color = MatchName)) + 
    geom_line(aes(group = MatchName)) +
    facet_wrap(. ~  Metabarcode) +
    geom_point() +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90))

ggsave(filename = paste0(plotdir, "dup_site_txlevel_", rare, "_", cutoff, ".pdf"), plot = pp1, width = 8, height = 6)
ggsave(filename = paste0(plotdir, "dup_site_txlevel_", rare, "_", cutoff, "_full.pdf"), plot = pp2, width = 8, height = 6)


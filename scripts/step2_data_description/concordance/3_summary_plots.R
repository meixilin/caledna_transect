# Title: generate summary plots for duplicated sites --------
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

indir <- "./derive_data/step2_data_description/concordance/raw_tables/"
outdir <- "./derive_data/step2_data_description/concordance/tables/"
plotdir <- "./plots/step2_data_description/concordance/"

dir.create(outdir, recursive = T)
dir.create(plotdir, recursive = T)

# get arguments --------
rare = "rare"; cutoff = 0
# rare = "deco"; cutoff = 2
# rare = "deco"; cutoff = 99

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

# for join long tables
full_full_union <- function(listxx, len) {
    data = listxx[[1]]
    for (ii in 2:len) {
        data = dplyr::union(data, listxx[[ii]])
    }
    return(data)
}

# read in data --------
# this is for only presence 
files <- list.files(path = indir, pattern = paste0(rare, "_", cutoff, ".csv"))
# this is for both presence and total absence 
fullfiles <- list.files(path = indir, pattern = paste0(rare, "_", cutoff, "_full.csv"))

# for writing the supplementary table ---------
table_data <- lapply(files, function(xx) {
    data <- read.csv(paste0(indir, xx), row.names = 1) 
    # data <- apply(data, c(1,2), filter_chr) # get rid of the annotation, convert to NA 
    data <- data %>%
        as.matrix() %>%
        reshape2::melt(varnames = c("MatchName", "Metabarcode"), value.name = "Overlap") %>%
        mutate(Taxonomy_level = strsplit(xx, "_")[[1]][1],
               Min_read_per_sample = cutoff, 
               Rarefied = (rare == "rare"), 
               presence_only = T) 
})

full_table_data <- lapply(fullfiles, function(xx) {
    data <- read.csv(paste0(indir, xx), row.names = 1) 
    # data <- apply(data, c(1,2), filter_chr) # get rid of the annotation, convert to NA 
    data <- data %>%
        as.matrix() %>%
        reshape2::melt(varnames = c("MatchName", "Metabarcode"), value.name = "Overlap") %>%
        mutate(Taxonomy_level = strsplit(xx, "_")[[1]][1],
               Min_read_per_sample = cutoff, 
               Rarefied = (rare == "rare"), 
               presence_only = F) 
})

table_data <- full_full_union(table_data, length(table_data))
full_table_data <- full_full_union(full_table_data, length(full_table_data))
table_data <- dplyr::union(table_data, full_table_data)
table_data$Taxonomy_level[table_data$Taxonomy_level == "Species"] = "Original_LCA"
write.csv(table_data, file = paste0(outdir, "percent_overlap_", rare, "_", cutoff, ".csv"))

# for plotting ---------
rawdata <- lapply(files, function(xx) {
    data <- read.csv(paste0(indir, xx), row.names = 1) 
    data <- apply(data, c(1,2), filter_chr) # get rid of the annotation, convert to NA
    data <- data %>%
        as.matrix() %>%
        reshape2::melt(varnames = c("MatchName", "Metabarcode"), value.name = "Overlap") %>%
        mutate(Taxonomy_level = strsplit(xx, "_")[[1]][1],
               Min_read_per_sample = cutoff, 
               Rarefied = (rare == "rare"), 
               presence_only = T) 
})

fullrawdata <- lapply(fullfiles, function(xx) {
    data <- read.csv(paste0(indir, xx), row.names = 1) 
    data <- apply(data, c(1,2), filter_chr) # get rid of the annotation, convert to NA
    data <- data %>%
        as.matrix() %>%
        reshape2::melt(varnames = c("MatchName", "Metabarcode"), value.name = "Overlap") %>%
        mutate(Taxonomy_level = strsplit(xx, "_")[[1]][1],
               Min_read_per_sample = cutoff, 
               Rarefied = (rare == "rare"), 
               presence_only = F) 
})

rawdata <- full_full_union(rawdata, length(rawdata))
rawdata$Taxonomy_level[rawdata$Taxonomy_level == "Species"] = "Original_LCA"
fullrawdata <- full_full_union(fullrawdata, length(fullrawdata))
fullrawdata$Taxonomy_level[fullrawdata$Taxonomy_level == "Species"] = "Original_LCA"

# plot the overlap by taxonomy level --------
pp1 <- ggplot(data = rawdata, aes(x = Taxonomy_level, y = Overlap, color = MatchName)) + 
    geom_line(aes(group = MatchName)) +
    facet_wrap(. ~  Metabarcode) +
    geom_point() +
    ylim(0,1) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90))

pp2 <- ggplot(data = fullrawdata, aes(x = Taxonomy_level, y = Overlap, color = MatchName)) + 
    geom_line(aes(group = MatchName)) +
    facet_wrap(. ~  Metabarcode) +
    geom_point() +
    ylim(0,1) + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90))

ggsave(filename = paste0(plotdir, "dup_site_txlevel_", rare, "_", cutoff, ".pdf"), plot = pp1, width = 8, height = 6)
ggsave(filename = paste0(plotdir, "dup_site_txlevel_", rare, "_", cutoff, "_full.pdf"), plot = pp2, width = 8, height = 6)



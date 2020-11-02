# Title: Generate figures for figure 4
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Tue Mar 31 10:29:59 2020

# preparation --------
rm(list = ls())
cat("\014")
setwd("~/Lab/caledna_transect/")

library(ggplot2)
library(dplyr)
library(gridExtra)
library(ggpubr)

source("./scripts/function_transect.R")

# def functions --------
# get ggplot2 object pp's title
extract_gg_title <- function(pp) {
    mytitle <- pp$labels$title
    return(mytitle)
}

get_plots <- function(myprimer, myvar, pplist) {
    id = paste(myprimer, myvar, sep = "|")
    temp = unlist(lapply(pplist, extract_gg_title))
    index = which(temp %in% id)
    if (length(index) != 1) {
        print("No match or multiple match, check!")
        return(-1)
    } else {
        pp = pplist[[index]]
        return(pp)
    }
}

get_leg <- function(pp) {
    myleg <- ggpubr::get_legend(pp)
    leg <- ggpubr::as_ggplot(myleg) +
        theme(legend.background = element_rect(fill = "transparent"))
    return(leg)
}

# final modification 
final_mod <- function(pp) {
    pp <- pp +
        theme(# text = element_text(size = 16),
              legend.position = "none")
    return(pp)
}

# def variables --------
indir = "./derive_data/step4_beta_diver/"
plotdir = "./plots/step4_beta_diver/manuscript/"
dir.create(plotdir)
# load data --------
# for A,D,E,F: PCoA 
pcoaplots <- loadRData(paste0(indir, "pcoa/pcoa_category_plot.RData"))
# for B: CAP scale 
capplots <- loadRData(paste0(indir, "cap_latlong/cap_category_plot_2020-03-31.RData"))
# for C: excluding coastal
pcoanocoastplots <- loadRData(paste0(indir, "nocoast/pcoa_nocoast_category_plot.RData"))

# main --------
ppA1 <- get_plots(myprimer = "16S", myvar = "majorhab", pplist = pcoaplots)
ppA2 <- get_plots(myprimer = "18S", myvar = "majorhab", pplist = pcoaplots)
ppB1 <- get_plots(myprimer = "16S", myvar = "majorhab", pplist = capplots)
ppB2 <- get_plots(myprimer = "18S", myvar = "majorhab", pplist = capplots)
ppC1 <- get_plots(myprimer = "16S", myvar = "majorhab", pplist = pcoanocoastplots)
ppC2 <- get_plots(myprimer = "18S", myvar = "majorhab", pplist = pcoanocoastplots)
ppD1 <- get_plots(myprimer = "16S", myvar = "transect", pplist = pcoaplots)
ppD2 <- get_plots(myprimer = "18S", myvar = "transect", pplist = pcoaplots)
ppE1 <- get_plots(myprimer = "CO1", myvar = "NLCD", pplist = pcoaplots)
ppE2 <- get_plots(myprimer = "PITS", myvar = "NLCD", pplist = pcoaplots)
ppF1 <- get_plots(myprimer = "FITS", myvar = "SoS", pplist = pcoaplots)
ppF2 <- get_plots(myprimer = "PITS", myvar = "SoS", pplist = pcoaplots)

# make a list of these 
pplist <- list(ppA1, ppA2,ppB1, ppB2,ppC1, ppC2,ppD1, ppD2,ppE1, ppE2,ppF1, ppF2)

# suppress legend and get legend 
leglist <- lapply(pplist, get_leg)
pplist <- lapply(pplist, final_mod)

allpp <- grid.arrange(pplist[[1]],pplist[[2]], pplist[[3]],pplist[[4]], 
                      pplist[[5]],pplist[[6]], pplist[[7]],pplist[[8]],
                      pplist[[9]],pplist[[10]], pplist[[11]],pplist[[12]],
                      nrow = 3, ncol = 4) 
ggsave(filename = paste0(plotdir, "beta_diver_combined_ordination.pdf"), plot = allpp, height = 9, width = 12)
# cleanup --------
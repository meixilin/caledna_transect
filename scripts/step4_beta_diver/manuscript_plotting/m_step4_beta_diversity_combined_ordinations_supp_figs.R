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

# supp figures: PCoA, NLCD --------
mypcoa <- lapply(primers_commeco, get_plots, myvar = "NLCD", pplist = pcoaplots)
mypcoa <- lapply(mypcoa, final_mod)
allpp <- grid.arrange(grobs = mypcoa, ncol = 3) 
ggsave(filename = paste0(plotdir, "beta_diver_ordination_NLCD.pdf"), plot = allpp, height = 8, width = 12)

# supp figures: PCoA, transect --------
mypcoa <- lapply(primers_commeco, get_plots, myvar = "transect", pplist = pcoaplots)
mypcoa <- lapply(mypcoa, final_mod)
allpp <- grid.arrange(grobs = mypcoa, ncol = 3) 
ggsave(filename = paste0(plotdir, "beta_diver_ordination_transect.pdf"), plot = allpp, height = 8, width = 12)

# supp figures: PCoA, majorhab --------
mypcoa <- lapply(primers_commeco, get_plots, myvar = "majorhab", pplist = pcoaplots)
mypcoa <- lapply(mypcoa, final_mod)
allpp <- grid.arrange(grobs = mypcoa, ncol = 3) 
ggsave(filename = paste0(plotdir, "beta_diver_ordination_majorhab.pdf"), plot = allpp, height = 8, width = 12)

# supp figures: PCoA, SoS --------
mypcoa <- lapply(primers_commeco, get_plots, myvar = "SoS", pplist = pcoaplots)
mypcoa <- lapply(mypcoa, final_mod)
allpp <- grid.arrange(grobs = mypcoa, ncol = 3) 
ggsave(filename = paste0(plotdir, "beta_diver_ordination_SoS.pdf"), plot = allpp, height = 8, width = 12)

# supp figures: PCoA, no coastal --------
mypcoa <- lapply(primers_commeco, get_plots, myvar = "majorhab", pplist = pcoanocoastplots)
mypcoa <- lapply(mypcoa, final_mod)
allpp <- grid.arrange(grobs = mypcoa, ncol = 3, 
                      top = textGrob("PCoA with coastal sites removed",gp=gpar(fontsize=20,font=3))) 
ggsave(filename = paste0(plotdir, "beta_diver_ordination_nocoast_majorhab.pdf"), plot = allpp, height = 8, width = 12)

mypcoa <- lapply(primers_commeco, get_plots, myvar = "transect", pplist = pcoanocoastplots)
mypcoa <- lapply(mypcoa, final_mod)
allpp <- grid.arrange(grobs = mypcoa, ncol = 3,
                      top = textGrob("PCoA with coastal sites removed",gp=gpar(fontsize=20,font=3))) 
ggsave(filename = paste0(plotdir, "beta_diver_ordination_nocoast_transect.pdf"), plot = allpp, height = 8, width = 12)

# supp figures: CAPscale majorhab -------
mypcoa <- lapply(primers_commeco, get_plots, myvar = "majorhab", pplist = capplots)
mypcoa <- lapply(mypcoa, final_mod)
allpp <- grid.arrange(grobs = mypcoa, ncol = 3,
                      top = textGrob("Constrained Ordinations",gp=gpar(fontsize=20,font=3))) 
ggsave(filename = paste0(plotdir, "beta_diver_ordination_capscale_majorhab.pdf"), plot = allpp, height = 8, width = 12)
# cleanup --------
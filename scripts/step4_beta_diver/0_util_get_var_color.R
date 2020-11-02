# Title: Utility scripts for generating colors for each variables 
# Used across beta diversity plotting 
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Mon Apr 20 00:13:25 2020

# preparation --------
rm(list = ls())
cat("\014")

setwd("~/Lab/caledna_transect/")
library(pals)
library(dplyr)
library(ggplot2)
date() # the execution date

source("./scripts/function_transect.R")

# set some variable --------
plotdir = "./plots/step4_beta_diver/util/"
outdir = "./derive_data/step4_beta_diver/util/"
dir.create(plotdir, recursive = T, showWarnings = F)
dir.create(outdir, recursive = T, showWarnings = F)

# define functions --------
# custom color scales for each variables based in biom 
get_palette_variable <- function(biom, groupv) {
    # use biom so variable has consistent names/colour pair 
    candidate_colors = c(pals::polychrome(n = 36)[3:36], pals::alphabet(n = 26), pals::alphabet2(n = 26))
    groupv_levels = levels(biom[,groupv])
    nnlevels = nlevels(biom[,groupv])
    if (nnlevels > length(candidate_colors)) {
        break;
    }
    groupv_pal = c(pals::polychrome(n = 2), candidate_colors[1:nnlevels])
    names(groupv_pal) = c("Others", "unknown", groupv_levels)
    pals::pal.bands(groupv_pal)
    return(groupv_pal)
}

# generate a dummy plot key for every possible categorical values 
plot_palette_variable <- function(mypal, myheader) {
    value = 1:length(mypal)
    xx = names(mypal)
    dummy = data.frame(xx, value) 
    dummy$xx = factor(dummy$xx, levels = dummy$xx)
    pp <- ggplot(data = dummy, aes(x = xx, y = value, fill = xx)) +
        geom_bar(stat = "identity") + 
        scale_fill_manual(values = mypal) + 
        labs(fill = myheader) + 
        theme_bw() + 
        theme(legend.title = element_text(size = 10),
              legend.text = element_text(size = 9),
              legend.key.size = unit(0.7, 'lines'), 
              legend.justification = c("left", "top"))
    myleg <- ggpubr::get_legend(pp)
    leg <- ggpubr::as_ggplot(myleg) +
        theme(legend.background = element_rect(fill = "transparent"))
    return(leg)
}

plot_majorhab_variable3 <- function(mypal, myheader) {
    value = 1:length(mypal)
    xx = names(mypal)
    dummy = data.frame(xx, value) 
    dummy$xx = factor(dummy$xx, levels = dummy$xx)
    pp <- ggplot(data = dummy, aes(x = xx, y = value, fill = xx)) +
        geom_bar(stat = "identity") + 
        scale_fill_manual(values = mypal,
                          labels = c("Non-vegetated", "unknown", "Aquatic", "Herbaceous-Dominated", "Shrub-Dominated", "Tree-Dominated")) + 
        labs(fill = myheader) + 
        theme_bw() + 
        theme(legend.title = element_text(size = 10),
              legend.text = element_text(size = 9),
              legend.key.size = unit(0.7, 'lines'), 
              legend.justification = c("left", "top"))
    myleg <- ggpubr::get_legend(pp)
    leg <- ggpubr::as_ggplot(myleg) +
        theme(legend.background = element_rect(fill = "transparent"))
    return(leg)
}

# load data -------
load(file = "./final_data/Final_metadata.RData")

# get colours --------
catlist_pal <- lapply(catlist, get_palette_variable, biom = biom)
names(catlist_pal) <- catlist

# edit the major habitat: 
majorhab_pal = catlist_pal[['majorhab']] 
majorhab_levels = names(majorhab_pal)[c(1:4,6:7,5)]
names(majorhab_pal) = majorhab_levels
catlist_pal[['majorhab']] = majorhab_pal

save(catlist_pal, file = paste0(outdir, "catlist_pal.RData"))

# generate plots for reference usages --------
for (ii in catlist) {
    leg = plot_palette_variable(catlist_pal[[ii]], ii)
    ggsave(filename = paste0(plotdir, ii, "_palette.pdf"), plot = leg, width = 4, height = 4)
}

# specific palettes 
load(file = "./derive_data/step4_beta_diver/util/catlist_pal.RData")
majorhab_pal2 = catlist_pal[['majorhab']][2:6]
leg = plot_palette_variable(majorhab_pal2, "majorhab") 
ggsave(filename = paste0(plotdir, "majorhab_palette2.pdf"), plot = leg, width = 3, height = 2)

majorhab_pal3 = catlist_pal[['majorhab']][1:6]
leg = plot_majorhab_variable3(majorhab_pal3, "majorhab") 
ggsave(filename = paste0(plotdir, "majorhab_palette3.pdf"), plot = leg, width = 2, height = 1.5)

majorhab_pal4 = catlist_pal[['majorhab']][3:6]
leg = plot_palette_variable(majorhab_pal4, "majorhab") 
ggsave(filename = paste0(plotdir, "majorhab_palette4.pdf"), plot = leg, width = 2, height = 1.5)

transect_pal2 = catlist_pal[['transect']][3:5]
leg = plot_palette_variable(transect_pal2, "transect")
ggsave(filename = paste0(plotdir, "transect_palette2.pdf"), plot = leg, width = 3, height = 2)

NLCD_pal2 = catlist_pal[['NLCD']][-2]
leg = plot_palette_variable(NLCD_pal2, "NLCD")
ggsave(filename = paste0(plotdir, "NLCD_palette2.pdf"), plot = leg, width = 4, height = 4)

SoS_pal2 = catlist_pal[['SoS']][-1]
leg = plot_palette_variable(SoS_pal2, "SoS")
ggsave(filename = paste0(plotdir, "SoS_palette2.pdf"), plot = leg, width = 3, height = 2)

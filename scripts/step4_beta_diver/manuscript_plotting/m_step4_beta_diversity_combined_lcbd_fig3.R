# Title: Basic beta diversity: LCBD --------
# Author: Meixi Lin
# Date: Tue Jun 26 15:12:01 2018
# Author: Meixi Lin
# Date: Sat Apr 18 15:37:42 2020

# preparation --------
rm(list = ls())
cat("\014")

setwd("~/Lab/caledna_transect")
library(tibble)
library(dplyr)
library(ggplot2)
library(ggedit)
library(phyloseq)
library(microbiomeSeq)
library(reshape2)
library(vegan) # microbiomeSeq calls a function "ordiellipse"
library(adespatial) # microbiomeSeq calls a function "beta.div"
library(grid) # microbiomeSeq calls a function "textGrob"
library(gridExtra) # microbiomeSeq calls a function "theme_minimal"
library(pals)
source("./scripts/function_transect.R")

date() # the execution date

# define functions --------
get_txlevel <- function(myprimer) {
    if (myprimer %in% c("FITS", "PITS")) {
        txlevel <- "Class"
    } else {
        txlevel <- "Phylum"
    }
    return(txlevel)
}

# get taxa names for palette generation
get_taxa_names_pal <- function(pplist) {
    # get taxonomy and link to colours 
    taxos1 <- unique(unlist(lapply(pplist, function(xx) {levels(xx$data$Taxa)})))
    taxos <- c(taxos1[match(c("Others", "unknown"), taxos1)],
               taxos1[-match(c("Others", "unknown"), taxos1)])
    return(taxos)
}

# edit plots, with an input palette
edit_lcbd_plot <- function(xx, taxo_pal, txlevel) {
    # get color palette 
    mytaxo <- as.character(unique(xx$data$Taxa))
    mytaxo_pal <- taxo_pal[match(mytaxo, names(taxo_pal))]
    xx1 <- ggedit::remove_geom(xx, geom = 'point')
    # change value labels
    xx1$data$Groups <- plyr::mapvalues(xx1$data$Groups, 
                                       from = c("Aquatic", "Herbaceous-Dominated", "Shrub-Dominated", "Tree-Dominated"), 
                                       to = c("Aquatic", "Herbaceous", "Shrub", "Tree"))
    xx1 <- xx1 + 
        scale_fill_manual(values = mytaxo_pal) +
        scale_y_continuous(labels = scales::percent) + 
        labs(fill = txlevel) +
        theme_bw() +
        theme(text = element_text(size = 10),
              title = element_text(size = 8), 
              legend.title = element_text(size = 10),
              legend.text = element_text(size = 8),
              legend.key.size = unit(0.5, 'lines'), 
              legend.justification = "left",
              strip.background = element_rect(colour = "black", fill = "white"),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title = element_blank()) 
    return(xx1)
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

# define variables --------
mycutoff <- 4
# (by major habitat)
groupv = "majorhab"
majorlist <- c("Aquatic", "Herbaceous-Dominated", "Shrub-Dominated", "Tree-Dominated")
plotdir = "./plots/step4_beta_diver/manuscript/"
# load --------
load(file = "./derive_data/phy_rare/phyrare.RData")

# start a list to store the lcbd plot --------
pplist <- vector(length = length(primers_commeco), mode = "list")
names(pplist) <- primers_commeco

# plot along each primers --------
for (ii in 1:length(primers_commeco)) {
    primer <- primers_commeco[ii]
    physeq1 <- phyrare[[ii]]
    txlevel <- get_txlevel(primer)
    
    # beta diversity: 
    phyglom <- glom_tax(phyloseq1 = physeq1, taxlevel = txlevel)
    # remove sites with zero occurrences 
    phyglom <-  phyloseq::prune_taxa(taxa_sums(phyglom) > 0, phyglom)
    phyglom <- phyloseq::subset_samples(phyglom, majorhab %in% majorlist)
    # normalize the data 
    mpg0 <- normalise_data(t(phyglom), norm.method = "relative") # calculate relative density
    # starts plotting 
    pp <- microbiomeSeq::plot_taxa(mpg0, grouping_column = groupv, number.taxa = 10) + 
        labs(title = paste(primer, txlevel, sep = "|"))
    pplist[[ii]] <- pp
}

phylaplots = pplist[1:3]
classplots = pplist[4:5]

# now do modifications to the plot, start with phyla --------
# for the taxonomy that used phylum level classifications, define color scheme
phyla_names = get_taxa_names_pal(phylaplots)
phyla_colors = c(pals::polychrome(n = 2), pals::alphabet(n = 21))
phyla_pal = phyla_colors
names(phyla_pal) = phyla_names
pals::pal.bands(phyla_pal)

# modify plots 
newphylapp = lapply(phylaplots, edit_lcbd_plot, taxo_pal = phyla_pal, txlevel = "Phylum")

# now do modifications to the class plot --------
# define colors
class_names = get_taxa_names_pal(classplots)
class_colors = c(pals::polychrome(n = 2), pals::alphabet2(n = 17))
class_pal = class_colors
names(class_pal) = class_names
pals::pal.bands(class_pal)

# modify plots 
newclasspp = lapply(classplots, edit_lcbd_plot, taxo_pal = class_pal, txlevel = "Class")

# arrange the five plots together --------
allpp = c(newphylapp, newclasspp)
leglist = lapply(allpp, get_leg)
nolegpp = lapply(allpp, final_mod)
finalpp = list(nolegpp[[1]], leglist[[1]], nolegpp[[2]], leglist[[2]], 
               nolegpp[[3]], leglist[[3]], nolegpp[[4]], leglist[[4]],
               nolegpp[[5]], leglist[[5]])
gridpp = grid.arrange(grobs = finalpp, ncol = 2, nrow = 5, widths = c(8,2))

ggsave(filename = "beta_diver_lcbd_combined.pdf", plot = gridpp, path = plotdir, height = 12, width = 8, device = "pdf")

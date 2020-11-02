# Title: PCoA plots for minor habitat within major habitat --------
# Author: Meixi Lin
# Date: Wed Sep 18 11:05:55 2019

# preparation --------
options(echo = TRUE)
rm(list = ls())
cat("\014")

setwd("~/Lab/caledna_transect/")

outdir = "./derive_data/step4_beta_diver/"
indir = outdir 
plotdir = "./plots/step4_beta_diver/manuscript/"

library(dplyr)
library(vegan)
library(phyloseq)
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)
date() # the execution date

source("./scripts/function_transect.R")
source("./scripts/step4_beta_diver/0_util_query_var_color.R")

# define function -------
plot_minor_pal <- function(biom, mycolours) {
    majorlist = c("Aquatic", "Herbaceous-Dominated", "Shrub-Dominated", "Tree-Dominated")
    leglist <- lapply(majorlist, function(mymajor, minor_pal = mycolours) {
        xx <- sort(unique(as.character(biom[biom$majorhab == mymajor, 'minorhab'])), na.last = T)
        xx[is.na(xx)] <- "unknown"
        # get color palette 
        myminor_pal <- minor_pal[match(xx, names(minor_pal))]
        value = 1:length(xx)
        dummy = data.frame(xx, value) 
        dummy$xx = factor(dummy$xx, levels = dummy$xx)
        pp <- ggplot(data = dummy, aes(x = xx, y = value, fill = xx)) +
            geom_bar(stat = "identity") + 
            scale_fill_manual(values = myminor_pal) + 
            labs(fill = mymajor) + 
            theme_bw() + 
            theme(legend.title = element_text(size = 9),
                  legend.text = element_text(size = 8),
                  legend.key.size = unit(0.7, 'lines'), 
                  legend.justification = c("left", "top"))
        myleg <- ggpubr::get_legend(pp)
        leg <- ggpubr::as_ggplot(myleg) +
            theme(legend.background = element_rect(fill = "transparent"))
        return(leg)
    })
    return(leglist)
}

# load data -------
load(file = "./derive_data/phy_rare/phyrare.RData")
load(file = paste0(indir, "jadiss.RData"))
load(file = "./final_data/Final_metadata.RData")

# define colour scheme for minorhab so same across figures --------
mycolours = query_var_color(groupv = 'minorhab')

# make a plot for any possible minor habitat within majorhabitat
leglist = plot_minor_pal(biom, mycolours)
legpp <- grid.arrange(grobs = leglist, 
                      nrow = 2, ncol = 2,
                      widths = c(5,4))
ggsave(filename = paste0(plotdir,"legend_beta_diver_minhab_mahab_pcoa.pdf"), plot = legpp, height = 5, width = 8)

# for each primers --------
for (ii in 1:length(primers_commeco)) {
    # get phyloseq objects 
    primer <- primers_commeco[ii]
    physeq1 <- phyrare[[primer]]
    mydiss <- jadiss[[primer]]
    
    # subset phyloseq by major habitat 
    aqua <- subset_samples(physeq1, majorhab == "Aquatic")
    herb <- subset_samples(physeq1, majorhab == "Herbaceous-Dominated")
    shru <- subset_samples(physeq1, majorhab == "Shrub-Dominated")
    tree <- subset_samples(physeq1, majorhab == "Tree-Dominated")
    
    # perform PCoA on all of them
    phy.majorhab = list(aqua, herb, shru, tree)
    
    # perform plotting 
    pplist <- lapply(phy.majorhab, function(xx) {
        # get majorhabitat type 
        mymajorhab = as.character(unlist(sample_data(xx)[1,'majorhab']))
        # start ordination 
        ord.res <- phyloseq::ordinate(xx, method = "PCoA", distance = "jaccard", binary = T)
        forplot <- plot_ordination(xx, ord.res, type = "samples", color = "minorhab", title = sample_data(xx)$majorhab[1])
        
        # make small changes to the plotting scheme 
        forplotdf <- forplot$data %>% 
            dplyr::select(Axis.1, Axis.2, MatchName, minorhab) 
        forplotdf[,'minorhab'] <- as.character(forplotdf[,'minorhab'])
        forplotdf[is.na(forplotdf[,'minorhab']),'minorhab'] <- "unknown"
        
        hull_forplot <- forplotdf %>%
            dplyr::group_by(minorhab) %>%
            dplyr::slice(chull(Axis.1, Axis.2))
        
        pp <- ggplot(forplotdf, aes(x = Axis.1, y = Axis.2, color = minorhab)) + 
            geom_point() + 
            scale_colour_manual(values = mycolours) +
            theme_bw() +
            labs(x = forplot$labels$x, y = forplot$labels$y,
                 title = mymajorhab) +
            theme(legend.position = "none")
        pp <- pp + 
            aes(fill = minorhab) + 
            geom_polygon(data = hull_forplot, alpha = 0.1) +
            scale_fill_manual(values = mycolours, guide = F) +
            theme(legend.position = "none")
        return(pp)
    })
    
    # arrange within minor habitats 
    allpp <- grid.arrange(grobs = pplist, 
                          nrow = 2, ncol = 2,
                          top = textGrob(primer,gp=gpar(fontsize=20,font=3)))
    ggsave(filename = paste0(plotdir, primer, "_beta_diver_minhab_mahab_pcoa.pdf"), plot = allpp, height = 8, width = 8)
}

# ending --------
date()
closeAllConnections()


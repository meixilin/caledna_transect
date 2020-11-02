# Title: alpha diversity and category variables plotting for figure 2 --------
# Author: Meixi Lin
# Date: Tue May  7 12:34:27 2019
# Author:
# Date:
# Modification:

# preparation --------
rm(list = ls())
cat("\014")
setwd("~/Lab/caledna_transect")

library(vegan)
library(tibble)
library(dplyr)
library(ggplot2)
library(phyloseq)
library(ggpubr)

source("./scripts/function_transect.R")
source("./scripts/step4_beta_diver/0_util_query_var_color.R")
date() # the execution date

# define functions --------
# get fortest
get_fortest <- function(myprimer, myvar, mymeasure) {
    psdata <- phyrare[[myprimer]]
    adivvar <- adivrare[[myprimer]][mymeasure] %>%
        tibble::rownames_to_column(var = "MatchName")
    category <- get_variable(psdata, varName = c("MatchName", myvar))
    print(sort(table(category[,myvar])))
    
    colnames(category) = c("MatchName", "category")
    sumby <- table(category[,"category"]) %>% as.data.frame() %>% 
        filter(Freq > mycutoff)
    # if all data are scatterred, 0 observation for sumby  
    if (nrow(sumby) < 3) {
        print(nrow(sumby))
        print(paste("for", myvar, "no data for stats test."))
        next
    }
    fortest <- merge(adivvar, category, by = "MatchName") %>%
        dplyr::filter(category %in% sumby$Var1)
    fortest$category <- droplevels(fortest$category)
    return(fortest)
}

# get post hoc dunn test 
get_dunn_test <- function(myprimer, myvar, mymeasure, dunndir) {
    mydunn = list.files(path = dunndir, pattern = paste0(paste(myprimer, myvar, mymeasure, sep = "_"), "_dunn"))
    dunnfile = loadRData(paste0(dunndir, mydunn))$res %>% 
        dplyr::filter(P.adj < mysig) 
    ypos = seq(mymax, by = (mymax - mymin)/20, length.out = nrow(dunnfile))
    dunnfile = dunnfile %>% 
        dplyr::mutate(term = myvar,
                      method = "Dunn-test",
                      p.signif = dplyr::case_when(P.adj < 0.001 ~ "***",
                                                  P.adj < 0.01 & P.adj >= 0.001 ~ "**",
                                                  P.adj < 0.05 & P.adj >= 0.01 ~ "*")) %>%
        tidyr::separate(col = Comparison, into = c("group1", "group2"), sep = " - ") %>%
        dplyr::arrange(group1, group2) %>%
        dplyr::mutate(y.position = ypos) %>%
        dplyr::rename(p.adj = P.adj) %>% 
        dplyr::select(term, group1, group2, p.adj, p.signif, method, y.position) %>%
        tidyr::as_tibble()
    return(dunnfile)
} 

# define variables --------
plotdir <- "./plots/step3_alpha_diver/boxplot/"
dunndir = "./derive_data/step3_alpha_diver/kruskalres/"
dir.create(plotdir, recursive = T)

# import phyloseq object --------
load("./derive_data/phy_rare/phyrare.RData")
load("./derive_data/step3_alpha_diver/adivrare_06012019.RData")

mymeasure = "Shannon"
myprimer = "FITS"
mysig <- 0.05 
mycutoff <- 4

# get min max --------
minmax <- adivrare[[myprimer]][mymeasure] %>%
    tibble::rownames_to_column(var = "MatchName")
mymin = min(minmax[, mymeasure])
mymax = max(minmax[, mymeasure])

# for transect --------
transect <- get_fortest(myprimer, myvar = "transect", mymeasure)
pp1 <- ggplot(data = transect, 
              aes(x = reorder(category, Shannon, FUN = median), y = Shannon, colour = category)) + 
    geom_boxplot() +
    geom_jitter(size = 0.8, alpha = 0.5, shape = 1) + 
    ggpubr::stat_compare_means(method = "kruskal.test", label.x.npc = "center", label.y.npc = "bottom") + 
    scale_x_discrete(labels=c("Coastal" = "Coast", "Forest" = "Forest", "ShrubScrub" = "Shrub")) +
    stat_pvalue_manual(data = get_dunn_test(myprimer, myvar = "transect", mymeasure, dunndir), label = "p.signif")  +
    scale_colour_manual(values = query_var_color(groupv = "transect")) +
    theme_bw() + 
    theme(text = element_text(size = 16), 
          legend.position = "none") +
    labs(x = "Transect",y = "Shannon Index")
ggsave(filename = paste0(plotdir, "fig2_FITS_transect.pdf"), 
       plot = pp1, width = 4, height = 4, device = "pdf")

# for substrate --------
substrate <- get_fortest(myprimer, myvar = "SoS", mymeasure)
pp2 <- ggplot(data = substrate, 
              aes(x = reorder(category, Shannon, FUN = median), y = Shannon, colour = category)) + 
    geom_boxplot() +
    geom_jitter(size = 0.8, alpha = 0.5, shape = 1) + 
    ggpubr::stat_compare_means(method = "kruskal.test", label.x.npc = "center", label.y.npc = "bottom") + 
    scale_x_discrete(labels=c("Coastal" = "Coast", "Forest" = "Forest", "ShrubScrub" = "Shrub")) +
    stat_pvalue_manual(data = get_dunn_test(myprimer, myvar = "SoS", mymeasure, dunndir), label = "p.signif")  +
    scale_colour_manual(values = query_var_color(groupv = "SoS")) +
    theme_bw() + 
    theme(text = element_text(size = 16), 
          legend.position = "none") +
    labs(x = "Substrate Type",y = "Shannon Index")
ggsave(filename = paste0(plotdir, "fig2_FITS_substrate.pdf"), 
       plot = pp2, width = 4, height = 4, device = "pdf")

# for majorhab --------
majorhab <- get_fortest(myprimer, myvar = "majorhab", mymeasure)
pp3 <- ggplot(data = majorhab, 
              aes(x = category, y = Shannon, colour = category)) + 
    geom_boxplot() +
    geom_jitter(size = 0.8, alpha = 0.5, shape = 1) + 
    ggpubr::stat_compare_means(method = "kruskal.test", label.x.npc = "center", label.y.npc = "bottom") + 
    scale_x_discrete(labels=c("Herbaceous-Dominated" = "Herbaceous", 
                              "Shrub-Dominated" = "Shrub", 
                              "Tree-Dominated" = "Tree")) +
    # stat_pvalue_manual(data = get_dunn_test(myprimer, myvar = "majorhab", mymeasure, dunndir), label = "p.signif")  +
    scale_colour_manual(values = query_var_color(groupv = "majorhab")) +
    theme_bw() + 
    ylim(mymin,mymax) + 
    theme(text = element_text(size = 16), 
          legend.position = "none") +
    labs(x = "Major Habitat",y = "Shannon Index")
ggsave(filename = paste0(plotdir, "fig2_FITS_majorhab.pdf"), 
       plot = pp3, width = 4, height = 4, device = "pdf")

# for location --------
loc <- get_fortest(myprimer, myvar = "loc", mymeasure)
pp4 <- ggplot(data = loc, 
              aes(x = reorder(category, Shannon, FUN = median), y = Shannon)) + 
    geom_boxplot() +
    geom_jitter(size = 0.8, alpha = 0.5, shape = 1) + 
    ggpubr::stat_compare_means(method = "kruskal.test", label.x.npc = "center", label.y.npc = "bottom") + 
    scale_x_discrete(labels=c("Coffee_Creek_Trinity_County" = "Coffee_Creek",
                              "San_Joaquin_Marsh_Reserve" = "San_Joaquin_Marsh",
                              "Mill_Creek_250_m_upstream_from_Kings_River" = "Mill_Creek",
                              "James_San_Jacinto_Mountains_Reserve" = "San_Jacinto_Mountains",
                              "North_Table_Mountain_Ecological_Reserve" = "North_Table_Mountain",
                              "Milpitas_Special_Interest_Area_Los_Padres_National_Forest" = "Milpitas_Los_Padres_Forest",
                              "Sweeney_Granite_Mountains_Desert_Research_Center" = "Sweeney_Granite_Mountains")) +
    stat_pvalue_manual(data = get_dunn_test(myprimer, myvar = "loc", mymeasure, dunndir), label = "p.signif")  +
    theme_bw() + 
    # ylim(mymin,mymax) + 
    theme(text = element_text(size = 16), 
          legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Location",y = "Shannon Index")

ggsave(filename = paste0(plotdir, "fig2_FITS_location.pdf"), 
       plot = pp4, width = 12, height = 7, device = "pdf")

# cleaning up --------
closeAllConnections()

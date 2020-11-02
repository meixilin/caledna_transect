# Title: alpha diversity and category variables plotting --------
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
library(reshape2)
library(dplyr)
library(ggplot2)
library(phyloseq)
library(ggpubr)

source("./scripts/function_transect.R")
source("./scripts/step4_beta_diver/0_util_query_var_color.R")
date() # the execution date

# import phyloseq object --------
load("./derive_data/phy_rare/phyrare.RData")
load("./derive_data/step3_alpha_diver/adivrare_06012019.RData")

# define variables 
plotdir <- "./plots/step3_alpha_diver/boxplot/"
mymeasures <- c("Observed", "Shannon")
mysig <- 0.05 
mycutoff <- 4

# define functions --------
# get fortest for each primer, all possible measures
get_fortest <- function(myprimer, myvar, mymeasures, mycutoff) {
    psdata <- phyrare[[myprimer]]
    adivvar <- adivrare[[myprimer]] %>%
        tibble::rownames_to_column(var = "MatchName")
    category <- get_variable(psdata, varName = c("MatchName", myvar))
    print(sort(table(category[,myvar])))
    
    colnames(category) = c("MatchName", "Category")
    sumby <- table(category[,"Category"]) %>% as.data.frame() %>% 
        filter(Freq > mycutoff)
    # if all data are scatterred, 0 observation for sumby  
    if (nrow(sumby) < 3) {
        print(nrow(sumby))
        print(paste("for", myvar, "no data for stats test."))
        next
    }
    # get the category and alpha diversity value 
    fortest <- base::merge(adivvar, category, by = "MatchName") %>%
        dplyr::filter(Category %in% sumby$Var1) 
    fortest$Category <- droplevels(fortest$Category)
    
    # melt the data set 
    fortest <- fortest %>%
        reshape2::melt(id.vars = c("MatchName", "Category"), measure.vars = mymeasures, value.name = "adiv.value") %>%
        dplyr::mutate(Metabarcode = myprimer) %>%
        dplyr::mutate_if(is.factor, as.character)
    
    colnames(fortest) = c("MatchName", "Category", "Measure", "adiv.value", "Metabarcode")
    return(fortest)
}

get_fortest_all_primer <- function(myprimers, myvar, mymeasures, mycutoff) {
    # automatically pass the myvar, mymeasures and mycutoff to the get_fortest function
    fortest.list <- lapply(myprimers, get_fortest, myvar = myvar, mymeasures = mymeasures, mycutoff = mycutoff)
    fortestall <- dplyr::bind_rows(fortest.list)
}

# main: perform plotting for all primers used and the variables of interest  --------
# for location ########
loc <- get_fortest_all_primer(myprimers = primers_commeco, myvar = "loc", mymeasures = mymeasures, mycutoff = mycutoff)
pp1 <- ggplot(data = loc, 
             aes(x = Category, y = adiv.value, color = Category)) + 
    facet_wrap(Metabarcode ~ Measure, scales = "free_y", ncol = 2) + 
    geom_boxplot() +
    geom_jitter(size = 0.8, alpha = 0.5, shape = 1) + 
    ggpubr::stat_compare_means(method = "kruskal.test", label.x.npc = "middle", label.y.npc = "bottom", size = 3) +
    labs(x = "Location",
         y = "Alpha diversity") +
    scale_colour_manual(values = query_var_color(groupv = "loc")) +
    theme_bw() + 
    theme(text = element_text(size = 10),
          axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "none")

ggsave(filename = paste0(plotdir, "supp_fig_alphadiv_location.pdf"),
       plot = pp1, width = 8, height = 12, device = "pdf")

# cleaning up --------
closeAllConnections()

# Title: PCoA plots for each category and primers_commeco --------
# Author: Meixi Lin
# Date: Mon Aug  5 19:37:31 2019
# Author:
# Date:
# Modification:

# preparation --------
rm(list = ls())
cat("\014")

setwd("~/Lab/caledna_transect/")
library(dplyr)
library(vegan)
library(phyloseq)
library(ggplot2)
date() # the execution date

source("./scripts/function_transect.R")
source("./scripts/step4_beta_diver/0_util_query_var_color.R")

# set some variable --------
outdir = "./derive_data/step4_beta_diver/"
sink(paste0(outdir, "logs/pcoa_plot_c4_", Sys.Date(), ".log"))
plotdir1 = "./plots/step4_beta_diver/"
myseed <- 17
mycutoff <- 4 # min items in one category 

# define functions --------
# get plotting df: # delete the sites that does not have at least 5 elements in that category
get_forplotdf <- function(forplot, physeq1, mycutoff, groupv) {
    # get the plotting df =========
    forplotdf <- forplot$data %>% 
        dplyr::select(Axis.1, Axis.2, MatchName, eval(groupv)) 
    forplotdf[,groupv] <- as.character(forplotdf[,groupv])
    # delete the sites that does not have at least 5 elements in that category
    category <- get_variable(physeq1, varName = groupv)
    sumby <- table(category) %>% as.data.frame() %>% 
        filter(Freq <= mycutoff) # here are the variables to change 
    
    otherid <- as.character(sampledf[(sampledf[,groupv] %in% sumby$category), 'MatchName'])
    forplotdf[otherid,groupv] <- "Others"
    forplotdf[is.na(forplotdf[,groupv]),groupv] <- "unknown"
    forplotdf[,groupv] <- as.factor(forplotdf[,groupv])
    return(forplotdf)
}

# custom plotting 
myplot_ordination <- function(forplot, forplotdf, mycolours, groupv, primer){
    plottitle = paste(primer, groupv, sep = "|")
    print(plottitle)
    # start plotting =========
    hull_forplot <- forplotdf %>%
        group_by(eval(as.name(groupv))) %>%
        slice(chull(Axis.1, Axis.2))
    
    pp <- ggplot(forplotdf, aes(x = Axis.1, y = Axis.2, color = eval(as.name(groupv)))) + 
        geom_point() + 
        scale_colour_manual(values = mycolours) +
        theme_bw() +
        labs(x = forplot$labels$x, y = forplot$labels$y,
             colour = groupv, title = plottitle)

    pp <- pp + 
        aes(fill = eval(as.name(groupv))) + 
        geom_polygon(data = hull_forplot, alpha = 0.1) +
        scale_fill_manual(values = mycolours, guide = F) 
    ggsave(paste0("beta_diver_pcoa_", groupv, ".pdf"), plot = pp, path = plotdir, width = 10, height = 6) # could be clear if needed 
    pp1 <- pp + 
        theme(legend.position = "none")
    ggsave(paste0("beta_diver_pcoa_", groupv, "_noleg.pdf"), plot = pp1, path = plotdir, width = 6, height = 6)
    return(pp)
}

# load data -------
load(file = "./derive_data/phy_rare/phyrare.RData")
load(file = "./derive_data/step4_beta_diver/jadiss.RData")
load(file = "./final_data/Final_metadata.RData")

# start a list to store the ordination result 
pcoalist <- vector(length = length(primers_commeco), mode = "list")
pplist <- vector(length = length(primers_commeco) * length(catlist), mode = "list")

# main --------
for (ii in 1:length(primers_commeco)) {
    # get phyloseq objects 
    primer <- primers_commeco[ii]
    physeq1 <- phyrare[[primer]]
    mydiss <- jadiss[[primer]]
    plotdir <- paste0(plotdir1, primer, "/pcoa/")
    dir.create(plotdir, recursive = T)
    
    set.seed(myseed)
    ord.res <- phyloseq::ordinate(physeq1, method = "PCoA", distance = "jaccard", binary = T)
    print(ord.res)
    
    pcoalist[[ii]] <- ord.res
    names(pcoalist)[ii] <- primer
    forplot <- plot_ordination(physeq1, ord.res, type = "samples") 
    sampledf <- data.frame(sample_data(physeq1))
    
    for (jj in 1:length(catlist)) {
        groupv = catlist[jj]
        forplotdf = get_forplotdf(forplot, physeq1, mycutoff, groupv)
        mycolours = query_var_color(groupv = groupv)
        pp = myplot_ordination(forplot, forplotdf, mycolours, groupv, primer)
        id = (ii -1) * length(catlist) + jj; print(id)
        pplist[[id]] <- pp
    }
}

# save the ordination plots 
pcoadir = paste0(outdir, "pcoa/")
dir.create(pcoadir, recursive = T)
save(pcoalist, file = paste0(pcoadir, "pcoa_result.RData"))
save(pplist, file = paste0(pcoadir, "pcoa_category_plot.RData"))

# end ========
sink()
closeAllConnections()

# note for PCoA output ========
# for all 6 primers_commeco, there were negative eigenvalues in PITS but not the others. No corrections were performed for PITS primer. 
# 
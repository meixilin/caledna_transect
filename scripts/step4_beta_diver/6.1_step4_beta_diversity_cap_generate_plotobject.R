# Title: cap plot controlling for the first order of geographical location --------
# Author: Meixi Lin
# Date: Fri Oct 11 15:41:14 2019
# Author:
# Date:
# Modification:

# preparation --------
rm(list = ls())
cat("\014")

setwd("~/Lab/caledna_transect/")

options(echo = T)
outdir = "./derive_data/step4_beta_diver/"
plotdir1 = "./plots/step4_beta_diver/"
sink(paste0(outdir, "logs/cap_plot_c4_latlong", Sys.Date(), ".log"), append = F)
library(dplyr)
library(vegan)
library(phyloseq)
library(ggplot2)
date() # the execution date

source("./scripts/function_transect.R")
source("./scripts/step4_beta_diver/0_util_query_var_color.R")

# set some variable --------
myseed <- 17
mycutoff <- 4

# define function --------
myplot_capscale <- function(forplot, mycolour, groupv, primer){
    plottitle = paste(primer, groupv, sep = "|")
    forplotdf <- forplot$data %>% 
        dplyr::select(CAP1, CAP2, MatchName, eval(groupv)) 
    
    hull_forplot <- forplotdf %>%
        group_by(eval(as.name(groupv))) %>%
        slice(chull(CAP1, CAP2))
    
    ct <- data.frame(sm.mycap$centroids) %>%
        tibble::rownames_to_column(var = "category") %>%
        dplyr::mutate(., category = 
                          sapply(strsplit(category, split = "eval(as.name(groupv))", fixed=TRUE),function(x) (x[2]))) # get rid of the name 
    
    pp <- ggplot() + 
        geom_point(data = forplotdf, aes(x = CAP1, y = CAP2, color = eval(as.name(groupv)))) + 
        scale_colour_manual(values = mycolour) +
        theme_bw() +
        labs(x = forplot$labels$x, y = forplot$labels$y,
             colour = groupv, title = paste0(primer, "|", groupv))
    pp <- pp +
        geom_polygon(data = hull_forplot, aes(x = CAP1, y = CAP2, fill = eval(as.name(groupv))), alpha = 0.1) +
        scale_fill_manual(values = mycolour, guide = F) 
    # show the centriods 
    pp1 <- pp +
        geom_text(data = ct, aes(x = CAP1, y = CAP2, colour = category), label = "X", size = 5, show.legend = F)
    return(pp1)
}

# load data -------
load(file = "./derive_data/phy_rare/phyrare.RData")
load(file = "./derive_data/step4_beta_diver/jadiss.RData")

# our category list of interest is shorter, only habitat and soil property 
catlist = catlist[c(2,3,4,5,6,8,9)]
print(catlist)

# start a list to store the ordination result 
caplist <- vector(length = length(primers_commeco) * length(catlist), mode = "list")
varparlist <- vector(length = length(primers_commeco) * length(catlist), mode = "list")
pplist <- vector(length = length(primers_commeco) * length(catlist), mode = "list")

# for each primers_commeco --------
for (ii in 1:length(primers_commeco)) {
    # get phyloseq objects --------
    primer <- primers_commeco[ii]
    physeq1 <- phyrare[[primer]]
    mydiss <- jadiss[[primer]]
    plotdir <- paste0(plotdir1, primer, "/cap_latlong/")
    dir.create(plotdir, recursive = T)
    dir.create(plotdir, recursive = T)
    sampledf <- data.frame(sample_data(physeq1))
    # # convert to euclidean, not affecting cap result 
    # biomxy <- SoDA::geoXY(sampledf[,"Latitude"], sampledf[,"Longitude"]) %>% 
    #     as.data.frame() 
    
    # now starts capscale ---------
    for (jj in 1:length(catlist)) {
        groupv <- catlist[jj]
        print(paste(primer, groupv, sep = "|"))
        # get the variables 
        capid <- (ii - 1) * length(catlist) + jj
        # delete the sites that does not have at least 5 elements in that category
        category <- get_variable(physeq1, varName = groupv)
        sumby <- table(category) %>% as.data.frame() %>% 
            filter(Freq > mycutoff)
        samid <- as.character(sampledf[(sampledf[,groupv] %in% sumby$category), 'MatchName'])
        diss1 <- usedist::dist_subset(mydiss, samid)
        attr(diss1, "method") = attr(mydiss, "method") 
        sampledf1 <- sampledf[samid, ]
        # subset physeq 
        physeq2 <- prune_samples(samid, x = physeq1)
        sample_data(physeq2)[,groupv][[1]] <- droplevels(sample_data(physeq2)[,groupv][[1]])
        
        # start the cap scale functions ========
        print("Perform default CAP.")
        mycap0 <- capscale(formula = diss1 ~ eval(as.name(groupv)), data = sampledf1)
        print(mycap0)
        print(RsquareAdj(mycap0))
        print(anova.cca(mycap0, by = "terms"))
        print("Perform capscale with location effect conditioned out.")
        mycap <- capscale(formula = diss1 ~ eval(as.name(groupv)) + Condition(Longitude + Latitude), data = sampledf1)
        print(mycap)
        (sm.mycap <- summary(mycap))
        print(RsquareAdj(mycap))
        print(anova.cca(mycap, by = "terms"))
        print("Perform variance partitioning test.")
        myvarpar <- varpart(diss1, ~ eval(as.name(groupv)), ~ (Longitude + Latitude), data = sampledf1)
        print(myvarpar)
        # save the output 
        caplist[[capid]] <- mycap 
        names(caplist)[capid] <- paste(primer, groupv, sep = "|")
        varparlist[[capid]] <- myvarpar 
        names(varparlist)[capid] <- paste(primer, groupv, sep = "|")
        # get colours 
        mycolour <- unname(pals::polychrome(n = length(levels(sample_data(physeq2)[,groupv][[1]])) + 2))
        mycolour <- mycolour[3:length(mycolour)]
        # start plotting =========
        forplot <- plot_ordination(physeq2, mycap, type = "samples", color = groupv) 
        pp1 <- myplot_capscale(forplot, mycolour, groupv, primer)
        ggsave(paste0("beta_diver_cap_", groupv, "_text.pdf"), plot = pp1, path = plotdir, width = 10, height = 6)
        pp2 <- pp1 + 
            theme(legend.position = "none")
        ggsave(paste0("beta_diver_cap_", groupv, "_noleg.pdf"), plot = pp2, path = plotdir, width = 6, height = 6)
        # save the plot 
        pplist[[capid]] <- pp1
        names(pplist)[capid] <- paste(primer, groupv, sep = "|")
    }
}

# save the ordination plots 
capdir = paste0(outdir, "cap_latlong/"); dir.create(capdir, recursive = T)
save(caplist, file = paste0(capdir, "cap_result_", Sys.Date(), ".RData"))
save(varparlist, file = paste0(capdir, "varpar_result_", Sys.Date(), ".RData"))
save(pplist, file = paste0(capdir, "cap_category_plot_", Sys.Date(), ".RData"))

# end ========
sink()
closeAllConnections()

# Note, CAP analysis will automatically remove sites with NA as response, so there were no "unknown" sites in the plots 

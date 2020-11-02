# Title: PCoA plots with ordisurf and envfit --------
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

# set some variable --------
myconf <- 0.001
outdir = "./derive_data/step4_beta_diver/"
indir = paste0(outdir, "envfit_ordisurf/")
plotdir = "./plots/step4_beta_diver/envfit_ordisurf/"; dir.create(plotdir, recursive = T)

# load data -------
load(file = "./derive_data/phy_rare/phyrare.RData")
load(file = "./final_data/Final_metadata.RData")
load(file = "./derive_data/step4_beta_diver/jadiss.RData")
load(file = paste0(indir, "envfit_pcoa_2020-03-30.RData"))

# here I will only use major habitat as the category --------
groupv = "majorhab"

# define point type 
# now set point type for each variables 
mypch = data.frame(1:(nlevels(biom[,groupv]) + 1), c(levels(biom[,groupv]), "unknown"))
colnames(mypch) <- c("pch", groupv)
mypch

# define a function for ii th primer and jj th variable ---------
# this only works with function style similar to R, that takes global environment in functions 
# prepare data for ordisurf 
surf_dtprep <- function(ii) {
    # get the baseline data for plotting 
    # get phyloseq objects 
    primer <- primers_commeco[ii]
    physeq <- phyrare[[primer]]
    mydiss <- jadiss[[primer]]
    myenvfit <- eflist[[primer]]
    # envfit does not take in phyloseq generated pcoa results 
    ord.res <- stats::cmdscale(d = mydiss, eig = T)
    # get sample dataframe 
    sampledf <- data.frame(sample_data(physeq)) %>%
        dplyr::select(contlist) 
    
    # get pch types 
    # now set point type for each variables 
    pchs = data.frame(sample_data(physeq)[, groupv]) %>%
        tibble::rownames_to_column(var = "MatchName")
    pchs[,groupv] = as.character(pchs[,groupv])
    pchs[is.na(pchs[, groupv]),groupv] = "unknown"
    # combine pch values 
    pchs = left_join(x = pchs, y = mypch, by = groupv)
    # check if result is the same for ord.res and pchs 
    if (!all(pchs[,'MatchName'] == rownames(ord.res$points))) {
        print("no matching majorhab and points, check!")
        break; 
    }
    to_return = list(sampledf, pchs, ord.res, myenvfit)
    names(to_return) = c("sampledf", "pchs", "ord.res", "myenvfit")
    return(to_return)
}

# plot for each ordination and results 
surf_plot <- function(to_return, primer, myvar, myconf, plotdir = NULL, legend = TRUE) {
    sampledf <- to_return[['sampledf']]
    pchs <- to_return[['pchs']]
    ord.res <- to_return[['ord.res']]
    ef <- to_return[['myenvfit']]
    # par(new = TRUE)
    plot(x = ord.res$points[,1], y = ord.res$points[,2], pch = pchs[,'pch'], main = paste(primer, myvar, sep = "|"), xlab = "Axis.1", ylab = "Axis.2")
    # this part needs manual changes 
    if (legend == TRUE) {
        legend("bottomright", legend=c("Aquatic", "Herbaceous-Dominated", "Shrub-Dominated", "Tree-Dominated", "Non-vegetated","unknown"), pch = c(1,2,4,5,3,6))
    }
    # plot the envfit result 
    plot(ef, p.max = myconf, add = TRUE)
    # add arrows for the variable of interest 
    # length are also scaled 
    myarr <- vegan::scores(ef, display = 'vectors')[myvar,] * ordiArrowMul(ef)
    # check if the result was significant 
    if (ef$vectors$pvals[myvar] <= myconf) {
        # significant 
        mycol = "orange"
    } else {
        # not significant
        mycol = "grey"
    }
    arrows(x0 = 0, y0 = 0, 
           x1 = myarr['Dim1'] , y1 = myarr['Dim2'], 
           col = mycol, lwd = 2)
    mp <- with(sampledf, ordisurf(ord.res, eval(as.name(myvar)), add = TRUE, cex = 2))
    return(0)
}

# main ========
for (jj in 1:length(contlist)) {
    pdf(file = paste0(plotdir, "beta_diver_envfit_surf_", contlist[jj],".pdf"), height = 9, width = 12)
    pp <- par(mfrow = c(2,3))
    for (ii in 1:length(primers_commeco)) {
        leg <- FALSE
        mydt <- surf_dtprep(ii = ii)
        if (ii == length(primers_commeco)) {
            leg <- TRUE
        }
        surf_plot(mydt, primer = primers_commeco[ii], myconf = myconf, myvar = contlist[jj], legend = leg)
    }
    par(pp)
    dev.off()
}


# end ========
closeAllConnections()

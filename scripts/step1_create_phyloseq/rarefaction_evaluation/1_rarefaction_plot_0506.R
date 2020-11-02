# Title: plot the rarefied phyloseq --------
# Author: Meixi Lin
# Date: Tue Jun 26 15:12:01 2018
# Author: Meixi Lin
# Date: Mon Jan 28 10:19:05 2019
# Modification: finish rarefaction plotting with new data

# preparation --------
rm(list = ls())
cat("\014")

setwd("~/Lab/caledna_transect/")
outdir = "./derive_data/step1_create_phyloseq/eval_rarefaction/"
plotdir = "./plots/step1_create_phyloseq/eval_rarefaction/"
dir.create(outdir, recursive = T)
dir.create(plotdir, recursive = T)

output <- file(paste0(outdir, "rarefaction_plot_log_05052019.log"))
sink(output, append = F)

library(tibble)
library(dplyr)
library(ggplot2)
library(ranacapa)
library(phyloseq)
library(microbiomeSeq)

source("./scripts/function_transect.R")

date() # the execution date

# define functions --------
# credit: rachel meyer 
calculate_rarefaction_curves <- function(psdata, measures, depths) {
    require('plyr') # ldply
    require('reshape2') # melt
    
    estimate_rarified_richness <- function(psdata, measures, depth, myseed) {
        if(max(sample_sums(psdata)) < depth) return()
        psdata <- prune_samples(sample_sums(psdata) >= depth, psdata)
        
        rarified_psdata <- phyloseq::rarefy_even_depth(psdata, depth)
        
        alpha_diversity <- phyloseq::estimate_richness(rarified_psdata, measures = measures)
        
        # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
        molten_alpha_diversity <- melt(as.matrix(alpha_diversity), varnames = c('Sample', 'Measure'), value.name = 'Alpha_diversity')
        
        molten_alpha_diversity
    }
    
    names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
    rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, psdata = psdata, measures = measures, .id = 'Depth', .progress = ifelse(interactive(), 'text', 'none'))
    
    # convert Depth from factor to numeric
    rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
    
    rarefaction_curve_data
}

calculate_rare_summary <- function(primer, rare_level = c(1, 2, 40, 100, 1000, 2000, 4000, 6000, 10000, 15000, 20000)) {
    psdata <- phydeco[[primer]]
    sample_sums(psdata)
    summary(sample_sums(psdata))

    rarefaction_curve_data <- calculate_rarefaction_curves(psdata, mymeasures, rep(rare_level, each = 10))
    summary(rarefaction_curve_data)
    save(rarefaction_curve_data, 
         file = paste0(outdir, primer, "_rarefaction_curve_data.RData"))
    
    rarefaction_curve_data_summary <- plyr::ddply(rarefaction_curve_data, c('Depth', 'Sample', 'Measure'), summarise, Alpha_diversity_mean = mean(Alpha_diversity), Alpha_diversity_sd = sd(Alpha_diversity))
    
    rarefaction_curve_data_summary_verbose <- base::merge(rarefaction_curve_data_summary, data.frame(sample_data(psdata)), by.x = 'Sample', by.y = 'row.names')
    save(rarefaction_curve_data_summary_verbose, 
         file = paste0(outdir, primer, "_rarefaction_summary_verbose.RData"))
    return(rarefaction_curve_data_summary_verbose)
}

get_abline <- function(primer) {
    psdata <- phydeco[[primer]]
    median = median(sample_sums(psdata))
    cutoff = rare_depth[primer]
    out = unname(c(median, cutoff))
    return(out)
}

# import phyloseq object -------- 
load("./derive_data/phy_deco/phydeco.RData")

# make a list and name variables -------
rarefaction_reps  <- 10
# myseeds 
# # seeds <- as.integer(runif(rarefaction_reps, min = 0, max = 1000)) 
# seeds <- c(825,598,744,276,705,872,895,901,759,544)
mymeasures <- c("Observed", "Shannon")
catlist <- c("transect")

# plot rarefaction --------
# note that this is for plotting and it was hard to set seed 
for (ii in primers) {
    if (ii == "all") {
        rare_level = c(1, 2, 40, 100, 1000, 2000, 4000, 6000, 10000, 15000, 20000, 30000, 40000, 60000)
        rarefaction_curve_data_summary_verbose = calculate_rare_summary(ii, rare_level)
    }
    rarefaction_curve_data_summary_verbose = calculate_rare_summary(ii)
    
    for (jj in catlist) {
        # get abline 
        pp1 <- ggplot(
            data = rarefaction_curve_data_summary_verbose,
            mapping = aes(
                x = Depth,
                y = Alpha_diversity_mean,
                ymin = Alpha_diversity_mean - Alpha_diversity_sd,
                ymax = Alpha_diversity_mean + Alpha_diversity_sd,
                colour = eval(as.name(jj)), ## CHANGE THIS ACCORDING TO MYCAT 
                group = Sample
            )) + 
            geom_vline(xintercept = get_abline(ii)[1], color = "blue") + 
            geom_vline(xintercept = get_abline(ii)[2], color = "red") + 
            geom_line() + 
            geom_pointrange(fatten = 2) + 
            facet_wrap(
                facets = ~ Measure,
                scales = 'free_y') + 
            labs(colour = jj, title = paste(ii, jj)) + 
            theme_bw()
        
        if (jj %in% c("loc", "clust")) {
            pp1 <- pp1 +
                theme(legend.position = "none")
        }
        
        ggsave(filename = paste0(plotdir, ii, "_", jj, "_rarefaction_plot.pdf"), 
               plot = pp1, device = "pdf",
               width = 8, height = 4)
    }
}

sink()
dev.off()
closeAllConnections()

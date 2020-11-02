# Title: data distribution in raster and variables --------
# Author: Meixi Lin
# Date: Wed May 15 20:49:18 2019
# Author:
# Date:
# Modification: parallel and hoffman version

# preparation --------
rm(list = ls())
cat("\014")

# setwd("/Users/linmeixi/UCLA/Lab/abiotic_transect")
setwd("/u/project/rwayne/meixilin/abiotic_transect")
library(dplyr)
library(tibble)
library(ggplot2)
library(ggcorrplot)
library(rgeos)
library(raster)
library(gstat)
library(rgdal)
library(doParallel)
library(foreach)
date() # the execution date

source("./r_codes/function_beta_0725.R")

# register cluster --------
# parallel::detectCores(all.tests = FALSE, logical = TRUE)
cl <- makeCluster(2)
registerDoParallel(cl)

# read data ---------
load("./final_data/Final_metadata.RData")
rasterdir <- read.csv(file = "./final_data/metadata/new_raster_dir.csv", stringsAsFactors = F)
myvar <- c("Longitude", "Latitude", "hfp", 
           "bio1", "bio2", "bio3", "bio4", "bio5", "bio6",
           "bio7", "bio8", "bio9", "bio10", "bio11", "bio12",
           "bio13", "bio14", "bio15", "bio16", "bio17","bio18", "bio19",
           "phihox", "orcdrc", "cecsol", "sndppt", "clyppt", "bldfie", "ntot",
           "elev", "Slope", "aspect","CTI", "Roughness", "Ruggedness", "DAH",
           "B1", "B2", "B3", "B4", "B5", "B6", "B7",
           "B8", "B8A", "B9", "B10", "B11", "B12",
           "NDVIS2", "NDVI32", "EVI", "NBRT", "greenness",
           "imprv", "ptrcv")
# delete the two factor 
rasterdir <- rasterdir[!(rasterdir[,'raster'] %in% catlist),]

# data distribution for raster --------

# plot distribution for the rest 
# cellstat <- as.data.frame(matrix(nrow = nrow(rasterdir), ncol = 4))
# colnames(cellstat) <- c("min", "max", "mean", "sd")
aa <- Sys.time()
cellstat <- foreach (ii = 1:2, .packages = "raster", .combine = "rbind") %dopar% {
    myraster <- raster(x = rasterdir[ii,'final_dir'])
    min <- cellStats(myraster, 'min')
    max <- cellStats(myraster, 'max')
    mean <- cellStats(myraster, 'mean')
    sd <- cellStats(myraster, 'sd')
    
    pdf(file = paste0("./plots_important/step0_prepare_data/raster_den_", rasterdir[ii, 'raster'], ".pdf"), width = 8, height = 8)
    par(cex = 1.5)
    pp <- raster::density(x = myraster, maxpixels = 10000, main = paste("CA raster", rasterdir[ii, 'raster']))
    mypar <- par(new = T)
    biom_den <- stats::density(biom[, rasterdir[ii, 'raster']], na.rm = T)
    plot(biom_den, col = "blue",
         xlim = c(min(pp$x), max(pp$x)), ylim = c(min(pp$y), max(pp$y)), xlab = "", main = "")
    mtext(side = 3, text = paste("N =", biom_den$n, "Bandwidth =", round(biom_den$bw, digits = 3)), col = "blue", cex = 1.5)
    legend("topright", legend=c("CA raster", "Transect data"),
           col=c("black", "blue"), lty=1, cex=1.2)
    par(mypar)
    dev.off()
    return(cbind(min, max, mean, sd))
}
Sys.time() - aa

# write out the cell stats
cellstat_r <- cbind(rasterdir[, 'raster'], cellstat)
write.csv(cellstat_r, file = "./final_data/metadata/raster_stats.csv")

# close the clusters --------
parallel::stopCluster(cl)
closeAllConnections()


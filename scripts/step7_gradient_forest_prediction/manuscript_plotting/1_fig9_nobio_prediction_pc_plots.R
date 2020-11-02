# Title: Plotting for figure insets
# Accompanying "predict_nogf_Scld_CAgrid_prcomp_RGB_20201025.pdf"
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Tue Apr 21 23:30:31 2020
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Sun Oct 25 18:39:42 2020

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE)
setwd("~/Lab/caledna_transect/")

library(gradientForest)
library(data.table)
library(ggplot2)

source("./scripts/function_transect.R")

# def functions --------

# def variables --------
# nobio/no
myseed = 17
nsites = 50449570
mysize = 50000 # number of sites to plot 
treatment = "nobio_CAgrid"
today = format(Sys.Date(), "%Y%m%d")
indir = "./derive_data/step7_gradient_forest_prediction/"
plotdir = paste0("./plots/step7_gradient_forest_prediction/gf_deco_all_Family_Presence_2000_2_0.05_FALSE_17/", treatment, "/")
dir.create(plotdir, recursive = T)
gfdir = "./derive_data/step6_gradient_forest/"
gfname = "gf_deco_all_Family_Presence_2000_2_0.05_FALSE_17"

# load data --------
# get the important variables 
gf = loadRData(file = paste0(gfdir, gfname, ".RData"))
imp.vars = names(importance(gf))
rm(gf)

# get the sampled Trns_CAgrb
set.seed(myseed)
id <- sample(x = 1:nsites, size = mysize)
load(file = paste0(indir, gfname, "/", treatment, "/", "predict_nogf_Scld_CAgrid_prcomp_RGB.RData"))
Scld_CArgbid = data.table::setDF(Scld_CArgb[id,])
rm(Scld_CArgb)

# load the PCs
load(file = paste0(indir, gfname, "/", treatment, "/", "predict_nogf_Scld_CAgrid_prcomp.RData"))

# load the sites 
load(file = paste0(indir, gfname, "/", treatment, "/", "predict_nogf_Scld_site.RData"))
# main --------
# get the loadings 
nvs <- dim(PCs$rotation)[1]
vec <- imp.vars[1:8]
lv <- length(vec)
vind <- rownames(PCs$rotation) %in% vec
scal <-0.025
xrng <- range(PCs$x[, 1], PCs$rotation[, 1]/scal) * 1.1
yrng <- range(PCs$x[, 2], PCs$rotation[, 2]/scal) * 1.1

# plot the back ground sample location #######
pdf(file = paste0(plotdir, "fig9_pcs_loading_nobio_CA_", today,".pdf"), height = 6, width = 6)
par(bg = NA)
plot(0, type = 'n', xlim = xrng, ylim = yrng,
        # suppress the boxes
     xaxt='n', yaxt='n', bty='n', xlab='', ylab='', asp = 1
)
# add plot axis 
##  Next add in your axis arrows:
arrows(min(xrng), 0, max(xrng), 0, lwd=1, col = "darkgrey", length=0)
arrows(0, min(yrng), 0, max(yrng), lwd=1, col = "darkgrey", length=0)
##  Add text annotations: 
# add axis labels 
xlab = paste0("PC1 (", scales::percent(summary(PCs)$importance['Proportion of Variance', 'PC1']), ")")
ylab = paste0("PC2 (", scales::percent(summary(PCs)$importance['Proportion of Variance', 'PC2']), ")")
text(min(xrng) - 1, 0, xlab, pos=3, col = "dimgrey")
text(0, min(yrng) + 10, ylab, pos=4, col = "dimgrey")

# plot the background samples 
points(PCs$x[id, 1:2], 
       pch = ".", 
       col = rgb(Scld_CArgbid[,'rr'], Scld_CArgbid[,'gg'], Scld_CArgbid[,'bb'], max = 255))

# now plot the 278 sites we analyzed
PCsites <- predict(PCs, Scld_site)
# combine them on the plots
points(PCsites[, 1:2], col = "grey", pch = "X", cex = 0.3)
# add arrow vector
arrows(rep(0, lv), rep(0, lv), PCs$rotation[vec, 1]/scal, PCs$rotation[vec, 2]/scal, lwd = 2, length = 0.0625)
jit <- 2
text(PCs$rotation[vec, 1]/scal + jit * sign(PCs$rotation[vec,1]), 
     PCs$rotation[vec, 2]/scal + jit * sign(PCs$rotation[vec,2]), 
     labels = vec, 
     cex = 1.2)
dev.off()

# cleanup --------
closeAllConnections()
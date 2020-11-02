# Title: alpha diversity and individual linear models --------
# Author: Meixi Lin
# Date: Sat Jan 25 16:04:31 2020
# Author:
# Date:
# Modification:

# preparation --------
rm(list = ls())
cat("\014")
setwd("~/Lab/caledna_transect")

# create output directory
plotdir <- "./plots/step3_alpha_diver/individual_linear_model/manuscript/"
dir.create(plotdir, recursive = T)

library(vegan)
library(tibble)
library(dplyr)
library(ggplot2)
library(phyloseq)

source("./scripts/function_transect.R")

date() # the execution date

# define function --------
do_lm <- function(phyrare, adivrare, myprimer, mymeasure, myvar) {
    # prepare the linear model data 
    print(paste0("Performing individual linear model for ", myprimer, "|", mymeasure, "alpha diversity:"))
    print(paste0("Model: ", myprimer, "|", mymeasure, " alpha_div ~ ", myvar))
    psdata <- phyrare[[myprimer]]
    adiv <- adivrare[[myprimer]]
    
    xxvar <- get_variable(psdata, varName = myvar)
    adivvar <- adiv[,mymeasure]
    
    fortest <- cbind(adivvar, xxvar) %>% 
        as.data.frame() %>% 
        tidyr::drop_na() 
    print(str(fortest))
    hh0 <- lm(adivvar ~ xxvar, data = fortest)
    print(summary(hh0))
    # get the variables 
    aa = hh0$coefficients[1]; aa <- unname(aa)
    bb = hh0$coefficients[2]; bb <- unname(bb)
    r2 = summary(hh0)$adj.r.squared; r2 <- unname(r2)
    padj = summary(hh0)$coefficients['xxvar', 'Pr(>|t|)']; padj <- unname(padj)
    res <- data.frame(myprimer, myvar, mymeasure, aa, bb, r2, padj)
    # first is the result summary, second is the model 
    output <- list(res, hh0)
    return(output)
}

# import phyloseq object and define variables --------
load("./derive_data/phy_rare/phyrare.RData")
load("./derive_data/step3_alpha_diver/adivrare_06012019.RData")

# main --------
# do lm for FITS, observed ========
myprimer = "FITS"
myvar = "greenness"
mymeasure = "Observed"

mydt = do_lm(phyrare, adivrare, myprimer, mymeasure, myvar)
hh0 = mydt[[2]]

# modelplotting ========
R2label = bquote(italic(R)^2 == .(format(mydt[[1]]$r2, digits = 3)))
pvallabel = bquote(italic(p-adj) == .(format(mydt[[1]]$padj, digits = 3)))
pdf(file = paste0(plotdir, "FITS_greenness_lm.pdf"), width = 8, height = 5)
op <- par(mfrow = c(1, 2)) 
plot(x = hh0$model$xxvar, y = hh0$model$adivvar, 
     xlab = myvar, ylab = mymeasure)
abline(hh0, col = "red")
plot(hh0, which = 2)
par(op)
dev.off()

# end --------
closeAllConnections()

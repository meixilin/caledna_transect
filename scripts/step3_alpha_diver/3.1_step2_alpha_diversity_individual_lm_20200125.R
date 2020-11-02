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
plotdir <- "./plots/step3_alpha_diver/individual_linear_model/"
dir.create(plotdir, recursive = T)
outdir <- "./derive_data/step3_alpha_diver/individual_linear_model/"
dir.create(outdir, recursive = T)

sink(file = paste0(outdir, "linear_model_logs_20200330.log"))
library(vegan)
library(tibble)
library(dplyr)
library(ggplot2)
library(phyloseq)

source("./scripts/function_transect.R")

date() # the execution date

# define function --------
do_lm <- function(phyrare, adivrare, myprimer, mymeasure, myvar, nobs.bonf, plotdir) {
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
    
    # model validation plots --------
    # if significant 
    sigif <- (padj * nobs.bonf) < mysig
    myplotdir = paste0(plotdir, sigif)
    dir.create(myplotdir, showWarnings = F, recursive = T)
    pdf(file = paste0(myplotdir,"/alpha_lm_individual_", myprimer, "_", myvar,"_", mymeasure, ".pdf"), width = 12, height = 6)
    op <- par(mfrow = c(2, 3), mar = c(5, 4, 1, 2), cex = 0.8) 
    plot(hh0, which = 1)
    plot(hh0, which = 2)
    plot(hh0, which = 3)
    plot(hh0, which = 5)
    E <- resid(hh0)
    hist(E, xlab = "Residuals", main = "")
    plot(x = fortest$xxvar, y = fortest$adivvar, xlab = myvar, ylab = mymeasure)
    abline(hh0, col = "red")
    mtext(side = 3, text = paste("primer =" , myprimer, "measure = ", mymeasure), cex = 0.6)
    mtext(side = 1, text = paste("a =" , round(hh0$coefficients[1], digits = 2), 
                                 "b =", round(hh0$coefficients[2], digits = 2), 
                                 "r2 =", round(summary(hh0)$adj.r.squared, digits = 2),
                                 "padj =", round(summary(hh0)$coefficients['xxvar', 'Pr(>|t|)'], digits = 2)),
          cex = 0.6)
    par(op)
    dev.off()
    return(res)
}

adjust_p_bonf <- function(pvalue, nobs.bonf) {
    palter = pvalue*nobs.bonf 
    if (palter > 1) {palter = 1}
    return(palter)
}

# import phyloseq object and define variables --------
load("./derive_data/phy_rare/phyrare.RData")
load("./derive_data/step3_alpha_diver/adivrare_06012019.RData")

mymeasures <- c("Observed", "Shannon")
mysig <- 0.05 
mycutoff <- 4
# bonferroni correction on significance 
nobs <- length(primers_commeco) * length(mymeasures) * length(contlist)
nobs.bonf <- length(mymeasures) * length(contlist)

# main --------
# initiate data frame
indi_lmdf <- data.frame(matrix(nrow = 0, ncol = 7))
colnames(indi_lmdf) <- c("myprimer", "myvar", "mymeasure", "aa", "bb", "r2", "padj")
for (ii in  1:length(primers_commeco)) {
    for (jj in 1:length(mymeasures)) {
        for (kk in 1:length(contlist)) {
            res <- do_lm(phyrare, adivrare, primers_commeco[ii], mymeasures[jj], contlist[kk], nobs.bonf, plotdir)
            indi_lmdf <- rbind(indi_lmdf, res)
        }
    }
}
colnames(indi_lmdf) <- c("Metabarcode", "Variable", "Measure", "Intercept_Estimate", "Variable_Estimate", "R2","P.value")

# test on the significance 
indi_lmdf <- indi_lmdf %>%
    dplyr::mutate(Adj_P.value = ifelse(P.value * nobs.bonf > 1, 1, P.value * nobs.bonf)) %>% 
    dplyr::mutate(Significant = (Adj_P.value < mysig)) %>%
    dplyr::arrange(P.value, desc(R2))

# write out the result --------
write.csv(indi_lmdf, file = paste0(outdir, "individual_lm_result.csv"))

# end --------
sink()
closeAllConnections()

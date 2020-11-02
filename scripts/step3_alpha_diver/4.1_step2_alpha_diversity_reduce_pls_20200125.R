# Title: alpha diversity and partial least square models --------
# Author: Meixi Lin
# Date: Tue Sep 17 14:35:58 2019

# preparation --------
rm(list = ls())
cat("\014")
setwd("~/Lab/caledna_transect")

# create output directory
plotdir <- "./plots/step3_alpha_diver/partial_least_square/"
dir.create(plotdir, recursive = T)
outdir <- "./derive_data/step3_alpha_diver/partial_least_square/"
dir.create(outdir, recursive = T)

sink(file = paste0(outdir, "partial_least_square_model_logs_20200330.log"))

library(vegan)
library(tibble)
library(dplyr)
library(ggplot2)
library(phyloseq)
library(pls)

source("./scripts/function_transect.R")
# source the VIP adds on: https://mevik.net/work/software/VIP.R 
source("https://mevik.net/work/software/VIP.R")

date() # the execution date

# define functions --------
do_pls <- function(phyrare, adivrare, myprimer, mymeasure, myvar, plotdir) {
    # prepare the linear model data 
    print(paste0("Performing partial least square model for ", myprimer, "|", mymeasure, " alpha diversity:"))
    psdata <- phyrare[[myprimer]]
    adiv <- adivrare[[myprimer]]
    
    xxvar <- get_variable(psdata, varName = myvar)
    adivvar <- adiv[,mymeasure]
    
    fortest <- cbind(adivvar, xxvar) %>% 
        as.data.frame() %>% 
        tidyr::drop_na() 
    # print(str(fortest))
    
    # start pls 
    ncomp.m <- 3 # explained ~ 95% variation in X variables
    plsmodel <- pls::plsr(adivvar ~ ., ncomp = ncomp.m, data = fortest, validation = "LOO", method = "oscorespls", jackknife = T)
    
    # output results 
    print(summary(plsmodel)) # root mean squared error of prediction RMSEP (dependent on the prediction)
    print(pls::R2(plsmodel)) # 1 - SSE/SST, where SST is the (corrected) total sum of squares of the response, and SSE is the sum of squared errors
    R2 = pls::R2(plsmodel)$val[1, 1, ncomp.m + 1]
    RMSEP = pls::RMSEP(plsmodel, estimate = "adjCV")$val[1, 1, ncomp.m + 1]
    
    # plot prediction and validation 
    pdf(file =paste0(plotdir, "PLS_predictions_", myprimer, "_",mymeasure, ".pdf"), width = 8, height = 6)
    pp <- par(mfrow = c(1,2))
    header = paste0("Alpha diversity PLS model ", myprimer, "|", mymeasure)
    plot(plsmodel, ncomp = ncomp.m, asp = 1, line = TRUE, plottype = "prediction",
         main = paste0("Prediction: ",header), cex.main = 0.8)
    
    text(0.95 * mean(adivvar), max(adivvar), labels = paste0("R2,", ncomp.m," comps = ", round(R2, digits = 3)))
    plot(plsmodel, plottype = "validation", main = "Validation plot", cex.main = 0.8, legendpos = "topright")
    par(pp)
    dev.off()
    
    # get coefficients; the variable importance (cutoff > 1) and significance 
    vip <- t(VIP(plsmodel)) %>% 
        data.frame %>%
        tibble::rownames_to_column(., var = "Var")
    colnames(vip)[-1] <- paste0("VIP.Comp", 1:3)
    # jack knife testing for the variance
    jk <- jack.test(object = plsmodel)
    print(jk)
    jk1 <- jk$pvalues %>%
        data.frame %>%
        tibble::rownames_to_column(., var = "Var")
    colnames(jk1)[-1] = "P.value"
    # combine values 
    # note: coef.mvr is used to extract the regression coefficients of a model, i.e. the B in y = XB (for the Q in y = TQ where T is the score). when comps is missing, this is cumulative 
    plscoef <- named_num2df(coef(plsmodel, ncomp = ncomp.m), yourcols = c("Var", "Coef")) %>%
        dplyr::arrange(desc(Coef)) %>%
        dplyr::left_join(., jk1, by = "Var") %>%
        dplyr::left_join(., vip, by = "Var") %>%
        as_tibble()
    print(plscoef)
    
    # plot coefficients and loading 
    pdf(file =paste0(plotdir, "PLS_coef_load_", myprimer, "_",mymeasure, ".pdf"), width = 8, height = 8)
    pp <- par(mfrow = c(2,1))
    plot(plsmodel, comps = 1:ncomp.m, plottype = "coefficients",main = paste0("Coefficients: ",header), cex.main = 0.8, legendpos = "topright", labels = contlist, pretty.xlabels = F)
    plot(plsmodel, comps = 1:ncomp.m, plottype = "loadings",main = paste0("Loadings: ",header), cex.main = 0.8, legendpos = "topright", labels = contlist, pretty.xlabels = F)
    par(pp)
    dev.off()
    
    res = data.frame(myprimer, mymeasure, ncomp.m, R2, RMSEP)
    to_return = list(res, plsmodel, plscoef)
    return(to_return)
}

# import phyloseq object and define variables --------
load("./derive_data/phy_rare/phyrare.RData")
load("./derive_data/step3_alpha_diver/adivrare_06012019.RData")

mymeasures <- c("Observed", "Shannon")
mysig <- 0.05 

# loop through measures and primers_commeco --------
# initiate data frame
res_plsdf <- data.frame(matrix(nrow = 0, ncol = 5))
colnames(res_plsdf) <- c("myprimer", "mymeasure", "ncomp.m", "R2", "RMSEP")
# initiate list for RData
plsmodel_list <- vector(mode = "list", length = length(primers_commeco) * length(mymeasures))
for (ii in  1:length(primers_commeco)) {
    for (jj in 1:length(mymeasures)) {
        res <- do_pls(phyrare, adivrare, primers_commeco[ii], mymeasures[jj], contlist, plotdir)
        res_plsdf <- rbind(res_plsdf, res[[1]])
        # store the plsmodel
        id <- (ii - 1) * length(mymeasures) + jj
        plsmodel_list[[id]] <- res[[2]]
        names(plsmodel_list)[id] <- paste(primers_commeco[ii], mymeasures[jj], sep = "|")
        # write the coefficient and VIP 
        coefdir <- paste0(outdir, "coef_vip_tables/")
        dir.create(coefdir, recursive = T, showWarnings = F)
        write.csv(res[[3]], file = paste0(coefdir, primers_commeco[ii], "_", mymeasures[jj], "_coef_vip_res.csv"))
    }
}

# write the data 
save(plsmodel_list, file = paste0(outdir, "plsmodel_list.RData"))
write.csv(res_plsdf, file = paste0(outdir, "pls_result_summary.csv"))

# ending --------
date()
sink(NULL)
closeAllConnections()

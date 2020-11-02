# Title: evaluate prediction uncertainty (Procrustes rotation) --------
# Author: Meixi Lin 
# Date: Fri Oct  2 16:41:07 2020
# Modification:

# preparation --------
rm(list = ls())
cat("\014")
options(echo = T)

setwd("~/Lab/caledna_transect/")
library(gradientForest)
library(vegan)
library(MASS)

library(dplyr)
library(ggplot2) 
library(ggpubr)
library(ggmap)
library(rgdal)
library(FSA) # for dunn test 

# set date
today = format(Sys.Date(), "%Y%m%d")

# some functions 
source("./scripts/function_transect.R")

# set variables -------
gfname <- "gf_deco_all_Family_Presence_2000_2_0.05_FALSE_17"
indir <- paste0("./derive_data/step7_gradient_forest_prediction/", gfname, "/")
sink(paste0(indir, "evaluate_prediction_uncertainty_", today, ".log"))
gfdir <- "./derive_data/step6_gradient_forest/"

plotdir0 <- "./plots/step7_gradient_forest_prediction/"
plotdir <- paste0(plotdir0, gfname,"/eval_prediction/")
dir.create(plotdir, recursive = T)

# define plotting functions --------
do_dunn_test <- function(forplot) {
    mymin = min(forplot$Procrustes_residuals)
    mymax = max(forplot$Procrustes_residuals)
    dunnfile <- FSA::dunnTest(Procrustes_residuals ~ transect, data = forplot, method = "bonferroni")$res 
    ypos = seq(mymax, by = (mymax - mymin)/15, length.out = nrow(dunnfile))
    dunnfile = dunnfile %>% 
        dplyr::mutate(.y. = "Procrustes_residuals",
                      method = "Dunn-test",
                      p.signif = dplyr::case_when(P.adj < 0.001 ~ "***",
                                                  P.adj < 0.01 & P.adj >= 0.001 ~ "**",
                                                  P.adj < 0.05 & P.adj >= 0.01 ~ "*",
                                                  TRUE ~ "ns")) %>%
        tidyr::separate(col = Comparison, into = c("group1", "group2"), sep = " - ") %>%
        dplyr::arrange(group1, group2) %>%
        dplyr::mutate(y.position = ypos) %>%
        dplyr::rename(p.adj = P.adj) %>% 
        dplyr::select(.y., group1, group2, p.adj, p.signif, method, y.position) %>%
        tidyr::as_tibble()
    return(dunnfile)
}

plot_residuals <- function(thispro, biom, cabase, calim_df, plotname) {
    # get the forplot 
    prores = residuals(thispro)
    longlat = biom %>% dplyr::select(MatchName, Longitude, Latitude, transect)
    forplot = named_num2df(prores, yourcols = c("MatchName", "Procrustes_residuals")) %>% 
        dplyr::left_join(longlat, by = "MatchName")
    levels(forplot$transect) = c("Coast", "Forest", "Shrub")
    
    # start plotting
    pp1 <- ggplot() + 
        geom_polygon(data = calim_df, aes(x = long, y = lat, group = group), colour = "snow4", fill = 'snow4', alpha = 0) +
        geom_point(data = forplot, aes(x = Longitude, y = Latitude, shape = transect, color = Procrustes_residuals), alpha = 0.6, size = 3) +
        scale_color_gradient(low = "blue", high = "red") + 
        labs(x = "Longitude", y = "Latitude", color = "Procrustes Residuals", title = plotname) +
        facet_grid(. ~ transect) +
        theme_bw() + 
        guides(shape = "none") +
        theme(legend.position = "bottom",
              legend.title = element_text(size = 9),
              legend.text = element_text(size = 8)) + 
        ggsn::scalebar(data = calim_df,
                       dist = 150, dist_unit = "km", border.size = 0.5, st.size = 2.5, st.dist = 0.05,
                       transform = TRUE, model = "WGS84", location = "bottomleft") 
    
    # get the residuals by site rotations
    forplot2 = forplot %>% dplyr::arrange(transect)
    forplot2$MatchName = factor(forplot2$MatchName, levels = forplot2$MatchName)
    pp2 <- ggplot(forplot2, aes(x = MatchName, y = Procrustes_residuals, fill = transect)) + 
        geom_bar(stat = "identity") + 
        theme_bw() + 
        scale_fill_discrete(name = "Transect") +
        labs(x = "Samples", y = "Procrustes Residuals", title = plotname) + 
        theme(axis.text.x = element_text(size = 3, angle = 90, hjust = 1),
              legend.position = "bottom",
              legend.title = element_text(size = 9),
              legend.text = element_text(size = 8)) 
    pp3 <- pp2 + 
        theme(axis.text.x = element_blank(), 
              axis.ticks.x = element_blank())
    # get dunn test results 
    print(kruskal.test(Procrustes_residuals ~ transect, data = forplot))
    print(FSA::dunnTest(Procrustes_residuals ~ transect, data = forplot, method = "bonferroni"))
    dunnres = do_dunn_test(forplot)
    pp4 <- ggplot(forplot, aes(x = transect, y = Procrustes_residuals)) + 
        geom_boxplot(aes(fill = transect)) + # make sure put in fill within geom_boxplot not in the ggplot section
        ggpubr::stat_pvalue_manual(data = dunnres, label = "p.signif", inherit.aes = F) + 
        ggpubr::stat_compare_means(method = "kruskal.test", label.x.npc = "left", label.y.npc = "top", size = 4) + 
        labs(x = "Transect", y = "Procrustes residuals") + 
        theme_bw() + 
        theme(legend.position = "none")
    pplist = list(pp1, pp2, pp3, pp4)
    return(pplist)
}

# run the procrustes for three ED objects
run_procrustes <- function(bio.diss_NMDS, ED.object, 
                           runname = c("extrap_FALSE", "extrap_TRUE", "nobio_CAgrid"), 
                           plotname = c("Gradient Forest Prediction", "Gradient Forest Prediction (extrap = T)", "Uninformed Environmental Matrix")) {
    print(paste0("INFO:", runname, " ----------------------------------------"))
    # Perform procrustes rotation from two ordinations --------
    # NOTE: X = target matrix, Y = matrix to be rotated. Arrows are from Y space to X matrix. (From environmental/transformed environmental to biological)
    # random, store the output 
    ED.object_NMDS <- metaMDS(ED.object, k = 3, try = 50, trymax = 500)
    pro_NMDS <- procrustes(X = bio.diss_NMDS, Y = ED.object_NMDS)
    print(summary(pro_NMDS))
    res_pro_NMDS <- residuals(pro_NMDS)
    
    pplist_NMDS = plot_residuals(thispro = pro_NMDS, biom, cabase, calim_df, plotname = paste0("eDNA Biological Matrix NMDS | ", plotname, " NMDS"))
    ggsave(filename = paste0(plotdir, "procrustes_NMDSk3_map_", runname, "_", today, ".pdf"), plot = pplist_NMDS[[1]], height = 4.5, width = 8)
    ggsave(filename = paste0(plotdir, "procrustes_NMDSk3_barplot_", runname, "_", today, ".pdf"), plot = pplist_NMDS[[2]], height = 4, width = 12)
    ggsave(filename = paste0(plotdir, "procrustes_NMDSk3_barplotNoName_", runname, "_", today, ".pdf"), plot = pplist_NMDS[[3]], height = 4, width = 8)
    ggsave(filename = paste0(plotdir, "procrustes_NMDSk3_boxplot_", runname, "_", today, ".pdf"), plot = pplist_NMDS[[4]], height = 5, width = 5)
    
    # plot the rotation 
    pdf(file = paste0(plotdir, "procrustes_rotation_", runname, "_", today, ".pdf"), height = 8, width = 8)
    plot(pro_NMDS, main = paste0("Procrustes errors: ", plotname, " NMDS (circles) to eDNA NMDS (arrows)"), cex.main = 0.8)
    dev.off()
    
    # Save the object used for procrustes
    saveRDS(ED.object_NMDS, file = paste0(indir, runname,"/EuclideanD_gf_predict_NMDSk3_", runname, "_", today, ".rds"))
    # perform the test 
    protest_NMDS <- protest(X = bio.diss_NMDS, Y = ED.object_NMDS)
    print(protest_NMDS)
    print(summary(protest_NMDS))
    
   
    return(pro_NMDS)
}

# load and prepare data --------
# full metadata ========
load("./final_data/Final_metadata.RData")
load("./derive_data/step0_prepare_data/sample_map/CA_ggplot_basemap_20201003.RData")
calim <- readOGR(dsn = "./maps_vectors_rasters/vectors/CA_boundary_TIGER/CA_boundary_TIGER.shp")
raster::crs(calim) 
calim_df <- ggplot2::fortify(calim)
# gradient forest data real ========
gf <- loadRData(file = paste0(gfdir, gfname, ".RData"))
bio.data <- apply(gf$Y, 2, function(xx) {yy = as.numeric(xx == "TRUE")})
rownames(bio.data) <- rownames(gf$Y)
bio.diss <- vegan::vegdist(bio.data, method = "jaccard")
# get the NMDS in one run
# random, store the output (but check the plots the residuals were very similar)
bio.diss_NMDS <- metaMDS(bio.diss, k = 3, try = 50, trymax = 500) # No convergence
saveRDS(bio.diss_NMDS, file = paste0(indir, "/eDNA_biological_jaccard_NMDSk3_", today, ".rds"))

# convert "ED" matrix of RF transformed sites ========
# gradient forest predictions, no extrapolations
load(file = paste0(indir, "/extrap_FALSE/predict_gf_Trns_site_noextrap.RData")) # Trns_site_noextrap
ED.bio.no <- stats::dist(Trns_site_noextrap, method = "euclidean")
# gradient forest predictions, with extrapolations
load(file = paste0(indir, "/extrap_TRUE/predict_gf_Trns_site_extrap.RData")) # Trns_site_extrap
ED.bio.ex <- stats::dist(Trns_site_extrap, method = "euclidean")
# Scaled environmental sites 
load(file = paste0(indir, "/nobio_CAgrid/predict_nogf_Scld_site.RData"))
ED.nobio <- stats::dist(Scld_site, method = "euclidean")

# check if the rowname is the same ========
table(rownames(Trns_site_extrap) == rownames(bio.data))
table(rownames(Trns_site_noextrap) == rownames(bio.data))
table(rownames(Scld_site) == rownames(bio.data))

# Perform procrustes rotation bio.diss VS ED.bio.no --------
pro_extrapF = run_procrustes(bio.diss_NMDS = bio.diss_NMDS, ED.object = ED.bio.no, 
               runname = "extrap_FALSE", plotname = "Gradient Forest Prediction")
saveRDS(pro_extrapF, file = paste0(indir, "extrap_FALSE/procrustes_extrap_FALSE_", today, ".rds"))

pro_extrapT = run_procrustes(bio.diss_NMDS = bio.diss_NMDS, ED.object = ED.bio.ex, 
               runname = "extrap_TRUE", plotname = "Gradient Forest Prediction (extrap = T)")
saveRDS(pro_extrapT, file = paste0(indir, "extrap_TRUE/procrustes_extrap_TRUE_", today, ".rds"))

pro_nobio = run_procrustes(bio.diss_NMDS = bio.diss_NMDS, ED.object = ED.nobio, 
               runname = "nobio_CAgrid", plotname = "Uninformed Environmental Matrix")
saveRDS(pro_nobio, file = paste0(indir, "nobio_CAgrid/procrustes_nobio_CAgrid_", today, ".rds"))

# clean up --------
sink()
closeAllConnections()

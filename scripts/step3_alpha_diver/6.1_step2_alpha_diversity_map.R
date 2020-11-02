# Title: alpha diversity and category variables plotting --------
# Author: Meixi Lin
# Date: Tue May  7 12:34:27 2019
# Author:
# Date:
# Modification:

# preparation --------
rm(list = ls())
cat("\014")
setwd("~/Lab/caledna_transect")

library(vegan)
library(tibble)
library(reshape2)
library(dplyr)
library(ggplot2)
library(phyloseq)
library(ggmap)
library(rgdal)
library(ggsn)
library(gridExtra)

source("./scripts/function_transect.R")
source("./scripts/step4_beta_diver/0_util_query_var_color.R")
date() # the execution date

# import phyloseq object --------
load("./derive_data/phy_rare/phyrare.RData")
load("./derive_data/step3_alpha_diver/adivrare_06012019.RData")

# define variables 
plotdir <- "./plots/step3_alpha_diver/map/"
dir.create(plotdir, recursive = T)
mymeasures <- c("Observed", "Shannon")
mysig <- 0.05 
mycutoff <- 4
myvar="transect"

# ggmap::register_google(key = "YOUR API KEY")

calim <- readOGR(dsn = "./maps_vectors_rasters/vectors/CA_boundary_TIGER/CA_boundary_TIGER.shp")
raster::crs(calim) 
calim_df <- ggplot2::fortify(calim)

# cabase <- get_map(location = "california", maptype = "roadmap", zoom = 6)
# load predefined cabase data
load("./derive_data/step0_prepare_data/sample_map/CA_ggplot_basemap_20201003.RData")

# define functions --------
# get fortest for each primer, all possible measures
get_fortest <- function(myprimer, myvar, mymeasures, mycutoff) {
    psdata <- phyrare[[myprimer]]
    adivvar <- adivrare[[myprimer]] %>%
        tibble::rownames_to_column(var = "MatchName")
    category <- get_variable(psdata, varName = c("MatchName", "Longitude", "Latitude", myvar))
    print(sort(table(category[,myvar])))
    
    colnames(category) = c("MatchName", "Longitude", "Latitude", "Category")
    sumby <- table(category[,"Category"]) %>% as.data.frame() %>% 
        filter(Freq > mycutoff)
    # if all data are scatterred, 0 observation for sumby  
    if (nrow(sumby) < 3) {
        print(nrow(sumby))
        print(paste("for", myvar, "no data for stats test."))
        next
    }
    # get the category and alpha diversity value 
    fortest <- base::merge(adivvar, category, by = "MatchName") %>%
        dplyr::filter(Category %in% sumby$Var1) 
    fortest$Category <- droplevels(fortest$Category)
    
    # melt the data set 
    fortest <- fortest %>%
        reshape2::melt(id.vars = c("MatchName", "Longitude", "Latitude", "Category"), measure.vars = mymeasures, value.name = "adiv.value") %>%
        dplyr::mutate(Metabarcode = myprimer) %>%
        dplyr::mutate_if(is.factor, as.character)
    
    colnames(fortest) = c("MatchName", "Longitude", "Latitude", "Category", "Measure", "adiv.value", "Metabarcode")
    return(fortest)
}

get_fortest_all_primer <- function(myprimers, myvar, mymeasures, mycutoff) {
    # automatically pass the myvar, mymeasures and mycutoff to the get_fortest function
    fortest.list <- lapply(myprimers, get_fortest, myvar = myvar, mymeasures = mymeasures, mycutoff = mycutoff)
    fortestall <- dplyr::bind_rows(fortest.list)
}

# main: perform plotting for all primers used and the variables of interest  --------
forplot <- get_fortest_all_primer(myprimers = primers_commeco, myvar = myvar, mymeasures = mymeasures, mycutoff = mycutoff)

pp1 <- ggmap(cabase) + 
    geom_polygon(data = calim_df, aes(x = long, y = lat, group = group), colour = "snow4", fill = 'snow4', alpha = 0.5) +
    geom_point(data = forplot[forplot$Measure == "Observed",], aes(x = Longitude, y = Latitude, color = adiv.value, shape = Category), alpha = 0.8) +
    facet_grid(. ~ Metabarcode) +
    scale_color_gradient2() + 
    labs(x = "Longitude", y = "Latitude", color = "Observed alpha diversity") +
    theme_bw() + 
    guides(shape = "none") + 
    theme(legend.position = "bottom",
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 8))

pp2 <- ggmap(cabase) + 
    geom_polygon(data = calim_df, aes(x = long, y = lat, group = group), colour = "snow4", fill = 'snow4', alpha = 0.5) +
    geom_point(data = forplot[forplot$Measure == "Shannon",], aes(x = Longitude, y = Latitude, color = adiv.value, shape = Category), alpha = 0.8) +
    facet_grid(. ~ Metabarcode) +
    scale_color_gradient2() + 
    labs(x = "Longitude", y = "Latitude", color = "Shannon alpha diversity", shape = "Transect") +
    theme_bw() + 
    guides(shape = "none") +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 8))
newpp = grid.arrange(grobs = list(pp1,pp2), nrow = 2, ncol = 1)

pp2.1 <- ggmap(cabase) + 
    geom_polygon(data = calim_df, aes(x = long, y = lat, group = group), colour = "snow4", fill = 'snow4', alpha = 0.5) +
    geom_point(data = forplot[forplot$Measure == "Shannon",], aes(x = Longitude, y = Latitude, color = adiv.value, shape = Category), alpha = 0.8) +
    facet_grid(. ~ Metabarcode) +
    scale_color_gradient2() + 
    labs(x = "Longitude", y = "Latitude", color = "Shannon alpha diversity", shape = "Transect") +
    theme_bw() + 
    # guides(shape = "none") + 
    scale_shape_discrete(labels = c("Coast", "Forest", "Shrub")) + 
    theme(legend.position = "left",
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 8))

ggsave(filename = "alpha_diversity_map.pdf", path = plotdir, plot = newpp, device = "pdf", width = 8, height = 6)
ggsave(filename = "alpha_diversity_legend.pdf", path = plotdir, plot = pp2.1, device = "pdf", width = 12, height = 4)
# cleaning up --------
closeAllConnections()

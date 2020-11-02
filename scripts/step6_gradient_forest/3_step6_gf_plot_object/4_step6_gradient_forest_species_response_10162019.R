# Title: evaluate species response --------
# Author: Meixi Lin
# Date: Wed Oct 16 14:36:25 2019
# Author:
# Date:
# Modification:

# preparation --------
options(echo = TRUE)
rm(list = ls())
cat("\014")

setwd("~/Lab/caledna_transect/")
# setwd("/u/home/m/meixilin/project-rwayne/abiotic_transect/")
# load packages
library(dplyr)
library(phyloseq)
library(gradientForest)
library(ggplot2) # 
library(ggmap)
library(rgdal)
library(ggsn)
library(gridExtra)
source("./scripts/function_transect.R")
functions.file <- list.files("./r_codes/step4_gradient_forest/functions/")
lapply(functions.file, function(xx) {
    source(file = paste0("./r_codes/step4_gradient_forest/functions/", xx))
    return(0)
    })

ggmap::register_google(key = "YOUR_API_KEY")

# define functions --------
# find matchname for a given taxa 
find_matchname <- function(gf, taxentry, imp) {
    # check if the entry is in the gf 
    if (!(taxentry %in% colnames(gf$Y))) {
        print("entry requested not in model, check!")
        break;
    }
    name.match <- gf$Y[,taxentry]
    id <- names(name.match)[name.match == "TRUE"]
    count <- length(id)
    to_return <- list(id, count, taxentry, as.numeric(imp))
    names(to_return) <- c("matchname", "n", "taxname", "imp")
    return(to_return)
}

# plot the corresponding matchname and phydeco 
plot_match <- function(biom, matchname, myvar = "majorhab") {
    # define colour scheme 
    mycolours <- pals::polychrome(n = nnlev + 1)
    names(mycolours) <- c("Non-vegetated", "unknown", "Aquatic", "Herbaceous-Dominated", "Shrub-Dominated", "Tree-Dominated")
    # get a subset of metadata 
    mybiom <- biom %>%
        dplyr::mutate(inmatch = MatchName %in% matchname[["matchname"]],
                      majorhab = as.character(majorhab))
    naid <- is.na(mybiom[,myvar])
    mybiom[naid,myvar] <- "unknown"
    pp <- ggmap(cabasemap) + 
        geom_polygon(data = calim_df, aes(x = long, y = lat, group = group), colour = "lightblue3", fill = 'lightblue3', alpha = 0.5, show.legend = FALSE) +
        geom_point(aes(x = Longitude, y = Latitude), data = mybiom[mybiom$inmatch == FALSE,], size = 1.5, colour = "slategrey", show.legend = FALSE) +
        geom_point(aes(x = Longitude, y = Latitude, fill = eval(as.name(myvar))), data = mybiom[mybiom$inmatch == TRUE,], shape = 23, colour = 'black', size = 2, show.legend = FALSE) +
        scale_fill_manual(values = mycolours) +
        labs(title = paste(matchname[['taxname']], "N =", matchname[['n']], "R^2 = ", round(matchname[["imp"]], digits = 2)),
             x = "Longitude", y = "Latitude") + 
        theme_bw() +
        theme(plot.title = element_text(size=12))
        
    return(pp)
}

# dig into the original sequences 

# get arguments --------
args <-  commandArgs(trailingOnly=TRUE)

# this script allow 3 inputs:
# args[1]: the absolute path to the original gradient forest object
# args[2]: filename
# args[3]: plotting directory

# load data --------
gf <- loadRData(as.character(args[1])) # RENAME WHAT EVER GF to "gf"
filename <- as.character(args[2]) # this should be the name of the "gf" object tested
plotdir <- as.character(args[3])

load("./final_data/Final_metadata.RData")

calim <- readOGR(dsn = "./maps_vectors_rasters_local/vectors/CA_boundary_TIGER/CA_boundary_TIGER.shp")
raster::crs(calim) 
calim_df <- ggplot2::fortify(calim)

# plot the value on map 
cabasemap <- get_map(location = "california", zoom = 6)

# main --------
ntax = 9
taxlist <- sort(gf$result, decreasing = T)[1:ntax] %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var = "taxentry")
colnames(taxlist) <- c("taxentry", "imp")
pplist <- apply(taxlist, 1,function(taxres) {
    match <- find_matchname(gf, taxres['taxentry'], taxres['imp'])
    pp <- plot_match(biom, match, "majorhab")
    return(pp)
})

# arrange the grid to plot it 
m1 <- gridExtra::grid.arrange(grobs = pplist, ncol = 3, nrow = 3)

# save the plot 
dir.create(plotdir, recursive = T)
ggplot2::ggsave(filename = paste0(plotdir, filename, "_top9_species.pdf"), plot = m1, device = "pdf", limitsize = F, height = 9, width = 9)

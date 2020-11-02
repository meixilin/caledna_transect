# Title: import original sequence reads files; build RData --------
# Author: Meixi Lin
# Date: Thu Oct 17 11:25:06 2019
# Author:
# Date:
# Modification:

# preparation --------
options(echo = TRUE)
rm(list = ls())
cat("\014")

setwd("~/Lab/caledna_transect/")
library(dplyr)
library(stringr)
date() # the execution date

# define functions ---------
# rename header 
# header: a vector containing original reads 
# transect: coastal/forest/shrub
rename_header <- function(matchname, transect) {
    # change names for blanks
    # print(matchname); print(length(matchname))
    blanks <- matchname %>%
        stringr::str_which(., pattern = "Blank|neg")
    if (length(blanks) > 0) {
        matchname[blanks] <- paste0("Blank", 1:length(blanks))
    }
    # for the not blank names 
    if (transect == "SCRUB") {
        mypattern = "K[0-9][0-9][0-9][0-9].[A-C][1-2]"
        myid <- matchname %>%
            stringr::str_which(., pattern = mypattern)
        matchname[myid] <- stringr::str_extract(matchname[myid], mypattern) %>%
            gsub("\\.", "", .)
    } else {
        mypattern = "K[0-9][0-9][0-9][0-9][A-C][1-2]"
        myid <- matchname %>%
            stringr::str_which(., pattern = mypattern)
        matchname[myid] <- stringr::str_extract(matchname[myid], mypattern)
    }
    
    # for shrub sites, check if it contains the duplicated samples and mark
    if (transect == "SCRUB") {
        dup <- c("K0058C2", "K0117C2", "K0142A2", "K0142C2", "K0177A1")
        matchname[matchname %in% dup] <- paste0(matchname[matchname %in% dup], "_S")
    }
    # print(matchname); print(length(matchname))
    return(matchname)
}

read_asv <- function(filepath, transect = c("COASTAL", "FOREST", "SCRUB"), primer) {
    # load file 
    filename <- paste0(filepath, transect, "/", primer, "_ASV_taxonomy_detailed.txt")
    raw <- read.table(file = filename, header = T, sep = "\t", stringsAsFactors = F, comment.char = "")
    # change matchnames 
    matchid <- str_which(colnames(raw), pattern = "L001")
    matchname <- colnames(raw)[matchid] %>% 
        rename_header(., transect)
    colnames(raw)[matchid] <- matchname
    # move asv to the last lines 
    neworder <- c(setdiff(1:ncol(raw), matchid), matchid)
    raw <- raw[,neworder]
    print(colnames(raw))
    # define a class 
    myasv <- list(asv = raw, primer = primer, transect = transect)
    class(myasv) <- "taxdetail"
    return(myasv)
}

# main --------
filepath <- "./raw_data/transect_ASV/"
primers <- c("16S", "18S", "CO1", "FITS", "PITS")
transects <- c("COASTAL", "FOREST", "SCRUB")
asvdb <- vector(mode = "list", length = length(primers) * length(transects))
for (ii in 1:length(primers)) {
    for (jj in 1:length(transects)) {
        id <- (ii - 1) * length(transects) + jj
        asvdb[[id]] <- read_asv(filepath, transect = transects[jj], primer = primers[ii])
        names(asvdb)[id] <- paste(primers[ii], transects[jj], sep = "|")
    }
}

save(asvdb, file = "./raw_data/transect_ASV/taxonomy_detail_all.RData")

# make it an lazyloading object 

# ending --------
date()
closeAllConnections()



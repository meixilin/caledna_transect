# Title: summary stats plotting --------
# This one is for use with using the interactive tree of life and ncbi taxonomy dump files 
# Author: Meixi Lin
# Date: Wed Jun 12 15:06:31 2019
# Author:
# Date:
# Modification:

# preparation --------
rm(list = ls())
cat("\014")

setwd("/Users/linmeixi/UCLA/Lab/abiotic_transect")
library(dplyr)
library(ggplot2)
library(ranacapa)
library(phyloseq)
library(reshape2)
library(ggsci)
library(ape)
library(stringr)

source("./scripts/function_transect.R")

date() # the execution date

# load taxonomy files --------
nodes <- read.table(file = "./final_data/TAXO/nodes.dmp", sep = "|", nrows = 1669688)
# get all phyla --------
phyla <- nodes[nodes[,3] == "\tphylum\t",1]
write.table(phyla, file = "./derive_data/step1_mk_phyloseq/phyla/taxid_all_phyla.txt", row.names = F, quote = F, col.names = F)
rm(nodes)

# get the common name -------
# phyla <- read.table(file = "./derive_data/step1_mk_phyloseq/phyla/taxid_all_phyla.txt")
# phyla <- phyla[,1]
names <- read.table(file = "./final_data/TAXO/names.dmp", sep = "|", comment.char = "", quote = "")
phynames <- names[(names[,1] %in% phyla & names[,4] == "\tscientific name\t"), 1:4]
colnames(phynames) <- c("tax_id", "name_txt", "unique_name", "name_class")
# strip "\t\t" 
phynames = apply(phynames, 2, function(xx) {
        xx <- stringr::str_remove_all(xx, pattern = "\t")
}) %>% 
    as.data.frame() 
write.csv(phynames, file = "./derive_data/step1_mk_phyloseq/phyla/taxname_all_phyla.txt")
rm(names)

# get the tree structure on NCBI common tree generator --------

# load the phylum summary --------
mydata <- read.csv(file = "./derive_data/step2_data_description/summary_stats/phy_deco_taxa_entry_count_Phylum.csv", row.names = 1)
mydata$Phylum[!(mydata$Phylum %in% phynames$name_txt)] # only empty and NA phylum were not matching 

colnames(phynames)[2] = "Phylum"
# generate data for tol --------
toldata <- dplyr::full_join(x = phynames, y = mydata, by = "Phylum") %>% 
    dplyr::select(Phylum, X16S, X18S, CO1, FITS, PITS, all) 
toldata$Phylum <- stringr::str_replace_na(toldata$Phylum, replacement = "unknown")
candiphy <- stringr::str_detect(toldata$Phylum, "^[c|C]andi")
toldata$Phylum[candiphy] <- paste0("'", toldata$Phylum[candiphy], "'")
toldata[is.na(toldata)] <- 0

write.csv(toldata, file = "./derive_data/step1_mk_phyloseq/phyla/summary_phycounts_itol.csv")

# merge the "unknown" and "" Phylum as "root" in excel 

# generate data for tol --------
toldata <- dplyr::right_join(x = phynames, y = mydata, by = "Phylum") %>% 
    dplyr::select(Phylum, X16S, X18S, CO1, FITS, PITS, all) 
toldata$Phylum <- stringr::str_replace_na(toldata$Phylum, replacement = "unknown")
candiphy <- stringr::str_detect(toldata$Phylum, "^[c|C]andi")
toldata$Phylum[candiphy] <- paste0("'", toldata$Phylum[candiphy], "'")
toldata[is.na(toldata)] <- 0

write.csv(toldata, file = "./derive_data/step1_mk_phyloseq/phyla/summary_phycounts_itol_v2.csv")

# delete the phyla that does not have reads --------
toldata2 <- dplyr::right_join(x = phynames, y = mydata, by = "Phylum") %>% 
    dplyr::select(tax_id, Phylum, X16S, X18S, CO1, FITS, PITS, all)
write.csv(toldata, file = "./derive_data/step1_mk_phyloseq/phyla/summary_phycounts_itol_no0.csv")
write.csv(toldata2$tax_id, file = "./derive_data/step1_mk_phyloseq/phyla/taxid_no0_phyla.txt", row.names = F, col.names = F, quote = F)
# DONE

# use ggplot to generate the circle plot --------
mytree <- read.tree(file = "./derive_data/step1_mk_phyloseq/phyla/phyliptree_no0_oldname.phy")
# reorder toldata2 to have the same order
toldata3 <- toldata2 
# new column "unknown" 
unknown <- toldata3[which(toldata3$Phylum == ""), 3:8] + toldata3[which(is.na(toldata3$Phylum)), 3:8]
toldata3[which(toldata3$Phylum == ""), 3:8] = unknown
toldata3[which(toldata3$Phylum == ""), 1:2] = c(NA, "unknown")
toldata3 = toldata3[-which(is.na(toldata3$Phylum)),-c(1, ncol(toldata3))]
# reorder toldata to 
# 1. match the phylums, NEED STR_REPLACE_ALL some had more than one space 
toldata3$Phylum = str_replace_all(toldata3$Phylum, " ", "_")
# 2. generate factor from Phylum using the order from my tree
toldata3$Phylum = factor(x = toldata3$Phylum, levels = c(mytree$tip.label, "unknown"))

# start plotting w/ ggplot2 --------
# use a transparent background 
MyTheme_transparent <- theme(
    panel.background = element_rect(fill = "transparent", colour = NA), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent", colour = NA), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent", colour = NA), # get rid of legend panel bg
    legend.key = element_rect(fill = "transparent", colour = NA), # get rid of key legend fill, and of the surrounding
    axis.line = element_line(colour = NA) # adding a black line for x and y axis
)

forplot <- melt(toldata3)
y_labels = levels(forplot$variable)
y_breaks = seq_along(y_labels) + 5
forplot$var2 <- as.numeric(factor(forplot$variable)) + 5

# create labels, first subset by variable level 
forlabs <- subset(forplot, variable==levels(forplot$variable)[nlevels(forplot$variable)])
forlabs <- forlabs[order(forlabs$Phylum),]
forlabs$ang <- seq(from=(360/nrow(forlabs))/1.5, to=(1.5*(360/nrow(forlabs)))-360, length.out=nrow(forlabs))+80
forlabs$hjust <- 0
forlabs$hjust[which(forlabs$ang < -90)] <- 1
forlabs$ang[which(forlabs$ang < -90)] <- (180+forlabs$ang)[which(forlabs$ang < -90)]
# change Phylum to integer 
# forplot$Phylum = as.integer(forplot$Phylum)

# use log10 transformation for better visualization result 
pp <- ggplot() +
    geom_tile(data = forplot, aes(x = Phylum, y = var2, fill = value), colour="white") +
    scale_fill_gradient(trans = "log10", name = "Number of ASVs") +
    geom_text(data = forlabs, aes(x=Phylum, y=var2 * 1.1, label = Phylum,
                                  angle = ang, hjust = hjust), size=3) +
    ylim(c(0, max(forplot$var2) + 5)) +
    coord_polar(theta="x") +
    theme(panel.background=element_blank(),
          axis.title=element_blank(),
          panel.grid=element_blank(),
          axis.text =element_blank(),
          axis.ticks=element_blank()
         ) +
    MyTheme_transparent 
    
    
ggsave(filename = "./plots_important/step1_mk_phyloseq/phylo_circ_asv.pdf", plot = pp, device = "pdf", width = 8, height = 6, bg = "transparent")
# save(pp, file = "./derive_data/step1_mk_phyloseq/phyla/pp.RData")
# Note: in pp, there are 87 levels, minus the "unknown" level, we have 86 phyla in total
# now we need to add polygon that shows the kingdom they are in
# polygon <- data.frame(
#     supk = rep(c('Archaea', 'Eukaryota', 'Bacteria'), each = 4),
#     Phylum = factor(x = rep(c("Thaumarchaeota", "Crenarchaeota", "Colponemidia", "Xanthophyceae", "Balneolaeota", "Bacteroidetes"), each = 2), levels = levels(forplot$Phylum)), 
#     var2 = rep(c(0,5,5,0), times = 3)
# )

polygon <- data.frame(
    supk = rep(c('Archaea', 'Eukaryota', 'Bacteria'), each = 4),
    Phylum = rep(c(0.5,3.5, 3.5, 48.5, 48.5,86.5), each = 2),
    var2 = rep(c(0,max(forplot$var2) + 5,max(forplot$var2) + 5,0), times = 3)
)


pp2 <- ggplot() +
    # geom_tile(data = forplot, aes(x = Phylum, y = var2, fill = value), colour="white") +
    geom_polygon(data = polygon, aes(x = Phylum, y = var2, fill = supk)) +
    xlim(c(0.5,87.5)) +
    ylim(c(0, max(forplot$var2) + 5)) +
    coord_polar(theta="x") +
    theme(panel.background=element_blank(),
          axis.title=element_blank(),
          panel.grid=element_blank(),
          axis.text =element_blank(),
          axis.ticks=element_blank()
    ) +
    MyTheme_transparent
ggsave(filename = "./plots_important/step1_mk_phyloseq/phylo_circ_asv_col.pdf", plot = pp2, device = "pdf", width = 8, height = 6, bg = "transparent")

# now add for eukaryo
as.integer(factor(c("Xenacoelomorpha", "Porifera", "Zoopagomycota", "Chytridiomycota", "Streptophyta", "Chlorophyta"), levels = levels(forplot$Phylum)))
polygon2 <- data.frame(
    king = rep(c('Metazoa', 'Fungi', 'Viridiplantae'), each = 4),
    Phylum = rep(c(8.5,32.5, 34.5, 42.5, 42.5,44.5), each = 2),
    var2 = rep(c(0,max(forplot$var2) + 5,max(forplot$var2) + 5,0), times = 3)
)


pp3 <- ggplot() +
    # geom_tile(data = forplot, aes(x = Phylum, y = var2, fill = value), colour="white") +
    geom_polygon(data = polygon2, aes(x = Phylum, y = var2, linetype = king), colour = "orange", fill = NA) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
    xlim(c(0.5,87.5)) +
    ylim(c(0, max(forplot$var2) + 5)) +
    coord_polar(theta="x") +
    theme(panel.background=element_blank(),
          axis.title=element_blank(),
          panel.grid=element_blank(),
          axis.text =element_blank(),
          axis.ticks=element_blank()
    ) +
    MyTheme_transparent

ggsave(filename = "./plots_important/step1_mk_phyloseq/phylo_circ_asv_king.pdf", plot = pp3, device = "pdf", width = 8, height = 6, bg = "transparent") 


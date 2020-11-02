# Title: Check TOS taxonomic entries 
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Mon Sep 28 02:27:44 2020

# preparation --------
rm(list = ls())
cat("\014")
setwd("~/Lab/caledna_transect/")

library(dplyr)
library(phyloseq)
library(stringr)

source("./scripts/function_transect.R")
today = format(Sys.Date(), "%Y%m%d")

# def functions --------
# check the TOS score
check_duptos <- function(duptos) {
    duptaxa = unique(duptos[duplicated(duptos$Phylum.Class.Order.Family.Genus.Species),'Phylum.Class.Order.Family.Genus.Species'])
    tos_check = duptos[duptos$Phylum.Class.Order.Family.Genus.Species %in% duptaxa,]
    result = 0
    for (taxa in duptaxa) {
        temp = tos_check[tos_check$Phylum.Class.Order.Family.Genus.Species == taxa,]
        # the shorter version of logical operators are vectorized
        if ((length(unique(temp$Adjusted_TOS)) != 1) | (length(unique(temp$TOS)) != 1)) {
            print(temp)
            result = 1
        }
    }
    return(result)
}

get_duptaxa <- function(df, taxacol) {
    duptaxa = unique(df[duplicated(df[,taxacol]),taxacol])
    outdf = df[df[,taxacol] %in% duptaxa,]
    return(outdf)
}

get_uniqtaxa <- function(df, taxacol) {
    uniqtaxa = as.matrix(table(df[,taxacol]))
    uniqtaxa = rownames(uniqtaxa)[uniqtaxa[,1] == 1]
    outdf = df[df[,taxacol] %in% uniqtaxa,]
    return(outdf)
}

# def variables --------

# load data --------
duptos = read.csv(file = "./derive_data/step2_data_description/gbif_tos/S6_TOS_GBIF.csv", stringsAsFactors = F)
load("./derive_data/phy_deco/phy_deco_all.RData")

# main --------
# change the name of "NA" to ";" for phydeco_all ========
# should be a total of 15848 unique reads that converted all `;NA;` to `;;` levels 
alltaxa = taxa_names(phy_deco_all)
table(duplicated(alltaxa))
alltaxa_table = as.data.frame(alltaxa) %>% 
    dplyr::mutate(alltaxa_nona = stringr::str_replace_all(
        stringr::str_replace_all(
            stringr::str_replace_all(alltaxa, pattern = ";NA;", replacement = ";;"), 
            pattern = ";NA", replacement = ";"), 
        pattern = "NA;", replacement = ";"))

table(duplicated(alltaxa_table$alltaxa_nona))

# generate the site summary ========
# a new otu table based on the 15848 taxonomic entries 
phy_deco_all_bool = tobool(phy_deco_all)
alltaxa_otu = otu_table(phy_deco_all_bool)@.Data %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var = "alltaxa") %>% 
    dplyr::full_join(., y = alltaxa_table, by = "alltaxa") 
# generate a remove duplicate table by for loop (really slow)
sites = sample_names(phy_deco_all_bool)
taxanames_nona  = unique(alltaxa_table$alltaxa_nona)
alltaxa_nona_otu = matrix(data = NA, nrow = length(taxanames_nona), ncol = length(sites),
                              dimnames = list(taxanames_nona, sites))
for (taxa in taxanames_nona) {
        temp = alltaxa_otu %>% 
            dplyr::filter(alltaxa_nona == taxa) %>% 
            dplyr::select(-alltaxa, -alltaxa_nona)
        temp = sign(colSums(temp))
        if(all(names(temp) == sites)) {
            alltaxa_nona_otu[taxa,] = temp
        } else {
            stop("Wrong Sites")
        }
}

# a site summary table get the frequency of occurences 
# https://dplyr.tidyverse.org/articles/rowwise.html
tos_meixi = data.frame(alltaxa_nona_otu) %>% 
    tibble::rownames_to_column(var = "alltaxa_nona") %>% 
    rowwise() %>% 
    dplyr::mutate(total = rowSums(across(where(is.numeric)))) 
tos_meixi2 = tos_meixi %>% dplyr::select(alltaxa_nona, total)

# remove duplicates in previous records ========
# check that the two names are the same ########
setequal(unique(duptos$Phylum.Class.Order.Family.Genus.Species), 
         taxanames_nona)
# YES! 
# check that duplicate records had the same TOS scores ########
check_duptos(duptos) # should return 0 
# YES! 
# compare the values (count of sites) ########
tos_compare = dplyr::left_join(duptos, tos_meixi2, by = c("Phylum.Class.Order.Family.Genus.Species" = "alltaxa_nona"))
tos_dup_comp = get_duptaxa(df = tos_compare, taxacol = "Phylum.Class.Order.Family.Genus.Species")
tos_uniq_comp = get_uniqtaxa(df = tos_compare, taxacol = "Phylum.Class.Order.Family.Genus.Species")
# check if the values for the uniq taxa is the same
table(tos_uniq_comp$Count_of_Samples_Observed_in_eDNA == tos_uniq_comp$total)

# generates stats used in the results (duptos) -------
table(duptos$Adjusted_TOS, useNA = "always")
# 0 0.142857143         0.2 0.285714286 0.333333333 0.428571429 0.666666667 0.714285714 0.857142857           1        <NA> 
# 974        1178          81         112         857        3567         134           5          60        8670        1700  
prop.table(table(duptos$Adjusted_TOS, useNA = "always"))
# 0  0.142857143          0.2  0.285714286  0.333333333  0.428571429  0.666666667  0.714285714  0.857142857            1         <NA> 
# 0.0561771831 0.0679432460 0.0046718191 0.0064597993 0.0494289999 0.2057330719 0.0077286884 0.0002883839 0.0034606068 0.5000576768 0.0980505249 
# 1. omitted entrie (NA entries): 1700 
# 2. proportion of TOS = 0 entries: 5.6%
# 3. proportion of TOS = 1 entries: 50.0%

# linear model 
thislm = lm(Adjusted_TOS ~ Count_of_Samples_Observed_in_eDNA, data = duptos)
summary(thislm)
thiscor = cor(x = duptos$Count_of_Samples_Observed_in_eDNA, y = duptos$Adjusted_TOS, use = "complete.obs", method="pearson")
thiscor

#######################################################################
# NOT USED IN FINAL OUTPUT
#######################################################################
# generates stats used in the results (remove duplicates; NOT USED IN THE FINAL OUTPUT) -------
tos_new = tos_compare %>% 
    dplyr::distinct(Phylum.Class.Order.Family.Genus.Species, .keep_all = TRUE) %>% 
    dplyr::select(-Count_of_Samples_Observed_in_eDNA) %>% 
    dplyr::rename(Count_of_Samples = total)
table(tos_new$Adjusted_TOS, useNA = "always")
# 0 0.142857143         0.2 0.285714286 0.333333333 0.428571429 0.666666667 0.714285714 0.857142857           1        <NA> 
# 878        1084          80         104         778        3313         130           5          55        7944        1477
prop.table(table(tos_new$Adjusted_TOS, useNA = "always"))
# 0  0.142857143          0.2  0.285714286  0.333333333  0.428571429  0.666666667  0.714285714  0.857142857            1         <NA> 
# 0.0554013125 0.0683997981 0.0050479556 0.0065623423 0.0490913680 0.2090484604 0.0082029278 0.0003154972 0.0034704695 0.5012619889 0.0931978799
# 1. omitted entrie (NA entries): 1477
# 2. proportion of TOS = 0 entries: 5.5%
# 3. proportion of TOS = 1 entries: 50.1%

thislm2 = lm(Adjusted_TOS ~ Count_of_Samples, data = tos_new)
summary(thislm2)
plot(thislm2)

write.csv(tos_new, file = paste0("./derive_data/step2_data_description/gbif_tos/S6_TOS_GBIF_removedup_", today, ".csv"))

# cleanup --------


# Title: UCNRS site record cleaning --------
# Author: Meixi Lin
# Date: Fri May 17 15:02:56 2019
# Author:
# Date:
# Modification:

# preparation --------
rm(list = ls())
cat("\014")

setwd("~/Lab/caledna_transect/")
library(dplyr)
date() # the execution date

# Notes before the raw data --------
# The `v0` file was directly from the master table, with the following modification: 
# line `3266` used to be: 
# Hastings Natural History Reserve ,Bird,,Gavia stellata,,Red-throated Loon

# load data -------
# the following data is not altered from the excel
fauna <- read.csv(file = "./raw_data/ucnrs/UC-NRS-Species-List_v0.csv", stringsAsFactors = F)
flora <- read.csv(file = "./raw_data/ucnrs/reserve_plant_list_v0.csv", stringsAsFactors = F)
rsnames <- read.csv(file = "./raw_data/ucnrs/reserve_name_binder.csv", stringsAsFactors = F)

# change the col name --------
colnames(fauna)[1] <- "Reserve"
flora <- flora[, 1:16]
colnames(flora)
colnames(flora) <- c("Reserve", "reserve_count_values", "Genus_Species", "Family", "Genus" , "Species", "ssp_var_designation", "ssp_var_name", "Division", "Reserve_count_formula", "Synonym", "Common_Name", "uncertainties_id","Native_Exotic" , "Habitat_approx", "Bloom_Time_approx")

sort(table(fauna$Reserve))
sort(table(flora$Reserve))

# only get the values for the reserve of interest --------
index <- fauna$Reserve %in% rsnames$rsfauna
fauna <- fauna[index, ]

index <- flora$Reserve %in% rsnames$rsflora
flora <- flora[index, ]

# write out the first version --------
write.csv(fauna, file = "./other_version_data/ucnrs_rec/UC-NRS-Species-List_v1.csv", row.names = F)
write.csv(flora, file = "./other_version_data/ucnrs_rec/reserve_plant_list_v1.csv", row.names = F)

# keep going on validating data --------
# replace the reserve name using edna loc character 
head(fauna$Reserve)
fauna$Reserve <- rsnames[match(fauna$Reserve, rsnames$rsfauna), 'locedna']
head(fauna$Reserve)

head(flora$Reserve)
flora$Reserve <- rsnames[match(flora$Reserve, rsnames$rsflora), 'locedna']
head(flora$Reserve)

# fauna specific columns --------
head(fauna)
table(fauna$Taxon)
fauna$Taxon[which(fauna$Taxon == 'Amphibiann ')] <- "Amphibian"
fauna$Taxon[which(fauna$Taxon == 'bird')] <- "Bird"
fauna$Taxon[which(fauna$Taxon == 'Bird ')] <- "Bird"
table(fauna$Taxon)

table(fauna$Family)
head(sort(table(fauna$Scientific.Name), decreasing = T))
table(fauna$Accepted.Name)
table(fauna$Common.Name)

fauna <- fauna[, -which(colnames(fauna) == "Accepted.Name")]

# change the taxon level a little clear
table(fauna$Taxon) 
fauna$Taxon[fauna$Taxon == 'Amphibian'] <- "Amphibia"
fauna$Taxon[fauna$Taxon == 'Bird'] <- "Aves"
fauna$Taxon[fauna$Taxon == 'Fish'] <- "Actinopterygii"
fauna$Taxon[fauna$Taxon == 'Insect'] <- "Insecta"
fauna$Taxon[fauna$Taxon == 'Mammal'] <- "Mammalia"
fauna$Taxon[fauna$Taxon == 'Reptile'] <- "Reptilia"

Phylum <- vector(length = nrow(fauna))
Phylum[fauna$Taxon == 'Insecta'] <- "Arthropoda"
Phylum[fauna$Taxon != 'Insecta'] <- "Chordata"
fauna <- cbind(fauna, Phylum)
colnames(fauna)
fauna <- fauna[,c(1,6,2:5)]
colnames(fauna) <- c("loc", "Phylum", "Class", "Family", "Species", "Common_Name")
write.csv(fauna, file = "./other_version_data/ucnrs_rec/UC-NRS-Species-List_v2.csv", row.names = F)

# Version 3 stem from visual inspection of the spaced sites etc. 

# flora specific columns --------
colnames(flora)
table(flora$Division)
flora$Division[flora$Division == "Angiosperms"] = "Angiosperm"

flora <- flora[,c(1, 9, 4, 5, 3, 12, 2, 6:8, 10, 11, 13:16)]
colnames(flora)
# drop the original species 
flora <- flora[, -which(colnames(flora) == "Species")]
colnames(flora)[1:6] <- c("loc", "Division", "Family", "Genus", "Species", "Common_Name")
# 
write.csv(flora, file = "./other_version_data/ucnrs_rec/reserve_plant_list_v2.csv", row.names = F)

# clean the white spaces --------
fauna <- apply(X = fauna, MARGIN = 2, FUN = trimws) %>% as.data.frame()
flora <- apply(X = flora, MARGIN = 2, FUN = trimws) %>% as.data.frame()

# find the bad characters --------
fauna <- data.frame(lapply(fauna, function(xx) {xx <- gsub("’s", "'s", xx)}))
flora <- data.frame(lapply(flora, function(xx) {xx <- gsub("’s", "'s", xx)}))

write.csv(fauna, file = "./other_version_data/ucnrs_rec/UC-NRS-Species-List_v3.csv", row.names = F)
write.csv(flora, file = "./other_version_data/ucnrs_rec/reserve_plant_list_v3.csv", row.names = F)

# fill in the empty specie by common name --------
common_nb <- fauna %>% 
    dplyr::select(Species, Common_Name) %>% 
    dplyr::distinct() %>%
    dplyr::filter(Species != "")
index <- match(fauna[fauna$Species == "", "Common_Name"], common_nb$Common_Name)
fauna[fauna$Species == "", "Species"] <- common_nb$Species[index]
write.csv(fauna, file = "./other_version_data/ucnrs_rec/UC-NRS-Species-List_v4.csv", row.names = F)

# load data -------
fauna <- read.csv(file = "./other_version_data/ucnrs_rec/UC-NRS-Species-List_v4.csv", stringsAsFactors = F)
flora <- read.csv(file = "./other_version_data/ucnrs_rec/reserve_plant_list_v3.csv", stringsAsFactors = F)

# fix fauna typo --------
colnames(fauna)
table(fauna$Phylum)
table(fauna$Class)
table(fauna$Family)

# fauna family ###### 
badnames <- c("Acciptridae", "Apida", "Charadriidea", "Columbridae", "Dioptidae", "Embrizidae", "Phalacrocoacidae", "Phalacrocoracide", "Plethodonitadae", "Plethodontide", "Procyoninae", "Sciurdae", "Trochildae", "Trogoldytidae")

goodnames <- c("Accipitridae", "Apidae", "Charadriidae", "Columbidae", "Dipodidae", "Emberizidae", "Phalacrocoracidae", "Phalacrocoracidae", "Plethodontidae", "Plethodontidae", "Procyonidae", "Sciuridae", "Trochilidae", "Troglodytidae")

for (ii in 1:length(badnames)) {
    index <- fauna$Family == badnames[ii]
    fauna$Family[index] <- goodnames[ii]
    fauna$Family[index]
}

# fauna species ########
table(fauna$Species)
# cat(names(table(fauna$Species)), sep = "\n")
# too much didn't dive 
length(table(fauna$Species))
# but fixed double space 
fauna$Species <- sub(pattern = "  ", replacement = " ", fauna$Species)

# fill in fauna Genus ########
tosplit <- fauna[,"Species"] %>% 
    reshape2::colsplit(., pattern = " ", names = c("Genus", "sp"))

# get the genus names 
length(table(tosplit$Genus))
# cat(names(table(tosplit$Genus)), sep = "\n")

# fix flora typo --------
# flora family #######
table(flora$Family)
# cat(names(table(flora$Family)), sep = "\n")

badnames <- c("Amaryllidacea","Amaryllidaecea", "Drypoteridaceae", "Myoperaceae", "Seleginellaceae")

goodnames <- c("Amaryllidaceae", "Amaryllidaceae", "Dryopteridaceae", "Myoporaceae", "Selaginellaceae")

for (ii in 1:length(badnames)) {
    index <- flora$Family == badnames[ii]
    flora$Family[index] <- goodnames[ii]
    flora$Family[index]
}

# flora genus ########
table(flora$Genus)
cat(names(table(flora$Genus)), sep = "\n")

badnames <- c("Anaphallis", "Anemposis", "Chamaeysce", "Custcuta", "Delaria", "Eriastrium", "Koelaria", "Lathryus", "Lepdium", "Linathus", "Lotis", "Mesembryathemum")

goodnames <- c("Anaphalis", "Anemopsis", "Chamaesyce", "Cuscuta", "Delairea", "Eriastrum", "Koeleria", "Lathyrus", "Lepidium", "Linanthus", "Lotus", "Mesembryanthemum")
for (ii in 1:length(badnames)) {
    index <- flora$Genus == badnames[ii]
    flora$Genus[index] <- goodnames[ii]
    flora$Genus[index]
}

# flora species ########
flora$Species <- sub(pattern = "  ", replacement = " ", flora$Species)

# write out the data --------
write.csv(fauna, file = "./derive_data/step2_data_description/ucnrs/UC-NRS-Species-List_cleaned.csv", row.names = F)
write.csv(flora, file = "./derive_data/step2_data_description/ucnrs/reserve_plant_list_cleaned.csv", row.names = F)
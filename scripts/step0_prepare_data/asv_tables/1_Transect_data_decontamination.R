# Title: performing decontamination --------
# Author: Rachel Meyer
asv16S_S<-read.csv("/Users/RachelSMeyer/Desktop/Box Sync/TRANSECT_RESULTS/SHRUB/16S_ASV_sum_by_taxonomy_60_DC_noB.txt", header=TRUE, sep="\t")
asv16S_F<-read.csv("/Users/RachelSMeyer/Desktop/Box Sync/TRANSECT_RESULTS/FOREST/16S_ASV_sum_by_taxonomy_60_DC_noB.txt", header=TRUE, sep="\t")
asv16S_C<-read.csv("/Users/RachelSMeyer/Desktop/Box Sync/TRANSECT_RESULTS/COASTAL/16S_ASV_sum_by_taxonomy_60_DC_noB.txt", header=TRUE, sep="\t")
asv18S_S<-read.csv("/Users/RachelSMeyer/Desktop/Box Sync/TRANSECT_RESULTS/SHRUB/18S_ASV_sum_by_taxonomy_60_DC_noB.txt", header=TRUE, sep="\t")
asv18S_F<-read.csv("/Users/RachelSMeyer/Desktop/Box Sync/TRANSECT_RESULTS/FOREST/18S_ASV_sum_by_taxonomy_60_DC_noB.txt", header=TRUE, sep="\t")
asv18S_C<-read.csv("/Users/RachelSMeyer/Desktop/Box Sync/TRANSECT_RESULTS/COASTAL/18S_ASV_sum_by_taxonomy_60_DC_noB.txt", header=TRUE, sep="\t")
asvFITS_C<-read.csv("/Users/RachelSMeyer/Desktop/Box Sync/TRANSECT_RESULTS/COASTAL/FITS_ASV_sum_by_taxonomy_60_DC_noB.txt", header=TRUE, sep="\t")
asvFITS_F<-read.csv("/Users/RachelSMeyer/Desktop/Box Sync/TRANSECT_RESULTS/FOREST/FITS_ASV_sum_by_taxonomy_60_DC_noB.txt", header=TRUE, sep="\t")
asvFITS_S<-read.csv("/Users/RachelSMeyer/Desktop/Box Sync/TRANSECT_RESULTS/SHRUB/FITS_ASV_sum_by_taxonomy_60_DC_noB.txt", header=TRUE, sep="\t")
asvPITS_S<-read.csv("/Users/RachelSMeyer/Desktop/Box Sync/TRANSECT_RESULTS/SHRUB/PITS_ASV_sum_by_taxonomy_60_DC_noB.txt", header=TRUE, sep="\t")
asvPITS_F<-read.csv("/Users/RachelSMeyer/Desktop/Box Sync/TRANSECT_RESULTS/FOREST/PITS_ASV_sum_by_taxonomy_60_DC_noB.txt", header=TRUE, sep="\t")
asvPITS_C<-read.csv("/Users/RachelSMeyer/Desktop/Box Sync/TRANSECT_RESULTS/COASTAL/PITS_ASV_sum_by_taxonomy_60_DC_noB.txt", header=TRUE, sep="\t")
asvCO1_C<-read.csv("/Users/RachelSMeyer/Desktop/Box Sync/TRANSECT_RESULTS/COASTAL/CO1_ASV_sum_by_taxonomy_60_DC_noB.txt", header=TRUE, sep="\t")
asvCO1_F<-read.csv("/Users/RachelSMeyer/Desktop/Box Sync/TRANSECT_RESULTS/FOREST/CO1_ASV_sum_by_taxonomy_60_DC_noB.txt", header=TRUE, sep="\t")
asvCO1_S<-read.csv("/Users/RachelSMeyer/Desktop/Box Sync/TRANSECT_RESULTS/SHRUB/CO1_ASV_sum_by_taxonomy_60_DC_noB.txt", header=TRUE, sep="\t")

BIOM_F<-read.csv("/Users/RachelSMeyer/Desktop/Box Sync/TRANSECT_RESULTS/FOREST/Forest_metadata_Oct25_2018.txt", header=TRUE, sep="\t")
BIOM_S<-read.csv("/Users/RachelSMeyer/Desktop/Box Sync/TRANSECT_RESULTS/SHRUB/Shrub_metadata_Oct25_2018.txt", header=TRUE, sep="\t")
BIOM_C<-read.csv("/Users/RachelSMeyer/Desktop/Box Sync/TRANSECT_RESULTS/COASTAL/Coastal_metadata_Oct25_2018.txt", header=TRUE, sep="\t")

PITS_S_physeq<-convert_anacapa_to_phyloseq(asvPITS_S,BIOM_S)
PITS_F_physeq<-convert_anacapa_to_phyloseq(asvPITS_F,BIOM_F)
PITS_C_physeq<-convert_anacapa_to_phyloseq(asvPITS_C,BIOM_C)
phy_PITS <- merge_phyloseq(PITS_C_physeq,PITS_F_physeq,PITS_S_physeq)
phy_PITS_3 <- filter_taxa(phy_PITS, function(x) sum(x) > 2, T)

FITS_S_physeq<-convert_anacapa_to_phyloseq(asvFITS_S,BIOM_S)
FITS_F_physeq<-convert_anacapa_to_phyloseq(asvFITS_F,BIOM_F)
FITS_C_physeq<-convert_anacapa_to_phyloseq(asvFITS_C,BIOM_C)
phy_FITS <- merge_phyloseq(FITS_C_physeq,FITS_F_physeq,FITS_S_physeq)
phy_FITS_3 <- filter_taxa(phy_FITS, function(x) sum(x) > 2, T)

a18S_S_physeq<-convert_anacapa_to_phyloseq(asv18S_S,BIOM_S)
a18S_F_physeq<-convert_anacapa_to_phyloseq(asv18S_F,BIOM_F)
a18S_C_physeq<-convert_anacapa_to_phyloseq(asv18S_C,BIOM_C)
phy_18S <- merge_phyloseq(a18S_C_physeq,a18S_F_physeq,a18S_S_physeq)
phy_18S_3 <- filter_taxa(phy_18S, function(x) sum(x) > 2, T)

a16S_S_physeq<-convert_anacapa_to_phyloseq(asv16S_S,BIOM_S)
a16S_F_physeq<-convert_anacapa_to_phyloseq(asv16S_F,BIOM_F)
a16S_C_physeq<-convert_anacapa_to_phyloseq(asv16S_C,BIOM_C)
phy_16S <- merge_phyloseq(a16S_C_physeq,a16S_F_physeq,a16S_S_physeq)
phy_16S_3 <- filter_taxa(phy_16S, function(x) sum(x) > 2, T)

CO1_S_physeq<-convert_anacapa_to_phyloseq(asvCO1_S,BIOM_S)
CO1_F_physeq<-convert_anacapa_to_phyloseq(asvCO1_F,BIOM_F)
CO1_C_physeq<-convert_anacapa_to_phyloseq(asvCO1_C,BIOM_C)
phy_CO1 <- merge_phyloseq(CO1_C_physeq,CO1_F_physeq,CO1_S_physeq)
phy_CO1_3 <- filter_taxa(phy_CO1, function(x) sum(x) > 2, T)

PITS_table<-phyloseq(otu_table(phy_PITS_3, taxa_are_rows=TRUE))
write.csv(PITS_table, file="PITS_table_3_Oct25_2018.csv")

FITS_table<-phyloseq(otu_table(phy_FITS_3, taxa_are_rows=TRUE))
write.csv(FITS_table, file="FITS_table_3_Oct25_2018.csv")

CO1_table<-phyloseq(otu_table(phy_CO1_3, taxa_are_rows=TRUE))
write.csv(CO1_table, file="CO1_table_3_Oct25_2018.csv")

a18S_table<-phyloseq(otu_table(phy_18S_3, taxa_are_rows=TRUE))
write.csv(a18S_table, file="18S_table_3_Oct25_2018.csv")

a16S_table<-phyloseq(otu_table(phy_16S_3, taxa_are_rows=TRUE))
write.csv(a16S_table, file="16S_table_3_Oct25_2018.csv")
###make merged OTU tables without the 3 read minimum too
PITS_table1<-phyloseq(otu_table(phy_PITS, taxa_are_rows=TRUE))
write.csv(PITS_table1, file="PITS_table_nomin_Oct25_2018.csv")

FITS_table1<-phyloseq(otu_table(phy_FITS, taxa_are_rows=TRUE))
write.csv(FITS_table1, file="FITS_table_nomin_Oct25_2018.csv")

CO1_table1<-phyloseq(otu_table(phy_CO1, taxa_are_rows=TRUE))
write.csv(CO1_table1, file="CO1_table_nomin_Oct25_2018.csv")

a18S_table1<-phyloseq(otu_table(phy_18S, taxa_are_rows=TRUE))
write.csv(a18S_table1, file="18S_table_nomin_Oct25_2018.csv")

a16S_table1<-phyloseq(otu_table(phy_16S, taxa_are_rows=TRUE))
write.csv(a16S_table1, file="16S_table_nomin_Oct25_2018.csv")
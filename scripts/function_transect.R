# Title: function used for CALeDNA transect analysis --------
# Author: Meixi Lin
# Credit:
# Date: Wed Jul 25 15:26:47 2018
# Author:
# Date:
# Modification:

# variable definitions --------
# 33 set of continous variable
contlist <- c("Longitude", "hfp", "bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio8", "bio14", "bio15", "phihox", "orcdrc", "cecsol", "sndppt", "bldfie", "ntot", "elev", "Slope", "aspect","CTI", "DAH", "B1", "B4", "B6", "B9", "B10", "B11", "NDVI32", "NBRT", "greenness", "imprv", "ptrcv")
# 56 full set of continuous variable
fullcontlist <- c("Longitude", "Latitude", "hfp",
                  "bio1", "bio2", "bio3", "bio4", "bio5", "bio6",
                  "bio7", "bio8", "bio9", "bio10", "bio11", "bio12",
                  "bio13", "bio14", "bio15", "bio16", "bio17","bio18", "bio19",
                  "phihox", "orcdrc", "cecsol", "sndppt", "clyppt", "bldfie", "ntot",
                  "elev", "Slope", "aspect","CTI", "Roughness", "Ruggedness", "DAH",
                  "B1", "B2", "B3", "B4", "B5", "B6", "B7",
                  "B8", "B8A", "B9", "B10", "B11", "B12",
                  "NDVIS2", "NDVI32", "EVI", "NBRT", "greenness",
                  "imprv", "ptrcv")

# categorical variable list
catlist <- c("loc","ecoregion","majorhab","minorhab","transect","SoS", "clust", "taxousda", "NLCD")
# a list of continuous variables in the metadata layer (used in mantel_test; and linear regression, etc.):
alllevels <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
primers <- c("16S", "18S", "CO1", "FITS", "PITS", "all")
# for community ecology analysis, the "all" metabarcode was not included
primers_commeco <- c("16S", "18S", "CO1", "FITS", "PITS")
rare_depth <- c(2000, 4000, 1000, 4000, 1000, 10000)
names(rare_depth) <- primers
# this value was derived from TIGER arcgis map
CAlimit <- matrix(data = c(-124.482, -114.1312, 32.52885, 42.0095), nrow = 2, ncol = 2,
                  dimnames = list(c("min", "max"), c("Longitude", "Latitude")))
# crs
mycrs = "+proj=longlat +datum=WGS84"
# ggplot2 theme definitions --------
require(ggplot2)

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

# function definitions --------

# custom rarefaction --------
custom_rarefaction_2 <- function(physeq_object, sample_size = 10000, replicates = 10, seeds = FALSE)
{
    reps <- vector("list", replicates)
    for (ii in 1:replicates) {
        myseed <- seeds[ii]
        reps[[ii]] <- phyloseq::rarefy_even_depth(physeq_object,
                                                sample.size = sample_size,
                                                rngseed = myseed)
    }
    dfs <- lapply(reps, function(x) as.data.frame(x@otu_table@.Data))
    dfs <- lapply(dfs, function(x) tibble::rownames_to_column(x,
                                                              var = "taxonomy"))
    dfs <- do.call(rbind.data.frame, dfs)
    otu <- dfs %>% dplyr::group_by(taxonomy) %>% dplyr::summarize_all(dplyr::funs(sum(.)/replicates)) %>%
        dplyr::mutate_if(is.numeric, dplyr::funs(round)) %>%
        data.frame %>% tibble::column_to_rownames("taxonomy") %>%
        as.matrix
    OTU <- phyloseq::otu_table(otu, taxa_are_rows = T)
    TAX <- physeq_object@tax_table
    physeq_to_return <- phyloseq::phyloseq(OTU, TAX)
    physeq_to_return <- phyloseq::merge_phyloseq(physeq_to_return,
                                                 physeq_object@sam_data)
    return(physeq_to_return)
}

# summarize reads --------
# a new way of doing get tax count
# date: 01/31/2019
get_tax_count <- function(phyloseq1, taxlevel = c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
    txcount <- tax_table(phyloseq1)@.Data %>%
        as.data.frame %>%
        # select by taxlevel
        group_by(eval(as.name(taxlevel))) %>%
        # get count for each value
        summarise(count = n())
    colnames(txcount) <- c(taxlevel, "count")
    return(txcount)
}

# construct phyloseq from anacapa asv table only
convert_asv_to_phyloseq <- function (taxon_table, cols_to_categorize = NULL)
{
    # validate_input_files(taxon_table, metadata_file)
    # metadata_file <- categorize_continuous_metadata(metadata_file,
                                                    # cols_to_categorize)
    taxon_table2 <- group_anacapa_by_taxonomy(taxon_table) %>%
        tibble::column_to_rownames("sum.taxonomy") %>% as.matrix
    taxon_table2 <- taxon_table2[, order(colnames(taxon_table2))]
    ana_taxon_table_physeq <- phyloseq::otu_table(taxon_table2,
                                                  taxa_are_rows = TRUE)
    taxon_names <- reshape2::colsplit(rownames(taxon_table2), ";", names = c("Phylum", "Class", "Order", "Family", "Genus", "Species")) %>% as.matrix
    rownames(taxon_names) <- rownames(taxon_table2)
    tax_physeq <- phyloseq::tax_table(taxon_names)
    colnames(tax_physeq) <- c("Phylum", "Class", "Order", "Family",
                              "Genus", "Species")
    physeq_object <- phyloseq::phyloseq(ana_taxon_table_physeq,
                                        tax_physeq)
    return(physeq_object)
}

# aggloremate asv records to specified taxonomy level and return new phyloseq object (Credit:Gaurav)

glom_tax <- function(phyloseq1, taxlevel) {
    glomphy <- data.frame(otu_table(phyloseq1))
    # split sum.taxonomy
    glomphy <- cbind(glomphy, reshape2::colsplit(rownames(glomphy), ";", names = c("Phylum", "Class", "Order", "Family", "Genus", "Species")))

    # change "NA" or "" to "unknown" for all cells
    glomphy <- glomphy %>%
        mutate(Phylum = ifelse(is.na(Phylum) | Phylum == "", "unknown", Phylum)) %>%
        mutate(Class = ifelse(is.na(Class) | Class == "", "unknown", Class)) %>%
        mutate(Order = ifelse(is.na(Order) | Order == "", "unknown", Order)) %>%
        mutate(Family = ifelse(is.na(Family) | Family == "", "unknown", Family)) %>%
        mutate(Genus = ifelse(is.na(Genus) | Genus == "", "unknown", Genus)) %>%
        mutate(Species = ifelse(is.na(Species)| Species == "", "unknown", Species))

    # sumup by standards
    glomphy <- glomphy %>%
        group_by_at(taxlevel) %>%
        # group_by(Phylum) %>%
        summarize_if(is.numeric, sum) %>%
        data.frame %>%
        tibble::column_to_rownames(taxlevel)
    glomphy <- glomphy[which(rowSums(glomphy) > 0),]
    # glomphy[glomphy == 0] <- NA

    # build tax_table for phyloseq
    taxmat <- as.matrix(row.names(glomphy))
    row.names(taxmat) <- row.names(glomphy)
    colnames(taxmat) <- taxlevel

    phyglom.tax <- tax_table(taxmat)
    phyglom.otu <- otu_table(glomphy, taxa_are_rows = T)

    # check if original phyloseq contained sample data slot
    if (is.null(sample_data(phyloseq1, errorIfNULL = F))) {
        phyglom <- phyloseq(phyglom.otu, phyglom.tax)
    } else {
        phyglom <- phyloseq(phyglom.otu, phyglom.tax, sample_data(phyloseq1))
    }
    return(phyglom)
}

glom_tax_df <- function(df, taxlevel = taxlevel) {
    # split sum.taxonomy
    glomphy <- cbind(df[,-1], reshape2::colsplit(df[,1], ";", names = c("Phylum", "Class", "Order", "Family", "Genus", "Species")))

    # change "NA" or "" to "unknown" for all cells
    glomphy <- glomphy %>%
        mutate(Phylum = ifelse(is.na(Phylum) | Phylum == "", "unknown", Phylum)) %>%
        mutate(Class = ifelse(is.na(Class) | Class == "", "unknown", Class)) %>%
        mutate(Order = ifelse(is.na(Order) | Order == "", "unknown", Order)) %>%
        mutate(Family = ifelse(is.na(Family) | Family == "", "unknown", Family)) %>%
        mutate(Genus = ifelse(is.na(Genus) | Genus == "", "unknown", Genus)) %>%
        mutate(Species = ifelse(is.na(Species)| Species == "", "unknown", Species))

    # sumup by standards
    glomphy <- glomphy %>%
        group_by_at(taxlevel) %>%
        # group_by(Phylum) %>%
        summarize_if(is.numeric, sum) %>%
        data.frame %>%
        tibble::column_to_rownames(taxlevel)
    glomphy <- glomphy[which(rowSums(glomphy) > 0),]

    return(glomphy)
}

# utility --------
tobool <- function(xx) {
    transform_sample_counts(xx, function(abund) 1*(abund>0))
}

loadRData <- function(fileName){
    #loads an RData file, and returns it
    load(fileName)
    output = get(ls()[ls() != "fileName"])
    return(output)
}

myveganifyOTU <- function (physeq)
{
    if (taxa_are_rows(physeq)) {
        physeq <- t(physeq)
    }
    return(as(otu_table(physeq), "matrix"))
}

named_num2df <- function(named_num, yourcols) {
    mydf <- named_num %>%
        as.data.frame() %>%
        tibble::rownames_to_column()
    if (ncol(mydf) != length(yourcols)) {
        print("colnames not match. check!")
        break;
    }
    colnames(mydf) <- yourcols
    return(mydf)
}

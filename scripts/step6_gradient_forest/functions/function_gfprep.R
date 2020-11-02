# Title: functions for preparing for gradient forest --------
# Author: Meixi Lin
# Date: Sat Jun  8 18:52:59 2019
# Author:
# Date:
# Modification:

# continuous variables to include ---------
# see 0_step4_gradient_forest_determine_var_06062019.R for description 

myvar <- c("Longitude", "hfp", "bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio8", "bio14", "bio15", "phihox", "orcdrc", "cecsol", "sndppt", "bldfie", "ntot", "elev", "Slope", "aspect","CTI", "DAH", "B1", "B4", "B6", "B9", "B10", "B11", "NDVI32", "NBRT", "greenness", "imprv", "ptrcv")

# clean up the data for gradient forest input --------
# Input: phyloseq object ########
gfprep <- function(physeq, glomtax = alllevels, metavar = myvar,
                   subsettax = NULL, subsettarget = NULL,
                   bothok = F, droptopo = F, abun = 4, cutoff = 0.05) 
{
    # subset taxa (e.g. only gf on arthropods) 
    if (!is.null(subsettax) & !is.null(subsettarget)) {
        physeq <- subset_taxa(physeq = physeq, eval(as.name(subsettax)) == subsettarget)
    }
    
    # prepare phyloseq 
    if (bothok == T) {
        physeq <- subset_samples(physeq, bothok == T)
    }
    if (glomtax %in% alllevels[1:5]) {
        physeq <- glom_tax(physeq, taxlevel = glomtax)
    } else {
        # keep the original physeq
    }
    
    # get Phys_sites
    # if drop topographical variable 
    if (droptopo == T) {
        topovar = c("aspect", "CTI", "DAH")
        metavar = metavar[!metavar %in% topovar]
    }
    
    Phys_site <- data.frame(sample_data(physeq))
    print(dim(Phys_site))
    # only keep certain variables for regression 
    # only keep the sites where all variables are available 
    Phys_site <- Phys_site[, eval(metavar)] %>%
        tidyr::drop_na()
    
    print(colnames(Phys_site))
    print(dim(Phys_site))
    
    # get Sp_mat 
    # transform to presence/absence
    physeq1 <- prune_samples(samples = rownames(Phys_site), x = physeq)
    physeq1 <- transform_sample_counts(physeq1, function(abund) {1*(abund > abun)})
    if (0 %in% taxa_sums(physeq1)) {
        physeq1 <- prune_taxa(taxa = taxa_sums(physeq1) > 0, x = physeq1)
    }
    if (0 %in% sample_sums(physeq1)) {
        physeq1 <- prune_samples(samples = sample_sums(physeq1) > 0, x = physeq1)
    }
    # get rid of the too abundant/sparse species 
    # get cutoff (should be both sides so "cutoff / 2" )
    # MODIFICATION FROM MAY 2019 COMMIT 
    lowtaxa <- floor(length(sample_names(physeq1)) * cutoff / 2)
    hightaxa <- length(sample_names(physeq1)) - lowtaxa
    physeq1 <- prune_taxa(taxa_sums(physeq1) > lowtaxa, physeq1)
    physeq1 <- prune_taxa(taxa_sums(physeq1) < hightaxa, physeq1)
    Sp_mat <-  t(as.matrix(otu_table(physeq1)@.Data))
    # convert the numeric value to factor, so that `randomForest` will treat the response as classification: 
    # y: A response vector. If a factor, classification is assumed, otherwise regression is assumed. If omitted, randomForest will run in unsupervised mode.
    Sp_mat1 <- apply(Sp_mat, 2, function(x) {factor(as.logical(x))})
    row.names(Sp_mat1) <- row.names(Sp_mat)
    # get rid of the unknown column in taxonomy classification
    unknown <- which(colnames(Sp_mat1) == "unknown")
    if (length(unknown) > 0) {
        Sp_mat <- Sp_mat1[,-unknown] # here reassign the y back to Sp_mat, not using Sp_mat2 anymore 
    } else {
        Sp_mat <- Sp_mat1
    }
    if (dim(Phys_site)[1] > dim(Sp_mat)[1]) {
        id <- rownames(Phys_site) %in% rownames(Sp_mat) 
        Phys_site = Phys_site[id,]
    }
    
    # get the lev needed 
    print(dim(Sp_mat))
    lev <- gf_lev(Sp_mat)
    
    # return the list 
    result <- list(Phys_site, Sp_mat, lev)
    return(result)
}

# (gradient forest get level) =================================
gf_lev <- function(Sp_mat) {
    # get some dimension stuff
    nSites <- dim(Sp_mat)[1]
    lev <- floor(log2(nSites * 0.368/2))
    return(lev)
}

# (gradient forest default wrapper) =================================
gfdefault <- function(Phys_site, Sp_mat, ntrees, lev, seed = NULL) {
    ## set the seed
    if (!is.null(seed)) {set.seed(seed)}
    gf <- gradientForest(cbind(Phys_site, Sp_mat), 
                         predictor.vars = colnames(Phys_site), response.vars = colnames(Sp_mat),
                         ntree = ntrees, transform = NULL, compact = F,
                         maxLevel = lev, corr.threshold = 0.5)
    return(gf)
}

# (gradient forest get two parameters) =================================
gf_res <- function(gf) {
    species.pos <- gf$species.pos.rsq
    ave.rsq <- sum(gf$result)/species.pos
    result <- c(species.pos, ave.rsq)
    names(result) <- c("species.pos", "ave.rsq")
    return(result)
}

# (gradient forest to phyloseq) ===================================
gf2phy <- function(gf, taxlevel = "Family") {
    otu <- gf$Y %>%
        as.data.frame() %>%
        dplyr::mutate_if(is.factor, as.numeric)
    otu <- otu - 1
    rownames(otu) <- rownames(gf$Y)
    gfotu <- otu_table(otu, taxa_are_rows = F)
    sample <- gf$X %>% as.data.frame() %>%
        tibble::rownames_to_column(var = "MatchName")
    rownames(sample) <- rownames(gf$X)
    gfsample <- sample_data(sample)
    tax <- matrix(colnames(gf$Y), 
                  nrow = ncol(gf$Y), ncol = length(taxlevel))
    rownames(tax) <- tax
    colnames(tax) <- taxlevel
    gftax <- tax_table(tax)
    gfphy <- phyloseq(gfotu, gftax, gfsample)
    return(gfphy)
}

gf2phydup <- function(gf, taxlevel = "Family", dup) {
    otu <- gf$Y %>%
        as.data.frame() %>%
        dplyr::mutate_if(is.factor, as.numeric)
    otu <- otu - 1
    rownames(otu) <- rownames(gf$Y)
    gfotu <- otu_table(otu, taxa_are_rows = F)
    sample <- gf$X %>% as.data.frame() %>%
        tibble::rownames_to_column(var = "MatchName") %>%
        dplyr::mutate(indup = ifelse(MatchName %in% dup, MatchName, NA))
    rownames(sample) <- rownames(gf$X)
    gfsample <- sample_data(sample)
    tax <- matrix(colnames(gf$Y), 
                  nrow = ncol(gf$Y), ncol = length(taxlevel))
    rownames(tax) <- tax
    colnames(tax) <- taxlevel
    gftax <- tax_table(tax)
    gfphy <- phyloseq(gfotu, gftax, gfsample)
    return(gfphy)
}

# more summary stats ---------
print_gf_summary <- function(gf) {
    print(gf)
    print(summary(gf))
    most_imp <- importance(gf)[1:10]
    print(most_imp)
    rsq <- gf_res(gf)
    print(rsq)
    imp_families <- sort(gf$result, decreasing = T)[1:10]
    print(imp_families)
    return("End of summary")
}

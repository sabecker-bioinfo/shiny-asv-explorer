library(shiny)
library(shinydashboard)
library(data.table)
library(DT)
library(magrittr)
library(ggplot2)
library(ggtree)
library(plotly)
library(genefilter)
library(dada2)
library(phyloseq)
library(vegan)
library(config)
library(DESeq2)
library(dplyr)

library(DECIPHER)
library(phangorn)

# converts empty strings, "none" and "NULL" to NULL
av <- function(x) {
    if(isTRUE(all.equal(x, "")) | isTRUE(all.equal(x, "none")) | isTRUE(all.equal(x, "NULL"))) {
        return(NULL)
    } else {
        return(x)
    }
}

# variable-to-facet-formula conversion function for facet_grid
get_facet_grid <- function(facetrow = NULL, facetcol = NULL) {
    if(is.null(av(facetrow)) & is.null(av(facetcol))) {
        return(NULL)
    } else if(is.null(av(facetcol))) {
        formstring <- paste(paste(facetrow, collapse = "+"), "~", ".")
    } else {
        formstring <- paste(
            paste(facetrow, collapse = "+"),
            "~",
            paste(facetcol, collapse = "+")
        )
    }
    return(as.formula(formstring))
}

compute_phylogenetic_tree <- function(physeq) {
    tt <- tax_table(physeq)
    sequences <- data.frame(row.names(tt@.Data), tt@.Data[,'SVhash'], tt@.Data[,'Phylum'])
    colnames(sequences) <- c('sequence', 'SVhash', "Phylum")
    row.names(sequences) <- paste0(sequences[,'Phylum'], '_', sequences[,'SVhash'])
    
    # it appears that the abundance is needed for the dada2 getSequences function, but is not actually used
    # seqs is just a named list and no data about abundance appears (just names and ASV sequences)
    # this could probably be anything (fixed value, NA, etc.)
    # currently defined as the total read count for each ASV across the selected samples
    sequences[,'abundance'] <- rowSums(otu_table(physeq)@.Data)
    
    seqs <- getSequences(sequences)
    names(seqs) <- row.names(sequences)
    
    alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)
    phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
    dm <- dist.ml(phang.align)
    treeNJ <- NJ(dm)
    
    return(treeNJ)
}

# supported ordination methods
ordlist <- sort(c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA"))

# supported distances
distlist <- c("Unifrac" = "unifrac", "Weighted Unifrac" = "wunifrac", "DPCoA" = "dpcoa", "Jensen-Shannon Divergence" = "JSD", "Jaccard" = "jaccard", "Gower" = "gower", "Gower Alt" = "altGower", "Morisita" = "morisita", "Horn" = "horn", "Bray-Curtis" = "bray", "Kulczynski" = "kulczynski", "Mountford" = "mountford", "Raup" = "raup", "Binomial" = "binomial", "Chao" = "chao", "Cao" = "cao", "Manhattan" = "manhattan", "Euclidean" = "euclidean", "Canberra" = "canberra", "w" = "w", "Minkowski" = "minkowski")
distlist <- distlist[order(names(distlist))]
distlist <- as.list(distlist)

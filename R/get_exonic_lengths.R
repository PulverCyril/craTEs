#' Returns the cumulated exonic length for each gene
#'
#' @param org organism
get_exonic_lengths <- function(org='hg19') {
    org <- GenomicFeatures::match.arg(org)
    ens <- GenomicFeatures::makeTxDbFromUCSC(genome=org,
                            tablename='ensGene',
                            url = 'http://genome-euro.ucsc.edu/cgi-bin/'
                           )
    exonic <- GenomicFeatures::exonsBy(ens, by='gene')
    red.exonic <- purrr::reduce(exonic)
    exon.lengths <- sum(width(red.exonic))
    # Removing ENSEMBL gene versions
    names(exon.lengths) <- unlist(lapply(strsplit(names(exon.lengths), '\\.'), function(x) x[[1]]))
    return(exon.lengths)
}

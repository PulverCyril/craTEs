#' Returns the cumulated exonic length for each gene
#'
#' @param org organism
get_exonic_lengths <- function(org='hg19') {
    org <- match.arg(org)
    exonic <- NULL
    if (org == 'hg19') {
        ens <- GenomicFeatures::makeTxDbFromUCSC(genome=org,
                                tablename='ensGene',
                                url = 'http://genome-euro.ucsc.edu/cgi-bin/'
                               )
        exonic <- GenomicFeatures::exonsBy(ens, by='gene')

        } else if (org == 'hg38') {
            exonic <- GenomicFeatures::exonsBy(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86, by='gene')
        }
    red.exonic <- GenomicRanges::reduce(exonic)
    exon.lengths <- sum(GenomicRanges::width(red.exonic))
    # Removing ENSEMBL gene versions
    names(exon.lengths) <- unlist(lapply(strsplit(names(exon.lengths), '\\.'), function(x) x[[1]]))
    return(exon.lengths)
}

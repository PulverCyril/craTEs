#' computes TPMs from a table of counts. Use get_tpm_local when possible to save time
#'
#' The count table should contain genes as rows and samples as columns
#'
#' @param x table with raw counts, ensembl id as rownames
#' @param org organism, one of c('hg19', 'mm9', 'rheMac8')
#' @return a TPM matrix with ensembl ID as rownames, samples as columns
get_tpm <- function(x, org=c('hg19', 'mm9', 'rheMac8')) {
    countTable = x
    org <- match.arg(org)
    ens <- GenomicFeatures::makeTxDbFromUCSC(genome=org,
                            tablename='ensGene',
                            url = 'http://genome-euro.ucsc.edu/cgi-bin/',
                            #goldenPath_url = 'http://hgdownload.cse.ucsc.edu/goldenPath',
                            )
    exonic <- GenomicFeatures::exonsBy(ens, by='gene')
    red.exonic <- purrr::reduce(exonic)
    exon.lengths <- sum(width(red.exonic))
    names(exon.lengths) <- unlist(lapply(strsplit(names(exon.lengths), '\\.'), function(x) x[[1]]))
    exon.lengths.o <- exon.lengths[rownames(countTable)]
    length_norm <- countTable / exon.lengths.o
    ans <- t(t(length_norm)/colSums(length_norm)) * 1e6 # see https://rpubs.com/bbolker/sweep_divide
    return(ans)
}

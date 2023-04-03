#' Computes TPM from a count table, stores a local copy of the exonic lengths to avoid extra computation
#'
#' @param countTable table with ensembl gene IDs as row names, sample IDs as columns
#' @param org organism (defaults to 'hg19')
#' @param exonic_lengths_path path to where computed exonic lengths are stored, to avoid recomputing them everytime
#' @param redownload_exons forces the download and computation of exon lengths from scratch
#' @return matrix of TPMs with ensembl id as rows and samples as columns 
get_tpm_local <- function(countTable, org='hg19', exonic_lengths_path=NULL, redownload_exons=FALSE) {
    # Defensive programming
    stopifnot(is.matrix(countTable))

    if(is.null(exonic_lengths_path)) {
        exonic_lengths_path = paste0('../data/temp/', org, '_exon_lengths.RDS')
    }
    
    # checking for the existence of a previously generated exon_length vector
    
    if(!dir.exists(dirname(exonic_lengths_path))){
        dir.create(dirname(exonic_lengths_path), showWarnings = FALSE)
    }
    
    exonic_lengths = NULL
    if(!file.exists(exonic_lengths_path) | redownload_exons) {
        file.remove(exonic_lengths_path, showWarnings = FALSE)
        exonic_lengths = get_exonic_lengths(org)
        saveRDS(exonic_lengths, exonic_lengths_path)
    } else {
        exonic_lengths = readRDS(exonic_lengths_path)
    }

    # normalizing to TPM
    
    # filtering of genes present in the countTable is done at this stage, therefore the full table should be stored somewhere.
    exonic_lengths_filtered <- exonic_lengths[rownames(countTable)]
    
    print(head(exonic_lengths_filtered))
    
    # normalizing for gene length
    length_norm <- countTable / exonic_lengths_filtered
    
    # normalizing for library size
    ans <- t(t(length_norm)/colSums(length_norm)) * 1e6 # see https://rpubs.com/bbolker/sweep_divide
    return(ans)
}

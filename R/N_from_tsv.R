#' fast loading of N from a .tsv, using data.table
#' 
#' loads N specified in a .tsv file, and returns it as a matrix ready for use in preprocess_E_N_for_activities()
#'
#' @param filename path to tsv file specifying N, with genes as rows and TE subfamilies as columns
#' @return matrix N
#' @export
N_from_tsv <- function(filename) {
	N = data.table::fread(filename, header=T, sep='\t', stringsAsFactors = F, data.table=FALSE)
	rownames(N) = N$V1
	N$V1 = NULL
	colnames(N) = gsub(pattern = '-', replacement = '.', x = colnames(N))
	return(as.matrix(N))
}

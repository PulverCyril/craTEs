#' Centers the matrix X by rows
#'
#' @param X matrix with rows to be centered.
#' @return row-centered matrix X
#' @export
row_center <- function(X) {
        stopifnot(is.matrix(X))

    return (t(scale(x = t(X), center = T, scale = F)))
}

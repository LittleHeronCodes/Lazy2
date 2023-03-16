#' removeNAs
#'
#' Remove NA values from a vector
#'
#' @param x input vector
#' @return NA removed vectors
#' @export
#' @examples
#' a <- c(1, 2, 3, NA, 6, NA, 7)
#' removeNAs(a)
removeNAs <- function(x) {
    x[which(!is.na(x))]
}


#' countUni
#'
#' count unique values in a vector
#'
#' @param x input vector
#' @return no of unique values in x
#' @export
#' @examples
#' a <- c(1, 2, 3, NA, 6, NA, 7, 7, 1)
#' countUni(a)
countUni <- function(x) {
    length(unique(x))
}

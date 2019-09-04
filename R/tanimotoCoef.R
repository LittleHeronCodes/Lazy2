#' Tanimoto coefficients
#' 
#' Calculates Tanimoto coefficients given two set vectors
#' @param A set input1 as vectors
#' @param B set input2 as vectors
#' @export
#' @examples
#' A = c(1:5)
#' B = c(3:7)
#' tanimotoCoef(A, B)

tanimotoCoef = function(A, B) {
	int = intersect(A, B)
	uni = union(A, B)
	return(length(int) / length(uni))
}

#' Enrichment factor
#' 
#' Calculates enrichment factor against random overlap number
#' @param A set input1 as vectors
#' @param B set input2 as vectors
#' @param total total set as vectors
#' @export
#' @examples
#' A = c(1:5)
#' B = c(3:7)
#' total = c(1:10)
#' enrichment(A, B)

enrichment <- function(A, B, total) {
	A = intersect(A, total)
	B = intersect(B, total)
	expected = length(A) * length(B) / length(total)
	observed = length(intersect(A, B))
	if(expected == 0) return(NA)
	return(observed / expected)
}

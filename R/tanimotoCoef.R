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
#' @param setA set input1 as vectors
#' @param setB set input2 as vectors
#' @param total total set as vectors
#' @param psc pseudocount for when expected count is too low.
#' @export
#' @examples
#' A = c(1:5)
#' B = c(3:7)
#' total = c(1:10)
#' getEnrichmentFactor(A, B, total)

getEnrichmentFactor <- function(setA,setB,setT, psc=0) {
	ef = NA
	A = length(setA)
	B = length(setB)
	T = length(setT)
	I = length(intersect(setA, setB))
	if((A / T * B) != 0) ef = (I+psc) / (A / T * B + psc)
	return(ef)
}

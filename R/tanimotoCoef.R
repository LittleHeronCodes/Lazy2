#' Tanimoto coefficients
#' 
#' Calculates Tanimoto coefficients given two set vectors
#' @param A set input1 as vectors
#' @param B set input2 as vectors
#' @param T background space (if given, will intersect set before calculating)
#' @export
#' @examples
#' A = c(1:5)
#' B = c(3:7)
#' tanimotoCoef(A, B)

tanimotoCoef = function(A, B, T=NULL) {
	if(!is.null(T)) {
		A = intersect(A, T)
		B = intersect(B, T)
	}
	int = intersect(A, B)
	uni = union(A, B)
	return(length(int) / length(uni))
}


#' Enrichment factor
#' 
#' Calculates enrichment factor against random overlap number
#' @param setA set input1 as vectors
#' @param setB set input2 as vectors
#' @param setT total set as vectors
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


#' hypergeoTest
#'
#' Hypergeometric test for set inputs (because I'm an idiot)
#' @param query     query set        (balls drawn)
#' @param reference reference set    (white balls)
#' @param bgspace   background space (balls in urn)
#' @export

hypergeoTest <- function(query, reference, bgspace) {
	query = intersect(query, bgspace)
	reference = intersect(reference, bgspace)

	N = length(bgspace)							# no of balls in urn
	k = length(query)							# no of balls drawn from urn (DEG no)
	m = length(reference) 						# no of white balls in urn
	q = length(intersect(reference, query))		# no of white balls drawn

	pVal = phyper(q-1, m, N-m, k, lower.tail = FALSE)
	odds = (q / k) / (m / N)
	return(data.frame(pVal = pVal, oddsRatio=odds, int=q, bg=m))
}


#' Harmonic Mean
#'
#' Calculate harmonic mean
#' @param v numeric vector
#' @export
#' @examples
#' v <- seq(0.01,.99, .01)
#' harmean(v)

harmean <- function(v, ...) { 1/mean(1/v, ..., na.rm = T) }



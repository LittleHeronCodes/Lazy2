#' Tanimoto coefficients
#'
#' Calculates Tanimoto coefficients (Jaccard Index) given two set vectors.
#' @param A set input1 as vectors
#' @param B set input2 as vectors
#' @param T background space (if given, A and B are intersected with T before calculating, causing the function to run slower)
#' @examples
#' A <- c(1:5)
#' B <- c(3:7)
#' tanimotoCoef(A, B)
#' overlapCoef(A, B)
#' @export

tanimotoCoef <- function(A, B, T = NULL) {
	if (!is.null(T)) {
		A <- intersect(A, T)
		B <- intersect(B, T)
	}

	length(intersect(A, B)) / length(union(A, B))
}

#' @describeIn tanimotoCoef
#' Calculate Overlap coefficients.
#' @export

overlapCoef <- function(A, B, T = NULL) {
	if (!is.null(T)) {
		A <- intersect(A, T)
		B <- intersect(B, T)
	}
	
	length(intersect(A, B)) / min(length(A), length(B))
}

#' Calculate fold enrichment (previously enrichment factor)
#'
#' Fold enrichment represents the "enrichment" of overlap between two sets in a given space.
#' This is calculated as (observed count) / (expected count).
#' pseudocount is added for when expected count is too low to add a 'weight'.
#' @param setA set input1 as vectors
#' @param setB set input2 as vectors
#' @param setT total set as vectors
#' @param psc pseudocount for when expected count is too low.
#' @examples
#' A <- c(1:5)
#' B <- c(3:7)
#' total <- c(1:10)
#' getEnrichmentFactor(A, B, total)
#' @export

get_fold_enrichment <- function(setA, setB, setT, psc = 0) {
	setA <- intersect(setA, setT)
	setB <- intersect(setB, setT)

	ef <- NA
	A <- length(setA)
	B <- length(setB)
	T <- length(setT)
	I <- length(intersect(setA, setB))
	if ((A / T * B) != 0) {
		ef <- (I + psc) / (A / T * B + psc)
	}
	ef
}

#' @describeIn get_fold_enrichment
#' Draw volcano plot (deprecated)
#' @export

getEnrichmentFactor <- function(...) {
	.Deprecated("get_fold_enrichment")
	get_fold_enrichment(...)
}

#' Hypergeometric test
#'
#' Run one-tailed hypergeometric test (fisher test) with set inputs.
#' 'hypergeoTest' uses hypergeometric distribution function phyper with lower.tail = FALSE fixed.
#' This is commonly used for testing significance in overlap between sets.
#'
#' @param query     query set        (balls drawn)
#' @param reference reference set    (white balls)
#' @param bgspace   background space (balls in urn)
#' @return Returns a dataframe of 4 columns. pVal is p value returned by phyper.
#' oddsRatio gives number of observation by number of expectation. int and bg refers to
#' the number of intersected items and background space respectively.
#' By nature of phyper, reference and query input should be interchangeable.
#' @export

hypergeoTest <- function(query, reference, bgspace) {
	query <- intersect(query, bgspace)
	reference <- intersect(reference, bgspace)

	N <- length(bgspace)                     # no of balls in urn (gene space)
	k <- length(query)                       # no of balls drawn from urn (DEG no)
	m <- length(reference)                   # no of white balls in urn (reference in gene space)
	q <- length(intersect(reference, query)) # no of white balls drawn (reference in DEG)

	pVal <- phyper(q - 1, m, N - m, k, lower.tail = FALSE)
	odds <- (q / k) / (m / N)
	
	data.frame(pVal = pVal, oddsRatio = odds, int = q, bg = N)
}


#' Harmonic Mean
#'
#' Calculate harmonic mean. Harmonic means are often used to calculate means of p-values.
#' @param v numeric vector
#' @param ... further arguments to be passed to mean
#' @export
#' @examples
#' v <- seq(0.01, .99, .01)
#' harmean(v)
harmean <- function(v, ...) {
	1 / mean(1 / v, ...)
}

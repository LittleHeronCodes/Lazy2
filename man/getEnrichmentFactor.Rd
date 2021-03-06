% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tanimotoCoef.R
\name{getEnrichmentFactor}
\alias{getEnrichmentFactor}
\title{Calculate Enrichment factor}
\usage{
getEnrichmentFactor(setA, setB, setT, psc = 0)
}
\arguments{
\item{setA}{set input1 as vectors}

\item{setB}{set input2 as vectors}

\item{setT}{total set as vectors}

\item{psc}{pseudocount for when expected count is too low.}
}
\description{
Enrichment factor represents the "enrichment" of overlap between two sets in a given space.
This is calculated as (observed count) / (expected count).
}
\examples{
A = c(1:5)
B = c(3:7)
total = c(1:10)
getEnrichmentFactor(A, B, total)
}

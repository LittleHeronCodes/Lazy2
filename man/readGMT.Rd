% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GMT_enrichments.R
\name{readGMT}
\alias{readGMT}
\alias{writeGMT}
\title{readGMT}
\usage{
readGMT(gmtfile, as.df = FALSE)

writeGMT(file, glist, geneset_desc = "")
}
\arguments{
\item{as.df}{Return as data frame?}

\item{file}{GMT file path}

\item{glist}{List of gene set. Should be un-nested level one named list.}

\item{geneset_desc}{Description meta information for gene set. Either length one or same length vector as glist.}
}
\value{
Gene set dataframe of 2 column or list
}
\description{
gmt file reader. A GMT file format is a tab delimited text file containing gene sets.
Each line should contain one geneset, delimited by tab. Genes should start from 3rd field.
}
\section{Functions}{
\itemize{
\item \code{writeGMT}: write gmt file for geneset list.
}}

\examples{
\dontrun{
readGMT('file/path.gmt')
}
}

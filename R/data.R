#' Human Gene Information
#'
#' Human gene ID dataframe containing Entrez, HGNC ID, Gene name, Gene symbol, Ensembl ID, gene type
#' Used for ent2sym and more gene mapper functions.
#' Last update : 2021-11-11 HGNC, build date 2021-11-25
#'
#' @docType data
#' @usage data(LazygeneInfo)
#' @format A dataframe with 41722 rows and 6 variables:
#' \describe{
#'   \item{entrez}{NCBI Entrez ID, Used in ent2sym}
#' 	 \item{hgnc_id}{HGNC ID}
#' 	 \item{hgnc_gene}{HGNC gene name}
#' 	 \item{hgnc_symbol}{HGNC approved gene symbol, Used in ent2sym}
#' 	 \item{ensembl}{Ensembl ID}
#' 	 \item{gene_type}{gene type}
#' }
#'
#' @keywords genes
#' @source \href{https://www.gencodegenes.org/human/}{Gencode Metadata}
#' @examples
#' data(LazygeneInfo)
"LazygeneInfo"

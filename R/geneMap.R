#' Gene Mapping dataframe for Lazy2
#'
#' Gene mapper dataframe containing Entrez, Symbol, Ensembl (optional)
#' used for ent2sym and more gene mapper functions.
#' Needs raw-data check for updates
#'
#' @docType data
#' @usage data(geneMap)
#' @format A dataframe of columns 'hgnc_symbol', 'entrez'
#' @keywords genes
#' @source \href{ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz}{NCBI human gene info}
#' @examples
#' data(geneMap)
"geneMap"



#' Gene symbol alias dataframe for Lazy2
#'
#' Gene mapper dataframe containing Entrez, Symbol, Ensembl (optional)
#' used for ent2sym and more gene mapper functions.
#' Needs raw-data check for updates
#'
#' @docType data
#' @usage data(geneAlias)
#' @format A dataframe of columns 'hgnc_symbol', 'entrez'
#' @keywords genes
#' @source \href{ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz}{NCBI human gene info}
#' @examples
#' data(geneAlias)
"geneAlias"


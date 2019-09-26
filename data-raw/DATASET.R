## code to prepare `DATASET` dataset goes here

## Gene Info processing from HGNC raw file

LazygeneInfo <- read.table('data-raw/gene_info.csv', sep = ',', header=TRUE)
LazygeneInfo$entrez = as.character(LazygeneInfo$entrez)
LazygeneInfo$hgnc_id = as.character(LazygeneInfo$hgnc_id)

usethis::use_data(LazygeneInfo)


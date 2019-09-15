## code to prepare `DATASET` dataset goes here

## Gene Info processing from HGNC raw file

LazygeneInfo <- read.table('data-raw/gene_info.csv', sep = ',', header=TRUE)

usethis::use_data(LazygeneInfo)


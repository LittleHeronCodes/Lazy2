## Gene Info processing from HGNC raw file

LazygeneInfo <- read.table('raw-data/gene_info.csv', sep = ',', header=TRUE)

usethis::use_data(LazygeneInfo)


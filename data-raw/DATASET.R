## code to prepare dataset goes here

library(data.table)
library(dplyr)
library(tidyr)

## Gene Info processing from HGNC raw file

LazygeneInfo <- read.table('data-raw/gene_info.csv', sep = ',', header=TRUE)
LazygeneInfo$entrez = as.character(LazygeneInfo$entrez)
LazygeneInfo$hgnc_id = as.character(LazygeneInfo$hgnc_id)

usethis::use_data(LazygeneInfo)


##  Gene info data processed from NCBI
## download link ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz

NCBIgeneInfo <- data.table::fread('data-raw/NCBIHomo_sapiens.gene_info.csv')
NCBIgeneInfo <- dplyr::rename(NCBIgeneInfo, 
	Entrez = GeneID, 
	chr = chromosome,
	gene_type = type_of_gene,
	HGNC_symbol = Symbol_from_nomenclature_authority,
	HGNC_name   = Full_name_from_nomenclature_authority)
NCBIgeneInfo <- NCBIgeneInfo[,c('Entrez', 'Symbol','Synonyms','dbXrefs', 'chr','gene_type','HGNC_symbol','HGNC_name')]


## gene alias map
geneAlias <- NCBIgeneInfo[,c('Entrez', 'Symbol', 'Synonyms')]
geneAlias <- tidyr::separate_rows(geneAlias, Synonyms, sep = '\\|')

# remove ambiguous synonyms
idx = with(geneAlias, !(duplicated(Synonyms) | duplicated(Synonyms,fromLast=TRUE)) | (Synonyms == '-'))
geneAlias = geneAlias[idx,]

## ent2sym map
geneMap <- NCBIgeneInfo[,c('Entrez','Symbol','HGNC_symbol','gene_type')]


usethis::use_data(geneMap)
usethis::use_data(geneAlias)




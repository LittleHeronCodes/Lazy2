## code to prepare dataset goes here

library(data.table)
library(tidyverse)

## Gene Info processing from HGNC raw file
## Downloaded from HGNC genenames.org/download/custom
hgnc = fread('data-raw/HGNC_Gene_Info.csv', quote='')
colnames(hgnc) = gsub(' ', '_', colnames(hgnc))
geneInfo = hgnc %>% 
  filter(Status == 'Approved') %>%
  select(NCBI_Gene_ID, HGNC_ID, Approved_name, Approved_symbol, Ensembl_gene_ID, Locus_group) %>%
  rename(entrez=NCBI_Gene_ID, hgnc_id = HGNC_ID, hgnc_gene = Approved_name, 
    hgnc_symbol = Approved_symbol, ensembl = Ensembl_gene_ID, gene_type = Locus_group) %>%
  mutate(hgnc_id = gsub('^HGNC:','',hgnc_id), entrez = as.character(entrez)) %>%
  arrange(entrez, hgnc_id)

LazygeneInfo = geneInfo

usethis::use_data(LazygeneInfo)

#   entrez hgnc_id                          hgnc_gene hgnc_symbol           gene_type
# 1      1       5             alpha-1-B glycoprotein        A1BG protein-coding gene
# 2      2       7              alpha-2-macroglobulin         A2M protein-coding gene
# 3      3       8 alpha-2-macroglobulin pseudogene 1       A2MP1          pseudogene
# 4     11      15     N-acetyltransferase pseudogene        NATP          pseudogene
# 5     12      16           serpin family A member 3    SERPINA3 protein-coding gene
# 6     13      17          arylacetamide deacetylase       AADAC protein-coding gene


##  Gene info data processed from NCBI
## download link ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz

# NCBIgeneInfo <- data.table::fread('data-raw/NCBIHomo_sapiens.gene_info.csv')
# NCBIgeneInfo <- dplyr::rename(NCBIgeneInfo, 
# 	Entrez = GeneID, 
# 	chr = chromosome,
# 	gene_type = type_of_gene,
# 	HGNC_symbol = Symbol_from_nomenclature_authority,
# 	HGNC_name   = Full_name_from_nomenclature_authority)
# NCBIgeneInfo <- NCBIgeneInfo[,c('Entrez', 'Symbol','Synonyms','dbXrefs', 'chr','gene_type','HGNC_symbol','HGNC_name')]


# ## gene alias map
# geneAlias <- NCBIgeneInfo[,c('Entrez', 'Symbol', 'Synonyms')]
# geneAlias <- tidyr::separate_rows(geneAlias, Synonyms, sep = '\\|')

# # remove ambiguous synonyms
# idx = with(geneAlias, !(duplicated(Synonyms) | duplicated(Synonyms,fromLast=TRUE)) | (Synonyms == '-'))
# geneAlias = geneAlias[idx,]

# ## ent2sym map
# geneMap <- NCBIgeneInfo[,c('Entrez','Symbol','HGNC_symbol','gene_type')]


# usethis::use_data(geneMap)
# usethis::use_data(geneAlias)





## code to prepare dataset goes here

library(data.table)
library(tidyverse)
library(usethis)

## Gene Info processing from HGNC raw file
## Downloaded from HGNC genenames.org/download/custom
hgnc <- fread('data-raw/HGNC_Gene_Info_20211112.csv', quote='')
colnames(hgnc) <- gsub(' ', '_', colnames(hgnc))
colnames(hgnc) <- gsub('\\(supplied_by.*\\)', '', colnames(hgnc))

geneInfo <- hgnc %>% 
  filter(Status == 'Approved') %>%
  select(NCBI_Gene_ID, HGNC_ID, Approved_name, Approved_symbol, Ensembl_ID, Locus_group) %>%
  rename(entrez=NCBI_Gene_ID, hgnc_id = HGNC_ID, hgnc_gene = Approved_name, 
    hgnc_symbol = Approved_symbol, ensembl = Ensembl_ID, gene_type = Locus_group) %>%
  mutate(hgnc_id = gsub('^HGNC:','',hgnc_id), entrez = as.character(entrez)) %>%
  arrange(entrez, hgnc_id)

LazygeneInfo <- as.data.frame(geneInfo)

usethis::use_data(LazygeneInfo, overwrite=TRUE)

#   entrez hgnc_id                          hgnc_gene hgnc_symbol           gene_type
# 1      1       5             alpha-1-B glycoprotein        A1BG protein-coding gene
# 2      2       7              alpha-2-macroglobulin         A2M protein-coding gene
# 3      3       8 alpha-2-macroglobulin pseudogene 1       A2MP1          pseudogene
# 4     11      15     N-acetyltransferase pseudogene        NATP          pseudogene
# 5     12      16           serpin family A member 3    SERPINA3 protein-coding gene
# 6     13      17          arylacetamide deacetylase       AADAC protein-coding gene



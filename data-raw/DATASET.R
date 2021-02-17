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

##===========================================================================================
## Mouse gene info

library(biomaRt)
library(rtracklayer)


## Parse GTF File for information
gtf_path <- "gencode.vM24.primary_assembly.annotation.gtf"

# GTF <- fread(gtf_path)
gtf <- readGFF(gtf_path)

# gtfDT <- as.data.table(gtf)

genes <- gtf[which(gtf$type == 'gene'),c('gene_id','seqid','gene_type','gene_name','mgi_id')]
table(genes$gene_type)

colnames(genes)[3] <- 'gene_type_raw'

genes$gene_type <- genes$gene_type_raw
genes$gene_type[grep('protein_coding', genes$gene_type_raw)] <- 'protein-coding gene'
genes$gene_type[grep('pseudogene$'   , genes$gene_type_raw)] <- 'pseudogene'

fwrite(genes, file='data-raw/mousegenes_from_gtf.csv')


# mouse gene IDs
mouse <- useMart('ensembl', dataset='mmusculus_gene_ensembl')
bmlist <- listFilters(mouse)
mouseGenes <- getBM(attributes = c('entrezgene_id','ensembl_gene_id', 'mgi_symbol', 'mgi_id'), mart = mouse)

mouseGeneTypes <- fread('data-raw/mousegenes_from_gtf.csv')
mouseGeneTypes$ensGene <- gsub('.[0-9]+$','', mouseGeneTypes$gene_id)

setdiff(mouseGeneTypes$mgi_id, mouseGenes$mgi_id)

mouseGenes <- mouseGenes[which(mouseGenes$mgi_id!=''),]
mouseGeneTypes <- mouseGeneTypes[which(mouseGeneTypes$mgi_id!=''),]

mouseGeneInfo <- merge(mouseGeneTypes, mouseGenes, by='mgi_id')
mouseGeneInfo <- mouseGeneInfo[,c('mgi_id','ensembl_gene_id','entrezgene_id','mgi_symbol','gene_type')]

save(mouseGeneInfo, file='data/others/mousegeneInfo.RData')




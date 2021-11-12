## GTF File parse

library(data.table)
library(rtracklayer)

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


##===========================================================================================
## Mouse gene info

library(biomaRt)
library(rtracklayer)


## Parse GTF File for information
gtf_path <- "data-raw/gencode.vM24.primary_assembly.annotation.gtf"

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
mouseGenes <- getBM(attributes=c('entrezgene_id','ensembl_gene_id', 'mgi_symbol', 'mgi_id'), mart = mouse, useCache=FALSE)

mouseGeneTypes <- fread('data-raw/mousegenes_from_gtf.csv')

mouseGeneTypes$ensGene <- gsub('.[0-9]+$','', mouseGeneTypes$gene_id)

setdiff(mouseGeneTypes$mgi_id, mouseGenes$mgi_id)

mouseGenes <- mouseGenes[which(mouseGenes$mgi_id!=''),]
# mouseGeneTypes <- mouseGeneTypes[which(mouseGeneTypes$mgi_id!=''),]

mouseGeneInfo <- merge(mouseGeneTypes, mouseGenes, by='mgi_id')
mouseGeneInfo <- mouseGeneInfo[,c('mgi_id','ensembl_gene_id','entrezgene_id','mgi_symbol','gene_type')]

save(mouseGeneInfo, file='data/others/mousegeneInfo.RData')




# library(org.Mm.eg.db)
# mapIds(org.Mm.eg.db, data, keytype="ENSEMBL", column="SYMBOL", multiVals = "first")




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





##==================================##
##  Pathway Analysis With Enricher  ##
##==================================##


# read data from msigdb

loadGMTs <- function(term=NULL, path = DIR_MSIG) {
	if(is.null(term)) {
		cat('loading c2.kegg, c5.gobp\n')
		c2.kegg <<- read.gmt(paste0(DIR_MSIG, '/c2.cp.kegg.v6.1.entrez.gmt'))
		c5.gobp <<- read.gmt(paste0(DIR_MSIG, '/c5.bp.v6.1.entrez.gmt'))		
	} else {
		return(read.gmt(sprintf('%s/%s.v6.1.entrez.gmt', DIR_MSIG, term)))
	}
}

##===================##
##  DEG List helper  ##
##===================##

makeGeneList <- function(resDF) {
	glist = list(
		upGene = with(df, geneEnt[which(t_stat >=  0.05 & t_pVal < .01)]),
		dnGene = with(df, geneEnt[which(t_stat <= -0.05 & t_pVal < .01)]),
		toGene = df$geneEnt )
	return(glist)
}

cellMap <- function(cells) {
	if(all(grepl('^ACH-', cells))) {
		out = cellInfo$CellName[match(cells, cellInfo$DepMap_ID)]
	} else {
		out = cellInfo$DepMap_ID[match(cells, cellInfo$CellName)]
	}
	return(out)
}


##======================================##
##  Gene info data processed from NCBI  ##
##======================================##
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





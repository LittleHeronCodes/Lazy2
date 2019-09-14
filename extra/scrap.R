##==================================##
##  Pathway Analysis With Enricher  ##
##==================================##


# read data from msigdb

# loadGMTs <- function(term=NULL, path = DIR_MSIG) {
# 	if(is.null(term)) {
# 		cat('loading c2.kegg, c5.gobp\n')
# 		c2.kegg <<- read.gmt(paste0(DIR_MSIG, '/c2.cp.kegg.v6.1.entrez.gmt'))
# 		c5.gobp <<- read.gmt(paste0(DIR_MSIG, '/c5.bp.v6.1.entrez.gmt'))		
# 	} else {
# 		return(read.gmt(sprintf('%s/%s.v6.1.entrez.gmt', DIR_MSIG, term)))
# 	}
# }

##===================##
##  DEG List helper  ##
##===================##

# makeGeneList <- function(resDF) {
# 	glist = list(
# 		upGene = with(df, geneEnt[which(t_stat >=  0.05 & t_pVal < .01)]),
# 		dnGene = with(df, geneEnt[which(t_stat <= -0.05 & t_pVal < .01)]),
# 		toGene = df$geneEnt )
# 	return(glist)
# }

# cellMap <- function(cells) {
# 	if(all(grepl('^ACH-', cells))) {
# 		out = cellInfo$CellName[match(cells, cellInfo$DepMap_ID)]
# 	} else {
# 		out = cellInfo$DepMap_ID[match(cells, cellInfo$CellName)]
# 	}
# 	return(out)
# }

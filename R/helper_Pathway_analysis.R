##==================================##
##  Pathway Analysis With Enricher  ##
##==================================##

library(clusterProfiler)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)


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

enricherForGeneListWrapper = function(glist, term, pcut=.1, qcut=.2, minGSSize=5, maxGSSize=100) {
	require(clusterProfiler)

	gspace = glist$toGene
	enrGeneLs = list(
		up   = glist$upGene,
		down = glist$dnGene
	)

	enr_res = lapply(enrGeneLs, function(gset) {
		output = enricher(gset, TERM2GENE = base::get(term), universe = gspace, 
			pvalueCutoff = pcut, pAdjustMethod = 'fdr', qvalueCutoff = qcut,
			minGSSize=minGSSize, maxGSSize=maxGSSize)
		return(output)
		})
	
	print(sapply(enr_res, nrow))
	return(enr_res)
}


enrObjectTransform_1set <- function(enr_obj, pco=.1, max_row=10) {
	require(tidyr)
	LS = lapply(names(enr_obj), function(dir) {
		df = data.frame(enr_obj[[dir]])
		df = df[order(df$pvalue),]
		df = df[which(df$pvalue < pco), c('ID', 'pvalue', 'p.adjust', 'Count')]
		if(nrow(df) >= max_row) df = df[order(df$pvalue)[1:max_row],]
		df$dir = factor(dir, levels=c('up','down'))
		df$logP = -log10(df$pvalue)
		return(df)
	})
	hmplot = do.call(rbind, LS)

	pathDF = tidyr::spread(hmplot[, c('dir','ID','logP')], dir, logP, fill = NA)
	pathDF$down = pathDF$down * -1
	plotMat = as.matrix(pathDF[,-1])
	rownames(plotMat) = gsub('^GO_|^KEGG_','',pathDF[,'ID'])
	plotMat = plotMat[order(apply(plotMat,1, sum, na.rm=T), decreasing=TRUE),]
	return(plotMat)
}


enrDrawHeatmap_1set <- function(enrLs, pco, max_row, mtitle) {
	colors = colorRamp2(c(-4,-2,0,2,4), 
		c('#134173','#46cbfc','#f7f7f7','#ff6249','#a63b39'))
	plotMat <- enrObjectTransform_1set(enrLs, pco = pco, max_row=max_row)
	hp = Heatmap(plotMat, cluster_rows = FALSE, cluster_columns = FALSE, col = colors,
		heatmap_width=unit(8,'cm'), heatmap_height = unit(16, 'cm'),
		rect_gp = gpar(col = "gray12", lty = 1, lwd = 0.2), 
		row_names_side='left', column_title = mtitle)
	return(hp)
}


enrObjectTransform_nset <- function(enr_obj, pco=.1, max_row=10) {
	LS = lapply(names(enr_obj), function(set) {
		df = data.frame(enr_obj[[set]])
		df = df[order(df$pvalue),]
		df$set = set
		df = df[which(df$pvalue < pco), c('set','ID', 'pvalue', 'p.adjust', 'Count')]
		if(nrow(df) >= max_row) df = df[order(df$pvalue)[1:max_row],]
		df$logP = -log10(df$pvalue)
		return(df)
	})
	hmplot = do.call(rbind, LS)

	pathDF = tidyr::spread(hmplot[, c('set','ID','logP')], set, logP, fill = NA)
	plotMat = as.matrix(pathDF[,-1])
	rownames(plotMat) = gsub('^GO_|^KEGG_','',pathDF[,'ID'])
	plotMat = plotMat[order(apply(plotMat,1, sum, na.rm=T), decreasing=TRUE),]
	return(plotMat)
}

enrDrawHeatmap_nset <- function(enrLs, pco, max_row, mtitle, colpal='up', clust = TRUE) {

	if(colpal == 'up') palette = c('#f7f7f7','#ff6249','#a63b39')
	if(colpal == 'dn') palette = c('#f7f7f7','#46cbfc','#134173')

	colors = colorRamp2(c(0,2,4), palette)
	plotMat <- enrObjectTransform_nset(enrLs, pco = pco, max_row=max_row)
	plotMat[which(is.na(plotMat))] = 0
	hp = Heatmap(plotMat, cluster_rows = TRUE, cluster_columns = clust, col = colors,
		heatmap_width=unit(12,'cm'), heatmap_height = unit(24, 'cm'),
		rect_gp = gpar(col = "gray12", lty = 1, lwd = 0.2), 
		show_row_dend = FALSE, show_column_dend = FALSE, name = 'pval',
		row_names_side='left', column_title = mtitle)
	return(hp)
}






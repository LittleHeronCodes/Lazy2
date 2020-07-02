#' enricherForGeneListWrapper
#' 
#' Wrapper function of enricher for differentially expressed gene list
#' @param glist List of genes (Entrez or Symbol) with names upGene, dnGene, toGene
#' @param term  Pathway gene set term in containing pathway list
#' @param pcut  p-value cutoff. Feeds to pvalueCutoff.
#' @param qcut  q-value cutoff (FDR adjusted). Feeds to qvalueCutoff.
#' @param minGSSize minimum gene set size for enricher
#' @param maxGSSize maximum gene set size for enricher
#' @return List of enricher Object: up, down
#' @export
#' @examples
#' \dontrun{
#' sapply(glist, length)
#' enricherForGeneListWrapper(glist, c2.kegg)
#' }


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



enricherForGLIST <- function(geneList, use) {
	require(clusterProfiler)
	enrLS <- list()
	for(dir in names(geneList[[1]])) {
		glist = lapply(geneList, function(ls) ls[[dir]])
		enrObj = lapply(glist, function(gset) {
			if(use == 'kegg') {
				output = enrichKEGG(gset, minGSSize=10, maxGSSize=200,
					pvalueCutoff = 1, pAdjustMethod = 'fdr', qvalueCutoff = 1)				
			} else if(use == 'gobp') {
				require(org.Hs.eg.db)
				output = enrichGO(gset, ont = 'BP', OrgDb = org.Hs.eg.db, minGSSize=50, maxGSSize=200,
					pvalueCutoff = 1, pAdjustMethod = 'fdr', qvalueCutoff = 1)
			} else {stop('Currently available use option : kegg, gobp')}
			})
		enrLS[[dir]] = enrObj
	}
	return(enrLS)
}



#' enrObjectTransform
#' 
#' enricher Object list format transformation for plotting
#' 
#' @param enr_obj enricher object list for single set (names: up, down)
#' @param pco     p-value cutoff Default .1
#' @param max_row maximum row of plotMat to return Default 10
#' @return Matrix to feed to enr_HeatmapOnly (plotMat)
#' @export
#' @examples
#' \dontrun{
#' plotMat <- enrObjectTransform_1set(enr_obj)
#' enrHeatmapOnly_1set(plotMat)
#' }

enrObjectTransform_1set <- function(enr_obj, pco=.1, max_row=10) {
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


#' @describeIn enrObjectTransform_1set For multiple dataset gene list
#' @export

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


#' @section Drawing Heatmaps
#' @describeIn enrObjectTransform_1set
#' Heatmap draw from plotMat from enrObjectTransform_1set 
#' @param plotMat Matrix for plotting from enrTransform functions
#' @param mtitle  Main title for pathway analysis heatmaps
#' @param colpal  Choosing color palette for nset pathway analysis heatmap. ('up' or 'dn')
#' @param clust   Cluster Heatmap or no?
#' @export

enrHeatmapOnly_1set <- function(plotMat, mtitle, clust = TRUE) {
	colors = colorRamp2(c(-4,-2,0,2,4), 
		c('#134173','#46cbfc','#f7f7f7','#ff6249','#a63b39'))
	hp = Heatmap(plotMat, cluster_rows = clust, cluster_columns = FALSE, col = colors,
		heatmap_width=unit(8,'cm'), heatmap_height = unit(16, 'cm'),
		rect_gp = gpar(col = "gray12", lty = 1, lwd = 0.2), 
		row_names_side='left', column_title = mtitle)
	return(hp)
}

#' @describeIn enrObjectTransform_1set
#' Heatmap draw from plotMat from enrObjectTransform_nset 
#' @export

enrHeatmapOnly_nset <- function(plotMat, max_row, mtitle, colpal='up', clust = TRUE) {

	if(colpal == 'up') palette = c('#f7f7f7','#ff6249','#a63b39')
	if(colpal == 'dn') palette = c('#f7f7f7','#46cbfc','#134173')

	colors = colorRamp2(c(0,2,4), palette)
	plotMat[which(is.na(plotMat))] = 0
	ord = apply(plotMat, 1, function(v) sum(v != 0))
	plotMat = plotMat[order(ord, decreasing=TRUE),]
	plotMat = plotMat[1:max_row,]

	hp = Heatmap(plotMat, cluster_rows = clust, cluster_columns = clust, col = colors,
		heatmap_width=unit(12,'cm'), heatmap_height = unit(24, 'cm'),
		rect_gp = gpar(col = "gray12", lty = 1, lwd = 0.2), 
		show_row_dend = FALSE, show_column_dend = FALSE, name = 'pval',
		row_names_side='left', column_title = mtitle)
	return(hp)
}


#' enrLS2plotMat
#'
#' Formerly enrObjectTransform nset, fixed to no filter, ordering function
#' 
#' @param enrLS enricher object list
#' @param ord.FUN function for ordering output plotMat
#' @export
#'

enrLS2plotMat <- function(enrLS, ord.FUN = NULL ) {
	hmplot = do.call(rbind, lapply(names(enrLS), function(set) {
		df = data.frame(enrLS[[set]])
		df = df[order(df$p.adjust),]
		df$set = set
		df$logP = -log10(df$p.adjust)
		return(df)
		}) )

	pathDF = tidyr::spread(hmplot[, c('set','ID','logP')], set, logP, fill = NA)
	plotMat = as.matrix(pathDF[,-1])
	rownames(plotMat) = pathDF[,'ID']
	plotMat = plotMat[order(apply(plotMat,1, mean, na.rm=T), decreasing=TRUE),]

	if(is.null(ord.FUN)) { 
		ord.FUN = function(v) sum(!is.na(v))
	}
	if(class(ord.FUN) == 'function') { 
		ord = apply(plotMat, 1, ord.FUN) 
	}
	if(class(ord.FUN) != 'function') { 
		ord = ord.FUN 
	}

	plotMat = plotMat[order(ord, decreasing=TRUE),]
	return(plotMat)
}



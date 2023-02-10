

## Functions to help draw plot
getSize <- function(pt) {
	# pt = draw(hp, heatmap_legend_side='left', merge=TRUE, padding = unit(c(2, 20, 2, 20), "mm"))
	W = ceiling(convertX(sum(component_width(pt)  + unit(1,'mm')), 'in',valueOnly=TRUE))
	H = ceiling(convertX(sum(component_height(pt) + unit(1,'mm')), 'in',valueOnly=TRUE))
	return(c(W, H))
}



drawHeatmap_hgeo2 <- function(plotMat, colpal, cell_size=4.9, asp=1.0,
	colsplit=NULL, rowsplit=NULL, name='-log10Q', legend_dir='vertical', ha=NULL, ...) {
	
	colgap <- ifelse(!is.null(colsplit), rep(1.5,length(unique(colsplit))-1), 1)
	rowgap <- ifelse(!is.null(rowsplit), rep(1.5,length(unique(rowsplit))-1), 1)
	wmm <- cell_size * 1.0 * ncol(plotMat) + ifelse(!is.null(colsplit), sum(colgap), 0)
	hmm <- cell_size * asp * nrow(plotMat) + ifelse(!is.null(rowsplit), sum(rowgap), 0)
	
	# Draw heatmap
	hp <- Heatmap(plotMat, name = name, col=colpal, bottom_annotation = ha,
			row_names_gp=gpar(col='grey2'), column_names_gp=gpar(col='grey2'), column_names_rot = 60,
			width = unit(wmm, 'mm'), height = unit(hmm, 'mm'),
			column_split=colsplit, column_gap= unit(colgap, 'mm'),
			row_split=rowsplit, row_gap= unit(rowgap, 'mm'),
			rect_gp = gpar(col = "gray12", lty = 1, lwd = 0.2),
			# cluster_rows = cluster, show_row_dend = cluster, show_row_names = TRUE,
			cluster_columns = FALSE, show_column_dend = FALSE, show_column_names = TRUE,
			border=TRUE, row_names_max_width = max_text_width( rownames(plotMat) ),
			column_title_gp = gpar(size=11),
			heatmap_legend_param=list(border=TRUE, direction=legend_dir),
	...)
	hp
}






tstatLS <- lapply(c(resultsLS.mm, resultsLS.ph), function(dff) {
	dff <- dff[which(!is.na(dff$entGene)),]
	dff$entGene = as.character(dff$entGene)
	dff <- merge(dff, mmu2hsa, by.x='entGene', by.y='MGI_ENTREZ', all.x=TRUE, all.y=FALSE)
	dff <- dff[which(!is.na(dff$HUMAN_ENTREZ)),]
	if(!'stat' %in% colnames(dff) & 't' %in% colnames(dff)) dff$stat = dff$t
	out <- unlist(tapply(dff$stat, dff$HUMAN_ENTREZ, function(v) v[which.max(abs(v))] ))
	# out <- unlist(tapply(dff$logFC, dff$entGene, mean))
	out <- sort(out, decreasing=TRUE)
	})


## GSEA (using fgsea)
set.seed(1234)
gseaLS <- mclapply(tstatLS, function(lfc) {
	gsea = fgsea(pathways = intgpath, stats = lfc, nperm=10000, minSize=20, maxSize=500)
	qv <- qvalue(gsea$pval)
	gsea$qVal <- qv$qvalues
	return(gsea)
	}, mc.cores=20)


###================================================================================

# (replaces geneListDistMat)
getPairmat <- function(pairdf, metric='ef') {
	pairm = reshape2::acast(pairdf, ix1~ix2, value.var=metric)
	diag(pairm) = NA
	return(pairm)
}




###  KEGG, GO, GSEA HELPER (ClusterProfiler)  ###
## KEGG wrapper
KEGG_enrichment <- function(glist, tglist=NULL, organism='hsa') {
	if(!organism %in% c('hsa','mmu')) stop('Supported organisms : hsa or mmu.')
	universe = NULL
	enrobj = lapply(names(glist), function(aid) {
		gset = glist[[aid]]
		if(!is.null(tglist)) universe = tglist[[aid]]
		enr = enrichKEGG(gset, organism = organism, universe=universe,
			pvalueCutoff=1.5, pAdjustMethod='fdr', qvalueCutoff=1.5, minGSSize=5, maxGSSize=500)
		return(enr)
		})
	names(enrobj) = names(glist)
	return(enrobj)
}


## GO wrapper
GO_enrichment <- function(glist, organism='hsa', ont='BP', ncore=1) {
	require(parallel)
	if(!organism %in% c('hsa','mmu')) stop('Supported organisms : hsa or mmu.')
	enrobj = mclapply(names(glist), function(aid) {
		gset = glist[[aid]]
		orgdb = NULL
		if(organism == 'hsa') orgdb = 'org.Hs.eg.db'
		if(organism == 'mmu') orgdb = 'org.Mm.eg.db'
		enr = enrichGO(gset, OrgDb = orgdb, ont=ont,
			pvalueCutoff=1.5, pAdjustMethod='fdr', qvalueCutoff=1.5, minGSSize=50, maxGSSize=150)
		return(enr)
		}, mc.cores=ncore)
	names(enrobj) = names(glist)
	return(enrobj)
}



## Draw Heatmap basic
enrSaveHeatmap <- function(plotMat2, mtitle, colpal_t, ha, name='logQ', fign=NULL, save=FALSE) {

	if(nrow(plotMat2) >= 40) {
		cat('Resort and plotting only top 40 mean values.\n')
		incidx = order(-apply(plotMat2, 1, function(v) mean(v)))[1:40]
		plotMat2 = plotMat2[incidx, ]
	}

	wmm = 5.5 * ncol(plotMat2)
	hmm = 5.0 * nrow(plotMat2)
	hp = Heatmap(plotMat2, col=colpal_t, name=name, column_title = mtitle, 
		width = unit(wmm, 'mm'), height = unit(hmm, 'mm'),
		row_names_gp=gpar(col='grey2'), column_names_gp=gpar(col='grey2'),
		bottom_annotation = ha, 
		rect_gp = gpar(col = "gray12", lty = 1, lwd = 0.2),
		cluster_rows = FALSE, show_row_dend = FALSE, show_row_names = TRUE,
		cluster_columns = FALSE, show_column_dend = FALSE, show_column_names = TRUE)

	if(save & !is.null(fign)) {
		W = ceiling(convertX(hp@matrix_param$width + unit(18,'cm'), 'in', valueOnly=T))
		H = ceiling(convertX(hp@matrix_param$height + unit(7,'cm'), 'in', valueOnly=T))
		png(file=fign, width=W, height=H, unit='in', res=150, type='cairo')
		draw(hp, heatmap_legend_side='left', padding=unit(c(2,10,5,70),'mm'),  merge=TRUE)
		dev.off()		
	} else {
		draw(hp, heatmap_legend_side='left', padding=unit(c(2,10,5,70),'mm'),  merge=TRUE)
	}

}






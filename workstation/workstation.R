## Supported formats : ebayes, DEseq2


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


## GSEA (from clusterProfiler)



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



## ggplot theme set

theme_transparent2 <- function(base_size=12, x.text.angle=0) {
	half_line = base_size/2
	xjust <- 0
	if(x.text.angle > 5) xjust <- 1
	.theme <- theme(
		text = element_text(face='bold', size=base_size),# family='Sans'),
		panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
		panel.border = element_blank(),
		plot.title = element_text(hjust = 0.5,vjust=1.5, size=rel(1), margin=unit(c(0.5,0,0.5,0), 'mm')),
		axis.line = element_line(colour='black', size=1),
		axis.text.x = element_text(size=rel(1), colour='black'),
		axis.text.y = element_text(size=rel(1), colour='black', face='plain'),
		legend.key = element_blank(), 
		panel.background = element_rect(fill = "transparent", colour = NA), panel.ontop=TRUE,
		plot.background = element_rect(fill = "transparent",colour = NA),
		legend.background = element_rect(fill = "transparent", colour = NA), 
		strip.background = element_rect(fill = "#F2F2F2", colour = "black", size = 0.7), 
		# plot.margin = margin(half_line, half_line, half_line, half_line), # to reduce space
		plot.margin = unit(c(1,1,1,1), 'mm'), #margin(half_line, half_line, half_line, half_line), # to reduce space
		complete = TRUE)
	if(x.text.angle != 0)
		.theme <- .theme + theme(axis.text.x = element_text(angle = x.text.angle, hjust = xjust, vjust=xjust))
	.theme
}

# theme_set(theme_transparent2())




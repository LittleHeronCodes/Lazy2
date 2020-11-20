## Updated pathway analysis set


drawHeatmap <- function(plotMat, col_pal, mtitle, name = 'logQ',ha = NULL) {
	colors = colorRamp2(c(0,.2,1,4), col_pal)

	hp = Heatmap(plotMat, col = colors, column_title = mtitle, name = name,
		cluster_rows = FALSE, cluster_columns = FALSE,
		row_names_gp = gpar(fontsize = 9), column_names_gp = gpar(fontsize = 9),
		show_row_dend = FALSE, show_column_dend = FALSE,
		row_names_max_width = unit(1, 'npc'), width = unit(3, 'cm'),
		row_names_side='right', rect_gp = gpar(col = "gray25", lty = 1, lwd = 0.4), 
		heatmap_legend_param = list(direction='vertical', title_position='topcenter'),
		right_annotation = ha)
	return(hp)
}


plotMatNamesClean <- function(plotMat) {
	rownames(plotMat) = gsub('^GO_|^KEGG_|^WIKI_','',rownames(plotMat))
	rownames(plotMat) = gsub('_',' ',rownames(plotMat))
	colnames(plotMat) = gsub('_vs_Control$', '/Con',colnames(plotMat))
	return(plotMat)
}


enrLS2plotMat <- function(enrLS, ord.FUN = NULL, use_adjP = TRUE, id.col='ID') {

	# detect P value column
	pvidx = grep('^p.?val', names(as.data.frame(enrLS[[1]])), ignore.case = TRUE)
	qvidx = grep('^p.?adj|^adj.?p', names(as.data.frame(enrLS[[1]])), ignore.case = TRUE)
	if( (length(pvidx) != 1) | ((length(qvidx) != 1) & use_adjP ) ) {
		stop("zero or multiple pvalue columns detected. set pvalue column name to pvalue, qvalue to adj.p")
	}
	if(length(intersect(pvidx, qvidx)) >= 1){
		stop("same columns for pvalue and p adjusted. I can't work in these conditions.")
	}

	if(use_adjP) ix = names(as.data.frame(enrLS[[1]]))[qvidx]
	if(!use_adjP) ix = names(as.data.frame(enrLS[[1]]))[pvidx]

	hmplot = do.call(rbind, lapply(names(enrLS), function(set) {
		df = data.frame(enrLS[[set]])
		df = df[order(df[,ix]),]
		df$set = set
		df$logP = -log10(df[,ix])
		return(df)
		}) )

	pathDF = tidyr::spread(hmplot[, c('set',id.col,'logP')], set, logP, fill = NA)
	plotMat = as.matrix(pathDF[,-1])
	rownames(plotMat) = pathDF[,id.col]
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



enrDrawPlotSave.go <- function(enrObj, use_adjP) {
	nm = ifelse(use_adjP, 'logQ','logP')
	plotMat <- enrLS2plotMat(enrObj, ord.FUN = ord.FUN, use_adjP = use_adjP)
	print(dim(plotMat))
	plotMat <- plotMat[1:30,c(1,5,3,2,6,4)]
	hp <- drawHeatmap(plotMat, col_pal$RdWh, sprintf('%s %s (fco %.1f)',dir, 'GO', fco), name=nm)
	# png(sprintf('%s/pathway/cell_death/CellDeathPathways_p_fco%.1f_%s_up.png',yml$fig, fco, 'GO'), 
	# 	width = 17, height = 22, unit = 'cm', res = 150)
	draw(hp, heatmap_legend_side = 'left')
	# dev.off()
}



enrDrawPlotSave.kg <- function(enrObj, use_adjP, col_pal) {
	nm = 'logP' #ifelse(use_adjP, 'logQ','logP')
	plotMat <- enrLS2plotMat(enrObj, ord.FUN = NULL, use_adjP = FALSE)
	# plotMat <- plotMat[grep(paste(cell_death_terms, collapse='|'), rownames(plotMat), ignore.case = TRUE), ]
	print(dim(plotMat))
	# plotMat <- plotMat[,c(1,5,3,2,6,4)]

	hacat = rep('',nrow(plotMat))
	for(x in cell_death_terms) {
		hacat = ifelse(grepl(x, rownames(plotMat), ignore.case=TRUE), x, hacat)
	}
	hacat = ifelse(hacat %in% c('apopto','necro','ferropto'), paste0(hacat, 'sis'),
		ifelse(hacat == 'autophag', paste0(hacat, 'y'), 'others'))
	hacol = structure(RColorBrewer::brewer.pal(8, 'Set2')[1:length(unique(hacat))], 
		names = sort(unique(hacat)))
	ha = HeatmapAnnotation(categories = hacat, which = 'row', 
		col = list(categories = hacol), width=unit(.5,'npc'))

	hp <- drawHeatmap(plotMat, col_pal, sprintf('%s %s (fco %.1f)', dir, 'kegg', fco), name=nm, ha=ha)
	png(sprintf('%s/pathway/cell_death/CellDeathPathways_%s_fco%.1f_%s_%s_.png', yml$fig,nm,fco,dir,'KEGG'), 
		width = 15, height = 10, unit = 'cm', res = 150)
	draw(hp, heatmap_legend_side = 'left')
	dev.off()
}

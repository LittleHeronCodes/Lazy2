## Gene list extraction from result
## Supported formats : ebayes, treat, confect, DEseq2

#' 
#' @param 
#' @param 
#' @param 
#' @param 
#' @param 
#' @param 
#' @return 
#' @export
#' @examples
#' \dontrun{
#' }

extractGeneList <- function(resultsLS, fco, qco, pco=NULL, cnt=NULL, remove_ambi=FALSE) {

	## Add logFC, p-value, adjusted p-value column detector ##
	## 
	## 
	## 

	if(length(qco) == 1) qco <- structure(rep(qco, length(resultsLS)), names=names(resultsLS))
	if(length(fco) == 1) fco <- structure(rep(fco, length(resultsLS)), names=names(resultsLS))

	geneList <- list(up = list(), dn = list(), to=list())
	for(aid in names(resultsLS)) {
		resultDF <- resultsLS[[aid]] %>% filter(adj.P.Val < qco[aid] & !is.na(entGene))

		geneList$up[[aid]] <- resultDF %>% filter(logFC >=  log2(fco[aid])) %>% arrange(-logFC) %>% .$entGene %>% unique()
		geneList$dn[[aid]] <- resultDF %>% filter(logFC <= -log2(fco[aid])) %>% arrange( logFC) %>% .$entGene %>% unique()
		geneList$to[[aid]] <- resultsLS[[aid]] %>% filter(!is.na(entGene)) %>% .$entGene %>% unique()

	}

	## Remove ambiguous option ##
	if(remove_ambi) geneList <- removeAmbigDEGs(geneList)
	## print how many were removed
	## 

	return(geneList)
}

## Gene list ambiguoity remove
## 

#' 
#' @param geneList gene list
#' @return 
#' @export
#' @examples
#' \dontrun{
#' }

removeAmbigDEGs <- function(geneList) {
	require(purrr)
	ambi <- map2(geneList$up, geneList$dn, function(x,y) intersect(x,y) )
	geneList$up <- map2(geneList$up, ambi, function(x,y) setdiff(x,y) )
	geneList$dn <- map2(geneList$dn, ambi, function(x,y) setdiff(x,y) )
	return(geneList)
}



## Overlap metrics
## calculate genelist set overlap
## compare with map2 methods

getOverlapDF <- function(gls, tgls) {
	pairdf <- expand.grid(names(gls), names(gls), stringsAsFactors=FALSE)
	LS <- apply(pairdf, 1, function(v) {
		ix1 = as.character(v[1])
		ix2 = as.character(v[2])
		gspace = intersect(tgls[[ix1]], tgls[[ix2]])
		setA = intersect(gls[[ix1]], gspace)
		setB = intersect(gls[[ix2]], gspace)
		hgeo = hypergeoTest(setA, setB, gspace)
		eff = getEnrichmentFactor(setA, setB, gspace, 2)
		tan = tanimotoCoef(setA, setB)
		data.frame(ix1 = ix1, ix2 = ix2, hgeo = hgeo, ef = eff, tan=tan)
		})
	pairdf <- do.call(rbind, LS)	
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


## General hypergeo test
Gen_enrichment <- function(glist, refgmt, tglist, ncore=1) {
	require(parallel)
	enrobj = mclapply(names(glist), function(aid) {
		hgeos = hypergeoTestForGeneset(glist[[aid]], refgmt, tglist[[aid]])
		hgeos$qVal = p.adjust(hgeos$pVal, method='fdr')
		hgeos$logQ = -log10(hgeos$qVal)
		return(hgeos)
		}, mc.cores = ncore)
	names(enrobj) = names(glist)
	return(enrobj)
}

## GSEA (from clusterProfiler)



## enrichment object list to matrix (genarilzed function compatible)
enrobj2Matrix <-function(enrobj, val.col='pvalue', log=TRUE) {
	LS = lapply(names(enrobj), function(set) {
		dff = data.frame(enrobj[[set]])
		if('Description' %in% names(dff)) dff = dff %>% dplyr::rename(termID=ID, ID=Description)
		dff$set = set
		dff = dff[order(dff$set),]
		return(dff)
	})
	hmplot = do.call(rbind, LS)

	# detect log values
	logdetect <- FALSE
	if(any(quantile(hmplot[,val.col], na.rm=TRUE) > 1)) logdetect <- TRUE
	if(log & logdetect) cat('Already in log values.'); log <- FALSE
	if(log & !logdetect) {
		hmplot$logV <- -log10(hmplot[,val.col])
		plotMat = reshape2::acast(hmplot, ID~set, value.var=logV, fill = NA)
	}

	if(!log) plotMat = reshape2::acast(hmplot, ID~set, value.var=val.col, fill = NA)
	plotMat = plotMat[order(apply(plotMat,1, sum, na.rm=TRUE), decreasing=TRUE),]
	return(plotMat)
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


##====================================================
## GSEA plot

plotEnrichment2 <- function(gset, stats, nes, qv, gseaParam = 1, mtitle=NULL, ylab='',
	ticksSize=0.4, base_size=7, line.col='green', lwd=2) {

	require(gtable)
	require(fgsea)
	rnk <- rank(-stats)
	ord <- order(rnk)
	statsAdj <- stats[ord]
	statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
	statsAdj <- statsAdj/max(abs(statsAdj))
	pathway <- unname(as.vector(na.omit(match(gset, names(statsAdj)))))
	pathway <- sort(pathway)
	gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, returnAllExtremes = TRUE)
	bottoms <- gseaRes$bottoms
	tops <- gseaRes$tops
	txt <- sprintf('NES : %.2f  \nq-value : %.2e  ', nes, qv)

	n <- length(statsAdj)
	xs <- as.vector(rbind(pathway - 1, pathway))
	ys <- as.vector(rbind(bottoms, tops))
	toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
	diff <- (max(tops) - min(bottoms))/8

	ln1=round(seq(min(bottoms), max(tops), 0.1), digits = 1)
	half_line = base_size/2

	g1 <- ggplot(toPlot, aes(x = x, y = y)) + geom_point(color = line.col, size = 0.1) + 
	  geom_hline(yintercept = ln1, colour = "grey85", linetype='dashed',size = lwd*0.8) +
	  geom_hline(yintercept = 0, colour = "black", linetype='dashed', size = lwd*0.8) +
	  geom_hline(yintercept = ifelse(nes>0, max(tops), min(bottoms)), colour = "red", linetype = "dashed", size = lwd*0.8) +
	  geom_line(color = line.col,size = lwd) +
	  annotate('text', x=max(toPlot$x), y=max(c(ln1,tops)), label=txt, hjust=1, vjust=1.5, fontface='plain', size=rel(3.0)) +
	  labs(y = ylab, title=mtitle) + 
	  theme_common(base_size=base_size) +
	  theme(
	    plot.title = element_text(hjust = 0.5, vjust=0.2, face='bold', margin=unit(c(0,0,1.5,0), 'mm')),
	    axis.title.y = element_text(face='bold', angle=90, margin=unit(c(0,1.5,0,0), 'mm'),size=rel(0.95)),
	    # plot.margin = unit(c(1.2,2.0,0,1.2), 'mm'),
	    plot.margin = margin(half_line, half_line*2.5, 0, half_line),
	  	axis.title.x=element_blank(), axis.text.x=element_blank()
	  	)
	g2 <- ggplot(data.frame(x=pathway),aes(x = x, y = -diff/2, xend = x, yend = diff/2)) +
	  geom_segment(size = ticksSize, colour='grey35') +
	  theme_common(base_size=base_size) +
	  theme(
		# plot.margin = unit(c(0,2.0,1.2,1.2), 'mm'),
		plot.margin = margin(0, half_line*2.5, half_line, half_line),
	  	axis.text.y=element_blank(), axis.title=element_blank()
	 	)

	gr1 <- ggplotGrob(g1)
	gr2 <- ggplotGrob(g2)
	gr <- rbind(gr1, gr2)
	gr$widths <- unit.pmax(gr1$widths, gr2$widths)

	# identify the position of the panels within the gtable
	panid <- gr$layout$t[grep(pattern="panel", gr$layout$name)]
	gr$heights[panid] <- unit(c(7,1), 'null')
	grid.newpage()
	grid.draw(gr)

	return(gr)
}


theme_gsea_common <- function(base_size=5) {
	half_line = base_size/2
	.theme <- theme(
		text = element_text(face='plain', size=base_size, colour='black', family='Arial'),
		plot.title = element_text(size=rel(1.1)),
		axis.ticks=element_blank(), 
		axis.text = element_text(size=rel(0.8), colour='black'),
		panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
		panel.background = element_rect(fill = "transparent", colour = NA),
		plot.background = element_rect(fill = "transparent",colour = NA),
		complete=TRUE)
	.theme
}





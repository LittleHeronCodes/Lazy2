## Gene list extraction from result
## Supported formats : ebayes, DEseq2

#' 
#' @param resultsLS list of result data frames. (limma results default)
#' @param fco fold change cut-offs (NOT LOG). Used to filter logFC column.
#' @param qco adjusted p value cut-offs. Used to filter adj.P.Val column.
#' @param cnt DEG count constrains (not used).
#' @param remove_ambi Remove genes both in up and down?
#' @return 
#' @export 
#' @examples
#' \dontrun{
#' }

extractGeneList <- function(resultsLS, fco, qco, cnt=NULL, remove_ambi=FALSE) {

	## Add logFC, p-value, adjusted p-value column detector ##
	## 
	## 
	## 

	if(length(qco) != 1 & length(qco) != length(resultsLS)) stop('qco length should be either 1 or same as resultsLS.')
	if(length(fco) != 1 & length(fco) != length(resultsLS)) stop('fco length should be either 1 or same as resultsLS.')

	if(length(qco) == 1) qco <- structure(rep(qco, length(resultsLS)), names=names(resultsLS))
	if(length(fco) == 1) fco <- structure(rep(fco, length(resultsLS)), names=names(resultsLS))

	geneList <- list(up = list(), dn = list(), to=list())
	for(aid in names(resultsLS)) {
		resultDF.f <- resultsLS[[aid]] %>% filter( !is.na(entGene) )

		geneList$up[[aid]] <- with(resultDF.f, unique(entGene[which(adj.P.Val < qco[aid] & logFC >=  log2(fco[aid]))]) )
		geneList$dn[[aid]] <- with(resultDF.f, unique(entGene[which(adj.P.Val < qco[aid] & logFC <= -log2(fco[aid]))]) )
		geneList$to[[aid]] <- unique(resultDF.f$entGene)

	}

	## Remove ambiguous option
	if(remove_ambi) geneList <- removeAmbigDEGs(geneList)

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




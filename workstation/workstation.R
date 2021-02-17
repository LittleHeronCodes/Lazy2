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

	## Add logFC, p-value, adjusted p-value detector ##
	## 
	## 
	## 

	if(length(qco) == 1) qco = structure(rep(qco, length(resultsLS)), names=names(resultsLS))
	if(length(fco) == 1) fco = structure(rep(fco, length(resultsLS)), names=names(resultsLS))

	geneList = list(up = list(), dn = list(), to=list())
	for(aid in names(resultsLS)) {
		resultDF = resultsLS[[aid]] %>% filter(adj.P.Val < qco[aid] & !is.na(entGene))

		geneList$up[[aid]] = resultDF %>% filter(logFC >=  log2(fco[aid])) %>% arrange(-logFC) %>% .$entGene %>% unique()
		geneList$dn[[aid]] = resultDF %>% filter(logFC <= -log2(fco[aid])) %>% arrange( logFC) %>% .$entGene %>% unique()
		geneList$to[[aid]] = resultsLS[[aid]] %>% filter(!is.na(entGene)) %>% .$entGene %>% unique()

		## Remove ambiguous option ##
		## if(remove_ambi) ambi = intersect(geneList$up[[aid]], geneList$dn[[aid]])
		## 
		## 

	}
	return(geneList)
}



## Overlap metrics
## calculate genelist set overlap
## compare with map2 methods

getOverlapDF <- function(gls, tgls) {
	pairdf = expand.grid(names(gls), names(gls), stringsAsFactors=FALSE)
	LS = apply(pairdf, 1, function(v) {
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
	pairdf = do.call(rbind, LS)	
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
	enrobj = mclapply(names(glist), function(aid) {
		hgeos = hypergeoTestForGeneset(glist[[aid]], refgmt, tglist[[aid]])
		hgeos$qVal = p.adjust(hgeos$pVal, method='fdr')
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
		dff = dff[order(dff[,val.col]),]
		dff$set = set
		if(log) dff$logP = -log10(dff[,val.col])
		return(dff)
	})
	hmplot = do.call(rbind, LS)
	logn = ifelse(grepl('q|fdr|Q|FDR', val.col), 'logQ', 'logP')
	if( log) plotMat = reshape2::acast(hmplot, ID~set, value.var=logn, fill = NA)
	if(!log) plotMat = reshape2::acast(hmplot, ID~set, value.var=val.col, fill = NA)
	plotMat = plotMat[order(apply(plotMat,1, sum, na.rm=T), decreasing=TRUE),]
	return(plotMat)
}


###================================================================================================
###================================================================================================
###================================================================================================
extractGeneList <- function(resultsLS, fco, qco, pcov=NULL, cnt=NULL, remove_ambi=FALSE) {

	# check inputs
	if(length(qcov) == 1) qco = structure(rep(qco, length(resultsLS)), names=names(resultsLS))
	if(length(fcov) == 1) fco = structure(rep(fco, length(resultsLS)), names=names(resultsLS))

	geneList = list(up = list(), dn = list(), to=list())
	for(aid in names(resultsLS)) {
		resultDF = resultsLS[[aid]] %>% filter(adj.P.Val < qco[aid] & !is.na(entGene))

		geneList$up[[aid]] = resultDF %>% filter(logFC >  log2(fco[aid])) %>% arrange(-logFC) %>% .$entGene %>% unique()
		geneList$dn[[aid]] = resultDF %>% filter(logFC < -log2(fco[aid])) %>% arrange( logFC) %>% .$entGene %>% unique()
		geneList$to[[aid]] = resultsLS[[aid]] %>% filter(!is.na(entGene)) %>% .$entGene %>% unique()
	}

	return(geneList)
}

# remove ambiguous DEGs
removeAmbigDEGs <- function(geneList) {
	require(purrr)
	ambi <- map2(geneList$up, geneList$dn, function(x,y) intersect(x,y) )
	geneList$up <- map2(geneList$up, ambi, function(x,y) setdiff(x,y) )
	geneList$dn <- map2(geneList$dn, ambi, function(x,y) setdiff(x,y) )
	return(geneList)
}


## mouse entrez to mgi
ent2mgi <- function(ent) {
	out = mouseGenes$mgi_symbol[match(ent, mouseGenes$entrezgene_id)]
	return(out)
}


## KEGG wrapper
KEGG_enrichment <- function(glist, tglist=NULL, organism='hsa') {
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
	require(Lazy2)
	require(parallel)
	enrobj = mclapply(names(glist), function(aid) {
		hgeos = hypergeoTestForGeneset(glist[[aid]], refgmt, tglist[[aid]])
		hgeos$qVal = p.adjust(hgeos$pVal, method='fdr')
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
		dff = dff[order(dff[,val.col]),]
		dff$set = set
		if(log) dff$logP = -log10(dff[,val.col])
		return(dff)
	})
	hmplot = do.call(rbind, LS)
	logn = ifelse(grepl('q|fdr|Q|FDR', val.col), 'logQ', 'logP')
	if( log) plotMat = reshape2::acast(hmplot, ID~set, value.var=logn, fill = NA)
	if(!log) plotMat = reshape2::acast(hmplot, ID~set, value.var=val.col, fill = NA)
	plotMat = plotMat[order(apply(plotMat,1, sum, na.rm=T), decreasing=TRUE),]
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






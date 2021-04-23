## Updated pathway analysis set


## General hypergeo test
Gen_enrichment <- function(glist, refgmt, tglist, ncore=1) {
	if(ncore <= 1) {
		enrobj <- lapply(names())
	}
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
	is_log <- FALSE
	if(any(quantile(hmplot[,val.col], na.rm=TRUE) > 1)) is_log <- TRUE
	if(log & is_log) cat('val.col is already in log values.'); log <- FALSE
	if(log & !is_log) {
		hmplot$logV <- -log10(hmplot[,val.col])
		plotMat = reshape2::acast(hmplot, ID~set, value.var=logV, fill = NA)
	}

	if(!log) plotMat = reshape2::acast(hmplot, ID~set, value.var=val.col, fill = NA)
	plotMat = plotMat[order(apply(plotMat,1, sum, na.rm=TRUE), decreasing=TRUE),]
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




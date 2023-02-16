#' readGMT
#' 
#' gmt file reader. A GMT file format is a tab delimited text file containing gene sets.
#' Each line should contain one geneset, delimited by tab. Genes should start from 3rd field.
#' 
#' @param gmtfile GMT file path
#' @param as.df Return as data frame?
#' @param genelist List of gene set. Should be un-nested level one named list.
#' @param geneset_desc Description meta information for gene set. Either length one or a named vector of same length as glist.
#' @return Gene set dataframe of 2 column or list
#' @export
#' @examples
#' \dontrun{
#' readGMT('file/path.gmt')
#' }

readGMT <- function(gmtfile, as.df=FALSE) {
	x <- readLines(gmtfile)
	res <- strsplit(x, "\t")
	names(res) <- vapply(res, function(y) y[1], character(1))
	out <- lapply(res, "[", -c(1:2))
	if(as.df) {
		ont2gene <- stack(out)
		ont2gene <- ont2gene[, c("ind", "values")]
		colnames(ont2gene) <- c("ont", "gene")
		out <- ont2gene
	}
	return(out)
}


#' @describeIn readGMT write gmt file for geneset list.
#' @export

writeGMT <- function(gmtfile, genelist, geneset_desc='') {
	if( !(is.list(genelist) & all(sapply(genelist, is.character))) ) {
		stop('Check genelist format. genelist should be a one-level list of genesets.')
	}
	if(length(names(geneset_desc)) != 0) {
		geneset_desc <- geneset_desc[names(genelist)]
	}

	concat <- sapply(genelist, function(v) paste(v, collapse='\t'))
	out <- paste0(names(concat), '\t',geneset_desc,'\t', concat)

	writeLines(out, con=gmtfile)
}


#' Hypergeometric test for gmt style list
#' 
#' Gene list hypergeometric test against gmt format list of gene set. Can be used for custom pathway analysis or CMAP style query.
#' 
#' @param query  gene set to query (eg. Differentially Expressed Genes)
#' @param refGMT list of reference gene set (eg. Pathways)
#' @param gspace background gene space. Should contain all genes in query.
#' @param minGeneSet minimum size of gene set. This is used to filter refGMT. Default 10
#' @param ef.psc pseudocount when calculating enrichment factor (oddsRatio). Default 0
#' @param ncore number of cores used in hypergeoTestForGeneset2 (soon to be removed)
#' @return Data frame of gene set analysis results
#' \describe{
#' 	  \item{pVal}{: hypergeometric test p values from phyper}
#' 	  \item{logP}{: -log10(p value)}
#' 	  \item{oddsRatio}{: odds ratio}
#' 	  \item{tanco}{: tanimoto coefficient (Jaccard index)}
#' 	  \item{int}{: intersected item count}
#' 	  \item{gsRatio}{: gene set ratio (selected genes in gene set / selected genes)}
#' 	  \item{bgRatio}{: background ratio (total genes in gene set / total gene space)}
#' }
#' @export
#' @examples
#' gset = c('A','B','C')
#' glist = list(ID1 = LETTERS[1:10], ID2 = LETTERS[4:25])
#' hypergeoTestForGeneset(gset, glist, LETTERS)

hypergeoTestForGeneset <- function(query, refGMT, gspace, minGeneSet=10, ef.psc=0) {

	# match gene space
	if(!all(query %in% gspace)) {
		query <- intersect(query, gspace)
		msg <- paste(length(setdiff(query, gspace)),'query items were found outside of background space.')
		warning(msg)
	}
	if(length(query) == 0) stop('Query length is zero.')

	if(!all(unlist(refGMT) %in% gspace)) {
		refGMT <- lapply(refGMT, function(g) intersect(g, gspace))		
	}

	# filter refGMT with less than minimum gene set
	exc <- which(sapply(refGMT, length) < minGeneSet)
	if(length(exc) != 0) {
		if(length(exc) <= 5) {
			msg <- paste('Ref set no.', paste(exc, collapse=', '), 'had less than', minGeneSet,'genes and were excluded.')
		} else {
			msg <- paste(length(exc), ' entries in refGMT had less than', minGeneSet,'genes and were excluded.')
		}
		warning(msg)
		refGMT <- refGMT[which(sapply(refGMT, length) >= minGeneSet)]
	}
	if(length(refGMT) == 0) stop('Length of refGMT after filtering is zero.')

	# hypergeometric test
	N <- length(gspace)								# no of balls in urn
	k <- length(query)								# no of balls drawn from urn (DEG no)

	intscts <- sapply(refGMT, function(x) intersect(x, query))
	qs <- sapply(intscts, length)	# no of white balls drawn
	ms <- sapply(refGMT, length) 	# no of white balls in urn

	pvals <- phyper(qs - 1, ms, N - ms, k, lower.tail = FALSE)
	odds <- (qs + ef.psc) / (ms / N * k + ef.psc)
	jacc <- qs / sapply(refGMT, function(x) length(union(x, query)))
	gs.ratio <- paste0(qs,'/',k)
	bg.ratio <- paste0(ms,'/',N)

	enrRes <- data.table(
		ID=names(refGMT), pVal=pvals, oddsRatio=odds, 
		tan=jacc, int=qs, gsRatio=gs.ratio, bgRatio=bg.ratio, intGenes=intscts
	)

	enrRes$logP <- -log10(enrRes$pVal)

	# Adjust p-value while filtering out 1 values
	pv <- ifelse(enrRes$int == 0, NA, enrRes$pVal)
	enrRes$qVal <- p.adjust(pv, method='fdr')
	enrRes$qVal <- ifelse(enrRes$int == 0, 1, enrRes$qVal)
	enrRes$logQ <- -log10(enrRes$qVal)

	enrRes <- enrRes[,c('ID','pVal','logP','qVal','logQ','oddsRatio','tan','int','gsRatio','bgRatio','intGenes')]
	return(enrRes)	
}


#' @describeIn hypergeoTestForGeneset
#' Using multiprocessing (Deprecated)
#' @importFrom data.table rbindlist
#' @export

hypergeoTestForGeneset2 <- function (query, refGMT, gspace, minGeneSet=10, ncore = 1, ef.psc=0) {
	.Deprecated("hypergeoTestForGeneset")

	if(!all(query %in% gspace)) {
		stop(paste(length(setdiff(query, gspace)),'Query items were found outside of background space. Check inputs.'))
	}
	# query = intersect(query, gspace)
	refGMT <- parallel::mclapply(refGMT, function(g) intersect(g,gspace), mc.cores=ncore)

	if(length(query) == 0) stop('Query length is zero.')

	exc <- which(sapply(refGMT, length) < minGeneSet)
	if(length(exc) != 0) {
		if(length(exc) <= 5) {
			mesg <- paste('Ref set no', paste(exc, collapse=', '), 'had less than 10 genes and were excluded.')
		} else {
			mesg <- paste(length(exc), 'entries in refGMT had less than 10 genes and were excluded.')
		}
		warning(mesg)
		refGMT <- refGMT[which(sapply(refGMT, length) >= minGeneSet)]
	}
	if(length(refGMT) == 0) stop('Length of refGMT after filtering is zero.')

	N <- length(gspace)
	k <- length(query)
	enrRes <- mclapply(refGMT, function(refgenes) {
		q <- length(intersect(refgenes, query))
		m <- length(intersect(gspace, refgenes))
		I <- intersect(refgenes, query)

		pVal <- phyper(q - 1, m, N - m, k, lower.tail = FALSE)
		odds <- (q + ef.psc) / (m / N * k + ef.psc)
		jacc <- q / length(union(query, refgenes))
		gs.ratio <- paste0(q,'/',k)
		bg.ratio <- paste0(m,'/',N)

		return(data.table(pVal = pVal, oddsRatio=odds, tan = jacc, int=q, gsRatio=gs.ratio, bgRatio=bg.ratio, intGenes=list(I)))
	}, mc.cores = ncore)

	enrRes <- rbindlist(enrRes, idcol='ID')
	# enrRes$ID <- names(refGMT)
	enrRes$logP <- -log10(enrRes$pVal)

	pv <- ifelse(enrRes$int == 0, NA, enrRes$pVal)
	enrRes$qVal <- p.adjust(pv, method='fdr')
	enrRes$qVal <- ifelse(enrRes$int == 0, 1, enrRes$qVal)
	enrRes$logQ <- -log10(enrRes$qVal)

	enrRes <- enrRes[,c('ID','pVal','logP','qVal','logQ','oddsRatio','tan','int','gsRatio','bgRatio','intGenes')]
	return(enrRes)
}


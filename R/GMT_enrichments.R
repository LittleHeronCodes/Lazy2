#' readGMT
#' 
#' gmt file reader. A GMT file format is a tab delimited text file containing gene sets.
#' Each line should contain one geneset, delimited by tab. Genes should start from 3rd field.
#' 
#' @param gmtfile GMT file path
#' @param as.df Return as data frame?
#' @param glist List of gene set. Should be un-nested level one named list.
#' @param geneset_desc Description meta information for gene set. Either length one or same length vector as glist.
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

writeGMT <- function(gmtfile, glist, geneset_desc='') {
	if( !(is.list(glist) & all(sapply(glist,is.character))) ) {
		stop('Check glist format. glist should be a one-level list of genesets.')
	}

	concat <- sapply(glist, function(v) paste(v, collapse='\t'))
	out <- paste0(names(concat), '\t',geneset_desc,'\t', concat)

	writeLines(out, con=gmtfile)
}


#' Hypergeometric test for gmt style list
#' 
#' Gene list hypergeometric test against gmt format list of gene set. Can be used for custom pathway analysis or CMAP style query.
#' 
#' @param query  gene set to query (ex. Differentially Expressed Genes)
#' @param refGMT list of reference gene set (ex. Pathways)
#' @param gspace background gene space
#' @param minGeneSet minimum size of gene set. (used to filter refGMT)
#' @return data frame
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
#' \dontrun{
#' gset = c('A','B','C')
#' glist = list(ID1 = LETTERS[1:10], ID2 = LETTERS[4:25])
#' hypergeoTestForGeneset(gset, glist, LETTERS)
#' }

hypergeoTestForGeneset <- function(query, refGMT, gspace, minGeneSet=10, ef.psc=0) {
	require(data.table)

	if(!all(query %in% gspace)) {
		stop(paste(length(setdiff(query, gspace)),'query items were found outside of background space. Check inputs.'))
	}
	# query <- intersect(query, gspace)
	refGMT <- lapply(refGMT, function(g) intersect(g,gspace))

	if(length(query) == 0) stop('Query length is zero.')

	exc <- which(sapply(refGMT, length) < minGeneSet)
	if(length(exc) != 0) {
		if(length(exc) <= 5) {
			mesg <- paste('Ref set no', paste(exc, collapse=', '), 'had less than 10 genes and were excluded.')
		} else {
			mesg <- paste(length(exc), ' entries in refGMT had less than 10 genes and were excluded.')
		}
		warning(mesg)
		refGMT <- refGMT[which(sapply(refGMT, length) >= minGeneSet)]
	}
	if(length(refGMT) == 0) stop('Length of refGMT after filtering is zero.')

	N <- length(gspace)							# no of balls in urn
	k <- length(query)							# no of balls drawn from urn (DEG no)
	enrRes <- lapply(refGMT, function(refgenes) {
		q <- length(intersect(refgenes, query))	# no of white balls drawn
		m <- length(intersect(gspace, refgenes)) # no of white balls in urn
		I <- intersect(refgenes, query)

		pVal <- phyper(q-1, m, N-m, k, lower.tail = FALSE)
		odds <- (q + ef.psc) / (m / N * k + ef.psc)
		jacc <- q / length(union(query, refgenes))
		gs.ratio <- paste0(q,'/',k)
		bg.ratio <- paste0(m,'/',N)

		return(data.table(pVal=pVal, oddsRatio=odds, tan=jacc, int=q, gsRatio=gs.ratio, bgRatio=bg.ratio, intGenes=list(I)))
		})

	enrRes <- rbindlist(enrRes)
	enrRes$ID <- names(refGMT)
	enrRes$logP <- -log10(enrRes$pVal)

	pv <- ifelse(enrRes$int == 0, NA, enrRes$pVal)
	enrRes$qVal <- p.adjust(pv, method='fdr')
	enrRes$qVal <- ifelse(enrRes$int == 0, 1, enrRes$qVal)
	enrRes$logQ <- -log10(enrRes$qVal)

	enrRes <- enrRes[,c('ID','pVal','logP','qVal','logQ','oddsRatio','tan','int','gsRatio','bgRatio','intGenes')]
	# enrRes <- enrRes[,c('ID', 'pVal', 'logP', 'oddsRatio', 'tan', 'int', 'gsRatio', 'bgRatio', 'intGenes')]
	return(enrRes)	
}



#' @describeIn hypergeoTestForGeneset
#' Using multiprocessing
#' @export

hypergeoTestForGeneset2 <- function (query, refGMT, gspace, minGeneSet=10, ncore = 1, ef.psc=0) {
	require(parallel)
	require(data.table)

	if(!all(query %in% gspace)) {
		stop(paste(length(setdiff(query, gspace)),'Query items were found outside of background space. Check inputs.'))
	}
	# query = intersect(query, gspace)
	refGMT <- parallel::mclapply(refGMT, function(g) intersect(g,gspace), mc.cores=ncore)

	if(length(query) == 0) stop('Query length is zero.')

	exc <- which(sapply(refGMT, length) <= minGeneSet)
	if(length(exc) != 0) {
		if(length(exc) <= 5) {
			mesg <- paste('Ref set no', paste(exc, collapse=', '), 'had less than 10 genes and were excluded.')
		} else {
			mesg <- paste(length(exc), 'entries in refGMT had less than 10 genes and were excluded.')
		}
		warning(mesg)
		refGMT <- refGMT[which(sapply(refGMT, length) > minGeneSet)]
	}
	if(length(refGMT) == 0) stop('Length of refGMT after filtering is zero.')

	N <- length(gspace)
	k <- length(query)
	enrRes <- parallel::mclapply(refGMT, function(refgenes) {
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

	# enrRes = do.call(rbind, enrRes)
	enrRes <- rbindlist(enrRes)
	enrRes$ID <- names(refGMT)
	enrRes$logP <- -log10(enrRes$pVal)

	pv <- ifelse(enrRes$int == 0, NA, enrRes$pVal)
	enrRes$qVal <- p.adjust(pv, method='fdr')
	enrRes$qVal <- ifelse(enrRes$int == 0, 1, enrRes$qVal)
	enrRes$logQ <- -log10(enrRes$qVal)

	enrRes <- enrRes[,c('ID','pVal','logP','qVal','logQ','oddsRatio','tan','int','gsRatio','bgRatio','intGenes')]
	
	# enrRes <- enrRes[,c('ID', 'pVal', 'logP', 'oddsRatio', 'tan', 'int', 'gsRatio', 'bgRatio', 'intGenes')]
	return(enrRes)
}

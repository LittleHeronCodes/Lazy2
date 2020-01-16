#' readGMT
#' 
#' gmt file reader. A GMT file format is a tab delimited text file containing gene sets.
#' Each line should contain one geneset, delimited by tab. Genes should start from 3rd field.
#' 
#' @param file GMT file path
#' @param skip Fields to skip. Usually 2.
#' @param as.df Return as data frame?
#' @param glist List of gene set. Should be un-nested level one named list.
#' @param geneset_desc Description meta information for gene set. Either length one or same length vector as glist.
#' @return Gene set dataframe of 2 column or list
#' @export
#' @examples
#' \dontrun{
#' readGMT('file/path.gmt')
#' }


readGMT <- function(file, skip=2, as.df = TRUE) {

	out = list()
	con = file(file, 'r')
	while( TRUE ) {
		line = readLines(con, n=1)
		if( length(line) == 0 ) break
		sp = unlist(strsplit(line, split='\t'))
		out[[sp[1]]] = sp[-skip]
	}

	if(as.df) {
		out = do.call(rbind,  lapply(names(out), function(i) data.frame(key=i, value=out[[i]]) ))
	}
	close(con)
	return(out)
}


#' @describeIn readGMT write gmt file for geneset list.
#' @export

writeGMT <- function(file, glist, geneset_desc='') {
	if( !(is.list(glist) & all(sapply(glist,is.character))) ) {
		stop('Check glist format. glist should be a one-level list of genesets.')
	}

	concat = sapply(glist, function(v) paste(v, collapse='\t'))
	out = paste0(names(concat), '\t',geneset_desc,'\t', concat)

	writeLines(out, con=file)
}


#' hypergeoTestForGeneset
#' 
#' gene list hypergeometric test against gmt format gene set. Can be used for pathway analysis or CMAP style query.
#' 
#' @param query  gene set to query (ex Differentially Expressed Genes)
#' @param refGMT reference gene set in list format
#' @param gspace gene space to query in
#' @return data frame
#' @export
#' @examples
#' \dontrun{
#' hypergeoTestForGeneset(glist, c2.kegg)
#' }

hypergeoTestForGeneset <- function(query, refGMT, gspace) {
	if(!all(query %in% gspace)) {
		stop(paste(length(setdiff(query, gspace)),'query items were found outside of background space. Check inputs.'))
	}
	if(length(query) == 0) stop('query length zero.')

	N = length(gspace)							# no of balls in urn
	k = length(query)							# no of balls drawn from urn (DEG no)
	enrRes = lapply(refGMT, function(refgenes) {
		q = length(intersect(refgenes, query))	# no of white balls drawn
		m = length(intersect(gspace, refgenes)) # no of white balls in urn

		pVal = phyper(q-1, m, N-m, k, lower.tail = FALSE)
		odds = (q / k) / (m / N)
		jacc = q / length(union(query, refgenes))

		return(data.frame(pVal = pVal, oddsRatio=odds, tanco=jacc, int=q, bg=m))
		})

	enrRes = do.call(rbind, enrRes)
	enrRes$ID = names(refGMT)
	enrRes$logP = -log10(enrRes$pVal)
	enrRes = enrRes[,c('ID','pVal','logP','oddsRatio','tanco','int','bg')]
	return(enrRes)	
}


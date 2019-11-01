#' readGMT
#' 
#' gmt file reader
#' 
#' @param fn gmt file path
#' @param skip fields to skip. Usually 2.
#' @param as.df return as data frame? 
#' @return gene set dataframe of 2 column or list
#' @export
#' @examples
#' \dontrun{
#' readGMT('file/path.gmt')
#' }


readGMT <- function(fn, skip=2, as.df = TRUE) {

	out = list()
	con = file(fn, 'r')
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



#' hypergeoTestForGeneset
#' 
#' gene hypergeometric test against gmt format gene set
#' 
#' @param query gene set to query (ex Differentially Expressed Genes)
#' @param refGMT reference gene set in list format
#' @param gspace gene space to query in
#' @return data frame
#' @export
#' @examples
#' \dontrun{
#' hypergeoTestForGeneset(glist, c2.kegg)
#' }

hypergeoTestForGeneset <- function(query, refGMT, gspace) {
	N = length(gspace)						# no of balls in urn
	k = length(query)						# no of balls drawn from urn (DEG no)
	enrRes = lapply(refGMT, function(refgenes) {
		q = length(intersect(refgenes, query))	# no of white balls drawn
		m = length(intersect(gspace, refgenes)) # no of white balls in urn

		pVal = phyper(q-1, m, N-m, k, lower.tail = FALSE)
		odds = (q / k) / (m / N)
		return(data.frame(pVal = pVal, oddsRatio=odds, int=q, bg=m))
		})
	enrRes = do.call(rbind, enrRes)
	enrRes$ID = names(refGMT)
	enrRes$logP = -log10(enrRes$pVal)
	enrRes = enrRes[,c('ID','pVal','logP','oddsRatio','int','bg')]
	return(enrRes)	
}


#' CMAP Style Enrichment Calculation
#'
#' Calculate enrichment factor from a dataframe by rank cutoff.
#' @param df input data frame
#' @param enrich_element category column from df to enrich
#' @param rank_by column used to sort rank
#' @param cutoff cut percentile. default 0.05
#' @param min_count minimum count to be valid. default 3
#' @param psc pseudocount when caclulating enrichment. default 1
#' @return dataframe of named vectors for each category
#' @export

## Calculate enrichment factor for drugs
enrichmentFactorForDataFrame = function(df, enrich_element, rank_by, cutoff = 0.05, min_count = 3, psc=1) {

	# index to use
	include = names(which(table(df[, enrich_element]) >= min_count))
	if( length(include) < length(unique(df[, enrich_element]))*.10 ) {
		stop("Included factors less than 10%. There's a problem with input or lower min_count.")
	}
	df = df[which(df[, enrich_element] %in% include),]

	# order df by rank, get only within cut off
	rowIdx = which(rank(-df[,rank_by], ties.method = 'max') < nrow(df)*cutoff)

	# actual coverage
	cover = length(rowIdx) / nrow(df)

	# counts
	cntAll = table(df[, enrich_element])		# total for each
	cntPct = table(df[rowIdx, enrich_element])	# observed
	idx = names(cntPct)
	N = nrow(df)
	k = floor(nrow(df) * cutoff)

	# enrichment factor
	ef = (cntPct[idx]+psc) / ( cntAll[idx] * cutoff + psc)
	ef = sort(ef, decreasing = TRUE)

	# hypergeo test
	hp = phyper(cntPct[idx]-1, cntAll[idx], N-cntAll[idx], k, lower.tail = FALSE)
	hp = hp[names(ef)]
	# lhp = -log10(hp)

	# ratio
	rat = structure(paste0(cntPct[idx], '/', cntAll[idx]), names = idx)
	rat = rat[names(ef)]


	# add zero values
	idx2= setdiff(names(cntAll), idx)
	ef2 = structure(rep(0, length(idx2)), names=idx2)
	hp2 = structure(rep(1, length(idx2)), names=idx2)
	# lhp2 = structure(rep(0, length(idx2)), names=idx2)
	rat2 = structure(paste0(0, '/', cntAll[idx2]), names=idx2)

	out = list(EF= c(ef, ef2), HyperP=c(hp, hp2), ratio=c(rat, rat2), cover = cover)
	out = data.frame(out)

	## More info
	out$logP = -log10(out$HyperP)
	out$qVal = p.adjust(out$HyperP, method='fdr')
	out$logQ = -log10(out$qVal)
	
	out = out[,c('EF','HyperP','logP','qVal','logQ','ratio','cover')]
	return(out)
}




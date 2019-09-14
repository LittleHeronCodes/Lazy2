#' ent2sym
#' 
#' gene mapper for entrez and symbol
#' @param genes genes either in Entrez or Symbol
#' @return genes
#' @export

ent2sym <- function(genes) {
	if(!exists('geneInfo')) {
		geneInfo = Lazy.geneInfo
	}
	if(all(grepl('^[0-9]+$', genes))) {
		out = geneInfo$hgnc_symbol[match(genes, geneInfo$entrez)]		
	} else {
		out = geneInfo$entrez[match(genes, geneInfo$hgnc_symbol)]
	}
	return(as.character(out))
}


#' TtestWithMat
#'
#' Row-wise t-Test for gene expression matrix
#' @param m 
#' @param idx1 index to feed to first vector in t.test
#' @param idx2 index to feed to second vector in t.test
#' @param alternative alternative argument for t.test
#' @return Dataframe with results
#' t_stat : T statistic
#' t_pVal : t-test p value
#' medDiff : median difference
#' menDiff : mean difference
#' @export
#' @examples
#' TtestWithMat(M)

TtestWithMat <- function(m, idx1, idx2, alternative ='two.sided') {
	resultsDF = apply(m, 1, function(v) {
		md = median(v[idx1]) - median(v[idx2])
		mn = mean(v[idx1]) - mean(v[idx2])
		t.res = t.test(v[idx1], v[idx2], alternative=alternative)
		df = data.frame(
			t_stat = t.res$statistic, t_pVal = t.res$p.value,
			medDiff = md, menDiff = mn)
		return(df)
		})
	resultsDF = do.call(rbind, resultsDF)
	columns = colnames(resultsDF)
	rownames(resultsDF) = rownames(m)
	# resultsDF$geneEnt = rownames(m)
	# resultsDF$geneSym = ent2sym(rownames(m))
	# resultsDF = resultsDF[,c('geneEnt','geneSym', columns)]
	return(resultsDF)
}



#' geneCount
#' 
#' Lazy function for gene number for geneList
#' 
#' @param geneList geneList of DEG
#' @export

geneCount <- function(geneList) { sapply(geneList, function(ls) sapply(ls, length)) }

#' convertDEGList2Matrix
#' 
#' DEG List to deg Count
#' 
#' @param geneList gene list
#' @param toSpace total gene space
#' @export

convertDEGList2Matrix <- function(geneList, toSpace = NULL) {

	if(is.null(toSpace)) {
		toSpace = Reduce('union', lapply(geneList, function(ls) ls$toGene))
		cat('Total space inferred from list.\n')
	}

	degMat <- matrix(nrow = length(toSpace), ncol = length(geneList),
		dimnames=list(toSpace, names(geneList)))

	for(can in names(geneList)) {
		degMat[geneList[[can]]$toGene, can] =  0
		degMat[geneList[[can]]$upGene, can] =  1
		degMat[geneList[[can]]$dnGene, can] = -1
	}
	return(degMat)
}

#' DEGMatSumm
#' 
#' DEG List to deg Count
#' 
#' @param degMat degMat from convertDEGList2Matrix results
#' @param sort   sort or no?
#' @export

DEGMatSumm <- function(degMat, sort = FALSE) {

	degMetaDF = data.frame(
		geneEnt = rownames(degMat),
		geneSym = ent2sym(rownames(degMat)),
		upRecord  = apply(degMat, 1, function(v) sum(v ==  1, na.rm=TRUE)),
		dnRecord  = apply(degMat, 1, function(v) sum(v == -1, na.rm=TRUE)),
		anyRecord = apply(degMat, 1, function(v) sum(v !=  0, na.rm=TRUE)),
		sigRecord = apply(degMat, 1, function(v) sum(v, na.rm = TRUE)),
		Record    = apply(degMat, 1, function(v) sum(!is.na(v))) )

	if(sort) degMetaDF = degMetaDF[order(degMetaDF$anyRecord, decreasing=TRUE),]

	return(degMetaDF)
}


#' geneListSetOverlap
#' 
#' DEG List overlap 
#' 
#' @param geneList 
#' @return dataframe of overlap measures
#' @export

geneListSetOverlap <- function(geneList) {
	pairs = expand.grid(names(geneList), names(geneList))
	pairs2 = apply(pairs, 1, function(v) {
		A = geneList[[v[1]]]
		B = geneList[[v[2]]]
		T = intersect(A$toGene, B$toGene)

		int = length(T)
		ef_up = getEnrichmentFactor(A$upGene, B$upGene, T)
		ef_dn = getEnrichmentFactor(A$dnGene, B$dnGene, T)
		tan_up = tanimotoCoef(A$upGene, B$upGene)
		tan_dn = tanimotoCoef(A$dnGene, B$dnGene)
		return(data.frame(int, ef_up, ef_dn, tan_up, tan_dn))
		})
	pairs2 = do.call(rbind,pairs2)
	pairs = cbind(pairs, pairs2)
	return(pairs)
}


#' @section Section
#' @describeIn geneListSetOverlap For parallel processing

geneListSetOverlap.parallel <- function(geneList, ncore) {
	pairs = expand.grid(names(geneList), names(geneList))
	pairs2 = mclapply(1:nrow(pairs), function(i) {
		v = as.matrix(pairs[i,])
		A = geneList[[v[1]]]
		B = geneList[[v[2]]]
		T = intersect(A$toGene, B$toGene)

		ef_up = getEnrichmentFactor(A$upGene, B$upGene, T)
		ef_dn = getEnrichmentFactor(A$dnGene, B$dnGene, T)
		tan_up = tanimotoCoef(A$upGene, B$upGene)
		tan_dn = tanimotoCoef(A$dnGene, B$dnGene)
		return(data.frame(ef_up, ef_dn, tan_up, tan_dn))
		}, mc.cores = ncore)
	pairs2 = do.call(rbind,pairs2)
	pairs = cbind(pairs, pairs2)
	return(pairs)
}


#' @section Section
#' @describeIn geneListSetOverlap 
#' Covert pairs dataframe from geneListSetOverlap to matrix for distance

geneListDistMat <- function(pairs, value.var) {
	dm = reshape2::acast(pairs, Var1~Var2, value.var = value.var)	
	return(dm)
}


#' @describeIn geneListSetOverlap
#' Draw heatmap for distance matrix 

geneListDistMat_HM <- function(distMat, name, mtitle = NULL, log = TRUE, km = NULL) {
	if(log) distMat = log2(distMat + .01)
	if(is.null(mtitle)) mtitle = paste('Overlap Enrichment')
	colors = colorRamp2(c(-3.5,-2,0,2,3.5), 
		c('#2166ac','#4393c3','#f7f7f7','#d6604d','#b2182b'))
	hp = Heatmap(distMat, col=colors, name=name, column_title=mtitle, 
		row_km = km, column_km = km, show_column_names = FALSE, show_row_names = FALSE)
	return(hp)
}



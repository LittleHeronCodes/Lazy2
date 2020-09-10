#' gene mapper for entrez and symbol
#' 
#' Entrez to Hugo symbol and vice versa
#' 
#' @param genes genes either in Entrez or Symbol (human)
#' @return genes in wanted IDs
#' @export
#' @examples
#' ent2sym('TP53')

ent2sym <- function(genes) {
	geneMap = Lazy2::LazygeneInfo
	genes = as.character(genes)
	if(all(grepl('^[0-9]+$', genes))) {
		out = geneMap$hgnc_symbol[match(genes, geneMap$entrez)]		
	} else {
		out = geneMap$entrez[match(genes, geneMap$hgnc_symbol)]
	}
	return(as.character(out))
}


#' ???
#' 
#' gene synonym mapper
#' 
#' @inheritParams ent2sym
#' @param return_ent genes either in Entrez or Symbol
#' @param verbose    print mapped percentage
#' @return genes
# ' @export
#' @examples
#' geneAliasMap('CHOP')

geneAliasMap <- function(genes, return_ent = FALSE, verbose=FALSE) {
	geneAlias = Lazy2::geneAlias
	genes = as.character(genes)
	out = ifelse(genes %in% geneAlias$Symbol, genes, geneAlias$Symbol[match(genes, geneAlias$Synonyms)])
	if(verbose) {
		cat('Mapped', length(which(!is.na(out))), 'out of', length(genes), 'input mapped.\n')
	}
	if(return_ent) out = ent2sym(out)

	return(as.character(out))
}


#' removeZeroVar
#'
#' Remove genes with zero variances before t-test matrix
#' @param m    Matrix (eg. Gene expression matrix)
#' @param idx1 Column index for group A
#' @param idx2 Column index for group B
#' @param verbose verbose
#' @return Matrix with zero variance samples removed (hopefully)
#' @export

removeZeroVar <- function(m, idx1, idx2, verbose=TRUE) {
	varsum = apply(m, 1, function(v) var(v[idx1]) + var(v[idx2]))
	m = m[which(varsum != 0),]
	if(verbose) cat('No of zero variances detected :',sum(varsum==0), '\n')
	return(m)
}


#' TtestWithMat
#'
#' Row-wise t-Test for gene expression matrix
#' @param m    Gene expression matrix
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
#' \dontrun{
#' TtestWithMat(M)
#' }

TtestWithMat <- function(m, idx1, idx2, alternative ='two.sided', na.rm=TRUE) {
	m1 = m[,idx1]
	m2 = m[,idx2]
	md = apply(m1, 1, median, na.rm=na.rm) - apply(m2, 1, median, na.rm=na.rm)
	mn = apply(m1, 1, mean, na.rm=na.rm) - apply(m2, 1, mean, na.rm=na.rm)
	t.res = t(apply(m,1, function(v) {
		tt = t.test(v[idx1], v[idx2], alternative=alternative)
		return(c(t_stat=tt$statistic, t_pVal=tt$p.value))
		}))
	resultsDF= data.frame(t.res, medDiff=md, menDiff=mn)
	colnames(resultsDF) = gsub('.t$','',colnames(resultsDF))
	rownames(resultsDF) = rownames(m)

	return(resultsDF)
}


#' geneCount
#' 
#' Lazy function for gene number for geneList
#' 
#' @param geneList nested gene list of DEG
#' @export

geneCount <- function(geneList) { sapply(geneList, function(ls) sapply(ls, length)) }


#' convertDEGList2Matrix
#' 
#' Take a gene list and tally counts for each set. Used to find consensus DEG signatures.
#' 
#' @param gls list of genes (unnested)
#' @param tgls list of gene space corresponding to each gene set
#' @export
#' @examples
#' glist = list(
#' 	a = c(1,2,3,4,5,6,7,8,9),
#' 	b = c(1,2,4,5,8,9,11,13),
#' 	c = c(3,5,7,8,10,11,15))
#' tglist = list(a = 1:15, b = 1:14, c = 1:16)
#' convertDEGList2Matrix(glist, tglist)

convertDEGList2Matrix <- function(gls, tgls) {

	toSpace = Reduce(union, tgls)
	degMat <- matrix(nrow = length(toSpace), ncol = length(gls), dimnames=list(toSpace, names(gls)))

	for(id in names(gls)) {
		degMat[tgls[[id]], id] =  0
		degMat[gls[[id]], id] =  1
	}
	hitcnt = apply(degMat, 1, function(v) sum(v == 1, na.rm = TRUE))
	cnt = apply(degMat, 1, function(v) sum(!is.na(v)))
	pct = hitcnt / cnt

	out = data.frame(cbind(degMat, hitcnt = hitcnt, cnt = cnt, pct = pct))
	return(out)
}


#' DEGMatSumm
#' 
#' DEG List to deg Count
#' 
#' @param degMat degMat from convertDEGList2Matrix results
#' @param sort   sort or no?
# ' @export

DEGMatSumm <- function(degMat, sort = FALSE) {

	degMetaDF = data.frame(
		gene = rownames(degMat),
		upRecord  = apply(degMat, 1, function(v) sum(v ==  1, na.rm=TRUE)),
		dnRecord  = apply(degMat, 1, function(v) sum(v == -1, na.rm=TRUE)),
		anyRecord = apply(degMat, 1, function(v) sum(v !=  0, na.rm=TRUE)),
		sigRecord = apply(degMat, 1, function(v) sum(v, na.rm = TRUE)),
		Record    = apply(degMat, 1, function(v) sum(!is.na(v))) )

	if(sort) degMetaDF = degMetaDF[order(degMetaDF$anyRecord, decreasing=TRUE),]

	return(degMetaDF)
}


#' pairsOverlap
#' Calculate overlap metric for pairs
#' This is an internal function used for geneListSetOverlap
#' 
#' @inheritParams geneCount
#' @param v length 2 index vector from pair created by expand.grid
#' @return 1-row dataframe

pairsOverlap <- function(v, geneList) {
		A = geneList[[v[1]]]
		B = geneList[[v[2]]]
		T = intersect(A$toGene, B$toGene)

		int_up = length(intersect(A$upGene, B$upGene))
		int_dn = length(intersect(A$dnGene, B$dnGene))
		ef_up = getEnrichmentFactor(A$upGene, B$upGene, T)
		ef_dn = getEnrichmentFactor(A$dnGene, B$dnGene, T)
		tan_up = tanimotoCoef(A$upGene, B$upGene)
		tan_dn = tanimotoCoef(A$dnGene, B$dnGene)

		return(data.frame(int_up, int_dn, ef_up, ef_dn, tan_up, tan_dn))
}


#' geneListSetOverlap
#' 
#' DEG List overlap 
#' 
#' @inheritParams geneCount
#' @return dataframe of overlap measures
# ' @export
#' @examples
#' \dontrun{
#'     pairs <- geneListSetOverlap(geneList)
#'     distMat <- geneListDistMat(pairs, 'ef_dn')
#'     geneListDistMat(distMat)
#' }

geneListSetOverlap <- function(geneList) {
	pairs = expand.grid(names(geneList), names(geneList))
	pairs2 = apply(pairs, 1, function(v) pairsOverlap(v, geneList))
	pairs2 = do.call(rbind,pairs2)
	pairs = cbind(pairs, pairs2)
	return(pairs)
}

#
#' @describeIn geneListSetOverlap For parallel processing
#' @param ncore number of cores to use
# ' @export

geneListSetOverlap.parallel <- function(geneList, ncore) {
	pairs = expand.grid(names(geneList), names(geneList))
	pairs2 = parallel::mclapply(1:nrow(pairs), function(i) {
		v = as.matrix(pairs[i,])
		return(pairsOverlap(v, geneList))
		}, mc.cores = ncore)
	pairs2 = do.call(rbind,pairs2)
	pairs = cbind(pairs, pairs2)
	return(pairs)
}


#' geneList pairs to distance matrix
#' 
#' @describeIn geneListSetOverlap 
#' Covert pairs dataframe from geneListSetOverlap to matrix for distance
#' @param pairs pairs dataframe from geneListSetOverlap
#' @param value.var variable in pairs to use as value in matrix. Feeds to reshape2::acast
#' @export

geneListDistMat <- function(pairs, value.var) {
	dm = reshape2::acast(pairs, Var1~Var2, value.var = value.var)	
	return(dm)
}


#' @describeIn geneListSetOverlap
#' 
#' Draw heatmap for distance matrix
#' 
#' @param distMat distance matrix from geneListDistMat
#' @param name name to feed to ComplexHeatmap::Heatmap
#' @param show_axis show column and row names?
#' @param mtitle main title
#' @param log log transform or no?
#' @param km kmeans cluster feed for ComplexHeatmap::Heatmap
#' @export

geneListDistMat_HM <- function(distMat, name, show_axis = TRUE, mtitle = NULL, log = TRUE, km = NULL) {
	require(circlize)
	require(ComplexHeatmap)

	if(log) distMat = log2(distMat + .01)
	if(is.null(mtitle)) mtitle = paste('Overlap Enrichment')
	colors = circlize::colorRamp2(c(-3.5,-2,0,2,3.5), c('#2166ac','#4393c3','#f7f7f7','#d6604d','#b2182b'))
	hp = ComplexHeatmap::Heatmap(distMat, col=colors, name=name, column_title=mtitle, row_km = km, 
		column_km = km, show_column_names = show_axis, show_row_names = show_axis)
	return(hp)
}


#' cleanGList
#' 
#' Cleanup concatenated gene list to remove na, blanks, filter protein coding genes
#' Intended to use with microarray probe map from AnnotationDbi.
#' 
#' @param g Gene vector
#' @param concat String used to concatenate multiple mapped genes. Default '///'
#' @param filter Genes to filter out. Default NULL
#' @export

cleanGList <- function(g, concat = '///', filter = NULL) {
	g = unlist(strsplit(g, concat))
	g = g[which(!is.na(g) & g != '')]
	if(!is.null(filter)) g = intersect(g, filter)
	g = unique(g)
	return(g)
}




# source('/R/set_overlap.R')


ent2sym <- function(genes) {
	if(!exists('geneInfo')) {
		geneInfo = read.table(paste0(DIR_RESOURCE, '/DB/gene_info.csv'), sep = ',', header=TRUE)
	}
	if(all(grepl('^[0-9]+$', genes))) {
		out = geneInfo$hgnc_symbol[match(genes, geneInfo$entrez)]		
	} else {
		out = as.character(geneInfo$entrez[match(genes, geneInfo$hgnc_symbol)])
	}
	return(out)
}

# cellMap <- function(cells) {
# 	if(all(grepl('^ACH-', cells))) {
# 		out = cellInfo$CellName[match(cells, cellInfo$DepMap_ID)]
# 	} else {
# 		out = cellInfo$DepMap_ID[match(cells, cellInfo$CellName)]
# 	}
# 	return(out)
# }



TtestWithMat <- function(m, idx1, idx2, alternative ='two.sided') {
	resultsDF = apply(m, 1, function(v) {
		md = median(v[id1]) - median(v[id2])
		mn = mean(v[id1]) - mean(v[id2])
		t.res = t.test(v[id1], v[id2], alternative=alternative)
		df = data.frame(
			t_stat = t.res$statistic, t_pVal = t.res$p.value,
			medDiff = md, menDiff = mn)
		return(df)
		})
	resultsDF = do.call(rbind, resultsDF)
	columns = colnames(resultsDF)
	resultsDF$geneEnt = rownames(m)
	resultsDF$geneSym = ent2sym(rownames(m))
	resultsDF = resultsDF[,c('geneEnt','geneSym', columns)]
	return(resultsDF)
}


makeGeneList <- function(resDF) {
	glist = list(
		upGene = with(df, geneEnt[which(t_stat >=  0.05 & t_pVal < .01)]),
		dnGene = with(df, geneEnt[which(t_stat <= -0.05 & t_pVal < .01)]),
		toGene = df$geneEnt )
	return(glist)
}



# DEG Count
geneCount <- function(geneList) { sapply(geneList, function(ls) sapply(ls, length)) }


convertDEGList2Matrix <- function(geneLs, toSpace = NULL) {

	if(is.null(toSpace)) {
		toSpace = Reduce('union', lapply(geneLs, function(ls) ls$toGene))
		cat('Total space inferred from list.\n')
	}

	degMat <- matrix(nrow = length(toSpace), ncol = length(geneLs),
		dimnames=list(toSpace, names(geneLs)))

	for(can in names(geneLs)) {
		degMat[geneLs[[can]]$toGene, can] =  0
		degMat[geneLs[[can]]$upGene, can] =  1
		degMat[geneLs[[can]]$dnGene, can] = -1
	}
	return(degMat)
}


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



geneListSetOverlap <- function(geneList) {
	pairs = expand.grid(names(geneList), names(geneList))
	pairs2 = apply(pairs, 1, function(v) {
		A = geneList[[v[1]]]
		B = geneList[[v[2]]]
		T = intersect(A$toGene, B$toGene)

		ef_up = getEnrichmentFactor(A$upGene, B$upGene, T)
		ef_dn = getEnrichmentFactor(A$dnGene, B$dnGene, T)
		tan_up = tanimotoCoef(A$upGene, B$upGene)
		tan_dn = tanimotoCoef(A$dnGene, B$dnGene)
		return(data.frame(ef_up, ef_dn, tan_up, tan_dn))
		})
	pairs2 = do.call(rbind,pairs2)
	pairs = cbind(pairs, pairs2)
	return(pairs)
}



geneListSetOverlap2 <- function(geneList) {
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
		}, mc.cores = 20)
	pairs2 = do.call(rbind,pairs2)
	pairs = cbind(pairs, pairs2)
	return(pairs)
}



geneListDistMat <- function(value.var) {
	require(reshape2)
	dm = reshape2::acast(pairs, Var1~Var2, value.var = value.var)	
	return(dm)
}


geneListDistMat_HM <- function(distMat, name, mtitle = NULL, log = TRUE, km = NULL) {
	if(log) distMat = log2(distMat + .01)
	if(is.null(mtitle)) mtitle = paste('Overlap Enrichment')
	colors = colorRamp2(c(-3.5,-2,0,2,3.5), 
		c('#2166ac','#4393c3','#f7f7f7','#d6604d','#b2182b'))
	hp = Heatmap(distMat, col=colors, name=name, column_title=mtitle, 
		row_km = km, column_km = km,
		show_column_names = FALSE, show_row_names = FALSE)
	return(hp)
}



# tanimoto/enrichment factor distance + plot


# Build report(?)


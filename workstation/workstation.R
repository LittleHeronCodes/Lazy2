## Gene list extraction from result
## Supported formats : ebayes, treat, confect, DEseq2

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




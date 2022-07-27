
.enrich.internal <- function(refgenes)

hypergeoTestForGeneset <- function(query, refGMT, gspace, minGeneSet=10, ef.psc=0, ncore=1) {
	require(data.table)

	if(!all(query %in% gspace)) {
		stop(paste(length(setdiff(query, gspace)),'query items were found outside of background space. Check inputs.'))
	}

	if(ncore <= 1) {
		refGMT <- lapply(refGMT, function(g) intersect(g,gspace))
	} else {
		refGMT <- parallel::mclapply(refGMT, function(g) intersect(g,gspace), mc.cores=ncore)
	}

	if(length(query) == 0) stop('Query length is zero.')

	exc <- which(sapply(refGMT, length) < minGeneSet)
	if(length(exc) != 0) {
		if(length(exc) <= 5) {
			mesg <- paste('Ref set no.', paste(exc, collapse=', '), 'had less than', minGeneSet,'genes and were excluded.')
		} else {
			mesg <- paste(length(exc), ' entries in refGMT had less than', minGeneSet,'genes and were excluded.')
		}
		warning(mesg)
		refGMT <- refGMT[which(sapply(refGMT, length) >= minGeneSet)]
	}
	if(length(refGMT) == 0) stop('Length of refGMT after filtering is zero.')

	N <- length(gspace)								# no of balls in urn
	k <- length(query)								# no of balls drawn from urn (DEG no)

	if(ncore <= 1) {

		enrRes <- lapply(refGMT, function(refgenes) {
			q <- length(intersect(refgenes, query))		# no of white balls drawn
			m <- length(intersect(gspace, refgenes)) 	# no of white balls in urn
			I <- intersect(refgenes, query)

			pVal <- phyper(q-1, m, N-m, k, lower.tail = FALSE)
			odds <- (q + ef.psc) / (m / N * k + ef.psc)
			jacc <- q / length(union(query, refgenes))
			gs.ratio <- paste0(q,'/',k)
			bg.ratio <- paste0(m,'/',N)
			# enrgenes <- list(intersect(refgenes, query))

			return(data.table(pVal=pVal, oddsRatio=odds, tan=jacc, int=q, gsRatio=gs.ratio, bgRatio=bg.ratio, intGenes=list(I)))
			})

	} else {

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
	}

	enrRes <- rbindlist(enrRes)
	enrRes$ID <- names(refGMT)
	enrRes$logP <- -log10(enrRes$pVal)

	# Adjust p-value while filtering out 1 values
	pv <- ifelse(enrRes$int == 0, NA, enrRes$pVal)
	enrRes$qVal <- p.adjust(pv, method='fdr')
	enrRes$qVal <- ifelse(enrRes$int == 0, 1, enrRes$qVal)
	enrRes$logQ <- -log10(enrRes$qVal)

	enrRes <- enrRes[,c('ID','pVal','logP','qVal','logQ','oddsRatio','tan','int','gsRatio','bgRatio','intGenes')]
	# enrRes <- enrRes[,c('ID', 'pVal', 'logP', 'oddsRatio', 'tan', 'int', 'gsRatio', 'bgRatio', 'intGenes')]
	return(enrRes)	
}




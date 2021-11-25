## Extract DEG from Matrix (for cmap, kmap, large expression change data)

getDEGlistFromMatrix <- function(lfcm, pvalm, fco, pco) {
	if(all(lfcm > 0)) warning('fcm may not be in log scale!')
	if(any(pvalm > 1, na.rm=TRUE)) stop('pvalm should not be in log scale.')
	if(!identical(dimnames(lfcm), dimnames(pvalm))) stop('lfcm and pvalm should be in same dimensions.')

	sids <- structure(colnames(lfcm), names = colnames(lfcm))

	glist_up <- lapply(sids, function(id) names(which(lfcm[,id] >  log2(fco) & pvalm[,id] < pco)) )
	glist_dn <- lapply(sids, function(id) names(which(lfcm[,id] < -log2(fco) & pvalm[,id] < pco)) )

	tglist <- structure(rep(list(rownames(lfcm)), ncol(lfcm)), names = colnames(lfcm))

	return(list(up = glist_up, dn = glist_dn, to = tglist))
}



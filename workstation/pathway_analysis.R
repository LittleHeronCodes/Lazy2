## Clusterprofiler

library(Clusterprofiler)

enrichKEGG(gset, 
	organism = 'mmu',
	pvalueCutoff=1, pAdjustMethod='fdr', qvalueCutoff=1, minGSSize=20, maxGSSize=500, universe = gspace)





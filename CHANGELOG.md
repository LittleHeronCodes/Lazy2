# Changelog

## Ver 0.1.5.9000 [2022.05.11]
- [v0.1.5.9000] Functions : `Gen_enrichment` : Fix background gene space.
- [v0.1.5.9000] Functions : `getEnrichmentFactor` : Fix total space for setA and setB.
- [v0.1.5.9001] Functions : `hypergeoTestForGeneset` : minimum gene set filtering  
  - `refGMT[which(sapply(refGMT, length) > minGeneSet)]` --> `refGMT[which(sapply(refGMT, length) >= minGeneSet)]`

## Ver 0.1.4.9000 [2022.04.01]
- [v0.1.4.9000] Style fix : Package startup message in `zzz.R`. Show version number among other things.
- [v0.1.4.9000] Datasets  : `Lazygeneinfo` HGNC data update, rebuild (2021-11-12). (for real)
- [v0.1.4.9000] Datasets  : Add last update info on data help page for `LazygeneInfo`

## Ver 0.1.3.9000 [2021.11.25]
- [v0.1.3.9000] Functions : logP, p adjusting part moved from `Gen_enrichment` to `hypergeoTestForGeneset.
- [v0.1.3.9000] Datasets  : `Lazygeneinfo` HGNC data update, rebuild (2021-11-11)
- [v0.1.3.9001] Functions : `hypergeoTestForGeneset`fixed intGenes list returning problems.

## Ver 0.1.2
- [v0.1.2.9001] Functions : `hypergeoTestForGeneset` : fixed warning message.
- [v0.1.2.9000] Style fix : Remove all unused functions (move to scrap)
- [v0.1.2.9000] Style fix : assign operation into arrows
- [v0.1.2.9000] Style fix : banish `df` from variables
- [v0.1.2.9000] Functions : **Hypergeometric test results background set return check.**
- [v0.1.2.9000] Functions : **GSEA plot function using fgsea.** function modified from `plotEnrichment`
- [v0.1.2.9000] Functions : Draw MA plot.
- [v0.1.2.9000] Functions : ggplot2 custom themes.
- [v0.1.2.9000] Functions : `hypergeoTestForGeneset` : returns intersected genes
- [v0.1.2.9000] Functions : `Gen_enrichment` : return dataframe column reorder

## Ver 0.1.1
- [v0.1.1.9002] `Gen_enrichment` : pseudocount for Odds Ratio option pass to `hypergeoTestForGeneset`.
- [v0.1.1.9002] `hypergeoTestForGeneset` : pseudocount for Odds Ratio.
- [v0.1.1.9001] `hypergeoTestForGeneset` : added bgRatio, gsRatio info on output
- [v0.1.1.9001] `Gen_enrichment` : multiprocessing option uses `hypergeoTestForGeneset2`

## Ver 0.1.1.9000 [2021.05.24]

* `hypergeoTestForGeneset` : warning message fixed for better readibility when too many gene sets are filtered. Stop if length of reference set is zero after filtering.
* `Gen_enrichment` : `hypergeoTestForGeneset` has been running twice in parallel mode for some reason??
* Added dependency import `parallel`.

## Ver 0.1.1.9000 [2021.05.18]

* `drawMA`, `drawVol` symmetrical axis for logFC

## Ver 0.0.4.1 [2021.04.25]

* `hypergeoTestSetForGenes` now uses data.table


## Ver 0.0.3.3 [2021.02.21]

* Added more descriptions for each function
* Resurrected `convertDEGList2Matrix` to be used in finding consensus genes
* `readGMT` fix
* Returned `ClusterProfiler` dependency for KEGG, GO enrichments
* Last commit before branching for ver 0.0.4

### Functions changed

* Added : `extractGeneList`, `GO_enrichment`, `KEGG_enrichment`, `enrobj2Matrix`, `readGTF`, `ens2sym`
* Changed : `Gen_enrichment`, `readGMT`
* Removed : `TtestWithMat`

## Ver 0.0.3.2 [2020.07.23]

* `enrichmentFactorForDataFrame`
  * returns dataframe
  * calculates actual coverage
  * Added results : logP, FDR, logQ
  * fixed potentially critical error with hgeo ordering

## Ver 0.0.3.1 [2020.07.16]

* patch 1 : Fix ent2sym (geneMap to LazygeneInfo)

## Ver 0.0.3.0 test [2020.07.08]

* Last testing build before minor deploy
* `TanimotoCoef` take universe argument `T`
* Removed `ClusterProfiler` dependency
* Lazygeneinfo rebuilt from HGNC dataset and documented.
* Removed genemap, genealias dataset
* Fixed documentations
* Added `hypergeoTestSetForGenes` multiprocessing


## Ver 0.0.2.0 [2020.01.29]

* Changelog moved to separate file.

Function changes :
 * Added `na.rm` option in `TtestWithMat`
 * Dependent function namespace control.

New functions :
 * Lazier functions for faster and lazier data survey. **(NOT TO BE USED IN ACTUAL SCRIPT)**
 * removeZeroVariances
 

## Version 0.0.1.2 patch [2019.11.05]

* Modifications in ...??


## Version 0.0.1.2 [2019.11.01]

* `gene_info` data file change to gene info from NCBI (`NCBIHomo_sapiens.gene_info.csv`)
  * added gene alias data, updataed entrez symbol map


## Version 0.0.1.1 [2019.10.07]

* geneSet Overlap fixed to include intersection counts


## Version 0.0.1.0 [2019.09.26]

* Pathway analysis function : new enrLs to plot possible matrix function
* LazygeneInfo entrez ID to character type


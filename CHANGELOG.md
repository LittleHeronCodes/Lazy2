# Changelog

## Ver 0.0.4.0 [2020.09.10]

* Resurrected `convertDEGList2Matrix` to be used in finding consensus genes



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


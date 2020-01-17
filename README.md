# Lazy2
Custom functions from Lazy for personal use, now as a package.

## Build Guide (local)

```r
library(devtools)
build('Lazy2')
```

## Installation Guide

```r
.libPaths('path/to/alt/RLib')
devtools::install_github('LittleHeronCodes/Lazy2')
library(Lazy2, lib.loc = 'path/to/alt/RLib')
```

## Update Guide

```r
library(roxygen2)
library(devtools)

getwd()
#'path/to/alt/RLib/Lazy2'
document()
build()
install()

# reload package
reload(inst('Lazy2'))
```

**NOTE** : Be sure to run `DATASET.R` when adding new data.

## To Do Notes (For ver 0.1.3)

* Add startup message for version
* Add Python based venn diagram function
* Add draw MA plot, log fc result get for microarray no replicate

## In Progress

### ~Removing ClusterProfiler dependency~

-[ ] read gmt file
-[ ] hypergeometric t-test against gmt list object
-[ ] draw overrepresentation plot
-[ ] test consistency against ClusterProfiler

---

# Changelog

## Ver 0.0.1.3

* Added na.rm option in `TtestWithMat`
* New functions :
 * removeZeroVariances
 * Added Lazier functions for faster data survey


## Version 0.0.1.2 patch

**2019.11.05**

* Modifications in ...??


## Version 0.0.1.2

**2019.11.01**

* `gene_info` data file change to gene info from NCBI (`NCBIHomo_sapiens.gene_info.csv`)
  * added gene alias data, updataed entrez symbol map

**2019.10.07**

* Version 0.0.1.1
* geneSet Overlap fixed to include intersection counts


**2019.09.26**

* Version 0.0.1.0
* Pathway analysis function : new enrLs to plot possible matrix function
* LazygeneInfo entrez ID to character type


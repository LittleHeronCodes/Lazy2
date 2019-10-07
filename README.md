# Lazy2
Custom functions from Lazy for personal use, now as a package.

## Build Guide (local)

```R
library(devtools)
build('Lazy2')
```

## Installation Guide

```R
.libLoc('path/to/alt/RLib')
devtools::install_github('LittleHeronCodes/Lazy2')
```

## Update Guide

```R
library(roxygen2)
library(devtools)

getwd()
#'path/to/alt/RLib/Lazy2'
document()
build()
install()
```


## To Do Notes

* Add startup message for version
* Add Python based venn diagram function
* Add draw MA plot, log fc result get for microarray no replicate


## Update Notes

**2019.10.07**

* Version 0.0.1.1
* geneSet Overlap fixed to include intersection counts


**2019.09.26**

* Version 0.0.1.0
* Pathway analysis function : new enrLs to plot possible matrix function
* LazygeneInfo entrez ID to character type
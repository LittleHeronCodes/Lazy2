# Lazy2

Custom functions from Lazy for personal use, now as a package!

## To Do -- Version 0.1

- [o] Dependency clean up (so package won't load everything automatically).
- [o] Update Gene info table, increase reliability
- [o] Gene list default input changed to unnested list
- [o] Remove overly specific functions (branch off for legacy)
- [o] `readGMT` : doesn't work, fix by adapting `read.gmt` from `clusterProfiler`
- [o] Version numbering soft : *major.minor.build.revision*
- [] Remove all unused functions (export commented out)
- [] Style fix : assign operation into arrows
- [] Style fix : banish `df` from variables
- [] Package startup message
- [] `hypergeoTestset` : multiprocessing into option
- [] GSEA plot function using fgsea. function modified from plotEnrichment

## To Do -- future ideas

- [] Think about utilizing cpp codes for simple calculation function. --> speed trade-off potentially not worth it
- [] Mathematical equations into help page. (if possible)
- [] Add mouse gene info
- [] ~`reshape2::acast` into tidyverse like~
- [] ~Add bioMart info???~



---

## Testing Build

```r
getwd()	#'path/to/alt/RLib'
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

## Note to Self

* Be sure to run `DATASET.R` when adding new data.
* Test build on local before deployment.
* Reduce dependency.

## To Do Notes (For ver 0.2.0)

* Startup message for version in `zzz.R`.
* Python based venn diagram function
* Draw MA plot, log fc result get for microarray no replicate.
* Tutorials
* Work to make Lazy2 as stand alone as possible.
* Remove deprecated functions.


---

## Useful Links and References

* [R Namespaces](http://r-pkgs.had.co.nz/namespace.html)
* [Package structure by CRAN](https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Package-structure)
* [Writing R Packages by Karl Broman](https://kbroman.org/Tools4RR/assets/lectures/08_rpack_withnotes.pdf)


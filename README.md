# Lazy2

Custom functions from Lazy for personal use, now as a package!

## To Do

* Dependency clean up (so package won't load everything automatically).
* `clusterprofiler` function wrapper removed.
* Update Gene info table, increase reliability
* Add bioMart info???
* Remove overly specific functions (branch off for legacy)
* Version numbering soft : *major.minor.build.revision*

## Ver 0.0.3.0 **[2020.07.31]**

* `TanimotoCoef` take universe argument `T`

---

## Testing Build

```r
getwd()	#'path/to/alt/RLib/Lazy2'
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
* Version number : 

## To Do Notes (For ver 0.2.0)

* Startup message for version in `zzz.R`.
* Python based venn diagram function
* Draw MA plot, log fc result get for microarray no replicate.
* Tutorials
* Work to make Lazy2 as stand alone as possible.
* Remove deprecated functions.

## Checklist

-[] Add function namespace
-[] Test build with minimal dependency

---

## Useful Links and References

* [R Namespaces](http://r-pkgs.had.co.nz/namespace.html)
* [Package structure by CRAN](https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Package-structure)
* [Writing R Packages by Karl Broman](https://kbroman.org/Tools4RR/assets/lectures/08_rpack_withnotes.pdf)


# Lazy2

Custom functions from Lazy for personal use, now as a package!

## Version 0.1.2 To Do

- [x] Style fix : Remove all unused functions (move to scrap)
- [ ] Style fix : assign operation into arrows
- [ ] Style fix : banish `df` from variables
- [ ] Style fix : Package startup message in `zzz.R`.
- [ ] Functions : `hypergeoTestForGeneset` : multiprocessing into option
- [x] Functions : **Hypergeometric test results background set return check.**
- [x] Functions : **GSEA plot function using fgsea.** function modified from `plotEnrichment`
- [x] Functions : Draw MA plot, ~log fc result get for microarray no replicate.~
- [x] Functions : ggplot2 custom themes.
- [ ] Datasets  : Add mouse gene info


## Checklist before moving to 0.2.0

Ver 0.0.4.1 is final backward incompatible change. (deploy number : 0.1.0.9000)

- [ ] Test all functions vigorously.
- [x] All functions used frequently in last two projects.
- [x] Descriptions in all functions.
- [ ] Read last reference article on releasing R package.
- [x] Read up on gitFlow. 


## Version 0.2 Ideas

- [ ] ens2sym ??
- [ ] Think about utilizing cpp codes for simple calculation function. --> speed trade-off potentially not worth it
- [ ] Mathematical equations into help page. (if possible)


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

### Version Numbers

**major.minor.patch.dev**

- minor : function add, modify, not interchangeable with previous versions
- patch : bug fix, tweaks
- dev : ???


---

## Useful Links and References

* [R Namespaces](http://r-pkgs.had.co.nz/namespace.html)
* [Package structure by CRAN](https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Package-structure)
* [Writing R Packages by Karl Broman](https://kbroman.org/Tools4RR/assets/lectures/08_rpack_withnotes.pdf)
* [Releasing a package](https://r-pkgs.org/release.html)
* [A successful Git branching model](https://nvie.com/posts/a-successful-git-branching-model/)

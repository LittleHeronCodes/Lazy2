# Lazy2

Custom functions from Lazy for personal use, now as a package!


## Version 0.1.6 To Do

- [ ] ~Functions : `hypergeoTestForGeneset` : multiprocessing into option.~
- [x] Functions : `hypergeoTestForGeneset` : Faster computing.
- [x] Functions : `getEnrichmentFactor` : Fix total space for setA and setB.
- [ ] Functions : `extractGeneList` : p value column into option. (Also check if column actually exists)
- [x] Functions : `Gen_enrichment` : Fix background gene 
- [ ] Issues    : Pathway analysis gene space problem. Need to update gene space after removing pathways less than minimum gene space.
- [x] Issues    :R version 4 compatibility check
- [ ] Functions :ens2sym ??
- [ ] Datasets  : Add mouse gene info
- [ ] Datasets  : Update LazygeneInfo


## Checklist before moving to 0.2.0

Ver 0.2.0 is final backward incompatible change.

- [x] Test all functions vigorously.
- [x] All functions used frequently in last two projects.
- [x] Descriptions in all functions.
- [x] Read last reference article on releasing R package.
- [x] Read up on gitFlow. 


## Version 0.2 Ideas

- [ ] Think about utilizing cpp codes for simple calculation function. --> speed trade-off potentially not worth it
- [ ] Mathematical equations into help page. (if possible)


## Note to Self

* Be sure to run `DATASET.R` when adding new data.
* Package must be rebuilt to update manuals.
* Test build on local before deployment.
* Reduce dependency.


### Version Numbers

**major.minor.patch.dev**

- minor : function add, modify, new builds (manual updates)
- patch : minor bug fix, function tweaks
- dev : 9000. Indicating in development. 


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

---

## Useful Links and References

* [R Namespaces](http://r-pkgs.had.co.nz/namespace.html)
* [Package structure by CRAN](https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Package-structure)
* [Writing R Packages by Karl Broman](https://kbroman.org/Tools4RR/assets/lectures/08_rpack_withnotes.pdf)
* [Releasing a package](https://r-pkgs.org/release.html)
* [A successful Git branching model](https://nvie.com/posts/a-successful-git-branching-model/)
* [Exploring Package Startup messages](https://www.rostrum.blog/2021/08/27/zzz/)

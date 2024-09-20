# Lazy2

Custom functions from Lazy for personal use, now as a package!

## :loudspeaker: New changes in Version 0.2.0

(TBA)


## :book: Installation Guide

```r
.libPaths('path/to/alt/RLib')
devtools::install_github('LittleHeronCodes/Lazy2')
library(Lazy2, lib.loc = 'path/to/alt/RLib')
```

### :moyai: Install other version (for R 3)

(TBA)

## :arrow_up: Update Guide

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


## :pushpin: Useful Links and References

* [R Namespaces](http://r-pkgs.had.co.nz/namespace.html)
* [Package structure by CRAN](https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Package-structure)
* [Writing R Packages by Karl Broman](https://kbroman.org/Tools4RR/assets/lectures/08_rpack_withnotes.pdf)
* [Releasing a package](https://r-pkgs.org/release.html)
* [A successful Git branching model](https://nvie.com/posts/a-successful-git-branching-model/)
* [Exploring Package Startup messages](https://www.rostrum.blog/2021/08/27/zzz/)

# Developer Notes

## :boom:Version 0.2.0 major changes

- [ ] Develop in R 4.3.0
- [ ] Re-do readme to be more user centric
- [ ] Function name changes to be consistent and in accordance to tidyverse style guide (keep old function names deprecated)


## :zap:Version 0.2.0 minor changes

### :heavy_exclamation_mark:High priority

- [ ] Functions : `plotMA`, `plotVol` -- turn into one function with option. (Deprecate `plotVol`)
- [x] Functions : change the variable `oddsRatio` to `FE` as fold enrichment is the more accurate term.
- [x] Functions : `getEnrichmentFactor` change name to `get_fold_enrichment` for more accurate term.
- [x] Functions : Add overlap coefficient calculation to `getOverlapDF`
- [ ] Functions : Remove `log` option from `enrobj2Matrix`
- [ ] Functions : `extractGeneList` : p value column into option. (Also check if column actually exists)
- [ ] Issues    : Pathway analysis gene space problem. Need to update gene space after removing pathways less than minimum gene space.
- [ ] Functions : Add option to `Gen_enrichment` on whether to silently control for background space using union set of refgmt, with default TRUE for backwards compatibility.

### :grey_exclamation:Low priority

- [x] Style fix : Remove `return` expression as per tidyverse style guide, change to just declaring the object to return for my benefit
- [ ] Style fix : Change script name of `theme_transparent2` to `theme_transparent3` because why is that still 2.
- [ ] Datasets  : Add mouse gene info
- [ ] Datasets  : Update `LazygeneInfo`
- [ ] Functions : `ent2sym` -- better type checking
- [ ] Functions : `ens2sym`


## :bangbang: Checklist before moving to 0.2.0

- [ ] Test all functions vigorously.
- [ ] Test build in **R 4.3.0**
- [ ] Backward compatibility for breaking changes
- [ ] Update DESCRIPTION
- [ ] Check if deprecation works the way I think it would
- [x] All functions used frequently in last two projects.
- [x] Descriptions in all functions.
- [x] Read last reference article on releasing R package.
- [x] Read up on gitFlow. 


## :bulb:Version 0.2 Ideas

- [ ] Mathematical equations into help page. (if possible)
- [ ] Commit/changelog ideas : start using icons for commit messages


## :memo:Note to Self

* Be sure to run `DATASET.R` when adding new data.
* Package must be rebuilt to update manuals.
* Test build on local before deployment.
* Reduce dependency.


### :1234:Version Numbers

*major.minor.patch.dev*

- **major** : major version, incompatible changes, probably will never be 1
- **minor** : function add, modify, new builds (manual updates)
- **patch** : minor bug fix, function tweaks
- **dev** : 9000. Indicating in development. 


## :rotating_light:Build test

```r
getwd()	#'path/to/alt/RLib/Lazy2'
library(devtools)
build()
document()
check()
```

## :video_game: Commit icon

| Commit type              | Emoji                                                     |
| :----------------------- | :-------------------------------------------------------- |
| Version tag              | :bookmark: `:bookmark:`                                   |
| README update            | :memo: `:memo:`                                           |
| General update           | :zap: `:zap:`                                             |
| Improve format/structure | :sparkles: `:sparkles:`                                   |
| New feature/function     | :tada: `:tada:`                                           |
| Documentation            | :books: `:books:`                                         |
| Bugfix                   | :bug: `:bug:`                                             |
| Metadata                 | :card_index: `:card_index:`                               |
| Documenting source code  | :bulb: `:bulb:`                                           |
| Note to self             | :clipboard: `:clipboard:`                                         |
| Work in progress         | :construction: `:construction:`                           |
| Breaking changes         | :boom: `:boom:`                                           |
| Critical hotfix          | :ambulance: `:ambulance:`                                 |
| Deprecation              | :warning: `:warning:`                                         |
| Tests                    | :rotating_light: `:rotating_light:`                       |
| Adding a test            | :white_check_mark: `:white_check_mark:`                   |
| Make a test pass         | :heavy_check_mark: `:heavy_check_mark:`                   |
| Removing a dependency    | :heavy_minus_sign: `:heavy_minus_sign:`                   |
| Adding a dependency      | :heavy_plus_sign: `:heavy_plus_sign:`                     |
| Upgrading dependencies   | :arrow_up: `:arrow_up:`                                   |
| Downgrading dependencies | :arrow_down: `:arrow_down:`                               |
| Performance              | :racehorse: `:racehorse:`                                 |
| Reverting changes        | :rewind: `:rewind:`                                       |
| Merging branches         | :twisted_rightwards_arrows: `:twisted_rightwards_arrows:` |
| Deploying stuff          | :rocket: `:rocket:`                                       |
| Funsies                  | :video_game: `:video_game:`                               |




<!-- README.md is generated from README.Rmd. Please edit that file -->

## Installation

You can download the development version from GitHub with:

``` r
if(!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("jmbadia/MS2ID", build_vignettes = TRUE, force=T)
```

## Introduction

This package annotates MS2 query spectra by using an in-house database.
The creation of such database from (mostly publicly) available resources
(HMDB, ChEBI, PubChem…) relies on the use of the
[CompoundDb](https://github.com/EuracBiomedicalResearch/CompoundDb)
package developed by J. Stanstrup and J. Rainer.

For more information visit the package
[website](https://jmbadia.github.io/MS2ID/articles/MS2ID_intro.html) or
[vignette.](https://github.com/jmbadia/MS2ID/tree/master/vignettes)

Discussions and suggestions are welcome:
<https://github.com/jmbadia/MS2ID/issues>

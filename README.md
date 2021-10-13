
<!-- README.md is generated from README.Rmd. Please edit that file -->

# semEff

<!-- badges: start -->

[![Repo
Status](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/)
[![Lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg?label=Lifecycle)](https://lifecycle.r-lib.org/articles/stages.html)
[![Licence](https://img.shields.io/badge/License-GPL3-green.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
![GitHub language
count](https://img.shields.io/github/languages/count/murphymv/semEff?label=Languages)
[![R-CMD-check](https://github.com/murphymv/semEff/workflows/R-CMD-check/badge.svg)](https://github.com/murphymv/semEff/actions)

[![CRAN](https://www.r-pkg.org/badges/version/semEff?color=blue)](https://CRAN.R-project.org/package=semEff)
![Downloads: Total](http://cranlogs.r-pkg.org/badges/grand-total/semEff)
![Downloads: Monthly](https://cranlogs.r-pkg.org/badges/semEff)

[![Donate](https://img.shields.io/badge/PayPal-Donate%20to%20Author-yellow.svg)](https://paypal.me/murphymv1)

<!-- badges: end -->

`semEff` provides functionality to automatically calculate direct,
indirect, and total effects for ‘piecewise’ structural equation models,
comprising lists of fitted models representing structured equations
(Shipley 2000, 2009; Lefcheck 2016). Confidence intervals are provided
via bootstrapping.

Currently supported model classes are `"lm"`, `"glm"`, `"lmerMod"`,
`"glmerMod"`, `"lmerModLmerTest"`, `"gls"`, and `"betareg"`.

## Installation

You can install the released version of `semEff` from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("semEff")
```

And the development version from [GitHub](https://github.com/) with:

``` r
devtools::install_github("murphymv/semEff@dev")
```

## Usage

The primary function is
[`semEff()`](https://murphymv.github.io/semEff/reference/semEff.html),
which returns an object of class `"semEff"` with associated `print()`
and `summary()` methods. Everything can technically be accomplished in a
single call to `semEff()`; however, since bootstrapping is employed to
generate resamples for confidence intervals (via
[`bootEff()`](https://murphymv.github.io/semEff/reference/bootEff.html)),
it is usually preferable to save these estimates separately prior to
calling `semEff()` – allowing more flexibility and saving time if/when
recalling the function.

## Examples

Package functions are well-documented and most include some short
examples. In addition, see the following vignettes for some longer
demonstrations:

-   [Analysing direct vs. indirect effects of landscape location on
    plant species
    richness](https://murphymv.github.io/semEff/articles/semEff.html)

-   [Predicting and plotting indirect effects of degree days to bud
    burst on tree
    growth](https://murphymv.github.io/semEff/articles/predicting-effects.html)

## References

Lefcheck, J. S. (2016). piecewiseSEM: Piecewise structural equation
modelling in R for ecology, evolution, and systematics. *Methods in
Ecology and Evolution*, **7**(5), 573–579. <https://doi.org/f8s8rb>

Shipley, B. (2000). A New Inferential Test for Path Models Based on
Directed Acyclic Graphs. *Structural Equation Modeling: A
Multidisciplinary Journal*, **7**(2), 206–218. <https://doi.org/cqm32d>

Shipley, B. (2009). Confirmatory path analysis in a generalized
multilevel context. *Ecology*, **90**(2), 363–368.
<https://doi.org/bqd43d>

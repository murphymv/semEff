
<!-- README.md is generated from README.Rmd. Please edit that file -->

# semEff

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/murphymv/semEff.svg?branch=master)](https://travis-ci.org/murphymv/semEff)
<!-- badges: end -->

semEff Provides functionality to automatically calculate direct,
indirect, and total effects from ‘piecewise’ structural equation models,
comprising lists of fitted models representing structured equations
(with local estimation). Confidence intervals are provided via
bootstrapping.

Currently supported model classes are “lm”, “glm”, and “merMod”.

## Installation

You can install the released version of semEff from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("semEff")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("murphymv/semEff")
```

## Example

``` r
# install.packages(c("semEff", "ggplot2"))
library(semEff)
library(ggplot2)

## Simulated data from Shipley (2009) on tree growth and survival (see ?Shipley)
head(Shipley)
#>   site tree      lat year     Date       DD   Growth  Survival Live
#> 1    1    1 40.38063 1970 115.4956 160.5703 61.36852 0.9996238    1
#> 2    1    2 40.38063 1970 118.4959 158.9896 43.77182 0.8433521    1
#> 3    1    3 40.38063 1970 115.8836 159.9262 44.74663 0.9441110    1
#> 4    1    4 40.38063 1970 110.9889 161.1282 48.20004 0.9568525    1
#> 5    1    5 40.38063 1970 120.9946 157.3778 50.02237 0.9759584    1
#> 6    1    1 40.38063 1972 114.2315 160.6120 56.29615 0.9983398    1

## Hypothesised SEM:
## latitude -> degree days to bud burst -> date of burst -> growth -> survival
lapply(Shipley.SEM, formula)
#> $DD
#> DD ~ lat + (1 | site) + (1 | tree)
#> 
#> $Date
#> Date ~ DD + (1 | site) + (1 | tree)
#> 
#> $Growth
#> Growth ~ Date + (1 | site) + (1 | tree)
#> 
#> $Live
#> Live ~ Growth + (1 | site) + (1 | tree)

# ## Bootstrap model effects (takes a while...)
# system.time(
#   Shipley.SEM.Boot <- bootEff(Shipley.SEM, ran.eff = "site", seed = 53908)
# )

## Calculate SEM effects and CI's (use bootstrapped SEM)
eff <- suppressWarnings(semEff(Shipley.SEM.Boot))

## Summary of effects for response "Growth"
eff$Summary$Growth
#> $Direct
#>           Date
#> Estimate 0.382
#> Lower CI 0.289
#> Upper CI 0.513
#>              *
#> 
#> $Indirect
#>            lat     DD
#> Estimate 0.165 -0.240
#> Lower CI 0.088 -0.351
#> Upper CI 0.290 -0.180
#>              *      *
#> 
#> $Total
#>            lat     DD  Date
#> Estimate 0.165 -0.240 0.382
#> Lower CI 0.088 -0.351 0.289
#> Upper CI 0.290 -0.180 0.513
#>              *      *     *
#> 
#> $Mediators
#>             DD   Date
#> Estimate 0.165 -0.075
#> Lower CI 0.088 -0.105
#> Upper CI 0.290 -0.048
#>              *      *

## Extract total effects for Growth
tot <- totEff(eff)[["Growth"]]
tot.b <- totEff(eff, type = "boot")[["Growth"]]

## Predict effects for "Date" (direct) and "DD" (indirect) on Growth
mod <- Shipley.SEM$Growth
dat <- na.omit(Shipley)
fit <- sapply(c("Date", "DD"), function(i) {
  x <- data.frame(seq(min(dat[i]), max(dat[i]), length = 100)); names(x) <- i
  c(x, predEff(mod, newdata = x, effects = tot[i], eff.boot = tot.b))
}, simplify = FALSE)

## Function to plot predictions
plotFit <- function(x, y, fit, x.lab = NULL, y.lab = NULL) {
  x2 <- fit[[1]]; f <- fit[[2]]; ci.l <- fit[[3]]; ci.u <- fit[[4]]
  ggplot () + 
    geom_point(aes(x, y)) +
    geom_ribbon(aes(x2, ymin = ci.l, ymax = ci.u), fill = "blue", alpha = "0.15") +
    geom_line(aes(x2, f), color = "blue", size = 1) +
    xlab(x.lab) + ylab(y.lab) +
    theme_bw() + theme(legend.position = "none")
}

## Direct effects of Date
plotFit(x = dat$Date, y = dat$Growth, fit = fit$Date, 
        x.lab = "Date of Bud Burst", y.lab = "Stem Growth")
```

<img src="man/figures/README-example-1.png" width="100%" />

``` r

## Indirect effects of DD (operating via Date)
plotFit(x = dat$DD, y = dat$Growth, fit = fit$DD, 
        x.lab = "Degree Days to Bud Burst", y.lab = "Stem Growth")
```

<img src="man/figures/README-example-2.png" width="100%" />

``` r

## Huge amount of scatter around each fit as random effects explain most
## variation in this model! Compare conditional vs. marginal R-squared:
round(c(R2c = R2(mod)[[1]], R2m = R2(mod, re.form = NA)[[1]]), 3)
#>   R2c   R2m 
#> 0.794 0.048
```

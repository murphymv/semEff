## semEff 0.3.0.9000

### New features:

* Added `offset` argument to functions `getY()` and `R2()`, to retain/remove the
offset variable in/from the response variable or fitted values (where
applicable).
* Added `incl.raw` argument to `stdCoeff()`, to append raw (unstandardised)
coefficients to the output. This allows simultaneous bootstrapping of
standardised and raw coefficients, and the latter can then also be used for
calculating (`semEff(..., use.raw = TRUE)`) or predicting (`predEff(..., use.raw
= TRUE)`) effects.

### Other changes:
* Added confidence interval attributes to `bootCI()`/`semEff()` output.
* `R2` no longer calculates predictive R-squared for GLMMs, as the
interpretation of the hat matrix is not reliable (see
https://rdrr.io/cran/lme4/man/hatvalues.merMod.html).
* Removed ability to pass arguments from `getY()` to `glt()`, allowing better
control over output of `getY(..., link = TRUE)`.
* Minor updates to function code and documentation, improvement and addition of
some new internal functions.

### Bugs fixed:

* `bootEff()` specified with correlated errors failed for mixed models of class
`"lmerModLmerTest"` (issue with `update()` inside the function).
* `predEff()` failed to evaluate some complex model terms (e.g. polynomials).
* `stdCoeff()` did not re-fit model properly to calculate correct VIFs for a
fully 'centred' model (i.e. did not account sufficiently for complex terms such
as polynomials or transformations, where centring should occur as the final
step).
* `xNam()` generated incorrect term names for categorical predictors under
certain circumstances (different contrast types, interactive effects without
main effects).
* `stdCoeff()` incorrectly calculated 'centred' intercept for models including
an offset.
* `predEff()` failed when a nested list of models and list of numeric weights
were supplied (i.e. a model averaging scenario).


## semEff 0.3.0

### New features:

* Support for mixed models of class `"lmerModLmerTest"`.
* New function `glt()`, for calculating 'generalised' link transformations for
non-gaussian variables.

### Other changes:

* Transfer of some functionality from `getY()` to `glt()`.
* Minor changes to arguments in `bootEff()` and `getY()`.
* Added ability in `stdCoeff()` to use variables not present in the model design
matrix (e.g. a 'missing' main effect for an interaction).
* Added ability to pass a boot object (from `bootEff()`) to the `effects`
argument of `predEff()`.
* Added a `refit.x` argument to `stdCoeff()`, allowing control over whether to
refit the model with centred predictors (for correct VIFs).
* Various updates to documentation.

### Bugs fixed:

* `xNam()` did not generate correct term names for categorical variables with
contrast types other than `contr.treatment()`.
* `stdCoeff()` did not correctly adjust for multicollinearity for a model
containing categorical variables when centring was specified (`cen.x = TRUE`).
* `getY()` failed to generate an estimated working response when a variable with
missing values (`NA`) was supplied (this functionality now in `glt()`).
* `predEff()` failed for models with categorical variables (did not access dummy
variables in model matrix).


## semEff 0.2.1

### Bugs fixed:

* Function `semEff()` did not output effects properly.


## semEff 0.2.0

### New features:

* Added support for generalised least squares models (class `"gls"`).
* Added support for beta regression models (class `"betareg"`).

### Bugs fixed:

* Function `xNam()` did not generate correct term names for interactions
involving multi-coefficient terms (e.g. factors).
* Function `xNam()` did not generate correct term names for factors when the
model intercept is suppressed.
* Function `R2()` with argument `pred = TRUE` threw an error for models where
any weights = 0.


## semEff 0.1.0

New package `semEff`, allowing the automatic calculation of effects for
'piecewise' structural equation models.


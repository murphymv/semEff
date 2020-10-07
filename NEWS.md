## semEff 0.4.0.9000


## semEff 0.4.0

### New features:

* Standardised vs raw effects: added `incl.raw` argument to `stdEff()`
(`stdCoeff()`), to append raw effects (unstandardised coefficients) to the
output. This facilitates simultaneous bootstrapping of both sets of effects,
allowing raw effects to be used alternatively for calculating (`semEff(...,
use.raw = TRUE)`) or predicting (`predEff(..., use.raw = TRUE)`) effects/CIs.

### Other changes:

* Renamed function `stdCoeff()` to `stdEff()`, to better reflect the concept of
standardised model coefficients as 'effects' (calling `stdCoeff()` will still
work - with a warning - until the next version at least).
* Added `offset` argument to `getY()` and `R2()`, to explicitly retain/remove an
offset (where present) in/from the response variable or fitted values. Offsets
are removed by default, which ensures, for example, that standardised effects
are scaled appropriately.
* Added `env` argument to multiple functions, for explicitly specifying the
location of data used to fit models (not necessary in most circumstances). This
replaces the `...` argument in many instances, which was previously used to pass
an environment to `eval()` (via `getData()`). `env` (and `data`) can also now be
passed (`...`) to `bootEff()` and `predEff()`.
* Added confidence interval attributes to `bootCI()`/`semEff()` output (i.e.
confidence level, type).
* `R2()` no longer calculates predictive R-squared for GLMMs, as the
interpretation of the hat matrix used in calculations is not reliable (see
<https://rdrr.io/cran/lme4/man/hatvalues.merMod.html>).
* Removed ability to pass arguments from `getY()` to `glt()`, allowing more
controlled output of `getY(..., link = TRUE)`.
* Various minor updates to function code and documentation, improvement and
addition of some new internal helper functions.

### Bugs fixed:

* `bootEff()` specified with correlated errors failed for mixed models of class
`"lmerModLmerTest"` (issue with re-fitting models using `update()`).
* `predEff()` failed to evaluate some complex model terms (e.g. polynomials).
* `stdEff()` (`stdCoeff()`) did not re-fit model properly to calculate correct
VIFs for a fully 'centred' model (i.e. did not account sufficiently for complex
terms such as polynomials or transformations, where mean-centring should occur
as the final step).
* `xNam()` generated incorrect term names for categorical predictors under
certain circumstances (different contrast types, interactive effects with no
'main' effects).
* `stdEff()` (`stdCoeff()`) incorrectly calculated 'centred' intercept for
models with an offset specified.
* `predEff()` failed when a nested list of models and list of numeric weights
were supplied (i.e. a model averaging scenario).
* `stdEff()` (`stdCoeff()`) did not return the 'phi' parameter(s) for beta
regression models.


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


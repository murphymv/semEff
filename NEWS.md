## semEff 0.3.0

New features:

* Support for mixed models of class "lmerModLmerTest"
* New function 'glt' for calculating link transformations for non-gaussian
variables

Other changes:

* Minor changes to argument order in 'bootEff' and 'getY'
* Added ability in 'stdCoeff' to use variables not present in the model design
matrix (e.g. a 'missing' main effect for an interaction)
* Added ability to pass a boot object (from 'bootEff') to the 'effects' argument
of 'predEff'
* Added a 'refit.x' argument to 'stdCoeff' - allowing control over whether to
refit the model with centred predictors (for correct VIFs)
* Various updates to documentation

Bugs fixed:

* 'xNam' did not generate correct term names for categorical variables with
contrast types other than 'contr.treatment'
* 'stdCoeff' did not correctly adjust for multicollinearity for a model
containing categorical variables when centring was specified ('cen.x = TRUE')
* 'getY' failed to generate an estimated working response when a variable with
missing values (NAs) was supplied to 'mod'
* 'predEff' failed for models with categorical variables (did not access dummy
variables in model matrix)

## semEff 0.2.1

Bugs fixed:

* Function 'semEff' does not output effects properly

## semEff 0.2.0

New features:

* Added support for generalised least squares models (class "gls")
* Added support for beta regression models (class "betareg")

Bugs fixed:

* Function 'xNam' does not generate correct term names for interactions
involving multi-coefficient terms (e.g. factors)
* Function 'xNam' does not generate correct term names for factors when the
model intercept is suppressed
* Function 'R2' with argument 'pred = TRUE' throws an error for models where any
weights = 0

## semEff 0.1.0

New package 'semEff', allowing the automatic calculation of effects for
'piecewise' structural equation models.

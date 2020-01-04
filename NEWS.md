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

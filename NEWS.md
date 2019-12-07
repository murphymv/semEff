## semEff 0.1.0.9000

New features:

* Added support for generalised least squares models (class "gls")

Bugs fixed:

* 'xNam' does not generate correct term names for interactions involving
multi-coefficient terms (e.g. factors)
* 'xNam' does not generate correct term names for factors in models where the
intercept is suppressed
* 'R2' with argument 'pred = TRUE' does not work for models with zero weights



## semEff 0.1.0
New package 'semEff', allowing the automatic calculation of effects for
'piecewise' structural equation models.

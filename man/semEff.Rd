% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/semEff-fun.R
\name{semEff}
\alias{semEff}
\title{SEM Effects}
\usage{
semEff(
  sem,
  predictors = NULL,
  mediators = NULL,
  excl.other.med = FALSE,
  use.raw = FALSE,
  ci.conf = 0.95,
  ci.type = "bca",
  digits = 3,
  bci.arg = NULL,
  ...
)
}
\arguments{
\item{sem}{A piecewise SEM, comprising a list of fitted model objects or of
boot objects (containing bootstrapped model effects). Alternatively, a
\code{"psem"} object from
\href{https://rdrr.io/cran/piecewiseSEM/man/psem.html}{\code{piecewiseSEM::psem()}}.
If list is unnamed, response variable names will be used.}

\item{predictors, mediators}{Names of variables for/through which to calculate
effects. If \code{NULL} (default), all predictors/mediators in the SEM will be
used.}

\item{excl.other.med}{Logical, whether to exclude other SEM mediators from
calculating indirect effects, i.e., those not specified in the \code{mediators}
argument. Useful for examining individual effect pathways with only the
specified mediators, rather than including all paths involving them
(default). Ignored if \code{mediators = NULL}.}

\item{use.raw}{Logical, whether to use 'raw' (unstandardised) effects for all
calculations (if present in \code{sem}).}

\item{ci.conf}{A numeric value specifying the confidence level for confidence
intervals on effects.}

\item{ci.type}{The type of confidence interval to return (defaults to \code{"bca"}
– see Details). See \code{\link[boot:boot.ci]{boot::boot.ci()}} for further specification details.}

\item{digits}{The number of decimal places to return for numeric values (for
summary tables).}

\item{bci.arg}{A named list of any additional arguments to \code{\link[boot:boot.ci]{boot::boot.ci()}},
excepting argument \code{index}.}

\item{...}{Arguments to \code{\link[=bootEff]{bootEff()}}.}
}
\value{
A list object of class \code{"semEff"} for which several methods and
extractor functions are available. Contains:
\enumerate{
\item Summary tables of effects and confidence intervals
\item All effects
\item All bootstrapped effects
\item All indirect effects (individual, not summed)
}
}
\description{
Automatically calculate direct, indirect, total, and mediator
effects for endogenous (response) variables in a 'piecewise' structural
equation model (SEM).
}
\details{
The eponymous function of this package calculates all direct,
indirect, total, and mediator effects for a 'piecewise' structural equation
model (SEM), that is, one where parameter estimation is local rather than
global (Lefcheck, 2016; Shipley, 2000, 2009). The SEM simply takes the form
of a list of fitted models, or bootstrapped estimates from such models,
describing hypothesised causal pathways from predictors to response
('endogenous') variables. These are either direct, or operate indirectly
via other response variables ('mediators'). This list should represent a
directed ('acyclic') causal model, which should be named exactly for each
response variable and ordered from 'upstream' or 'causal' variables through
to 'downstream' (i.e. those at the end of the pathway). If \code{sem} is a list
of fitted models, effects will first be bootstrapped using \code{\link[=bootEff]{bootEff()}}
(this may take a while!).

Direct effects are calculated as fully standardised model coefficients for
each response variable (see \code{\link[=stdEff]{stdEff()}} for details), while indirect effects
are the product of these direct effects operating along causal pathways in
the SEM. The total effects of any given predictor on a response are then
the sum of its direct and (all) its indirect effects. 'Mediator effects'
are also calculated, as the sum of all indirect paths which operate through
each individual mediator – useful to assess the relative importance of
different mediators in affecting the response. All of these effect types
can be calculated automatically for all (default) or for a specified subset
of predictors and/or mediators in the SEM. As indirect, total, and mediator
effects are not directly bootstrapped using the fitted models for response
variables (i.e. via \code{\link[=bootEff]{bootEff()}}), their equivalent 'bootstrapped' estimates
are calculated instead using each bootstrapped direct effect.

Confidence intervals for all effects are returned in summary tables for
each response (see \code{\link[=bootCI]{bootCI()}}), with BC\emph{a} intervals calculated by default
using the bootstrapped estimates for each effect type (Cheung, 2009; Hayes
& Scharkow, 2013; MacKinnon et al., 2004). Effects for which the confidence
intervals do not contain zero are highlighted with a star (i.e.
'significant' at the \code{ci.conf} level). Bootstrap standard errors (standard
deviations of the samples) and biases (sample means minus original
estimates) are also included. Correlated errors (and confidence intervals)
are also returned if their bootstrapped values are present in \code{sem}, or if
they are specified to argument \code{cor.err} or as part of a \code{"psem"} object
(see \code{\link[=bootEff]{bootEff()}}). These represent residual relationships among response
variables, unaccounted for by the hypothesised SEM paths. Use \code{summary()}
for effect summary tables and \code{print()} to return a table of variable names
and associated details.

All calculated effects and bootstrapped effects are also returned in lists
for each response variable, with all except mediator effects also including
the model intercept(s) – required for prediction (these will be zero for
ordinary linear models with fully standardised effects). Effects can be
conveniently extracted with \code{\link[=getEff]{getEff()}} and related functions.
}
\examples{
# SEM effects
(shipley.sem.eff <- semEff(shipley.sem.boot))
summary(shipley.sem.eff)

# Effects for selected variables
summary(shipley.sem.eff, response = "Live")
# summary(semEff(shipley.sem.boot, predictor = "lat"))
# summary(semEff(shipley.sem.boot, mediator = "DD"))

# Effects calculated using original SEM (models)
# (not typically recommended – better to use saved boot objects)
# system.time(
#  shipley.sem.eff <- semEff(shipley.sem, R = 1000, seed = 13,
#                            ran.eff = "site")
# )
}
\references{
Cheung, M. W. L. (2009). Comparison of methods for constructing
confidence intervals of standardized indirect effects. \emph{Behavior Research
Methods}, \emph{41}(2), 425-438. \doi{10/fnx7xk}

Hayes, A. F., & Scharkow, M. (2013). The Relative Trustworthiness of
Inferential Tests of the Indirect Effect in Statistical Mediation Analysis:
Does Method Really Matter? \emph{Psychological Science}, \emph{24}(10), 1918-1927.
\doi{10/bbhr}

Lefcheck, J. S. (2016). piecewiseSEM: Piecewise structural equation
modelling in \code{R} for ecology, evolution, and systematics. \emph{Methods in
Ecology and Evolution}, \emph{7}(5), 573-579. \doi{10/f8s8rb}

MacKinnon, D. P., Lockwood, C. M., & Williams, J. (2004). Confidence Limits
for the Indirect Effect: Distribution of the Product and Resampling
Methods. \emph{Multivariate Behavioral Research}, \emph{39}(1), 99. \doi{10/chqcnx}

Shipley, B. (2000). A New Inferential Test for Path Models Based on
Directed Acyclic Graphs. \emph{Structural Equation Modeling: A Multidisciplinary
Journal}, \emph{7}(2), 206-218. \doi{10/cqm32d}

Shipley, B. (2009). Confirmatory path analysis in a generalized multilevel
context. \emph{Ecology}, \emph{90}(2), 363-368. \doi{10/bqd43d}
}

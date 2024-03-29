% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stdEff-fun.R
\name{getY}
\alias{getY}
\title{Get Model Response Variable}
\usage{
getY(mod, data = NULL, link = FALSE, offset = FALSE, env = NULL)
}
\arguments{
\item{mod}{A fitted model object, or a list or nested list of such objects.}

\item{data}{An optional dataset, used to first refit the model(s).}

\item{link}{Logical. If \code{TRUE}, return the GLM response variable on the link
scale (see Details).}

\item{offset}{Logical. If \code{TRUE}, include model offset(s) in the response.}

\item{env}{Environment in which to look for model data (if none supplied).
Defaults to the \code{\link[=formula]{formula()}} environment.}
}
\value{
A numeric vector comprising the response variable on the original or
link scale, or an array, list of vectors/arrays, or nested list.
}
\description{
Extract the response variable from a fitted model on the
original or link scale.
}
\details{
\code{getY()} will return the response variable from a model by summing
the fitted values and the response residuals. If \code{link = TRUE} and the
model is a GLM, the response is transformed to the link scale. If this
results in undefined values, an estimate based on the 'working' response
variable of the GLM is returned instead (see \code{\link[=glt]{glt()}}).

Any offset variables are subtracted from the response by default. This
means that, for example, rates rather than raw counts will be returned for
poisson GLMs (where applicable).
}
\examples{
# All SEM responses (original scale)
head(getY(shipley.sem))

# Estimated response in link scale from binomial model
head(getY(shipley.sem$Live, link = TRUE))
}



#' @title Weighted Variance
#' @description Calculate the weighted variance of \code{x}.
#' @param x A numeric vector from which to calculate the weighted variance.
#' @param w A numeric vector of weights of the same length as \code{x}.
#' @param na.rm Logical, should NA's in \code{x} be removed?
#' @param ... Not currently used.
#' @details This function calculates the weighted variance of \code{x} via the
#'   weighted covariance matrix (\code{cov.wt}). If no weights are supplied, the
#'   simple variance is returned instead. As in \code{weighted.mean},
#'   \code{NA}'s in \code{w} are not handled specially and will return \code{NA}
#'   as result.
#' @return A numeric value, the weighted variance of \code{x}.
#' @seealso \code{\link[stats]{var}}, \code{\link[stats]{cov.wt}},
#'   \code{\link[stats]{weighted.mean}}
#' @examples
#' ## Weighted variance
#' x <- rnorm(30)
#' w <- runif(30, 0, 1)
#' varW(x, w)
#'
#' ## Simple variance
#' stopifnot(varW(x) == var(x))
#'
#' ## NA handling
#' varW(c(x[1:29], NA), w, na.rm = TRUE)  # NA in x (removed)
#' varW(c(x[1:29], NA), w, na.rm = FALSE)  # NA in x (NA returned)
#' varW(x[1:29], w = c(w[1:29], NA))  # NA in w (NA returned)
#' @export
varW <- function(x, w = NULL, na.rm = FALSE, ...) {
  if (!is.null(w)) {
    x.na <- is.na(x)
    w.na <- is.na(w)
    if (!(any(x.na) && !na.rm || any(w.na))) {
      x <- cbind(x[!x.na])
      w <- as.vector(w)[!x.na]
      cov.wt(x, w)$cov[1, ]
    } else NA
  } else var(x, na.rm = na.rm)
}


#' @title Weighted Standard Deviation
#' @description Calculate the weighted standard deviation of \code{x}.
#' @param ... Arguments to \code{varW}.
#' @details This is a simple wrapper for \code{varW}, which applies the square
#'   root to the output.
#' @return A numeric value, the weighted standard deviation of \code{x}.
#' @seealso \code{\link[stats]{sd}}, \code{\link[semEff]{varW}}
#' @export
sdW <- function(...) {
  sqrt(varW(...))
}


#' @title Get Model Term Names
#' @description Extract term names from a fitted model object as a character
#'   vector or list.
#' @param m A fitted model object of class \code{lm}, \code{glm}, or
#'   \code{merMod}, or a list or nested list of such objects.
#' @param intercept Logical, whether the intercept term should be included.
#' @param aliased Logical, whether names of aliased terms should be included
#'   (see Details).
#' @param list Logical, whether names should be returned as a (named) list, with
#'   all multi-coefficient terms grouped under their term (variable) name.
#' @param ... Not currently used.
#' @details This function can be used to extract term names from a fitted model.
#'   Names of terms for which coefficients cannot be estimated are also included
#'   if \code{aliased = TRUE} (default). These may be terms which are perfectly
#'   correlated with other terms in the model, so that the model design matrix
#'   is rank deficient.
#' @return A character vector or named list of term names.
#' @examples
#' ## Term names from Shipley SEM
#' m <- Shipley.SEM
#' xNam(m)
#' xNam(m, intercept = FALSE)  # only predictors
#'
#' ## Model with different types of predictor (some multi-coefficient terms)
#' x1 <- poly(rnorm(100), 2)  # polynomial
#' x2 <- as.factor(rep(c("a", "b", "c", "d"), each = 25))  # categorical
#' x3 <- rep(1, 100)  # no variation
#' m <- lm(rnorm(100) ~ x1 + x2 + x3)
#' xNam(m)
#' xNam(m, aliased = FALSE)  # drop term that cannot be estimated (x3)
#' xNam(m, aliased = FALSE, list = TRUE)  # as named list
#' @export
xNam <- function(m, intercept = TRUE, aliased = TRUE, list = FALSE, ...) {

  ## Function
  xNam <- function(m) {

    ## All names as list
    xn <- labels(terms(m))
    xn2 <- rownames(summary(m)$coef)
    xn <- c(xn2[isInt(xn2)], xn)
    mf <- model.frame(m)
    XN <- sapply(xn, function(i) {
      if (i %in% names(mf)) {
        xi <- mf[, i]
        xic <- class(xi)
        if ("factor" %in% xic) {
          paste0(i, levels(xi)[-1])
        } else {
          if ("matrix" %in% xic) {
            paste0(i, colnames(xi))
          } else i
        }
      } else i
    }, simplify = FALSE)

    ## Drop intercept?
    if (!intercept) XN <- XN[!isInt(names(XN))]

    ## Drop aliased terms?
    if (!aliased) {
      XN <- XN[sapply(XN, function(i) all(i %in% xn2))]
    }

    ## Return as list?
    if (!list) unlist(unname(XN)) else XN

  }

  ## Apply recursively
  rMapply(xNam, m, SIMPLIFY = FALSE)

}


#' @title Get Model Data
#' @description Extract the data used to fit a model.
#' @param m A fitted model object of class \code{lm}, \code{glm}, or
#'   \code{merMod}, or a list or nested list of such objects.
#' @param subset Logical. If \code{TRUE}, rows with missing observations
#'   (\code{NA}) are removed - if those variable(s) were used to fit model(s).
#' @param merge Logical. If \code{TRUE}, and \code{m} is a list or nested list,
#'   all unique data columns are merged into one data frame, in the order in
#'   which they are encountered.
#' @param ... Arguments to \code{eval}.
#' @details This is a simple function to extract the data used to fit a model,
#'   by evaluating the data slot of the model call object. If the model was fit
#'   without using the \code{data} argument, data is returned via
#'   \code{model.frame} instead (with a warning). If \code{m} is a list and
#'   \code{merge = TRUE}, a single dataset containing all (unique) variables
#'   used to fit models is returned. This will return an error if \code{subset =
#'   TRUE} results in datasets with different numbers of observations (rows).
#' @return A data frame of the variables used to fit the model(s).
#' @seealso \code{\link[base]{eval}}, \code{\link[stats]{getCall}},
#'   \code{\link[stats]{model.frame}}
#' @examples
#' ## Get data used to fit SEM from Shipley (2009)
#' getData(Shipley.SEM[[1]])  # from single model
#' getData(Shipley.SEM)  # from SEM (list of datasets)
#' getData(Shipley.SEM, merge = TRUE)  # from SEM (combined dataset)
#' @export
getData <- function(m, subset = FALSE, merge = FALSE, ...) {

  ## Function
  getData <- function(m) {
    d <- eval(getCall(m)$data, ...)
    mf <- model.frame(m)
    if (is.null(d)) {
      d <- mf
      wo <- as.character(getCall(m))[c("weights", "offset")]
      names(d)[names(d) %in% c("(weights)", "(offset)")] <- wo
      names(d) <- gsub("^offset\\((.*)\\)$", "\\1", names(d))
      warning("Model frame returned. Variable names may not match originals.")
    }
    if (subset) {
      w <- weights(m)
      if (!is.null(w)) mf <- mf[w > 0, ]
      d[rownames(mf), ]
    } else d
  }

  ## Apply recursively
  rgetData <- function(m) {
    if (isList(m)) {
      d <- lapply(m, rgetData)
      if (merge) {
        d <- do.call(cbind, unname(d))
        d[unique(names(d))]
      } else d
    } else getData(m)
  }
  rgetData(m)

}


#' @title Get Model Response Variable
#' @description Extract the response variable from a fitted model in the
#'   original or link scale.
#' @param m A fitted model object of class \code{lm}, \code{glm}, or
#'   \code{merMod}, or alternatively a numeric vector corresponding to a
#'   variable to be transformed. \code{m} can also be a list or nested list of
#'   such objects.
#' @param family Optional, the error distribution family containing the link
#'   function which will be used to transform the response (see
#'   \code{\link[stats]{family}} for specification details).
#' @param data An optional dataset used to first re-fit the model(s) (if
#'   \code{m} is a model object/list).
#' @param link Logical. If \code{TRUE}, return the GLM response variable on the
#'   link scale (see Details).
#' @param ... Not currently used.
#' @details This function will return the response variable from a model by
#'   summing the fitted values and the response residuals. If \code{link = TRUE}
#'   and the model is a GLM, the response is then transformed using the link
#'   function. However, if any values are undefined, this transformation is
#'   replaced by an estimate based on the 'working' response variable of the GLM
#'   (see below). The function can also be used to transform a variable
#'   (supplied to \code{m}) using the link function from the specified
#'   \code{family} - in which case the \code{link} argument is ignored.
#'
#'   \strong{Estimating the link-transformed response}
#'
#'   A key challenge in generating fully standardised model coefficients for a
#'   generalised linear model (GLM) with a non-gaussian link function is the
#'   frequent inability to calculate appropriate standardised ranges (typically
#'   the standard deviation) for the response variable in the link scale. This
#'   is because directly transforming the response will often produce undefined
#'   values. Although methods for circumventing this issue by indirectly
#'   estimating the standard deviation of the link-transformed response have
#'   been proposed - including a latent-theoretic approach for binomial models
#'   (McKelvey & Zavoina 1975) and a more general variance-based method using
#'   pseudo R-squared (Menard 2011) - here an alternative approach is used.
#'   Where transformed values are undefined, the function will instead return
#'   the synthetic 'working' response arising from the iteratively reweighted
#'   least squares (IRLS) algorithm of the GLM (McCullagh & Nelder 1989). This
#'   is reconstructed as the sum of the linear predictor and the working
#'   residuals - with the latter comprising the error of the model in the link
#'   scale. The advantage of this approach is that a relatively straightforward
#'   'transformation' of any non-gaussian response is readily attainable in all
#'   cases. The standard deviation (or other relevant range) can then be
#'   calculated using values of the transformed response and used to scale the
#'   coefficients. An additional benefit for piecewise SEM's is that the
#'   transformed rather than original response can then be specified as a
#'   predictor in other models, ensuring that standardised indirect and total
#'   effects are calculated correctly (i.e. using the same standard deviation).
#'
#'   To ensure a high level of 'accuracy' in the working response - in the sense
#'   that the inverse-transformed values are practically indistinguishable from
#'   the original response - the function uses the following iterative fitting
#'   procedure to calculate a 'final' working response:
#'
#'   \enumerate{\item A new GLM is fit with the original response as the sole
#'   predictor, with the same error family, and using a single IWLS iteration.
#'   \item The working response is calculated from this model \item The inverse
#'   transformation of the working response is then calculated \item If the
#'   inverse transformation is effectively equal to the original response
#'   (testing using \code{all.equal} with default tolerance), the working
#'   response is returned; otherwise, the GLM is re-fit with the working
#'   response now as the sole predictor, and steps 2-4 are repeated - each time
#'   with an additional IWLS iteration in the model}
#'
#'   This approach will generate a very reasonable transformation of the
#'   response variable, which will closely resemble the direct
#'   link-transformation where this can be compared - see Examples. It also
#'   ensures that the transformed values, and hence the standard deviation, are
#'   the same for any GLM fitting the same response - provided it uses the same
#'   link function - and so facilitates model comparisons, selection, and
#'   averaging.
#' @note As we often cannot directly observe the response variable on the link
#'   scale, any method estimating its values or statistics will be wrong to some
#'   degree. The aim should be to try to minimise this error as far as
#'   (reasonably) possible, while also generating standardised coefficients
#'   whose interpretation most closely resembles those of the ordinary linear
#'   model - something which the current method achieves. The solution of using
#'   the working response from the GLM to scale coefficients is a purely
#'   practical, but reasonable one, and one that takes advantage of modern
#'   computing power to minimise error through iterative model fitting. An added
#'   bonus is that the estimated variance is constant across models, which
#'   cannot be said of previous methods (Menard 2011). The overall approach
#'   would be classed as 'observed-empirical' by Grace \emph{et al.} (2018), as
#'   it utilises model error variance (the working residuals) rather than
#'   theoretical distribution-specific variance.
#' @return A numeric vector comprising the response variable in the original or
#'   link scale, or an array, list or nested list of such vectors.
#' @references Grace, J.B., Johnson, D.J., Lefcheck, J.S. and Byrnes, J.E.K.
#'  (2018) Quantifying relative importance: computing standardized effects in
#'  models with binary outcomes. \emph{Ecosphere} \strong{9}, e02283.
#'  \url{https://doi.org/10.1002/ecs2.2283}
#'
#'  McCullagh P. and Nelder, J. A. (1989) \emph{Generalized Linear Models} (2nd
#'  Edition). London: Chapman and Hall.
#'
#'  McKelvey, R. D., & Zavoina, W. (1975). A statistical model for the analysis
#'  of ordinal level dependent variables. \emph{The Journal of Mathematical
#'  Sociology}, \strong{4}(1), 103-120.
#'  \url{https://doi.org/10.1080/0022250x.1975.9989847}
#'
#'  Menard, S. (2011) Standards for Standardized Logistic Regression
#'  Coefficients. \emph{Social Forces} \strong{89}, 1409-1428.
#'  \url{https://doi.org/10.1093/sf/89.4.1409}
#' @seealso \code{\link[stats]{glm.fit}}, \code{\link[base]{all.equal}}
#' @examples
#' ## SEM responses (original scale)
#' getY(Shipley.SEM)
#'
#' ## Estimated response in link scale from binomial model
#' m <- Shipley.SEM$Live
#' getY(m, link = TRUE)
#' getY(m, link = TRUE, family = binomial("probit"))  # diff. link function
#'
#' ## Same estimate calculated using variable instead of model
#' y <- Shipley$Live
#' getY(y, binomial)
#'
#' ## Compare working response with a direct link transformation
#' ## (test with a poisson model, log link)
#' set.seed(1)
#' y <- rpois(30, lambda = 10)
#' y2 <- y
#' m <- suppressWarnings(
#'   glm(y ~ y2, poisson, control = list(maxit = 1))
#' )
#' i <- 0
#' repeat {
#'   yl <- predict(m) + resid(m, "working")
#'   yli <- family(m)$linkinv(yl)
#'   eql <- isTRUE(all.equal(yli, y, check.names = FALSE))
#'   if (eql) return(yl) else {
#'     i <- i + 1
#'     m <- suppressWarnings(
#'       update(m, . ~ yl, control = list(maxit = i))
#'     )
#'   }
#' }
#'
#' ## Effectively equal?
#' all.equal(yl, log(y), check.names = FALSE)
#' # TRUE
#'
#' ## Actual difference?
#' all.equal(yl, log(y), check.names = FALSE, tolerance = .Machine$double.eps)
#' # "Mean relative difference: 1.05954e-12"
#' @export
getY <- function(m, family = NULL, data = NULL, link = FALSE, ...) {

  f <- family; d <- data

  ## Function
  getY <- function(m) {

    ## Update model with any supplied data
    mod <- isMod(m)
    if (mod && !is.null(d)) {
      m <- eval(update(m, data = d, evaluate = FALSE))
    }

    ## Model response
    y <- if (mod) {
      y <- fitted(m) + resid(m, "response")
      w <- weights(m)
      if (!is.null(w)) y[w > 0] else y
    } else {
      y <- as.numeric(m)
      setNames(y, names(m))
    }

    ## Return in original or link scale
    if (isGlm(m) && link || !mod && !is.null(f)) {

      ## Error family
      if (is.character(f)) {
        f <- get(f, mode = "function", envir = parent.frame())
      }
      if (is.function(f)) f <- f()
      if (is.null(f)) f <- family(m)

      ## Transform response to link scale
      yl <- f$linkfun(y)

      ## Return the transformed (or working) response
      if (any(is.infinite(yl))) {
        y2 <- y
        suppressWarnings(
          m <- do.call(glm, list(y ~ y2, f, control = list(maxit = 1)))
        )
        i <- 0
        repeat {
          yl <- predict(m) + resid(m, "working")
          yli <- f$linkinv(yl)
          eql <- isTRUE(all.equal(yli, y, check.names = FALSE))
          if (eql) return(yl) else {
            i <- i + 1
            suppressWarnings(
              m <- update(m, . ~ yl, control = list(maxit = i))
            )
          }
        }
      } else yl

    } else y

  }

  ## Apply recursively
  rMapply(getY, m)

}


#' @title Generalised Variance Inflation Factors
#' @description Calculate (generalised) variance inflation factors for terms in
#'   a fitted model via the variance-covariance matrix of coefficients.
#' @param m A fitted model object of class \code{lm}, \code{glm}, or
#'   \code{merMod}, or a list or nested list of such objects.
#' @param data An optional dataset used to first re-fit the model(s).
#' @param ... Not currently used.
#' @details This function calculates generalised variance inflation factors
#'   (GVIF) as described in Fox & Monette (1992), and also implemented in the
#'   \code{vif} function in the \pkg{car} package. However, whereas \code{vif}
#'   returns both GVIF and GVIF^(1/(2*Df)) values, function \code{VIF} simply
#'   returns the squared result of the latter measure, which equals the standard
#'   VIF for single-coefficient terms and is the equivalent measure for
#'   multi-coefficient terms (e.g. categorical or polynomial). Also, while
#'   \code{vif} returns values per model term (i.e. predictor variable),
#'   \code{VIF} returns values per coefficient, meaning that the same VIF will
#'   be returned per coefficient for multi-coefficient terms. Finally, \code{NA}
#'   is returned for any coefficients which could not be estimated in the model
#'   (e.g. aliased terms).
#' @return A numeric vector of the VIF's, or an array, list or nested list of
#'   such vectors.
#' @references Fox, J. and Monette, G. (1992) Generalized Collinearity
#'   Diagnostics. \emph{Journal of the American Statistical Association}
#'   \strong{87}, 178-183. \url{https://doi.org/10.2307/2290467}
#' @seealso
#' \href{https://www.rdocumentation.org/packages/car/versions/3.0-3/topics/vif}{
#' vif (web)}
#' @examples
#' ## Model with two correlated terms
#' m <- Shipley.Growth[[3]]
#' VIF(m)  # Date & DD somewhat correlated
#' VIF(update(m, . ~ . - DD))  # drop DD
#'
#' ## Model with different types of predictor (some multi-coefficient terms)
#' x1 <- poly(rnorm(100), 2)  # polynomial
#' x2 <- as.factor(rep(c("a", "b", "c", "d"), each = 25))  # categorical
#' x3 <- rep(1, 100)  # no variation
#' m <- lm(rnorm(100) ~ x1 + x2 + x3)
#' VIF(m)
#' @export
VIF <- function(m, data = NULL, ...) {

  d <- data

  ## Function
  VIF <- function(m) {

    ## Update model with any supplied data
    if (!is.null(d)) m <- eval(update(m, data = d, evaluate = FALSE))

    ## Term names (all)
    XN <- xNam(m, intercept = FALSE, list = TRUE)
    xn <- unlist(unname(XN))

    ## Term names (drop aliased)
    XN2 <- xNam(m, intercept = FALSE, list = TRUE, aliased = FALSE)
    xn2 <- unlist(unname(XN2))

    ## VIF's
    if (length(xn2) > 1) {

      ## T/F for terms as matrices?
      mf <- model.frame(m)
      mat <- sapply(names(XN), function(i) {
        if (i %in% names(mf)) class(mf[, i])[1] == "matrix" else FALSE
      })

      ## var-cov & cor matrix
      V <- as.matrix(vcov(m))[xn2, xn2]
      R <- cov2cor(V)
      det.R <- det(R)

      ## Function
      VIF <- function(i) {
        if (all(i %in% xn2)) {
          df <- length(i)
          ni <- !xn2 %in% i
          Ri <- R[i, i, drop = FALSE]
          Rni <- R[ni, ni, drop = FALSE]
          vif <- det(Ri) * det(Rni) / det.R
          (vif^(1 / (2 * df)))^2
        } else NA
      }

      ## Calculate
      vif <- Map(function(i, j) {
        if (j) sapply(i, VIF) else {
          vif <- VIF(i)
          sapply(i, function(x) vif)
        }
      }, XN, mat)
      unlist(unname(vif))

    } else {
      sapply(xn, function(i) if (i %in% xn2) 1 else NA)
    }

  }

  ## Apply recursively
  rMapply(VIF, m, SIMPLIFY = FALSE)

}


#' @title R-squared/Pseudo R-squared
#' @description Calculate R-squared or pseudo R-squared for a fitted model, as
#'   the squared multiple correlation between the observed and fitted values for
#'   the response variable. 'adjusted' and 'predicted' R-squared values are also
#'   calculated (see Details).
#' @param m A fitted model object of class \code{lm}, \code{glm}, or
#'   \code{merMod}, or a list or nested list of such objects.
#' @param data An optional dataset used to first re-fit the model(s).
#' @param adj,pred Logical. If \code{TRUE} (default), adjusted and/or predicted
#'   R squared are also returned.
#' @param re.form For mixed models of class \code{merMod}, the formula for
#'   random effects to condition on when generating fitted values used in the
#'   calculation of R-squared. Defaults to \code{NULL}, meaning all random
#'   effects are included. See \code{\link[lme4]{predict.merMod}} for further
#'   specification details.
#' @param ... Not currently used.
#' @details Various approaches to the calculation of a goodness-of-fit measure
#'   for GLM's analogous to R-squared in the ordinary linear model have been
#'   proposed. Generally termed 'pseudo R-squared' measures, they include
#'   variance-based, likelihood-based, and distribution-specific approaches.
#'   Here however, a straightforward definition is used, which can be applied to
#'   any model for which fitted values of the response variable are generated:
#'   R-squared is calculated as the squared (weighted) correlation between the
#'   observed and fitted values of the response (in the original units). This is
#'   simply the squared version of the correlation measure advocated by Zheng &
#'   Agresti (2000), itself an intuitive measure of goodness-of-fit describing
#'   the predictive power of a model. As the measure does not depend on any
#'   specific error distibution or model estimating procedure, it is also
#'   generally comparable across many different types of model (Kvalseth 1985).
#'   In the case of the ordinary linear model, the measure equals the more
#'   traditional R-squared based on sums of squares.
#'
#'   If argument \code{adj} is \code{TRUE} (default), the 'adjusted' R-squared
#'   value is also returned, which provides an estimate of the population - as
#'   opposed to sample - R-squared, via an analytical formula which adjusts
#'   R-squared for the 'degrees of freedom' of the model (i.e. the ratio of
#'   observations to parameters). Here, this is calculated via the 'Pratt'
#'   rather than standard 'Ezekiel/Wherry' formula, as this was shown in a
#'   previous simulation to be the most effective of a range of formulas at
#'   estimating the population R-squared, across a range of model specification
#'   scenarios (Yin & Fan 2001).
#'
#'   If \code{pred = TRUE} (default), then a 'predicted' R-squared is also
#'   returned, which is calculated via the same formula as for R-squared but
#'   using cross-validated rather than standard model predictions. These are
#'   obtained by dividing model response residuals by the complement of the
#'   observation leverages (diagonals of the hat matrix), then subtracting these
#'   inflated 'predicted' residuals from the response variable. This is
#'   essentially a short cut to obtaining out-of-sample predictions, normally
#'   arising via a leave-one-out cross validation procedure (in a GLM however
#'   they are not exactly equal to such predictions). The resulting R-squared is
#'   an estimate of the R-squared that would occur were the model to be fitted
#'   to new data, and will be lower than the original R-squared, and likely also
#'   the adjusted R-squared - highlighting the degree of noise in the original
#'   sample. This measure is a variant of an existing one (see
#'   \url{http://bit.ly/2IovBaP}), calculated by substituting the 'PRESS'
#'   statistic, i.e. the sum of squares of the predictive residuals (Allen
#'   1974), for the residual sum of squares in the classic R-squared formula.
#'
#'   For mixed models, the function will, by default, calculate all R-squared
#'   metrics using fitted values incorporating both the fixed and random effects
#'   - equivalent to the 'conditional' R-squared of Nakagawa \emph{et al.}
#'   (2017). To include only selected or no random effects, simply set the
#'   appropriate formula using the argument \code{re.form}, which is passed
#'   directly to \code{predict.merMod}. If \code{re.form = NA}, the measure is
#'   equivalent to the 'marginal' R-squared of Nakagawa \emph{et al.} (2017).
#'
#'   R-squared values produced by this function will always be bounded between
#'   zero (no fit) and one (perfect fit), meaning that any negative values
#'   arising from calculations will be rounded up to zero. Negative values
#'   typically mean that the fit is 'worse' than the null expectation of no
#'   relationship between the variables, which is difficult to interpret in
#'   practice and in any case usually only occurs in rare situations, such as
#'   where the intercept is suppressed. Hence, for simplicity and ease of
#'   interpretation, values <= 0 are presented as a complete lack of model fit.
#'
#' @note Caution must be exercised in interpreting the values of any pseudo
#'   R-squared measure calculated for a GLM or mixed model (including those
#'   produced by this function), as such measures do not hold all the properties
#'   of R-squared in the ordinary linear model and as such may not always behave
#'   as expected. They are, at best, approximations. Care must also be taken in
#'   comparing the measures to their equivalents from ordinary linear models.
#'   This is particularly the case for the adjusted and predicted versions,
#'   which have previously only been defined for ordinary linear models, and
#'   which could be described as 'approximations of approximations' of what they
#'   intend to measure. For example, for the adjusted R-squared for mixed
#'   models, its not entirely clear what the sample size (n) in the formula
#'   should represent - no. of observations? groups? something else? (the
#'   default interpretation of no. of observations is used). With all that being
#'   said, the value of standardised R-squared measures for even 'rough' model
#'   fit assessment and comparison may outweigh such reservations, and the
#'   adjusted and predicted versions in particular may aid the user in
#'   diagnosing and preventing overfitting. They should NOT, however, replace
#'   measures such as AIC or BIC for comparing and/or ranking competing models
#'   fit to the same response variable.
#' @return A numeric vector of the R-squared value(s), or an array, list or
#'   nested list of such vectors.
#' @references Allen, D. M. (1974). The Relationship Between Variable Selection
#'   and Data Agumentation and a Method for Prediction. \emph{Technometrics},
#'   \strong{16}(1), 125-127.
#'   \url{https://doi.org/10.1080/00401706.1974.10489157}
#'
#'   Kvalseth, T. O. (1985) Cautionary Note about R2. \emph{The American
#'   Statistician}, \strong{39}(4), 279-285.
#'   \url{https://doi.org/10.2307/2683704}
#'
#'   Nakagawa, S., Johnson, P.C.D. and Schielzeth, H. (2017) The coefficient of
#'   determination R2 and intra-class correlation coefficient from generalized
#'   linear mixed-effects models revisited and expanded. \emph{Journal of the
#'   Royal Society Interface} \strong{14}(134).
#'   \url{https://doi.org/10.1098/rsif.2017.0213}
#'
#'   Yin, P. and Fan, X. (2001) Estimating R2 Shrinkage in Multiple Regression:
#'   A Comparison of Different Analytical Methods. \emph{The Journal of
#'   Experimental Education} \strong{69}(2), 203-224.
#'   \url{https://doi.org/10.1080/00220970109600656}
#'
#'   Zheng, B. and Agresti, A. (2000) Summarizing the predictive power of a
#'   generalized linear model. \emph{Statistics in Medicine} \strong{19}(13),
#'   1771-1781.
#'   \url{https://doi.org/10.1002/1097-0258(20000715)19:13<1771::aid-sim485>3.0
#'   .co;2-p}
#' @examples
#' ## Pseudo R-squared for mixed models
#' R2(Shipley.SEM)  # fixed + random
#' R2(Shipley.SEM, re.form = ~ (1 | tree))  # fixed + 'tree'
#' R2(Shipley.SEM, re.form = ~ (1 | site))  # fixed + 'site'
#' R2(Shipley.SEM, re.form = NA)  # fixed only ('marginal')
#'
#' ## Predicted R-squared: compare cross-validated predictions calculated/
#' ## approximated via the hat matrix to more standard method (leave-one-out)
#'
#' \dontrun{
#'
#' ## Fit test models using Shipley data - compare lm vs glm
#' d <- na.omit(Shipley)
#' # m <- lm(Live ~ Date + DD + lat, d)
#' m <- glm(Live ~ Date + DD + lat, binomial, d)
#' ## Manual cv predictions (leave-one-out)
#' cvf1 <- sapply(1:nrow(d), function(i) {
#'   m.ni <- update(m, data = d[-i, ])
#'   predict(m.ni, d[i, ], type = "response")
#' })
#'
#' ## Short-cut via the hat matrix
#' y <- getY(m)
#' f <- fitted(m)
#' cvf2 <- y - (y - f) / (1 - hatvalues(m))
#'
#' ## Compare predictions (not exactly equal for GLM's)
#' all.equal(cvf1, cvf2)
#' # lm: TRUE; glm: "Mean relative difference: 1.977725e-06"
#' cor(cvf1, cvf2)
#' # lm: 1; glm: 0.9999987
#'
#' }
#'
#' # NOTE: comparison not tested here for mixed models, as hierarchical data can
#' # complicate the choice of an appropriate leave-one-out procedure. However,
#' # there is no reason why use of the leverage values (diagonals of the hat
#' # matrix) to calculate/estimate CV predictions shouldn't generalise
#' # (roughly?) to the mixed model case. In any case, users should exercise
#' # caution in interpretation of the predicted R-squared for mixed models,
#' # especially GLMM's.
#' @export
R2 <- function(m, data = NULL, adj = TRUE, pred = TRUE, re.form = NULL,
               ...) {

  d <- data

  ## Function
  R2 <- function(m) {

    ## Update model with any supplied data
    if (!is.null(d)) m <- eval(update(m, data = d, evaluate = FALSE))

    ## R squared
    i <- attr(terms(m), "intercept")
    k <- nrow(summary(m)$coef) - i
    if (isMerMod(m)) k <- k + length(m@theta)
    R2 <- if (k > 0) {
      y <- getY(m)
      n <- length(y)
      w <- weights(m)
      if (is.null(w)) w <- rep(1, n)
      s <- w > 0; w <- w[s]
      f <- predict(m, type = "response", re.form = re.form)[s]
      R <- cov.wt(cbind(y, f), w, cor = TRUE)$cor[1, 2]
      if (is.na(R)) R <- 0
      if (R > 0) R^2 else 0
    } else 0

    ## Adjusted R squared (Pratt formula)
    R2a <- if (adj) {
      if (R2 > 0) {
        R2a <- 1 - ((n - 3) * (1 - R2) / (n - k - i)) *
          (1 + (2 * (1 - R2) / (n - k - 2.3)))
        # R2a <- 1 - (1 - R2) * ((n - i) / (n - k - i))  # standard formula (Ezekiel/Wherry)
        if (R2a > 0) R2a else 0
      } else 0
    }

    ## Predictive R squared
    R2p <- if (pred) {
      if (R2 > 0) {
        Hii <- suppressWarnings(hatvalues(m))
        s <- Hii < 1
        f <- y - (y - f) / (1 - Hii)
        Rp <- cov.wt(cbind(y, f)[s, ], w[s], cor = TRUE)$cor[1, 2]
        if (Rp > 0) Rp^2 else 0
      } else 0
    }

    ## Return values
    c("(r.squared)" = R2, "(adj.r.squared)" = R2a, "(pred.r.squared)" = R2p)

  }

  ## Apply recursively
  rMapply(R2, m)

}


#' @title Weighted Average of Model Estimates
#' @description Calculate a weighted average of model estimates (e.g.
#'   coefficients, fitted values, residuals) for a set of models.
#' @param e A list or nested list of numeric vectors, comprising the model
#'   estimates. In the latter case, these should correspond to estimates for
#'   candidate models for each of a set of different response variables.
#' @param weights An optional numeric vector of weights to use for model
#'   averaging, or a named list of such vectors. The former should be supplied
#'   when \code{m} is a list, and the latter when it is a nested list (with
#'   matching list names). If \code{weights = "equal"} (default), a simple
#'   average is calculated instead.
#' @param est.names An optional character vector of names used to extract and/or
#'   sort estimates from the output.
#' @param ... Not currently used.
#' @details This function can be used to calculate a weighted average of model
#'   estimates such as coefficients, fitted values, or residuals, where models
#'   are typically competing candidate models fit to the same response variable
#'   and data. Weights are typically a 'weight of evidence' type metric such as
#'   Akaike model weights (Burnham & Anderson 2002, Burnham \emph{et al.} 2011),
#'   which can be calculated in \code{R} using packages such as \pkg{MuMIn} or
#'   \pkg{AICcmodavg}. However, weights of any sort can be used, provided they
#'   are numeric. If none are supplied, the simple average is calculated
#'   instead.
#'
#' @note Averaging is performed via the 'full/zero' rather than
#'   'subset/conditional/natural' method, meaning that zero is substituted for
#'   any missing estimates (coefficients) prior to calculations - providing a
#'   form of 'shrinkage' (Burnham & Anderson 2002, Grueber \emph{et al.} 2011).
#' @return A numeric vector of the model-averaged estimates, or a list of such
#'   vectors.
#' @references Burnham, K. P., & Anderson, D. R. (2002). \emph{Model Selection
#'   and Multimodel Inference: A Practical Information-Theoretic Approach} (2nd
#'   ed.). New York: Springer-Verlag. Retrieved from
#'   \url{https://www.springer.com/gb/book/9780387953649}
#'
#'   Burnham, K. P., Anderson, D. R., & Huyvaert, K. P. (2011). AIC model
#'   selection and multimodel inference in behavioral ecology: some background,
#'   observations, and comparisons. \emph{Behavioral Ecology and Sociobiology},
#'   \strong{65}(1), 23-35. \url{https://doi.org/10.1007/s00265-010-1029-6}
#'
#'   Dormann, C. F., Calabrese, J. M., Guillera-Arroita, G., Matechou, E., Bahn,
#'   V., Barton, K., ... Hartig, F. (2018). Model averaging in ecology: a review
#'   of Bayesian, information-theoretic, and tactical approaches for predictive
#'   inference. \emph{Ecological Monographs}, \strong{88}(4), 485-504.
#'   \url{https://doi.org/10.1002/ecm.1309}
#'
#'   Grueber, C. E., Nakagawa, S., Laws, R. J., & Jamieson, I. G. (2011).
#'   Multimodel inference in ecology and evolution: challenges and solutions.
#'   \emph{Journal of Evolutionary Biology}, \strong{24}(4), 699-711.
#'   \url{https://doi.org/10.1111/j.1420-9101.2010.02210.x}
#'
#'   Walker, J. A. (2019). Model-averaged regression coefficients have a
#'   straightforward interpretation using causal conditioning. \emph{BioRxiv},
#'   133785. \url{https://doi.org/10.1101/133785}
#' @seealso \code{\link[stats]{weighted.mean}}
#' @examples
#' ## Model-averaged coefficients
#' m <- Shipley.Growth  # candidate models
#' b <- lapply(m, function(i) coef(summary(i))[, 1])
#' avgEst(b)
#'
#' ## Using weights
#' w <- runif(length(b), 0, 1)
#' avgEst(b, w)
#'
#' ## Model-averaged predictions
#' f <- lapply(m, predict)
#' avgEst(f, w)
#' @export
avgEst <-  function(e, weights = "equal", est.names = NULL, ...) {

  w <- weights; en <- est.names

  ## Weights
  if (all(w == "equal")) {
    f <- function(i) rep(1, length(i))
    w <- if (any(sapply(e, isList))) lapply(e, f) else f(e)
  }

  ## Function
  avgEst <- function(e, w) {

    ## Sort names (for coefs)
    en2 <- unique(unlist(lapply(e, names)))
    num <- suppressWarnings(as.numeric(en2))
    if (all(is.na(num))) {
      i <- isInt(en2)
      r2 <- isR2(en2)
      en3 <- sort(en2[!i & !r2])
      en3 <- names(sort(sapply(en3, function(i) {
        lengths(regmatches(i, gregexpr(":", i)))
      })))
      en2 <- c(en2[i], en3, en2[r2])
    }
    en <- if (!is.null(en)) en[en %in% en2] else en2

    ## Combine estimates into table (missing -> zero)
    e <- do.call(cbind, lapply(e, function(i) {
      sapply(en, function(j) {
        if (j %in% names(i)) i[[j]] else 0
      })
    }))

    ## Weighted average
    apply(e, 1, weighted.mean, w)

  }

  ## Apply recursively
  rMapply(avgEst, e, w, SIMPLIFY = FALSE)

}


#' @title Standardised Coefficients
#' @description Calculate fully standardised model coefficients in standard
#'   deviation units, and adjusted for multicollinearity among predictors.
#' @param m A fitted model object of class \code{lm}, \code{glm}, or
#'   \code{merMod}, or a list or nested list of such objects.
#' @param weights An optional numeric vector of weights to use for model
#'   averaging, or a named list of such vectors. The former should be supplied
#'   when \code{m} is a list, and the latter when it is a nested list (with
#'   matching list names). If set to \code{"equal"}, a simple average is
#'   calculated instead.
#' @param data An optional dataset used to first re-fit the model(s).
#' @param term.names An optional vector of term names used to extract and/or
#'   sort coefficients from the output.
#' @param cen.x,cen.y Logical, whether the intercept and coefficients should be
#'   calculated from mean-centred variables.
#' @param std.x,std.y Logical, whether coefficients should be scaled by the
#'   standard deviations of variables.
#' @param unique.x Logical, whether coefficients should be adjusted for
#'   multicollinearity among predictors.
#' @param r.squared Logical, whether R-squared values should also be returned.
#' @param ... Arguments to \code{R2}.
#' @details This function will calculate fully standardised coefficients in
#'   standard deviation units for linear, generalised linear, and mixed models
#'   (of class \code{merMod}). It achieves this via adjusting the 'raw' model
#'   coefficients, so no standardisation of input variabes is required
#'   beforehand. Users can simply specify the model with all variables in their
#'   original units and the function will do the rest. However, the user is free
#'   to scale and/or centre any input variables should they choose, which will
#'   not affect the outcome of standardisation (provided any scaling is by
#'   standard deviations). This may be desirable in some cases, such as to
#'   increase numerical stability during model fitting when variables are on
#'   widely different scales.
#'
#'   If arguments \code{cen.x} or \code{cen.y} are \code{TRUE}, model estimates
#'   will be calculated as if all predictors (x) and/or the response variable
#'   (y) were mean-centred prior to model-fitting. Thus, for an ordinary linear
#'   model where centring of x and y is specified, the intercept will be zero -
#'   the mean (or weighted mean) of y. In addition, if \code{cen.x = TRUE} and
#'   there are interacting terms in the model, all coefficients for lower order
#'   terms of the interation are adjusted using an expression which ensures that
#'   each main effect or lower order term is estimated at the mean values of the
#'   terms they interact with (zero in a 'centred' model) - typically improving
#'   the interpretation of coefficients. The expression used comprises a
#'   weighted sum of all the coefficients that contain the lower order term,
#'   with the weight for the term itself being zero and those for 'containing'
#'   terms being the product of the means of the other variables involved in
#'   that term (i.e. those not in the lower order term itself). For example, for
#'   a three-way interaction (x1 * x2 * x3), the expression for main effect
#'   \eqn{\beta1} would be:
#'
#'   \deqn{\beta_{1} + \beta_{12} \bar{x}_{2} + \beta_{13} \bar{x}_{3} +
#'   \beta_{123} \bar{x}_{2} \bar{x}_{3}}{\beta1 + (\beta12 * Mx2) + (\beta13 *
#'   Mx3) + (\beta123 * Mx2 * Mx3)} (adapted from here:
#'   \url{http://bit.ly/2FOJPk8})
#'
#'   In addition, if \code{std.x = TRUE} or \code{unique.x = TRUE} (see below),
#'   product terms for interactive effects will be recalculated using
#'   mean-centred variables, to ensure that standard deviations and variance
#'   inflation factors for predictors are calculated correctly (the model must
#'   be re-fit for this latter purpose, to recalculate the variance-covariance
#'   matrix).
#'
#'   If \code{std.x = TRUE}, coefficients are standardised by multiplying by the
#'   standard deviations of predictor variables (or terms), while if \code{std.y
#'   = TRUE} they are divided by the standard deviation of the response. If the
#'   model is a GLM, this latter is calculated from the link-transformed
#'   response (or its estimate) generated using the function \code{getY}. If
#'   both arguments are true, the coefficients are regarded as 'fully'
#'   standardised in the traditional sense, often referred to as 'betas'.
#'
#'   If \code{unique.x = TRUE}, coefficients are adjusted for multicollinearity
#'   among predictors by dividing by the square root of the variance inflation
#'   factors (Dudgeon 2016, Thompson \emph{et al.} 2017). If they have also been
#'   standardised by the standard deviations of x and y, this converts them to
#'   semipartial correlations, i.e. the correlation between the unique
#'   components of predictors (residualised on other predictors) and the
#'   response variable. This measure of effect size is arguably much more
#'   interpretable and useful than the traditional standardised coefficient, as
#'   it is always estimated independent of the confounding effects of other
#'   predictors, and so can more readily be compared both within and across
#'   models. Values range from zero (no effect) to +/-1 (perfect relationship),
#'   rather than from zero to +/- infinity (as in the case of betas) - putting
#'   them on the same scale as the bivariate correlation between predictor and
#'   response. In the case of GLM's however, the measure is not exactly equal to
#'   the semipartial correlation, so its value may not be always be bound
#'   between +/-1 (such cases are likely rare). Crucially, for ordinary linear
#'   models, the square of the semipartial correlation equals the increase in
#'   R-squared when that variable is added last in the model - directly linking
#'   the measure to model fit and 'variance explained'. See
#'   \url{http://bit.ly/2GmyrMA} for additional arguments in favour of
#'   semipartial correlations.
#'
#'   If \code{r.squared = TRUE}, R-squared values are also returned via the
#'   \code{R2} function.
#'
#'   Finally, if \code{weights} are specified, the function calculates a
#'   weighted average of the standardised coefficients across models (Burnham &
#'   Anderson 2002).
#' @return A numeric vector of the standardised coefficients, or a list or
#'   nested list of such vectors.
#' @references Burnham, K. P., & Anderson, D. R. (2002). \emph{Model Selection
#'   and Multimodel Inference: A Practical Information-Theoretic Approach} (2nd
#'   ed.). New York: Springer-Verlag. Retrieved from
#'   \url{https://www.springer.com/gb/book/9780387953649}
#'
#'   Dudgeon, P. (2016). A Comparative Investigation of Confidence Intervals for
#'   Independent Variables in Linear Regression. \emph{Multivariate Behavioral
#'   Research}, \strong{51}(2-3), 139-153.
#'   \url{https://doi.org/10.1080/00273171.2015.1121372}
#'
#'   Thompson, C. G., Kim, R. S., Aloe, A. M., & Becker, B. J. (2017).
#'   Extracting the Variance Inflation Factor and Other Multicollinearity
#'   Diagnostics from Typical Regression Results. \emph{Basic and Applied Social
#'   Psychology}, \strong{39}(2), 81-90.
#'   \url{https://doi.org/10.1080/01973533.2016.1277529}
#' @seealso \code{\link[stats]{coef}}, \code{\link[semEff]{VIF}},
#'   \code{\link[semEff]{getY}}, \code{\link[semEff]{R2}},
#'   \code{\link[semEff]{avgEst}}
#' @examples
#' library(lme4)
#'
#' ## Standardised coefficients for SEM (i.e. direct effects)
#' m <- Shipley.SEM
#' stdCoeff(m)
#' stdCoeff(m, std.y = FALSE)  # x-only
#' stdCoeff(m, std.x = FALSE, std.y = FALSE)  # centred only
#' stdCoeff(m, cen.x = FALSE, cen.y = FALSE)  # standardised only
#' stdCoeff(m, r.squared = TRUE)  # add R-squared
#'
#' ## Demonstrate equality with manually-standardised variables (gaussian)
#' m <- Shipley.Growth[[3]]
#' d <- data.frame(scale(na.omit(Shipley)))
#' b1 <- stdCoeff(m, unique.x = FALSE)
#' b2 <- coef(summary(update(m, data = d)))[, 1]
#' stopifnot(all.equal(b1, b2))
#'
#' ## Demonstrate equality with increment in R-squared (ordinary linear model)
#' m <- lm(Growth ~ Date + DD + lat, data = Shipley)
#' r2 <- summary(m)$r.squared
#' b1 <- stdCoeff(m)[-1]
#' bn <- names(b1)
#' b2 <- sqrt(sapply(bn, function(i) {
#'   f <- reformulate(bn[!bn %in% i])
#'   r2i <- summary(update(m, f))$r.squared
#'   r2 - r2i
#' }))
#' stopifnot(all.equal(b1, b2))
#'
#' ## Model-averaged standardised coefficients
#' m <- Shipley.Growth  # candidate models
#' w <- runif(length(m), 0, 1)  # weights
#' stdCoeff(m, w)
#' @export
stdCoeff <- function(m, weights = NULL, data = NULL, term.names = NULL,
                     cen.x = TRUE, cen.y = TRUE, std.x = TRUE, std.y = TRUE,
                     unique.x = TRUE, r.squared = FALSE, ...) {

  w <- weights; d <- data; bn <- term.names

  ## Function
  stdCoeff <- function(m) {

    ## Update model with any supplied data
    if (!is.null(d)) m <- eval(update(m, data = d, evaluate = FALSE))

    ## Coefficients
    b <- summary(m)$coef
    b <- setNames(b[, 1], rownames(b))
    xn <- names(b)

    ## Intercept?
    i <- isInt(xn)
    int <- if (any(i)) {xn <- xn[!i]; TRUE} else FALSE

    ## Model weights
    n <- nobs(m)
    w <- weights(m)
    if (is.null(w)) w <- rep(1, n)
    s <- w > 0; w <- w[s]

    ## Centre/standardise x
    k <- length(xn)
    if (k > 0) {

      ## Predictors
      x <- model.matrix(m)[s, , drop = FALSE]
      x <- data.frame(x, check.names = FALSE)[xn]
      obs <- rownames(x)

      ## Interactions?
      inx <- any(isInx(xn))

      ## Centre predictors and adjust coefs/intercept
      if (cen.x) {

        ## For interactions, adjust coefs and centre predictors
        if (inx) {

          ## List of main effect names for all terms
          XN <- sapply(xn, function(i) unlist(strsplit(i, ":")))

          ## Predictor means
          xm <- colMeans(x)

          ## Adjust lower-order coefs
          b[xn] <- sapply(xn, function(i) {
            bi <- b[[i]]
            XNi <- XN[[i]]
            ti <- xn[sapply(xn, function(j) all(XNi %in% XN[[j]]))]
            ti <- ti[ti != i]
            if (length(ti) > 0) {
              nim <- sapply(ti, function(j) {
                ni <- XN[[j]][!XN[[j]] %in% XNi]
                prod(xm[ni])
              })
              bi + sum(b[ti] * nim)
            } else bi
          })

          ## Centre predictors (for correct SD's/VIF's)
          if (std.x || unique.x) {
            x <- sapply(XN, function(i) {
              xi <- sweep(x[i], 2, xm[i])
              apply(xi, 1, prod)
            })
            x <- data.frame(x, row.names = obs, check.names = FALSE)
          }

        }

        ## Adjust intercept (weighted mean of predicted y)
        if (int) b[1] <- weighted.mean(predict(m, re.form = NA)[s], w)

      }

      ## Standardise by x
      if (std.x) b[xn] <- b[xn] * sapply(x, sdW, w)

      ## Calculate unique effects of predictors (adjust for multicollinearity)
      if (unique.x && k > 1) {

        ## Re-fit model with centred predictors
        ## (to calculate correct VIF's for interactions)
        m2 <- if (cen.x && inx) {

          ## Add centred predictors to data (list)
          if (is.null(d)) d <- getData(m)
          d <- c(list(x = as.matrix(x)), as.list(d[obs, ]))

          ## Add any offset(s)
          o <- model.offset(model.frame(m)[s, ])
          if (is.null(o)) o <- rep(0, n)
          d <- c(list(o = o), d)

          ## New model formula
          ran <- if (isMerMod(m)) {
            bars <- lme4::findbars(formula(m))
            sapply(bars, function(i) {
              paste0("(", deparse(i), ")")
            })
          }  # ran. effects (http://bit.ly/2V1yDeu)
          f <- reformulate(c("x", ran), response = ".")

          ## Update model
          update(m, f, data = d, offset = o)

        } else m

        ## Divide coefs by square root of VIF's
        vif <- na.omit(VIF(m2))
        b[xn] <- b[xn] / sqrt(vif)

      }

    }

    ## Centre/standardise y
    if (cen.y && int) {
      ym <- weighted.mean(getY(m), w)
      b[1] <- b[1] - family(m)$linkfun(ym)
    }
    if (std.y) b <- b / sdW(getY(m, link = TRUE), w)

    ## Return standardised coefficients
    sapply(xNam(m), function(i) unname(b[i]))

  }

  ## Add R-squared?
  stdCoeff2 <- if (r.squared) {
    function(m) c(stdCoeff(m), R2(m, d, ...))
  } else stdCoeff

  ## Apply recursively
  b <- rMapply(stdCoeff2, m, SIMPLIFY = FALSE)

  ## Output coefs or weighted average
  if (!is.null(w) && isList(b)) avgEst(b, w, bn)
  else {
    if (!is.null(bn)) {
      f <- function(i) i[bn[bn %in% names(i)]]
      rMapply(f, b, SIMPLIFY = FALSE)
    } else b
  }

}


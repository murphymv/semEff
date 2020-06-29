

#' @title Weighted Variance
#' @description Calculate the weighted variance of \code{x}.
#' @param x A numeric vector.
#' @param w A numeric vector of weights of the same length as \code{x}.
#' @param na.rm Logical, whether NAs in \code{x} should be removed.
#' @param ... Not currently used.
#' @details Calculate the weighted variance of \code{x} via the weighted
#'   covariance matrix (\code{cov.wt}). If no weights are supplied, the simple
#'   variance is returned instead. As in \code{weighted.mean}, \code{NA}s in
#'   \code{w} are not handled specially and will return \code{NA} as result.
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
#' varW(x)
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
#' @details This is simply a wrapper for \code{varW}, applying the square root
#'   to the output.
#' @return A numeric value, the weighted standard deviation of \code{x}.
#' @seealso \code{\link[stats]{sd}}, \code{\link[semEff]{varW}}
#' @export
sdW <- function(...) {
  sqrt(varW(...))
}


#' @title Get Model Data
#' @description Extract the data used to fit a model.
#' @param mod A fitted model object, or a list or nested list of such objects.
#' @param subset Logical. If \code{TRUE}, only observations used to fit the
#'   model(s) are returned (i.e. missing observations (\code{NA}) are removed).
#' @param merge Logical. If \code{TRUE}, and \code{mod} is a list or nested
#'   list, a single dataset containing all variables used to fit models is
#'   returned.
#' @param ... Arguments to \code{eval}.
#' @details This is a simple convenience function to return the data used to fit
#'   a model, by evaluating the 'data' slot of the model call object. If the
#'   'data' argument of the model call was not specified, or is not a data frame
#'   (or coercible to such) containing all variables referenced in the model
#'   formula, an error will be thrown - this restriction is largely to ensure
#'   that a single coherent dataset of all model variables can be made available
#'   for resampling purposes.
#'
#'   If \code{mod} is a list of models and \code{merge = TRUE}, all (unique)
#'   variables used to fit models are merged into a single data frame. This will
#'   return an error if \code{subset = TRUE} results in datasets with different
#'   numbers of observations (rows).
#' @return A data frame of the variables used to fit the model(s), or a list or
#'   nested list of same.
#' @seealso \code{\link[stats]{getCall}}, \code{\link[base]{eval}}
#' @examples
#' ## Get data used to fit SEM from Shipley (2009)
#' getData(Shipley.SEM[[1]])  # from single model
#' getData(Shipley.SEM)  # from SEM (list)
#' getData(Shipley.SEM, merge = TRUE)  # from SEM (single dataset)
#' @export
getData <- function(mod, subset = FALSE, merge = FALSE, ...) {

  m <- mod

  ## Function
  getData <- function(m) {

    ## Data from 'data' argument of model call
    mc <- getCall(m)
    d <- eval(mc$data, ...)
    if (!is.null(d)) d <- data.frame(d) else
      stop("'data' argument of model call not specified.")

    ## All var names from model call
    f <- c(formula(m), mc$subset, mc$weights, mc$offset, mc$correlation)
    vn <- unlist(lapply(f, all.vars))
    if (!all(vn %in% names(d)))
      stop("'data' does not contain all model variables.")

    ## Subset data for model observations?
    if (subset) {
      obs <- names(fitted(m)); w <- weights(m)
      if (!is.null(w)) obs <- obs[w > 0 & !is.na(w)]
      d[obs, ]
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


#' @title Get Model Term Names
#' @description Extract term names from a fitted model object.
#' @param mod A fitted model object, or a list or nested list of such objects.
#' @param data An optional dataset used to construct the model frame.
#' @param intercept Logical, whether the intercept should be included.
#' @param aliased Logical, whether names of aliased terms should be included
#'   (see Details).
#' @param list Logical, whether names should be returned as a list, with all
#'   multi-coefficient terms grouped under their main term names.
#' @param ... Arguments to \code{eval} (for evaluating model data).
#' @details Extract term names from a fitted model. Names of terms for which
#'   coefficients cannot be estimated are also included if \code{aliased = TRUE}
#'   (default). These may be terms which are perfectly correlated with other
#'   terms in the model, so that the model design matrix is rank deficient.
#' @return A character vector or list/nested list of term names.
#' @examples
#' ## Term names from Shipley SEM
#' m <- Shipley.SEM
#' xNam(m)
#' xNam(m, intercept = FALSE)
#'
#' ## Model with different types of predictor (some multi-coefficient terms)
#' d <- data.frame(
#'   y = rnorm(100),
#'   x1 = rnorm(100),
#'   x2 = as.factor(rep(c("a", "b", "c", "d"), each = 25)),
#'   x3 = rep(1, 100)
#' )
#' m <- lm(y ~ poly(x1, 2) + x2 + x3, data = d)
#' xNam(m)
#' xNam(m, aliased = FALSE)  # drop term that cannot be estimated (x3)
#' xNam(m, aliased = FALSE, list = TRUE)  # names as list
#' @export
xNam <- function(mod, data = NULL, intercept = TRUE, aliased = TRUE,
                 list = FALSE, ...) {

  m <- mod; d <- data

  ## Function
  xNam <- function(m) {

    ## Term names
    tt <- terms(m)
    xn <- labels(tt)
    int <- attr(tt, "intercept")
    if (int) xn <- c("(Intercept)", xn)

    ## Expand interaction terms (list)
    XN <- sapply(xn, function(i) {
      unlist(strsplit(i, "(?<!:):(?!:)", perl = TRUE))
    })

    ## Predictor values
    if (is.null(d)) d <- getData(m, ...)
    mf <- suppressWarnings(model.frame(m, data = d))
    x <- mf[names(mf) %in% unlist(unname(XN))]
    x <- sapply(x, "[", simplify = FALSE)

    ## Expand factor/matrix terms (list)
    f <- sapply(x, is.factor); f1 <- names(f[f][1])
    xn2 <- unique(c(xn, names(x)))
    XN2 <- sapply(xn2, function(i) {
      if (i %in% names(x)) {
        xi <- x[[i]]
        if (f[i]) {
          ci <- contrasts(d[[i]])
          ct <- getCall(m)$contrasts
          xi <- if (!is.null(ct)) {
            ci <- ct[names(ct) %in% i][[1]]; chr <- is.character(ci)
            ci <- eval(if (chr) parse(text = ci) else ci)
            if (is.function(ci)) {
              l <- levels(xi); ci(if (!chr) length(l) else l)
            } else ci
          } else ci
        }
        j <- colnames(xi); n <- ncol(xi)
        if (is.null(j) && (f[i] || isTRUE(n > 1))) j <- 1:n
        paste0(i, j)
      } else i
    }, simplify = FALSE)

    ## Expand interaction terms involving factors/matrices (list)
    XN <- sapply(xn, function(i) {
      if (isInx(i)) {
        j <- expand.grid(XN2[XN[[i]]])
        apply(j, 1, paste, collapse = ":")
      } else XN2[[i]]
    }, simplify = FALSE)

    ## If no intercept, generate all levels for first factor
    if (!int && any(f)) XN[[f1]] <- paste0(f1, levels(x[[f1]]))

    ## Drop some names? (aliased, subsetted factor levels, etc.)
    b <- if (isMer(m)) lme4::fixef(m, add.dropped = TRUE) else coef(m)
    bn <- names(b); if (!aliased) bn <- bn[!is.na(b)]
    XN <- lapply(XN, function(i) i[i %in% bn])
    XN <- XN[sapply(XN, length) > 0]

    ## Drop intercept?
    if (!intercept && int) XN <- XN[-1]

    ## Return as list?
    if (!list) unlist(unname(XN)) else XN

  }

  ## Apply recursively
  rMapply(xNam, m, SIMPLIFY = FALSE)

}


#' @title Generalised Link Transformation
#' @description Transform a numeric variable using a GLM link function, or
#'   return an estimate of same.
#' @param x a positive numeric vector, corresponding to a variable to be
#'   transformed. Can also be a list or nested list of such vectors.
#' @param family Optional, the error distribution family containing the link
#'   function which will be used to transform \code{x} (see
#'   \code{\link[stats]{family}} for specification details). If not supplied, it
#'   is determined from \code{x} (see Details).
#' @param force.est Logical, whether to force the return of the estimated rather
#'   than direct transformation, where the latter is available (i.e. does not
#'   contain undefined values).
#' @param ... Not currently used.
#' @details \code{glt} can be used to provide a 'generalised' transformation of
#'   a numeric variable using the link function from a generalised linear model
#'   (GLM) fit to the variable. The transformation is generalised in the sense
#'   that it can always be generated, even where a standard link transformation
#'   would produce undefined values. It achieves this via an estimate based on
#'   the 'working' response variable of the GLM (see below). If the error
#'   distribution \code{family} is not specified (default), then it is
#'   determined (roughly) from \code{x}, with \code{binomial(link = "logit")}
#'   used when all x <= 1 and \code{poisson(link = "log")} otherwise. Although
#'   the function is generally intended for binomial or poisson variables, any
#'   variable which can be fit using \code{glm} can be supplied. One of the key
#'   purposes of \code{glt} is to allow the calculation of fully standardised
#'   model coefficients for GLMs (in which case \code{x} = the response
#'   variable), while it can also facilitate the proper calculation of SEM
#'   indirect effects (see below).
#'
#'   \strong{Estimating the link transformation}
#'
#'   A key challenge in generating fully standardised model coefficients for a
#'   GLM with a non-gaussian link function is the difficulty in calculating
#'   appropriate standardised ranges (typically the standard deviation) for the
#'   response variable in the link scale. This is because directly transforming
#'   the response will often produce undefined values. Although methods for
#'   circumventing this issue by indirectly estimating the variance of the
#'   response on the link scale have been proposed - including a
#'   latent-theoretic approach for binomial models (McKelvey & Zavoina 1975) and
#'   a more general variance-based method using a pseudo R-squared (Menard 2011)
#'   - here an alternative approach is used. Where transformed values are
#'   undefined, the function will instead return the synthetic 'working'
#'   response from the iteratively reweighted least squares algorithm (IRLS) of
#'   the GLM (McCullagh & Nelder 1989). This is reconstructed as the sum of the
#'   linear predictor and the working residuals - with the latter comprising the
#'   error of the model on the link scale. The advantage of this approach is
#'   that a relatively straightforward 'transformation' of any non-gaussian
#'   response is readily attainable in all cases. The standard deviation (or
#'   other relevant range) can then be calculated using values of the
#'   transformed response and used to scale the coefficients. An additional
#'   benefit for piecewise SEMs is that the transformed rather than original
#'   response can be specified as a predictor in other models, ensuring that
#'   standardised indirect and total effects are calculated correctly (i.e.
#'   using the same units).
#'
#'   To ensure a high level of 'accuracy' in the working response - in the sense
#'   that the inverse-transformation is practically indistinguishable from the
#'   original response variable - the function uses the following iterative
#'   fitting procedure to calculate a 'final' working response:
#'
#'   \enumerate{\item A new GLM of the same error family is fit with the
#'   original response variable as both predictor and response, and using a
#'   single IWLS iteration. \item The working response is extracted from this
#'   model. \item The inverse transformation of the working response is then
#'   calculated. \item If the inverse transformation is 'effectively equal' to
#'   the original response (tested using \code{all.equal}), the working response
#'   is returned; otherwise, the GLM is re-fit with the working response now as
#'   the predictor, and steps 2-4 are repeated - each time with an additional
#'   IWLS iteration.}
#'
#'   This approach will generate a very reasonable transformation of the
#'   response variable, which will also be practically indistinguishable from
#'   the direct transformation where this can be compared (see Examples). It
#'   also ensures that the transformed values, and hence the standard deviation,
#'   are the same for any GLM fitting the same response - provided it uses the
#'   same link function - facilitating model comparisons, selection, and
#'   averaging.
#'
#' @note As we often cannot directly observe the GLM response variable on the
#'   link scale, any method estimating its values or statistics will be 'wrong'
#'   to a greater or lesser degree. The aim should be to try to minimise this
#'   error as far as (reasonably) possible, while also generating standardised
#'   coefficients whose interpretation most closely resembles those of the
#'   ordinary linear model - something which the current method achieves. The
#'   solution of using the working response from the GLM to scale coefficients
#'   is a purely practical, but reasonable one, and one that takes advantage of
#'   modern computing power to minimise error through iterative model fitting.
#'   An added bonus is that the estimated variance is constant across models fit
#'   to the same response variable, which cannot be said of previous methods
#'   (Menard 2011). The overall approach would be classed as
#'   'observed-empirical' by Grace \emph{et al.} (2018), as it utilises model
#'   error variance (the working residuals) rather than theoretical
#'   distribution-specific variance.
#' @return A numeric vector of the transformed values, or an array, list of
#'   vectors/arrays, or nested list.
#' @references Grace, J.B., Johnson, D.J., Lefcheck, J.S. and Byrnes, J.E.K.
#'   (2018) Quantifying relative importance: computing standardized effects in
#'   models with binary outcomes. \emph{Ecosphere} \strong{9}, e02283.
#'   \url{https://doi.org/gdm5bj}
#'
#'   McCullagh P. and Nelder, J. A. (1989) \emph{Generalized Linear Models} (2nd
#'   Edition). London: Chapman and Hall.
#'
#'   McKelvey, R. D., & Zavoina, W. (1975). A statistical model for the analysis
#'   of ordinal level dependent variables. \emph{The Journal of Mathematical
#'   Sociology}, \strong{4}(1), 103-120. \url{https://doi.org/dqfhpp}
#'
#'   Menard, S. (2011) Standards for Standardized Logistic Regression
#'   Coefficients. \emph{Social Forces} \strong{89}, 1409-1428.
#'   \url{https://doi.org/bvxb6s}
#' @seealso \code{\link[stats]{glm}}, \code{\link[base]{all.equal}}
#' @examples
#' ## Compare estimate with a direct link transformation
#' ## (test with a poisson variable, log link)
#' set.seed(1)
#' y <- rpois(30, lambda = 10)
#' yl <- glt(y, force.est = TRUE)
#'
#' ## Effectively equal?
#' all.equal(log(y), yl, check.names = FALSE)
#' # TRUE
#'
#' ## Actual difference...
#' all.equal(log(y), yl, check.names = FALSE, tolerance = .Machine$double.eps)
#' # "Mean relative difference: 1.05954e-12"
#' @export
glt <- function(x, family = NULL, force.est = FALSE, ...) {

  f <- family

  ## Function
  glt <- function(x) {

    ## Error family
    if (is.character(f)) {
      f <- get(f, mode = "function", envir = parent.frame())
    }
    if (is.function(f)) f <- f()
    if (is.null(f)) {
      f <- if (all(x <= 1)) binomial() else poisson()
    }

    ## Transform variable to link scale
    xl <- f$linkfun(x)

    ## Return transformation or estimation (GLM working response)
    if (any(is.infinite(xl)) || force.est) {
      x2 <- x
      suppressWarnings(
        m <- do.call(glm, list(x ~ x2, family = f, control = list(maxit = 1),
                               na.action = na.exclude))
      )
      i <- 0
      repeat {
        xl <- predict(m) + resid(m, "working")
        xli <- f$linkinv(xl)
        eql <- isTRUE(all.equal(x, xli, check.names = FALSE))
        if (eql) return(xl) else {
          i <- i + 1
          suppressWarnings(
            m <- update(m, . ~ xl, control = list(maxit = i))
          )
        }
      }
    } else xl

  }

  ## Apply recursively
  rMapply(glt, x)

}


#' @title Get Model Response Variable
#' @description Extract the response variable from a fitted model on the
#'   original or link scale.
#' @param mod A fitted model object, or a list or nested list of such objects.
#' @param data An optional dataset used to first re-fit the model(s).
#' @param link Logical. If \code{TRUE}, return the GLM response variable on the
#'   link scale (see Details).
#' @param offset Logical. If \code{TRUE}, include model offset(s) in the
#'   response.
#' @param ... Arguments to \code{glt} (not including \code{family}, which is
#'   determined from \code{mod}).
#' @details \code{getY} will return the response variable from a model by
#'   summing the fitted values and the response residuals. If \code{link = TRUE}
#'   and the model is a GLM, the response is returned on the link scale. If this
#'   results in undefined values, it is replaced by an estimate based on the
#'   'working' response variable of the GLM (see \code{\link[semEff]{glt}}).
#'
#'   Any offset variables are removed from the response by default. This means
#'   that, for example, rates rather than raw counts will be returned for
#'   poisson GLMs (where applicable).
#' @return A numeric vector comprising the response variable on the original or
#'   link scale, or an array, list of vectors/arrays, or nested list.
#' @examples
#' ## All SEM responses (original scale)
#' getY(Shipley.SEM)
#'
#' ## Estimated response in link scale from binomial model
#' getY(Shipley.SEM$Live, link = TRUE)
#' @export
getY <- function(mod, data = NULL, link = FALSE, offset = FALSE, ...) {

  m <- mod; d <- data

  ## Function
  getY <- function(m) {

    ## Update model with any supplied data
    if (!is.null(d)) m <- eval(update(m, data = d, evaluate = FALSE))

    ## Model error family
    f <- if (isBet(m)) m$link$mean else family(m)

    ## Model weights and offset
    w <- weights(m)
    o <- if (!offset) {
      if (is.null(d)) d <- getData(m, subset = TRUE)
      mf <- suppressWarnings(model.frame(m, data = d))
      model.offset(mf)
    }

    ## Model response
    y <- fitted(m) + resid(m, "response")
    if (!is.null(w)) y <- y[w > 0 & !is.na(w)]
    if (!is.null(o)) y <- f$linkinv(f$linkfun(y) - o)
    a <- names(attributes(y))
    attributes(y)[a != "names"] <- NULL

    ## Return response on original or link scale
    if (isGlm(m) && link) glt(y, family = f, ...) else y

  }

  ## Apply recursively
  rMapply(getY, m)

}


#' @title Generalised Variance Inflation Factors
#' @description Calculate generalised variance inflation factors for terms in a
#'   fitted model via the variance-covariance matrix of coefficients.
#' @param mod A fitted model object, or a list or nested list of such objects.
#' @param data An optional dataset used to first re-fit the model(s).
#' @param ... Arguments to \code{eval} (for evaluating model data).
#' @details \code{VIF} calculates generalised variance inflation factors (GVIF)
#'   as described in Fox & Monette (1992), and also implemented in the
#'   \code{vif} function in the \pkg{car} package. However, whereas \code{vif}
#'   returns both GVIF and GVIF^(1/(2*Df)) values, \code{VIF} simply returns the
#'   squared result of the latter measure, which equals the standard VIF for
#'   single-coefficient terms and is the equivalent measure for
#'   multi-coefficient terms (e.g. categorical or polynomial). Also, while
#'   \code{vif} returns values per model term (i.e. predictor variable),
#'   \code{VIF} returns values per coefficient, meaning that the same VIF will
#'   be returned per coefficient for multi-coefficient terms. Finally, \code{NA}
#'   is returned for any coefficients which could not be estimated in the model
#'   (e.g. aliased terms).
#' @return A numeric vector of the VIFs, or an array, list of vectors/arrays,
#'   or nested list.
#' @references Fox, J. and Monette, G. (1992) Generalized Collinearity
#'   Diagnostics. \emph{Journal of the American Statistical Association}
#'   \strong{87}, 178-183. \url{https://doi.org/dm9wbw}
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
#' d <- data.frame(
#'   y = rnorm(100),
#'   x1 = rnorm(100),
#'   x2 = as.factor(rep(c("a", "b", "c", "d"), each = 25)),
#'   x3 = rep(1, 100)
#' )
#' m <- lm(y ~ poly(x1, 2) + x2 + x3, data = d)
#' VIF(m)
#' @export
VIF <- function(mod, data = NULL, ...) {

  m <- mod; d <- data

  ## Function
  VIF <- function(m) {

    ## Update model with any supplied data
    if (!is.null(d)) m <- eval(update(m, data = d, evaluate = FALSE))

    ## Term names
    XN <- xNam(m, intercept = FALSE, list = TRUE, ...)
    xn <- xNam(m, intercept = FALSE, aliased = FALSE, ...)

    ## VIFs
    if (length(xn) > 1) {

      ## T/F for terms as matrices
      if (is.null(d)) d <- getData(m, ...)
      mf <- suppressWarnings(model.frame(m, data = d))
      mat <- sapply(names(XN), function(i) {
        if (i %in% names(mf)) class(mf[, i])[1] == "matrix" else FALSE
      })

      ## Var-cov/cor matrix
      V <- as.matrix(vcov(m))[xn, xn]
      R <- cov2cor(V)
      det.R <- det(R)

      ## Function
      VIF <- function(i) {
        if (all(i %in% xn)) {
          ni <- !xn %in% i
          Ri <- R[i, i, drop = FALSE]
          Rni <- R[ni, ni, drop = FALSE]
          vif <- det(Ri) * det(Rni) / det.R
          (vif^(1 / (2 * length(i))))^2
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
      sapply(unlist(unname(XN)), function(i) {
        if (i %in% xn) 1 else NA
      })
    }

  }

  ## Apply recursively
  rMapply(VIF, m, SIMPLIFY = FALSE)

}


#' @title R-squared/Pseudo R-squared
#' @description Calculate R-squared or pseudo R-squared for a fitted model,
#'   defined as the squared multiple correlation between the observed and fitted
#'   values for the response variable. 'Adjusted' and 'Predicted' versions are
#'   also calculated (see Details).
#' @param mod A fitted model object, or a list or nested list of such objects.
#' @param data An optional dataset used to first re-fit the model(s).
#' @param adj,pred Logical. If \code{TRUE} (default), adjusted and/or predicted
#'   R-squared are also returned (the latter is not available for all model
#'   types).
#' @param offset Logical. If \code{TRUE}, include model offset(s) in the
#'   calculations (i.e. in the response and fitted values).
#' @param re.form For mixed models of class \code{"merMod"}, the formula for
#'   random effects to condition on when generating fitted values used in the
#'   calculation of R-squared. Defaults to \code{NULL}, meaning all random
#'   effects are included. See \code{\link[lme4]{predict.merMod}} for further
#'   specification details.
#' @param ... Not currently used.
#' @details Various approaches to the calculation of a goodness-of-fit measure
#'   for GLMs analogous to R-squared in the ordinary linear model have been
#'   proposed. Generally termed 'pseudo R-squared' measures, they include
#'   variance-based, likelihood-based, and distribution-specific approaches.
#'   Here however, a more straightforward definition is used, which can be
#'   applied to any model for which fitted values of the response variable are
#'   generated: R-squared is calculated as the squared (weighted) correlation
#'   between the observed and fitted values of the response (in the original
#'   units). This is simply the squared version of the correlation measure
#'   advocated by Zheng & Agresti (2000), itself an intuitive measure of
#'   goodness-of-fit describing the predictive power of a model. As the measure
#'   does not depend on any specific error distribution or model estimating
#'   procedure, it is also generally comparable across many different types of
#'   model (Kvalseth 1985). In the case of the ordinary linear model, the
#'   measure equals the more traditional R-squared based on sums of squares.
#'
#'   If \code{adj = TRUE} (default), the 'adjusted' R-squared value is also
#'   returned, which provides an estimate of the population - as opposed to
#'   sample - R-squared. This is achieved via an analytical formula which
#'   adjusts R-squared for the 'degrees of freedom' of the model (i.e. the ratio
#'   of observations to parameters). Here, this is calculated via the 'Pratt'
#'   rather than standard 'Ezekiel/Wherry' formula, shown in a previous
#'   simulation to be the most effective of a range of existing formulas at
#'   estimating the population R-squared across a range of model specification
#'   scenarios (Yin & Fan 2001). Adjusted R-squared can be used to safeguard
#'   against overfitting of the model to the original sample.
#'
#'   If \code{pred = TRUE} (default), a 'predicted' R-squared is also returned,
#'   which is calculated via the same formula as for R-squared but using
#'   cross-validated rather than standard fitted values. These are obtained by
#'   dividing model response residuals by the complement of the observation
#'   leverages (diagonals of the hat matrix), then subtracting these inflated
#'   'predicted' residuals from the response variable. This is essentially a
#'   short cut to obtaining out-of-sample predictions, normally arising via a
#'   leave-one-out cross-validation procedure (in a GLM they will not be exactly
#'   equal to such predictions). The resulting R-squared is an estimate of the
#'   R-squared that would occur were the model to be fitted to new data, and
#'   will be lower than the original - and likely also the adjusted - R-squared,
#'   highlighting the loss of explanatory power when predicting to new data.
#'   This measure is a variant of an
#'   \href{https://www.r-bloggers.com/can-we-do-better-than-r-squared/}{existing
#'   one}, calculated by substituting the 'PRESS' statistic, i.e. the sum of
#'   squares of the predicted residuals (Allen 1974), for the residual sum of
#'   squares in the classic R-squared formula.
#'
#'   For mixed models, the function will, by default, calculate all R-squared
#'   measures using fitted values incorporating both the fixed and random
#'   effects, thus encompassing all variation captured by the model. This is
#'   equivalent to the 'conditional' R-squared of Nakagawa \emph{et al.} (2017).
#'   To include only some or no random effects, simply set the appropriate
#'   formula using the argument \code{re.form}, which is passed directly to
#'   \code{predict.merMod}. If \code{re.form = NA}, R-squared is calculated for
#'   the fixed effects only - equivalent to the 'marginal' R-squared of Nakagawa
#'   \emph{et al.} (2017).
#'
#'   R-squared values produced by this function will always be bounded between
#'   zero (no fit) and one (perfect fit), meaning that any negative values
#'   arising from calculations will be rounded up to zero. Negative values
#'   typically mean that the fit is 'worse' than the null expectation of no
#'   relationship between the variables, which can be difficult to interpret in
#'   practice and in any case usually only occurs in rare situations, such as
#'   where the intercept is suppressed. Hence, for simplicity and ease of
#'   interpretation, values less than zero are presented here as a complete lack
#'   of model fit.
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
#'   models, it's not entirely clear what the sample size (n) in the formula
#'   should represent - the no. of observations? independent groups? something
#'   else? (the default interpretation of no. of observations is used). With all
#'   that being said, the value of standardised R-squared measures for even
#'   'rough' model fit assessment and comparison may outweigh such reservations,
#'   and the adjusted and predicted versions in particular may aid the user in
#'   diagnosing and preventing overfitting. They should NOT, however, replace
#'   other measures such as AIC or BIC for comparing and/or ranking competing
#'   models fit to the same data.
#' @return A numeric vector of the R-squared value(s), or an array, list of
#'   vectors/arrays, or nested list.
#' @references Allen, D. M. (1974). The Relationship Between Variable Selection
#'   and Data Augmentation and a Method for Prediction. \emph{Technometrics},
#'   \strong{16}(1), 125-127. \url{https://doi.org/gfgv57}
#'
#'   Kvalseth, T. O. (1985) Cautionary Note about R2. \emph{The American
#'   Statistician}, \strong{39}(4), 279-285. \url{https://doi.org/b8b782}
#'
#'   Nakagawa, S., Johnson, P.C.D. and Schielzeth, H. (2017) The coefficient of
#'   determination R2 and intra-class correlation coefficient from generalized
#'   linear mixed-effects models revisited and expanded. \emph{Journal of the
#'   Royal Society Interface} \strong{14}(134). \url{https://doi.org/gddpnq}
#'
#'   Yin, P. and Fan, X. (2001) Estimating R2 Shrinkage in Multiple Regression:
#'   A Comparison of Different Analytical Methods. \emph{The Journal of
#'   Experimental Education} \strong{69}(2), 203-224.
#'   \url{https://doi.org/fbdq5g}
#'
#'   Zheng, B. and Agresti, A. (2000) Summarizing the predictive power of a
#'   generalized linear model. \emph{Statistics in Medicine} \strong{19}(13),
#'   1771-1781. \url{https://doi.org/db7rfv}
#' @examples
#' ## Pseudo R-squared for mixed models
#' R2(Shipley.SEM)  # fixed + random ('conditional')
#' R2(Shipley.SEM, re.form = ~ (1 | tree))  # fixed + 'tree'
#' R2(Shipley.SEM, re.form = ~ (1 | site))  # fixed + 'site'
#' R2(Shipley.SEM, re.form = NA)  # fixed only ('marginal')
#'
#' ## Predicted R-squared: compare cross-validated predictions calculated/
#' ## approximated via the hat matrix to standard method (leave-one-out)
#'
#' ## Fit test models using Shipley data - compare lm vs glm
#' d <- na.omit(Shipley)
#' m <- lm(Live ~ Date + DD + lat, d)
#' # m <- glm(Live ~ Date + DD + lat, binomial, d)
#'
#' ## Manual CV predictions (leave-one-out)
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
#' ## Compare predictions (not exactly equal for GLMs)
#' all.equal(cvf1, cvf2)
#' # lm: TRUE; glm: "Mean relative difference: 1.977725e-06"
#' cor(cvf1, cvf2)
#' # lm: 1; glm: 0.9999987
#'
#' # NOTE: comparison not tested here for mixed models, as hierarchical data can
#' # complicate the choice of an appropriate leave-one-out procedure. However,
#' # there is no obvious reason why use of the leverage values (diagonals of the
#' # hat matrix) to estimate CV predictions shouldn't generalise to the mixed
#' # model case (at least for LMMs). In any case, users should exercise the
#' # appropriate caution in interpretation.
#' @export
R2 <- function(mod, data = NULL, adj = TRUE, pred = TRUE, offset = FALSE,
               re.form = NULL, ...) {

  m <- mod; d <- data; rf <- re.form

  ## Function
  R2 <- function(m) {

    ## Update model with any supplied data
    if (!is.null(d)) m <- eval(update(m, data = d, evaluate = FALSE))

    ## No. observations/parameters
    n <- nobs(m)
    k <- length(na.omit(if (isMer(m)) lme4::fixef(m) else coef(m)))
    i <- attr(terms(m), "intercept"); k <- k - i
    if (isMer(m)) k <- k + length(m@theta)

    ## R-squared
    R2 <- if (k > 0) {

      ## Model link function
      f <- if (isBet(m)) m$link$mean else family(m)
      lF <- f$linkfun; lI <- f$linkinv

      ## Model weights
      w <- weights(m)
      if (is.null(w)) w <- rep(1, n)
      w <- w[w > 0 & !is.na(w)]

      ## Model offset
      o <- if (!offset) {
        if (is.null(d)) d <- getData(m, subset = TRUE)
        mf <- suppressWarnings(model.frame(m, data = d))
        model.offset(mf)
      }
      if (is.null(o)) o <- 0

      ## Response and fitted values
      y <- getY(m, offset = offset); obs <- names(y)
      f <- lI(predict(m, re.form = rf)[obs] - o)

      ## R-squared
      R <- cov.wt(cbind(y, f), w, cor = TRUE)$cor[1, 2]
      if (is.na(R)) R <- 0
      if (R > 0) R^2 else 0

    } else 0

    ## Adjusted R-squared
    R2a <- if (adj) {
      if (R2 > 0) {

        ## Pratt formula
        R2a <- 1 - ((n - 3) * (1 - R2) / (n - k - i)) *
          (1 + (2 * (1 - R2) / (n - k - 2.3)))

        ## 'Standard' formula
        # R2a <- 1 - (1 - R2) * (n - i) / (n - k - i)

        if (R2a > 0) R2a else 0

      } else 0
    }

    ## Predictive R-squared
    R2p <- if (pred && !(isGls(m) || isMer(m) && isGlm(m))) {
      if (R2 > 0) {

        ## Leverage values (diagonals of the hat matrix, hii)
        hii <- hatvalues(m)[obs]
        s <- hii < 1

        ## CV fitted values (response minus 'predictive' residuals)
        pr <- (y - f) / (1 - hii)
        f <- y - pr

        ## Predictive R-squared
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
#' @param est A list or nested list of numeric vectors, comprising the model
#'   estimates. In the latter case, these should correspond to estimates for
#'   candidate models for each of a set of different response variables.
#' @param weights An optional numeric vector of weights to use for model
#'   averaging, or a named list of such vectors. The former should be supplied
#'   when \code{est} is a list, and the latter when it is a nested list (with
#'   matching list names). If \code{weights = "equal"} (default), a simple
#'   average is calculated instead.
#' @param est.names An optional vector of names used to extract and/or sort
#'   estimates from the output.
#' @param ... Not currently used.
#' @details This function can be used to calculate a weighted average of model
#'   estimates such as coefficients, fitted values, or residuals, where models
#'   are typically competing candidate models fit to the same response variable.
#'   Weights are typically a 'weight of evidence' type metric such as Akaike
#'   model weights (Burnham & Anderson 2002, Burnham \emph{et al.} 2011), which
#'   can be conveniently calculated in \emph{R} using packages such as
#'   \pkg{MuMIn} or \pkg{AICcmodavg}. However, numeric weights of any sort can
#'   be used. If none are supplied, the simple average is calculated instead.
#'
#'   Averaging is performed via the 'full'/'zero' rather than
#'   'subset'/'conditional'/'natural' method, meaning that zero is substituted
#'   for estimates for any 'missing' parameters (e.g. coefficients) prior to
#'   calculations. This provides a form of shrinkage and thus reduces
#'   \href{https://stackoverflow.com/questions/53055050/predicted-values-with-mumin-throwing-error-when-full-false}{estimate
#'   bias} (Burnham & Anderson 2002, Grueber \emph{et al.} 2011).
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
#'   \strong{65}(1), 23-35. \url{https://doi.org/c4mrns}
#'
#'   Dormann, C. F., Calabrese, J. M., Guillera-Arroita, G., Matechou, E., Bahn,
#'   V., Barton, K., ... Hartig, F. (2018). Model averaging in ecology: a review
#'   of Bayesian, information-theoretic, and tactical approaches for predictive
#'   inference. \emph{Ecological Monographs}, \strong{88}(4), 485-504.
#'   \url{https://doi.org/gfgwrv}
#'
#'   Grueber, C. E., Nakagawa, S., Laws, R. J., & Jamieson, I. G. (2011).
#'   Multimodel inference in ecology and evolution: challenges and solutions.
#'   \emph{Journal of Evolutionary Biology}, \strong{24}(4), 699-711.
#'   \url{https://doi.org/b7b5d4}
#'
#'   Walker, J. A. (2019). Model-averaged regression coefficients have a
#'   straightforward interpretation using causal conditioning. \emph{BioRxiv},
#'   133785. \url{https://doi.org/c8zt}
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
avgEst <-  function(est, weights = "equal", est.names = NULL, ...) {

  e <- est; w <- weights; en <- est.names

  ## Weights
  if (all(w == "equal")) {
    eqW <- function(i) rep(1, length(i))
    w <- if (any(sapply(e, isList))) lapply(e, eqW) else eqW(e)
  }

  ## Function
  avgEst <- function(e, w) {

    ## Sort names (for coefs)
    en2 <- unique(unlist(lapply(e, names)))
    num <- suppressWarnings(as.numeric(en2))
    if (all(is.na(num))) {
      i <- isInt(en2); r2 <- isR2(en2)
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
#'   deviation units, adjusted for multicollinearity among predictors.
#' @param mod A fitted model object, or a list or nested list of such objects.
#' @param weights An optional numeric vector of weights to use for model
#'   averaging, or a named list of such vectors. The former should be supplied
#'   when \code{mod} is a list, and the latter when it is a nested list (with
#'   matching list names). If set to \code{"equal"}, a simple average is
#'   calculated instead.
#' @param data An optional dataset used to first re-fit the model(s).
#' @param term.names An optional vector of term names used to extract and/or
#'   sort coefficients from the output.
#' @param cen.x,cen.y Logical, whether the intercept and coefficients should be
#'   calculated using mean-centred variables.
#' @param std.x,std.y Logical, whether coefficients should be scaled by the
#'   standard deviations of variables.
#' @param unique.x Logical, whether coefficients should be adjusted for
#'   multicollinearity among predictors.
#' @param refit.x Logical, whether the model should be re-fit with centred
#'   predictors.
#' @param r.squared Logical, whether R-squared values should also be returned.
#' @param ... Arguments to \code{R2}.
#' @details \code{stdCoeff} will calculate fully standardised coefficients in
#'   standard deviation units for a fitted model or list of models. It achieves
#'   this via adjusting the 'raw' model coefficients, so no standardisation of
#'   input variables is required beforehand. Users can simply specify the model
#'   with all variables in their original units and the function will do the
#'   rest. However, the user is free to scale and/or centre any input variables
#'   should they choose, which should not affect the outcome of standardisation
#'   (provided any scaling is by standard deviations). This may be desirable in
#'   some cases, such as to increase numerical stability during model fitting
#'   when variables are on widely different scales.
#'
#'   If arguments \code{cen.x} or \code{cen.y} are \code{TRUE}, model estimates
#'   will be calculated as if all predictors (x) and/or the response variable
#'   (y) were mean-centred prior to model-fitting (including any dummy variables
#'   arising from categorical predictors). Thus, for an ordinary linear model
#'   where centring of x and y is specified, the intercept will be zero - the
#'   mean (or weighted mean) of y. In addition, if \code{cen.x = TRUE} and there
#'   are interacting terms in the model, all coefficients for lower order terms
#'   of the interaction are adjusted using an expression which ensures that each
#'   main effect or lower order term is estimated at the mean values of the
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
#'   Mx3) + (\beta123 * Mx2 * Mx3)} (adapted from
#'   \href{https://stats.stackexchange.com/questions/65898/why-could-centering-independent-variables-change-the-main-effects-with-moderatio}{here})
#'
#'   In addition, if \code{std.x = TRUE} or \code{unique.x = TRUE} (see below),
#'   product terms for interactive effects will be recalculated using
#'   mean-centred variables, to ensure that standard deviations and variance
#'   inflation factors (VIF) for predictors are calculated correctly (the model
#'   must be re-fit for this latter purpose, to recalculate the
#'   variance-covariance matrix).
#'
#'   If \code{std.x = TRUE}, coefficients are standardised by multiplying by the
#'   standard deviations of predictor variables (or terms), while if \code{std.y
#'   = TRUE} they are divided by the standard deviation of the response. If the
#'   model is a GLM, this latter is calculated using the link-transformed
#'   response (or an estimate of same) generated using the function \code{getY}.
#'   If both arguments are true, the coefficients are regarded as 'fully'
#'   standardised in the traditional sense, often referred to as 'betas'.
#'
#'   If \code{unique.x = TRUE} (default), coefficients are adjusted for
#'   multicollinearity among predictors by dividing by the square root of the
#'   VIFs (Dudgeon 2016, Thompson \emph{et al.} 2017). If they have also been
#'   standardised by the standard deviations of x and y, this converts them to
#'   semipartial correlations, i.e. the correlation between the unique
#'   components of predictors (residualised on other predictors) and the
#'   response variable. This measure of effect size is arguably much more
#'   interpretable and useful than the traditional standardised coefficient, as
#'   it is always estimated independent of other predictors and so can more
#'   readily be compared both within and across models. Values range from zero
#'   to +/-1 rather than +/- infinity (as in the case of betas) - putting them
#'   on the same scale as the bivariate correlation between predictor and
#'   response. In the case of GLMs however, the measure is analogous but not
#'   exactly equal to the semipartial correlation, so its values may not always
#'   be bound between +/-1 (such cases are likely rare). Crucially, for ordinary
#'   linear models, the square of the semipartial correlation equals the
#'   increase in R-squared when that variable is added last in the model -
#'   directly linking the measure to model fit and 'variance explained'. See
#'   \href{https://www.daviddisabato.com/blog/2016/4/8/on-effect-sizes-in-multiple-regression}{here}
#'   for additional arguments in favour of the use of semipartial correlations.
#'
#'   If \code{refit.x = TRUE}, the model will be re-fit with any (newly-)centred
#'   continuous predictors. This will occur (and will normally be desired) when
#'   \code{cen.x} and \code{unique.x} are \code{TRUE} and there are interaction
#'   terms in the model, in order to calculate correct VIFs from the var-cov
#'   matrix. However, re-fitting may not be necessary in some cases - for
#'   example where predictors have already been centred (and whose values will
#'   not subsequently be resampled during bootstrapping) - and disabling this
#'   option may save time with larger models and/or bootstrap runs.
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
#'   Research}, \strong{51}(2-3), 139-153. \url{https://doi.org/gfww3f}
#'
#'   Thompson, C. G., Kim, R. S., Aloe, A. M., & Becker, B. J. (2017).
#'   Extracting the Variance Inflation Factor and Other Multicollinearity
#'   Diagnostics from Typical Regression Results. \emph{Basic and Applied Social
#'   Psychology}, \strong{39}(2), 81-90. \url{https://doi.org/gfww2w}
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
stdCoeff <- function(mod, weights = NULL, data = NULL, term.names = NULL,
                     cen.x = TRUE, cen.y = TRUE, std.x = TRUE, std.y = TRUE,
                     unique.x = TRUE, refit.x = TRUE, r.squared = FALSE, ...) {

  m <- mod; w <- weights; d <- data; bn <- term.names

  ## Function
  stdCoeff <- function(m) {

    ## Update model with any supplied data
    if (!is.null(d)) m <- eval(update(m, data = d, evaluate = FALSE))

    ## Coefficients
    b <- if (isMer(m)) lme4::fixef(m, add.dropped = TRUE) else coef(m)
    bn <- names(b); b <- na.omit(b)
    xn <- names(b)

    ## Intercept/no. parameters
    int <- isInt(xn[1]); if (int) xn <- xn[-1]
    k <- length(xn)

    ## Model weights
    w <- weights(m)
    if (is.null(w)) w <- rep(1, nobs(m))
    w <- w[w > 0 & !is.na(w)]

    ## Response
    y <- getY(m); obs <- names(y)

    ## Centre/standardise x
    if (k > 0) {

      ## Predictors
      dF <- function(...) data.frame(..., check.names = FALSE)
      if (is.null(d)) d <- getData(m, subset = TRUE)
      x <- dF(model.matrix(m, data = d))[xn]

      ## Interactions?
      inx <- any(isInx(xn))

      ## Centre predictors and adjust coefs/intercept
      if (cen.x) {

        ## Model offset
        mf <- suppressWarnings(model.frame(m, data = d))
        o <- model.offset(mf); if (is.null(o)) o <- 0

        ## For interactions, adjust coefs and centre predictors
        if (inx) {

          ## Predictor means
          sI <- function(x) {
            unlist(strsplit(x, "(?<!:):(?!:)", perl = TRUE))
          }
          XN <- lapply(labels(terms(m)), sI)
          x2 <- dF(model.matrix(reformulate(unlist(XN)), data = d))
          xm <- colMeans(x2)

          ## Adjust lower-order terms
          ## (ti = terms containing term i; ni = non-i components of ti)
          b[xn] <- sapply(xn, function(i) {
            bi <- b[[i]]; ii <- sI(i)
            ti <- xn[xn != i & sapply(xn, function(j) all(ii %in% sI(j)))]
            if (length(ti) > 0) {
              bi + sum(sapply(ti, function(j) {
                jj <- sI(j)
                ni <- jj[!jj %in% ii]
                prod(b[j], xm[ni])
              }))
            } else bi
          })

          ## Centre predictors (for correct product term SDs)
          if (std.x) {
            x <- dF(sapply(xn, function(i) {
              ii <- sI(i)
              xi <- sweep(x2[ii], 2, xm[ii])
              apply(xi, 1, prod)
            }))
          }

        }

        ## Adjust intercept (set to mean of predicted y)
        if (int) {
          f <- predict(m, re.form = NA)[obs] - o
          b[1] <- weighted.mean(f, w)
        }

      }

      ## Standardise by x
      if (std.x) {
        xs <- sapply(x, sdW, w)
        b[xn] <- b[xn] * xs
      }

      ## Calculate unique effects of predictors (adjust for multicollinearity)
      if (unique.x && k > 1) {

        ## Re-fit model with centred predictors
        ## (to calculate correct VIFs for interacting terms)
        m2 <- if (cen.x && inx && refit.x) {

          ## Update dataset
          ## (add response, weights, offset; set sum contrasts for factors)
          d2 <- data.frame(y, w, o, d)
          d2 <- dF(sapply(d2, function(i) {
            if (!is.numeric(i)) {
              i <- factor(i); contrasts(i) <- contr.sum(levels(i)); i
            } else i
          }))

          ## Update term names (scale numeric predictors)
          xn2 <- unlist(lapply(XN, function(i) {
            paste(sapply(i, function(j) {
              if (is.numeric(mf[, j])) paste0("scale(", j, ")") else j
            }), collapse = ":")
          }))

          ## Add random effects
          if (isMer(m)) {
            re <- lme4::findbars(formula(m))
            re <- sapply(re, function(i) paste0("(", deparse(i), ")"))
            xn2 <- c(xn2, re)
          }

          ## Rename any "y", "w", or "o" in terms
          if (any(c("y", "w", "o") %in% names(d))) {
            xn2 <- sapply(xn2, function(i) {
              i <- gsub("([^\\w.])", "~\\1~", i, perl = TRUE)
              i <- unlist(strsplit(i, "~"))
              s <- i %in% c("y", "w", "o")
              i[s] <- paste0(i[s], ".1")
              paste(i, collapse = "")
            })
          }

          ## Re-fit model
          eval(update(m, reformulate(xn2, "y", int), data = d2, weights = w,
                      offset = o, contrasts = NULL, evaluate = FALSE))

        } else m

        ## Divide coefs by square root of VIFs
        vif <- na.omit(VIF(m2, envir = environment()))
        b[xn] <- b[xn] / sqrt(vif)

      }

    }

    ## Centre/standardise y
    if (cen.y && int) {
      ym <- weighted.mean(y, w)
      if (isGlm(m)) {
        f <- if (isBet(m)) m$link$mean else family(m)
        ym <- f$linkfun(ym)
      }
      b[1] <- b[1] - ym
    }
    if (std.y) {
      if (isGlm(m)) y <- getY(m, link = TRUE)
      ys <- sdW(y, w)
      b <- b / ys
    }

    ## Return standardised coefficients
    b <- sapply(bn, function(i) unname(b[i]))
    if (r.squared) c(b, R2(m, ...)) else b

  }

  ## Apply recursively
  b <- rMapply(stdCoeff, m, SIMPLIFY = FALSE)

  ## Output coefs or weighted average
  if (!is.null(w) && isList(b)) avgEst(b, w, bn)
  else {
    if (!is.null(bn)) {
      f <- function(i) i[bn[bn %in% names(i)]]
      rMapply(f, b, SIMPLIFY = FALSE)
    } else b
  }

}


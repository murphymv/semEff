

#' @title Weighted Variance
#' @description Calculate the weighted variance of `x`.
#' @param x A numeric vector.
#' @param w A numeric vector of weights of the same length as `x`.
#' @param na.rm Logical, whether NAs in `x` should be removed.
#' @details Calculate the weighted variance of `x` via the weighted covariance
#'   matrix ([cov.wt()]). If no weights are supplied, the simple variance is
#'   returned instead. As in [weighted.mean()], `NA`s in `w` are not handled
#'   specially and will return `NA` as result.
#' @return A numeric value, the weighted variance of `x`.
#' @seealso [var()]
#' @examples
#' # Weighted variance
#' x <- rnorm(30)
#' w <- runif(30, 0, 1)
#' varW(x, w)
#'
#' # Simple variance
#' varW(x)
#' stopifnot(varW(x) == var(x))
#'
#' # NA handling
#' varW(c(x[1:29], NA), w, na.rm = TRUE)  # NA in x (removed)
#' varW(c(x[1:29], NA), w, na.rm = FALSE)  # NA in x (NA returned)
#' varW(x[1:29], w = c(w[1:29], NA))  # NA in w (NA returned)
#' @export
varW <- function(x, w = NULL, na.rm = FALSE) {
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
#' @description Calculate the weighted standard deviation of `x`.
#' @param ... Arguments to [varW()].
#' @details This is simply a wrapper for [varW()], applying the square root to
#'   the output.
#' @return A numeric value, the weighted standard deviation of `x`.
#' @seealso [sd()]
#' @export
sdW <- function(...) {
  sqrt(varW(...))
}


#' @title Get Model Data
#' @description Extract the data used to fit a model.
#' @param mod A fitted model object, or a list or nested list of such objects.
#' @param subset Logical. If `TRUE`, only observations used to fit the model(s)
#'   are returned (i.e. missing observations (`NA`) or those with zero weight
#'   are removed).
#' @param merge Logical. If `TRUE`, and `mod` is a list or nested list, a single
#'   dataset containing all variables used to fit models is returned (variables
#'   must be the same length).
#' @param env Environment in which to look for data (passed to [eval()]).
#'   Defaults to the [formula()] environment.
#' @details This is a simple convenience function to return the data used to fit
#'   a model, by evaluating the 'data' slot of the model call object. If the
#'   'data' argument of the model call was not specified, or is not a data frame
#'   (or coercible to such) containing all variables referenced in the model
#'   formula, an error will be thrown – this restriction is largely to ensure
#'   that a single coherent dataset of all model variables can be made available
#'   for resampling purposes.
#'
#'   If `mod` is a list of models and `merge = TRUE`, all (unique) variables
#'   used to fit models are merged into a single data frame. This will return an
#'   error if `subset = TRUE` results in datasets with different numbers of
#'   observations (rows).
#' @return A data frame of the variables used to fit the model(s), or a list or
#'   nested list of same.
#' @seealso [getCall()]
#' @examples
#' # Get data used to fit SEM from Shipley (2009)
#' head(getData(shipley.sem[[1]]))  # from single model
#' lapply(getData(shipley.sem), head)  # from SEM (list)
#' head(getData(shipley.sem, merge = TRUE))  # from SEM (single dataset)
#' @export
getData <- function(mod, subset = FALSE, merge = FALSE, env = NULL) {

  m <- mod

  # Function
  getData <- function(m) {

    # Data from 'data' argument of model call
    mc <- getCall(m)
    f <- formula(m)
    if (is.null(env)) env <- environment(f)
    d <- eval(mc$data, env)
    if (!is.null(d)) d <- data.frame(d) else
      stop("'data' argument of model call not specified, or data is NULL.")

    # All var names from model call
    f <- c(f, mc$subset, mc$weights, mc$offset, mc$correlation)
    vn <- unlist(lapply(f, all.vars))
    if (!all(vn %in% names(d)))  #  if (!all(!vn[grepl("tree.name",vn)] %in% names(d))) ??
      stop("'data' does not contain all variables used to fit model.")

    # Subset data for model observations?
    if (subset) {
      obs <- names(fitted(m))
      w <- weights(m)
      if (!is.null(w)) obs <- obs[w > 0 & !is.na(w)]
      d[obs, ]
    } else d

  }

  # Apply recursively
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


#' @title Get Model Design Matrix
#' @description Return the design matrix for a fitted model, with some
#'   additional options.
#' @param mod A fitted model object, or a list or nested list of such objects.
#'   Can also be a model formula(s) or character vector(s) of term names (in
#'   which case `data` must be supplied).
#' @param data An optional dataset, used to refit the model(s) and/or construct
#'   the design matrix.
#' @param contrasts Optional, a named list of contrasts to apply to factors (see
#'   the `contrasts.arg` argument of [model.matrix()] for specification). These
#'   will override any existing contrasts in the data or model call.
#' @param add.data Logical, whether to append data not specified in the model
#'   formula (with factors converted to dummy variables).
#' @param centre,scale Logical, whether to mean-centre and/or scale terms by
#'   standard deviations (for interactions, this is carried out prior to
#'   construction of product terms). Alternatively, a numeric vector of
#'   means/standard deviations (or other statistics) can be supplied, whose
#'   names must match term names.
#' @param as.df Logical, whether to return the matrix as a data frame (without
#'   modifying names).
#' @param merge Logical. If `TRUE`, and `mod` is a list or nested list, a single
#'   matrix containing all terms is returned (variables must be the same
#'   length).
#' @param env Environment in which to look for model data (if none supplied).
#'   Defaults to the [formula()] environment.
#' @details This is primarily a convenience function to enable more flexible
#'   construction of design matrices, usually for internal use and for further
#'   processing. Use cases include processing and/or return of terms which may
#'   not be present in a typical design matrix (e.g. constituents of product
#'   terms, dummy variables).
#' @return A matrix or data frame of model(s) terms, or a list or nested list of
#'   same.
#' @seealso [model.matrix()]
#' @examples
#' # Model design matrix (original)
#' m <- shipley.growth[[3]]
#' x1 <- model.matrix(m)
#' x2 <- getX(m)
#' stopifnot(all.equal(x1, x2, check.attributes = FALSE))
#'
#' # Using formula or term names (supply data)
#' d <- shipley
#' x1 <- getX(formula(m), data = d)
#' x2 <- getX(names(lme4::fixef(m)), data = d)
#' stopifnot(all.equal(x1, x2))
#'
#' # Scaled terms
#' head(getX(m, centre = TRUE, scale = TRUE))
#'
#' # Combined matrix for SEM
#' head(getX(shipley.sem, merge = TRUE))
#' head(getX(shipley.sem, merge = TRUE, add.data = TRUE))  # add other variables
#' @export
getX <- function(mod, data = NULL, contrasts = NULL, add.data = FALSE,
                 centre = FALSE, scale = FALSE, as.df = FALSE, merge = FALSE,
                 env = NULL) {

  m <- mod; d <- data; ct <- contrasts; cn <- centre; sc <- scale

  # Function
  getX <- function(m) {

    # Function to convert to data frame without modifying names
    dF <- function(...) {
      data.frame(..., check.names = FALSE)
    }

    # Function to evaluate a term name (x) in data (d)
    eT <- function(x, d) {
      eT <- function(x) {
        if (!x %in% names(d)) {
          eval(parse(text = x), d)
        } else d[[x]]
      }
      if (grepl("poly\\(.*[0-9]$", x)) {
        n <- nchar(x)
        xd <- eT(substr(x, 1, n - 1))
        xd[, substr(x, n, n)]
      } else eT(x)
    }

    # Refit model with any supplied data
    if (isMod(m) && !is.null(d)) {
      m <- eval(update(m, data = d, evaluate = FALSE))
      env <- environment()
    }

    # Data
    if (isMod(m)) d <- getData(m, subset = TRUE, env = env)
    if (is.null(d))
      stop("'data' must be supplied.")

    # Factor contrasts
    if (isMod(m) && is.null(ct)) {
      ct <- suppressWarnings(
        attr(model.matrix(m, data = d), "contrasts")
      )
    }

    # Formulas
    i <- if (!is.character(m)) attr(terms(m), "intercept") else any(isInt(m))
    xn <- if (!is.character(m)) labels(terms(m)) else m
    xn <- xn[!isInt(xn) & !grepl("[|]", xn)]
    fd <- reformulate(c(xn, names(d)), intercept = i)
    fx <- reformulate(xn, intercept = i)

    # Design matrix (x)
    d <- model.frame(fd, data = d, na.action = na.pass)
    x <- suppressWarnings(
      dF(model.matrix(fx, data = d, contrasts.arg = ct))
    )
    obs <- rownames(x)
    
    # Term/variable names
    xn <- names(x)
    xn2 <- names(model.frame(fx, data = d))
    
    # Add other variables from data (convert factors to dummy vars)
    d <- sapply(names(d), function(i) {
      di <- d[[i]]
      if (!is.numeric(di)) {
        if (i %in% xn2) {
          di <- as.factor(di)
          ci <- if (!is.null(ct)) list(di = ct[[i]])
          di <- cbind(
            model.matrix( ~ 0 + di),
            model.matrix( ~ di, contrasts.arg = ci)
          )
          colnames(di) <- gsub("^di", i, colnames(di))
          di
        } else 0
      } else di
    }, simplify = FALSE)
    d <- dF(do.call(cbind, d))
    d <- d[unique(names(d))]
    x <- dF(x, d[!names(d) %in% xn])

    # Means/SDs (for centring/scaling)
    xm <- sapply(names(x), function(i) {
      if (!isFALSE(cn)) {
        if (i %in% names(cn)) cn[[i]] else mean(x[[i]])
      } else 0
    })
    xs <- sapply(names(x), function(i) {
      if (!isFALSE(sc)) {
        if (i %in% names(sc)) sc[[i]] else sd(x[[i]])
      } else 1
    })

    # Term names to return
    if (add.data) {
      rn <- if (isMer(m)) {
        unlist(lapply(lme4::VarCorr(m), rownames))
      }
      xn <- unique(c(names(x), rn))
    }

    # Construct final matrix
    x <- sapply(xn, function(i) {
      if (!isInt(i)) {
        ii <- unlist(strsplit(i, "(?<!:):(?!:)", perl = TRUE))
        xi <- sapply(ii, eT, x)
        xn <- colnames(xi)
        xi <- sweep(xi, 2, xm[xn])
        xi <- sweep(xi, 2, xs[xn], "/")
        apply(xi, 1, prod)
      } else x[, i]
    })
    rownames(x) <- obs

    # Return
    if (as.df) dF(x) else x

  }

  # Apply recursively
  rgetX <- function(m) {
    if (isList(m)) {
      x <- lapply(m, rgetX)
      if (merge) {
        x <- do.call(cbind, unname(x))
        x[, unique(colnames(x)), drop = FALSE]
      } else x
    } else getX(m)
  }
  rgetX(m)

}


#' @title Get Model Term Names
#' @description Extract term names from a fitted model object.
#' @param mod A fitted model object, or a list or nested list of such objects.
#' @param intercept Logical, whether the intercept should be included.
#' @param aliased Logical, whether names of aliased terms should be included
#'   (see Details).
#' @param list Logical, whether names should be returned as a list, with all
#'   multi-coefficient terms grouped under their main term names.
#' @param env Environment in which to look for model data (used to construct the
#'   model frame). Defaults to the [formula()] environment.
#' @details Extract term names from a fitted model. Names of terms for which
#'   coefficients cannot be estimated are also included if `aliased = TRUE`
#'   (default). These may be terms which are perfectly correlated with other
#'   terms in the model, so that the model design matrix is rank deficient.
#' @return A character vector or list/nested list of term names.
#' @examples
#' # Term names from Shipley SEM
#' m <- shipley.sem
#' xNam(m)
#' xNam(m, intercept = FALSE)
#'
#' # Model with different types of predictor (some multi-coefficient terms)
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
xNam <- function(mod, intercept = TRUE, aliased = TRUE, list = FALSE,
                 env = NULL) {

  m <- mod

  # Function
  xNam <- function(m) {

    # Term names
    tt <- terms(m)
    xn <- labels(tt)
    int <- attr(tt, "intercept")
    if (int) xn <- c("(Intercept)", xn)

    # Expand interaction terms (list)
    XN <- sapply(xn, function(i) {
      unlist(strsplit(i, "(?<!:):(?!:)", perl = TRUE))
    }, simplify = FALSE)

    # Predictors
    d <- getData(m, env = env)
    x <- suppressWarnings(model.frame(m, data = d))
    x <- x[names(x) %in% unlist(unname(XN))]
    xn2 <- unique(c(xn, names(x)))

    # Factor contrasts
    ct <- suppressWarnings(
      attr(model.matrix(m, data = d), "contrasts")
    )

    # Expand factor/matrix terms (list)
    XN2 <- sapply(xn2, function(i) {
      if (i %in% names(x)) {
        xi <- x[[i]]
        j <- colnames(xi)
        n <- ncol(xi)
        if (!is.numeric(xi)) {
          xi <- as.factor(xi)
          ci <- list(xi = ct[[i]])
          xi <- model.matrix( ~ xi, contrasts.arg = ci)
          j <- gsub("^xi", "", colnames(xi)[-1])
          j <- c(levels(xi), j)
          n <- nlevels(xi)
        }
        j <- unique(c("", seq(n), j))
        paste0(i, j)
      } else i
    }, simplify = FALSE)

    # Expand interaction terms involving factors/matrices (list)
    XN <- sapply(xn, function(i) {
      if (isInx(i)) {
        j <- expand.grid(XN2[XN[[i]]])
        apply(j, 1, paste, collapse = ":")
      } else XN2[[i]]
    }, simplify = FALSE)

    # Drop some names? (aliased, subsetted factor levels, etc.)
    b <- if (isMer(m)) lme4::fixef(m, add.dropped = TRUE) else coef(m)
    bn <- names(b)
    if (!aliased) bn <- bn[!is.na(b)]
    XN <- lapply(XN, function(i) bn[bn %in% i])
    XN <- XN[sapply(XN, length) > 0]

    # Drop intercept?
    if (!intercept && int) XN <- XN[-1]

    # Return names
    if (!list) unlist(unname(XN)) else XN

  }

  # Apply recursively
  rMapply(xNam, m, SIMPLIFY = FALSE)

}


#' @title Generalised Link Transformation
#' @description Transform a numeric variable using a GLM link function, or
#'   return an estimate of same.
#' @param x a positive numeric vector, corresponding to a variable to be
#'   transformed. Can also be a list or nested list of such vectors.
#' @param family Optional, the error distribution family containing the link
#'   function which will be used to transform `x` (see [family()] for
#'   specification details). If not supplied, it is determined from `x` (see
#'   Details).
#' @param force.est Logical, whether to force the return of the estimated rather
#'   than direct transformation, where the latter is available (i.e. does not
#'   contain undefined values).
#' @details `glt()` can be used to provide a 'generalised' transformation of a
#'   numeric variable, using the link function from a generalised linear model
#'   (GLM) fit to the variable. The transformation is generalised in the sense
#'   that it can extend even to cases where a standard link transformation would
#'   produce undefined values. It achieves this by using an estimate based on
#'   the 'working' response variable of the GLM (see below). If the error
#'   distribution `family` is not specified (default), then it is determined
#'   (roughly) from `x`, with `binomial(link = "logit")` used when all x <= 1
#'   and `poisson(link = "log")` otherwise. Although the function is generally
#'   intended for variables with a binomial or Poisson distribution, any
#'   variable which can be fit using [glm()] can be supplied. One of the key
#'   purposes of `glt()` is to allow the calculation of fully standardised
#'   effects (coefficients) for GLMs (in which case `x` = the response
#'   variable), while it can also facilitate the proper calculation of SEM
#'   indirect effects (see below).
#'
#'   **Estimating the direct link transformation**
#'
#'   A key challenge in generating fully standardised effects for a GLM with a
#'   non-gaussian link function is the difficulty in calculating appropriate
#'   standardised ranges (typically the standard deviation) for the response
#'   variable in the link scale. This is because a direct transformation of the
#'   response will often produce undefined values. Although methods for
#'   circumventing this issue by indirectly estimating the variance of the
#'   response on the link scale have been proposed – including a
#'   latent-theoretic approach for binomial models (McKelvey & Zavoina, 1975)
#'   and a more general variance-based method using pseudo R-squared (Menard,
#'   2011) – here an alternative approach is used. Where transformed values are
#'   undefined, the function will instead return the synthetic 'working'
#'   response from the iteratively reweighted least squares (IRLS) algorithm of
#'   the GLM (McCullagh & Nelder, 1989). This is reconstructed as the sum of the
#'   linear predictor and the working residuals – with the latter comprising the
#'   error of the model on the link scale. The advantage of this approach is
#'   that a relatively straightforward 'transformation' of any non-gaussian
#'   response is readily attainable in all cases. The standard deviation (or
#'   other relevant range) can then be calculated using values of the
#'   transformed response and used to scale the effects. An additional benefit
#'   for piecewise SEMs is that the transformed rather than original response
#'   can be specified as a predictor in other models, ensuring that standardised
#'   indirect and total effects are calculated correctly (i.e. using the same
#'   units).
#'
#'   To ensure a high level of 'accuracy' in the working response – in the sense
#'   that the inverse-transformation is practically indistinguishable from the
#'   original response variable – the function uses the following iterative
#'   fitting procedure to calculate a 'final' working response:
#'
#'   1. A new GLM of the same error family is fit with the original response
#'   variable as both predictor and response, and using a single IRLS iteration.
#'
#'   2. The working response is reconstructed from this model.
#'
#'   3. The inverse transformation of the working response is calculated.
#'
#'   4. If the inverse transformation is 'effectively equal' to the original
#'   response (tested using [all.equal()] with the default tolerance of
#'   `1.5e-8`), the working response is returned; otherwise, the GLM is refit
#'   with the working response now as the predictor, and steps 2-4 are repeated
#'   – each time with an additional IRLS iteration.
#'
#'   This approach will generate a very reasonable transformation of the
#'   response variable, which will also be practically indistinguishable from
#'   the direct transformation, where this can be compared (see Examples). It
#'   also ensures that the transformed values, and hence the standard deviation,
#'   are the same for any GLM fitting the same response (provided it uses the
#'   same link function) – facilitating model comparisons, selection, and
#'   averaging.
#'
#' @note As we often cannot directly observe the GLM response variable on the
#'   link scale, any method estimating its values or statistics will be wrong to
#'   some degree. The heuristic approach described here aims to reduce this
#'   error as far as (reasonably) possible, while also generating standardised
#'   effects whose interpretation most closely resembles those of the ordinary
#'   linear model. The solution of using the working response from the GLM to
#'   scale effects is a practical, but reasonable one, and one that takes
#'   advantage of modern computing power to minimise error through iterative
#'   model fitting. An added bonus is that the estimated variance is constant
#'   across models fit to the same response variable, which cannot be said of
#'   previous methods (Menard, 2011). The overall approach would be classed as
#'   'observed-empirical' by Grace et al. (2018), as it utilises model error
#'   variance (the estimated working residuals) rather than theoretical
#'   distribution-specific variance.
#' @return A numeric vector of the transformed values, or an array, list of
#'   vectors/arrays, or nested list.
#' @references Grace, J. B., Johnson, D. J., Lefcheck, J. S., & Byrnes, J. E. K.
#'   (2018). Quantifying relative importance: computing standardized effects in
#'   models with binary outcomes. *Ecosphere*, *9*, e02283. \doi{10/gdm5bj}
#'
#'   McCullagh P., & Nelder, J. A. (1989). *Generalized Linear Models* (2nd
#'   Edition). Chapman and Hall.
#'
#'   McKelvey, R. D., & Zavoina, W. (1975). A statistical model for the analysis
#'   of ordinal level dependent variables. *The Journal of Mathematical
#'   Sociology*, *4*(1), 103-120. \doi{10/dqfhpp}
#'
#'   Menard, S. (2011). Standards for Standardized Logistic Regression
#'   Coefficients. *Social Forces*, *89*, 1409-1428. \doi{10/bvxb6s}
#' @examples
#' # Compare estimate with a direct link transformation
#' # (test with a poisson variable, log link)
#' set.seed(13)
#' y <- rpois(30, lambda = 10)
#' yl <- unname(glt(y, force.est = TRUE))
#'
#' # Effectively equal?
#' all.equal(log(y), yl)
#' # TRUE
#'
#' # Actual difference...
#' all.equal(log(y), yl, tolerance = .Machine$double.eps)
#' # "Mean relative difference: 2.489317e-10"
#' @export
glt <- function(x, family = NULL, force.est = FALSE) {

  f <- family

  # Function
  glt <- function(x) {

    # Error family
    if (is.character(f)) {
      f <- get(f, mode = "function", envir = parent.frame())
    }
    if (is.function(f)) f <- f()
    if (is.null(f)) {
      f <- if (all(x <= 1)) binomial() else poisson()
    }

    # Transform variable to link scale
    xl <- f$linkfun(as.numeric(x))

    # Return transformation or estimate (GLM working response)
    if (any(is.infinite(xl)) || force.est) {
      x2 <- x
      m <- suppressWarnings(
        do.call(glm, list(x ~ x2, family = f, control = list(maxit = 1),
                          na.action = na.exclude))
      )
      i <- 0
      repeat {
        xl <- predict(m) + resid(m, type = "working")
        xli <- f$linkinv(xl)
        eql <- isTRUE(all.equal(x, xli, tolerance = sqrt(.Machine$double.eps),
                                check.names = FALSE))
        if (eql) {
          return(xl)
        }
        i <- i + 1
        m <- suppressWarnings(
          update(m, . ~ xl, control = list(maxit = i))
        )
      }
    } else xl

  }

  # Apply recursively
  rMapply(glt, x)

}


#' @title Get Model Error Distribution Family
#' @description Extract the error distribution family for a fitted model as a
#'   [family()] object.
#' @param mod A fitted model object, or a list or nested list of such objects.
#' @details `getFamily()` will return a(n appropriate) family object for a range
#'   of different model classes, even those without an existing family method.
#' @return A model [family()] object, or a list or nested list of such objects.
#' @examples
#' # SEM model error distributions
#' getFamily(shipley.sem)
#' @export
getFamily <- function(mod) {
  
  m <- mod
  
  # Function
  getFamily <- function(m) {
    if (isGls(m)) gaussian() else {
      if (isBet(m)) m$link$mean else family(m)
    }
  }

  # Apply recursively
  rMapply(getFamily, m)
  
}


#' @title Get Model Response Variable
#' @description Extract the response variable from a fitted model on the
#'   original or link scale.
#' @param mod A fitted model object, or a list or nested list of such objects.
#' @param data An optional dataset, used to first refit the model(s).
#' @param link Logical. If `TRUE`, return the GLM response variable on the link
#'   scale (see Details).
#' @param offset Logical. If `TRUE`, include model offset(s) in the response.
#' @param env Environment in which to look for model data (if none supplied).
#'   Defaults to the [formula()] environment.
#' @details `getY()` will return the response variable from a model by summing
#'   the fitted values and the response residuals. If `link = TRUE` and the
#'   model is a GLM, the response is transformed to the link scale. If this
#'   results in undefined values, an estimate based on the 'working' response
#'   variable of the GLM is returned instead (see [glt()]).
#'
#'   Any offset variables are subtracted from the response by default. This
#'   means that, for example, rates rather than raw counts will be returned for
#'   poisson GLMs (where applicable).
#' @return A numeric vector comprising the response variable on the original or
#'   link scale, or an array, list of vectors/arrays, or nested list.
#' @examples
#' # All SEM responses (original scale)
#' head(getY(shipley.sem))
#'
#' # Estimated response in link scale from binomial model
#' head(getY(shipley.sem$Live, link = TRUE))
#' @export
getY <- function(mod, data = NULL, link = FALSE, offset = FALSE, env = NULL) {

  m <- mod; d <- data

  # Function
  getY <- function(m) {

    # Refit model with any supplied data
    if (!is.null(d)) {
      m <- eval(update(m, data = d, evaluate = FALSE))
      env <- environment()
    }

    # Model error family
    f <- getFamily(m)

    # Model weights and offset
    w <- weights(m)
    o <- if (!offset) {
      d <- getData(m, subset = TRUE, env = env)
      mf <- suppressWarnings(model.frame(m, data = d))
      model.offset(mf)
    }

    # Model response
    y <- fitted(m) + resid(m, type = "response")
    if (!is.null(w)) y <- y[w > 0 & !is.na(w)]
    if (!is.null(o)) y <- f$linkinv(f$linkfun(y) - o)
    y[zapsmall(y) == 0] <- 0
    an <- names(attributes(y))
    attributes(y)[an != "names"] <- NULL

    # Return response on original or link scale
    if (isGlm(m) && link) glt(y, family = f) else y

  }

  # Apply recursively
  rMapply(getY, m)

}


#' @title Generalised Variance Inflation Factors
#' @description Calculate generalised variance inflation factors for terms in a
#'   fitted model(s) via the variance-covariance matrix of coefficients.
#' @param mod A fitted model object, or a list or nested list of such objects.
#' @param data An optional dataset, used to first refit the model(s).
#' @param env Environment in which to look for model data (if none supplied).
#'   Defaults to the [formula()] environment.
#' @details `VIF()` calculates generalised variance inflation factors (GVIF) as
#'   described in Fox & Monette (1992), and also implemented in
#'   [`car::vif()`](https://rdrr.io/cran/car/man/vif.html). However, whereas
#'   `vif()` returns both GVIF and GVIF^(1/(2*Df)) values, `VIF()` simply
#'   returns the squared result of the latter measure, which equals the standard
#'   VIF for single-coefficient terms and is the equivalent measure for
#'   multi-coefficient terms (e.g. categorical or polynomial). Also, while
#'   `vif()` returns values per model term (i.e. predictor variable), `VIF()`
#'   returns values per coefficient, meaning that the same value will be
#'   returned per coefficient for multi-coefficient terms. Finally, `NA` is
#'   returned for any terms which could not be estimated in the model (e.g.
#'   aliased).
#' @return A numeric vector of the VIFs, or an array, list of vectors/arrays, or
#'   nested list.
#' @references Fox, J., & Monette, G. (1992). Generalized Collinearity
#'   Diagnostics. *Journal of the American Statistical Association*, *87*,
#'   178-183. \doi{10/dm9wbw}
#' @examples
#' # Model with two correlated terms
#' m <- shipley.growth[[3]]
#' VIF(m)  # Date & DD somewhat correlated
#' VIF(update(m, . ~ . - DD))  # drop DD
#'
#' # Model with different types of predictor (some multi-coefficient terms)
#' d <- data.frame(
#'   y = rnorm(100),
#'   x1 = rnorm(100),
#'   x2 = as.factor(rep(c("a", "b", "c", "d"), each = 25)),
#'   x3 = rep(1, 100)
#' )
#' m <- lm(y ~ poly(x1, 2) + x2 + x3, data = d)
#' VIF(m)
#' @export
VIF <- function(mod, data = NULL, env = NULL) {

  m <- mod; d <- data

  # Function
  VIF <- function(m) {

    # Refit model with any supplied data
    if (!is.null(d)) {
      m <- eval(update(m, data = d, evaluate = FALSE))
      env <- environment()
    }

    # Term names
    XN <- xNam(m, intercept = FALSE, list = TRUE, env = env)
    xn <- xNam(m, intercept = FALSE, aliased = FALSE, env = env)

    # VIFs
    if (length(xn) > 1) {

      # T/F for terms as matrices
      if (is.null(d)) d <- getData(m, env = env)
      mf <- suppressWarnings(model.frame(m, data = d))
      mat <- sapply(names(XN), function(i) {
        if (i %in% names(mf)) class(mf[, i])[1] == "matrix" else FALSE
      })

      # Var-cov & cor matrix
      V <- as.matrix(vcov(m))[xn, xn]
      R <- cov2cor(V)
      det.R <- det(R)

      # Function
      VIF <- function(i) {
        if (all(i %in% xn)) {
          j <- !xn %in% i
          Ri <- R[i, i, drop = FALSE]
          Rj <- R[j, j, drop = FALSE]
          vif <- det(Ri) * det(Rj) / det.R
          p <- 1 / (2 * length(i))
          (vif^p)^2
        } else NA
      }

      # Calculate
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

  # Apply recursively
  rMapply(VIF, m, SIMPLIFY = FALSE)

}


#' @title Root Variance Inflation Factors
#' @description Calculate root variance inflation factors (RVIF) for terms in a
#'   fitted model(s), i.e. the square root of the (generalised) VIFs.
#' @param ... Arguments to [VIF()].
#' @details RVIFs quantify the inflation of estimate standard errors due to
#'   multicollinearity among predictors, and also of estimates themselves
#'   compared to the 'unique' (residualised) effects. RVIFs may often be more
#'   practical than VIFs for assessing multicollinearity, relating more directly
#'   to the parameters of interest.
#' @return A numeric vector of the RVIFs, or an array, list of vectors/arrays,
#'   or nested list.
#' @export
RVIF <- function(...) {
  rMapply(sqrt, VIF(...), SIMPLIFY = FALSE)
}


#' @title R-squared
#' @description Calculate (Pseudo) R-squared for a fitted model, defined here as
#'   the squared multiple correlation between the observed and fitted values for
#'   the response variable. 'Adjusted' and 'Predicted' versions are also
#'   calculated (see Details).
#' @param mod A fitted model object, or a list or nested list of such objects.
#' @param data An optional dataset, used to first refit the model(s).
#' @param adj,pred Logical. If `TRUE` (default), adjusted and/or predicted
#'   R-squared are also returned (the latter is not available for all model
#'   types).
#' @param offset Logical. If `TRUE`, include model offset(s) in the calculations
#'   (i.e. in the response variable and fitted values).
#' @param re.form For mixed models of class `"merMod"`, the formula for random
#'   effects to condition on when generating fitted values used in the
#'   calculation of R-squared. Defaults to `NULL`, meaning all random effects
#'   are included. See [predict.merMod()] for further specification details.
#' @param type The type of correlation coefficient to use. Can be `"pearson"`
#'   (default) or `"spearman"`.
#' @param adj.type The type of adjusted R-squared estimator to use. Can be
#'   `"olkin-pratt"` (default) or `"ezekiel"`. See Details.
#' @param positive.only Logical, whether to return only positive values for
#'   R-squared (negative values replaced with zero).
#' @param env Environment in which to look for model data (if none supplied).
#'   Defaults to the [formula()] environment.
#' @details Various approaches to the calculation of a goodness of fit measure
#'   for GLMs analogous to R-squared in the ordinary linear model have been
#'   proposed. Generally termed 'pseudo R-squared' measures, they include
#'   variance-based, likelihood-based, and distribution-specific approaches.
#'   Here however, a more straightforward definition is used, which can be
#'   applied to any model for which fitted values of the response variable are
#'   generated: R-squared is calculated as the squared (weighted) correlation
#'   between the observed and fitted values of the response (in the original
#'   units). This is simply the squared version of the correlation measure
#'   advocated by Zheng & Agresti (2000), itself an intuitive measure of
#'   goodness of fit describing the predictive power of a model. As the measure
#'   does not depend on any specific error distribution or model estimating
#'   procedure, it is also generally comparable across many different types of
#'   model (Kvalseth, 1985). In the case of the ordinary linear model, the
#'   measure is exactly equal to the traditional R-squared based on sums of
#'   squares.
#'
#'   If `adj = TRUE` (default), the 'adjusted' R-squared value is also returned,
#'   which provides an estimate of the population – as opposed to sample –
#'   R-squared. This is achieved via an analytical formula which adjusts
#'   R-squared using the 'degrees of freedom' of the model (i.e. the ratio of
#'   observations to parameters), helping to counter multiple R-squared's
#'   positive bias and guard against overfitting of the model to noise in the
#'   original sample. By default, this is calculated via the exact 'Olkin-Pratt'
#'   estimator, shown in recent simulations to be the optimal unbiased
#'   population R-squared estimator across a range of estimators and
#'   specification scenarios (Karch, 2020), and thus a good general first
#'   choice, even for smaller sample sizes. Setting `adj.type = "ezekiel"`
#'   however will use the simpler and more common 'Ezekiel' formula, which can
#'   be more appropriate where minimising the mean squared error (MSE) of the
#'   estimate is more important than strict unbiasedness (Hittner, 2019; Karch,
#'   2020).
#'
#'   If `pred = TRUE` (default), a 'predicted' R-squared is also returned, which
#'   is calculated via the same formula as for R-squared but using
#'   cross-validated, rather than original, fitted values. These are obtained by
#'   dividing the model residuals (in the response scale) by the complement of
#'   the observation leverages (diagonals of the hat matrix), then subtracting
#'   these inflated 'predicted' residuals from the response variable. This is
#'   essentially a short cut to obtaining 'out-of-sample' predictions, normally
#'   arising via a 'leave-one-out' cross-validation procedure (they are not
#'   exactly equal for GLMs). The resulting R-squared is an estimate of the
#'   R-squared that would result were the model fit to new data, and will be
#'   lower than the original – and likely also the adjusted – R-squared,
#'   highlighting the loss of explanatory power due to sample noise. Predicted
#'   R-squared [may be a more powerful and general indicator of overfitting than
#'   adjusted
#'   R-squared](https://stats.stackexchange.com/questions/242770/difference-between-adjusted-r-squared-and-predicted-r-squared),
#'   as it provides a true out-of-sample test. This measure is a variant of an
#'   [existing
#'   one](https://www.r-bloggers.com/2014/05/can-we-do-better-than-r-squared/),
#'   calculated by substituting the 'PRESS' statistic, i.e. the sum of squares
#'   of the predicted residuals (Allen, 1974), for the residual sum of squares
#'   in the classic R-squared formula. It is not calculated here for GLMMs, as
#'   the interpretation of the hat matrix is not reliable (see
#'   [hatvalues.merMod()]).
#'
#'   For models fitted with one or more offsets, these will be removed by
#'   default from the response variable and fitted values prior to calculations.
#'   Thus R-squared will measure goodness of fit only for the predictors with
#'   estimated, rather than fixed, coefficients. This is likely to be the
#'   intended behaviour in most circumstances, though if users wish to include
#'   variation due to the offset(s) they can set `offset = TRUE`.
#'
#'   For mixed models, the function will, by default, calculate all R-squared
#'   measures using fitted values incorporating both the fixed and random
#'   effects, thus encompassing all variation captured by the model. This is
#'   equivalent to the 'conditional' R-squared of Nakagawa et al. (2017) (though
#'   see that reference for a more advanced approach to R-squared for mixed
#'   models). To include only some or no random effects, simply set the
#'   appropriate formula using the argument `re.form`, which is passed directly
#'   to [predict.merMod()]. If `re.form = NA`, R-squared is calculated for the
#'   fixed effects only, i.e. the 'marginal' R-squared of Nakagawa et al.
#'   (2017).
#'
#'   As R-squared is calculated here as a squared correlation, the `type` of
#'   correlation coefficient can also be specified. Setting this to `"spearman"`
#'   may be desirable in some cases, such as where the relationship between
#'   response and fitted values is not necessarily bivariate normal or linear,
#'   and a correlation of the ranks may be more informative and/or general. This
#'   purely monotonic R-squared can also be considered a [useful goodness of fit
#'   measure](https://stats.stackexchange.com/questions/44268/reporting-coefficient-of-determination-using-spearmans-rho),
#'   and may be more appropriate for comparisons between GLMs and ordinary
#'   linear models in some scenarios.
#'
#'   R-squared values produced by this function will by default be in the range
#'   0-1, meaning that any negative values arising from calculations will be
#'   converted to zero. Negative values essentially mean that the fit is 'worse'
#'   than the null expectation of no relationship between the variables, which
#'   can be difficult to interpret in practice and in any case usually only
#'   occurs in rare situations, such as where the intercept is suppressed or
#'   where a low value of R-squared is adjusted downwards via an analytic
#'   estimator. Such values are also 'impossible' in practice, given that
#'   R-squared is a strictly positive measure (as generally known). Hence, for
#'   simplicity and ease of interpretation, values less than zero are presented
#'   as a complete lack of model fit. This is also recommended by Shieh (2008),
#'   who shows for adjusted R-squared that such 'positive-part' estimators have
#'   lower MSE in estimating the population R-squared (though higher bias). To
#'   allow return of negative values however, set `positive.only = FALSE`. This
#'   may be desirable for simulation purposes, and/or where strict unbiasedness
#'   is prioritised.
#'
#' @note Caution must be exercised in interpreting the values of any pseudo
#'   R-squared measure calculated for a GLM or mixed model (including those
#'   produced by this function), as such measures do not hold all the properties
#'   of R-squared in the ordinary linear model and as such may not always behave
#'   as expected. Care must also be taken in comparing the measures to their
#'   equivalents from ordinary linear models, particularly the adjusted and
#'   predicted versions, as assumptions and/or calculations may not generalise
#'   well. With that being said, the value of standardised R-squared measures
#'   for even 'rough' model fit assessment and comparison may outweigh such
#'   reservations, and the adjusted and predicted versions in particular may aid
#'   the user in diagnosing and preventing overfitting. They should NOT,
#'   however, replace other measures such as AIC or BIC for formally comparing
#'   and/or ranking competing models fit to the same response variable.
#' @return A numeric vector of the R-squared value(s), or an array, list of
#'   vectors/arrays, or nested list.
#' @references Allen, D. M. (1974). The Relationship Between Variable Selection
#'   and Data Augmentation and a Method for Prediction. *Technometrics*,
#'   *16*(1), 125-127. \doi{10/gfgv57}
#'
#'   Hittner, J. B. (2019). Ezekiel’s classic estimator of the population
#'   squared multiple correlation coefficient: Monte Carlo-based extensions and
#'   refinements. *The Journal of General Psychology*, *147*(3), 213–227.
#'   \doi{10/gk53wb}
#'
#'   Karch, J. (2020). Improving on Adjusted R-Squared. *Collabra: Psychology*,
#'   *6*(1). \doi{10/gkgk2v}
#'
#'   Kvalseth, T. O. (1985). Cautionary Note about R2. *The American
#'   Statistician*, *39*(4), 279-285. \doi{10/b8b782}
#'
#'   Nakagawa, S., Johnson, P. C. D., & Schielzeth, H. (2017). The coefficient
#'   of determination R2 and intra-class correlation coefficient from
#'   generalized linear mixed-effects models revisited and expanded. *Journal of
#'   the Royal Society Interface*, *14*(134). \doi{10/gddpnq}
#'
#'   Shieh, G. (2008). Improved Shrinkage Estimation of Squared Multiple
#'   Correlation Coefficient and Squared Cross-Validity Coefficient.
#'   *Organizational Research Methods*, *11*(2), 387–407. \doi{10/bcwqf3}
#'
#'   Zheng, B., & Agresti, A. (2000). Summarizing the predictive power of a
#'   generalized linear model. *Statistics in Medicine*, *19*(13), 1771-1781.
#'   \doi{10/db7rfv}
#' @examples
#' # Pseudo R-squared for mixed models
#' R2(shipley.sem)  # fixed + random ('conditional')
#' R2(shipley.sem, re.form = ~ (1 | tree))  # fixed + 'tree'
#' R2(shipley.sem, re.form = ~ (1 | site))  # fixed + 'site'
#' R2(shipley.sem, re.form = NA)  # fixed only ('marginal')
#' R2(shipley.sem, re.form = NA, type = "spearman")  # using Spearman's Rho
#'
#' # Predicted R-squared: compare cross-validated predictions calculated/
#' # approximated via the hat matrix to standard method (leave-one-out)
#'
#' # Fit test models using Shipley data – compare lm vs glm
#' d <- na.omit(shipley)
#' m <- lm(Live ~ Date + DD + lat, d)
#' # m <- glm(Live ~ Date + DD + lat, binomial, d)
#'
#' # Manual CV predictions (leave-one-out)
#' cvf1 <- sapply(1:nrow(d), function(i) {
#'   m.ni <- update(m, data = d[-i, ])
#'   predict(m.ni, d[i, ], type = "response")
#' })
#'
#' # Short-cut via the hat matrix
#' y <- getY(m)
#' f <- fitted(m)
#' cvf2 <- y - (y - f) / (1 - hatvalues(m))
#'
#' # Compare predictions (not exactly equal for GLMs)
#' all.equal(cvf1, cvf2)
#' # lm: TRUE; glm: "Mean relative difference: 1.977725e-06"
#' cor(cvf1, cvf2)
#' # lm: 1; glm: 0.9999987
#'
#' # NOTE: comparison not tested here for mixed models, as hierarchical data can
#' # complicate the choice of an appropriate leave-one-out procedure. However,
#' # there is no obvious reason why use of the leverage values (diagonals of the
#' # hat matrix) to estimate CV predictions shouldn't generalise, roughly, to
#' # the mixed model case (at least for LMMs). In any case, users should
#' # exercise the appropriate caution in interpretation.
#' @export
R2 <- function(mod, data = NULL, adj = TRUE, pred = TRUE, offset = FALSE,
               re.form = NULL, type = c("pearson", "spearman"),
               adj.type = c("olkin-pratt", "ezekiel"), positive.only = TRUE,
               env = NULL) {

  m <- mod; d <- data; rf <- re.form; type <- match.arg(type);
  adj.type <- match.arg(adj.type)

  # Function
  R2 <- function(m) {

    # Refit model with any supplied data
    if (!is.null(d)) {
      m <- eval(update(m, data = d, evaluate = FALSE))
      env <- environment()
    }

    # Degrees of freedom
    n <- nobs(m)
    i <- attr(terms(m), "intercept")
    ndf <- n - i
    rdf <- df.residual(m)

    # No. predictors
    k <- ndf - rdf

    # R-squared
    R2 <- if (k > 0) {

      # Model link function
      f <- getFamily(m)
      lF <- f$linkfun
      lI <- f$linkinv

      # Model weights
      w <- weights(m)
      if (is.null(w)) w <- rep(1, n)
      w <- w[w > 0 & !is.na(w)]

      # Model offset
      d <- getData(m, subset = TRUE, env = env)
      o <- if (!offset) {
        mf <- suppressWarnings(model.frame(m, data = d))
        model.offset(mf)
      }
      if (is.null(o)) o <- 0

      # Response
      y <- getY(m, offset = offset, env = env)
      obs <- names(y)

      # Fitted values 
      # (need to supply data to predict() to avoid issues w/ re.form argument)
      f <- predict(m, d, re.form = rf)[obs]
      f <- lI(f - o)

      # Correlation
      yf <- cbind(y, f)
      if (type == "spearman") yf <- apply(yf, 2, rank)
      cm <- cov.wt(yf, w, cor = TRUE, center = as.logical(i))
      R <- cm$cor[1, 2]
      if (is.na(R)) R <- 0

      # Return R-squared
      if (positive.only) R <- max(R, 0)
      if (R < 0) -R^2 else R^2

    } else 0

    # Adjusted R-squared
    R2a <- if (adj) {

      if (k > 0) {

        # Variance unexplained (model 'error')
        e <- 1 - R2

        # Estimator (some code from altR2:::OPExactEstimator())
        R2a <- if (adj.type == "olkin-pratt") {
          1 - e * (ndf - 2) / rdf * gsl::hyperg_2F1(1, 1, (rdf + 2) / 2, e)
        } else {
          1 - e * ndf / rdf
        }

        # Return adjusted R-squared
        if (positive.only) R2a <- max(R2a, 0)
        R2a

      } else 0

    }

    # Predicted R-squared
    R2p <- if (pred) {

      if (!isGls(m) && !(isMer(m) && isGlm(m))) {

        if (k > 0) {

          # Leverage values (diagonals of the hat matrix, hii)
          hii <- hatvalues(m)[obs]
          s <- hii < 1

          # CV fitted values (response minus 'predicted' residuals)
          rp <- (y - f) / (1 - hii)
          f <- y - rp

          # Correlation
          yf[, 2] <- if (type == "spearman") rank(f) else f
          cm <- cov.wt(yf[s, ], w[s], cor = TRUE, center = as.logical(i))
          Rp <- cm$cor[1, 2]

          # Return predicted R-squared
          if (positive.only) Rp <- max(Rp, 0)
          if (Rp < 0) -Rp^2 else Rp^2

        } else 0

      } else NA

    }

    # Return values
    c("R.squared" = R2, "R.squared.adj" = R2a, "R.squared.pred" = R2p)

  }

  # Apply recursively
  rMapply(R2, m)

}


#' @title Weighted Average of Model Estimates
#' @description Calculate a weighted average of model estimates (e.g. effects,
#'   fitted values, residuals) for a set of models.
#' @param est A list or nested list of numeric vectors, comprising the model
#'   estimates. In the latter case, these should correspond to estimates for
#'   candidate models for each of a set of different response variables.
#' @param weights An optional numeric vector of weights to use for model
#'   averaging, or a named list of such vectors. The former should be supplied
#'   when `est` is a list, and the latter when it is a nested list (with
#'   matching list names). If `weights = "equal"` (default), a simple average is
#'   calculated instead.
#' @param est.names An optional vector of names used to extract and/or sort
#'   estimates from the output.
#' @details This function can be used to calculate a weighted average of model
#'   estimates such as effects, fitted values, or residuals, where models are
#'   typically competing candidate models fit to the same response variable.
#'   Weights are typically a 'weight of evidence' type metric such as Akaike
#'   model weights (Burnham & Anderson, 2002; Burnham et al., 2011), which can
#'   be conveniently calculated in *R* using packages such as
#'   [MuMIn](https://cran.r-project.org/package=MuMIn) or
#'   [AICcmodavg](https://cran.r-project.org/package=AICcmodavg). However,
#'   numeric weights of any sort can be used. If none are supplied, a simple
#'   average is calculated instead.
#'
#'   Averaging is performed via the 'full'/'zero' rather than
#'   'subset'/'conditional'/'natural' method, meaning that zero is substituted
#'   for estimates for any 'missing' parameters (e.g. effects) prior to
#'   calculations. This provides a form of shrinkage and thus reduces [estimate
#'   bias](https://stackoverflow.com/questions/53055050/predicted-values-with-mumin-throwing-error-when-full-false)
#'   (Burnham & Anderson, 2002; Grueber et al., 2011).
#' @return A numeric vector of the model-averaged estimates, or a list of such
#'   vectors.
#' @references Burnham, K. P., & Anderson, D. R. (2002). *Model Selection and
#'   Multimodel Inference: A Practical Information-Theoretic Approach* (2nd
#'   ed.). Springer-Verlag. <https://link.springer.com/book/10.1007/b97636>
#'
#'   Burnham, K. P., Anderson, D. R., & Huyvaert, K. P. (2011). AIC model
#'   selection and multimodel inference in behavioral ecology: some background,
#'   observations, and comparisons. *Behavioral Ecology and Sociobiology*,
#'   *65*(1), 23-35. \doi{10/c4mrns}
#'
#'   Dormann, C. F., Calabrese, J. M., Guillera‐Arroita, G., Matechou, E., Bahn,
#'   V., Bartoń, K., Beale, C. M., Ciuti, S., Elith, J., Gerstner, K., Guelat,
#'   J., Keil, P., Lahoz‐Monfort, J. J., Pollock, L. J., Reineking, B., Roberts,
#'   D. R., Schröder, B., Thuiller, W., Warton, D. I., … Hartig, F. (2018).
#'   Model averaging in ecology: A review of Bayesian, information-theoretic,
#'   and tactical approaches for predictive inference. *Ecological Monographs*,
#'   *88*(4), 485–504. \doi{10/gfgwrv}
#'
#'   Grueber, C. E., Nakagawa, S., Laws, R. J., & Jamieson, I. G. (2011).
#'   Multimodel inference in ecology and evolution: challenges and solutions.
#'   *Journal of Evolutionary Biology*, *24*(4), 699-711. \doi{10/b7b5d4}
#'
#'   Walker, J. A. (2019). Model-averaged regression coefficients have a
#'   straightforward interpretation using causal conditioning. *BioRxiv*,
#'   133785. \doi{10/c8zt}
#' @examples
#' # Model-averaged effects (coefficients)
#' m <- shipley.growth  # candidate models
#' e <- lapply(m, function(i) coef(summary(i))[, 1])
#' avgEst(e)
#'
#' # Using weights
#' w <- runif(length(e), 0, 1)
#' avgEst(e, w)
#'
#' # Model-averaged predictions
#' f <- lapply(m, predict)
#' head(avgEst(f, w))
#' @export
avgEst <-  function(est, weights = "equal", est.names = NULL) {

  e <- est; w <- weights; en <- est.names

  # Weights
  if (all(w == "equal")) {
    eqW <- function(x) {
      rep(1, length(x))
    }
    w <- if (any(sapply(e, isList))) lapply(e, eqW) else eqW(e)
  }

  # Function
  avgEst <- function(e, w) {

    # Sort names (for effects)
    en2 <- unique(unlist(lapply(e, names)))
    num <- suppressWarnings(as.numeric(en2))
    if (all(is.na(num))) {
      s1 <- isInt(en2)
      s2 <- isPhi(en2) | isR2(en2)
      en3 <- sort(en2[!s1 & !s2])
      en3 <- names(sort(sapply(en3, function(i) {
        lengths(regmatches(i, gregexpr(":", i)))
      })))
      en2 <- c(en2[s1], en3, en2[s2])
    }
    en <- if (!is.null(en)) en[en %in% en2] else en2

    # Combine estimates into table (missing = zero)
    e <- do.call(cbind, lapply(e, function(i) {
      sapply(en, function(j) {
        if (j %in% names(i)) i[[j]] else 0
      })
    }))

    # Weighted average
    apply(e, 1, weighted.mean, w)

  }

  # Apply recursively
  rMapply(avgEst, e, w, SIMPLIFY = FALSE)

}


#' @title Standardised Effects
#' @description Calculate fully standardised effects (model coefficients) in
#'   standard deviation units, adjusted for multicollinearity.
#' @param mod A fitted model object, or a list or nested list of such objects.
#' @param weights An optional numeric vector of weights to use for model
#'   averaging, or a named list of such vectors. The former should be supplied
#'   when `mod` is a list, and the latter when it is a nested list (with
#'   matching list names). If set to `"equal"`, a simple average is calculated
#'   instead.
#' @param data An optional dataset, used to first refit the model(s).
#' @param term.names An optional vector of names used to extract and/or sort
#'   effects from the output.
#' @param unique.eff Logical, whether unique effects should be calculated
#'   (adjusted for multicollinearity among predictors).
#' @param cen.x,cen.y Logical, whether effects should be calculated as if from
#'   mean-centred variables.
#' @param std.x,std.y Logical, whether effects should be scaled by the standard
#'   deviations of variables.
#' @param refit.x Logical, whether the model should be refit with mean-centred
#'   predictor variables (see Details).
#' @param incl.raw Logical, whether to append the raw (unstandardised) effects
#'   to the output.
#' @param R.squared Logical, whether R-squared values should also be calculated
#'   (via [R2()]).
#' @param R2.arg A named list of additional arguments to [R2()] (where
#'   applicable), excepting argument `env`. Ignored if `R.squared = FALSE`.
#' @param env Environment in which to look for model data (if none supplied).
#'   Defaults to the [formula()] environment.
#' @details `stdEff()` will calculate fully standardised effects (coefficients)
#'   in standard deviation units for a fitted model or list of models. It
#'   achieves this via adjusting the 'raw' model coefficients, so no
#'   standardisation of input variables is required beforehand. Users can simply
#'   specify the model with all variables in their original units and the
#'   function will do the rest. However, the user is free to scale and/or centre
#'   any input variables should they choose, which should not affect the outcome
#'   of standardisation (provided any scaling is by standard deviations). This
#'   may be desirable in some cases, such as to increase numerical stability
#'   during model fitting when variables are on widely different scales.
#'
#'   If arguments `cen.x` or `cen.y` are `TRUE`, effects will be calculated as
#'   if all predictors (x) and/or the response variable (y) were mean-centred
#'   prior to model-fitting (including any dummy variables arising from
#'   categorical predictors). Thus, for an ordinary linear model where centring
#'   of x and y is specified, the intercept will be zero – the mean (or weighted
#'   mean) of y. In addition, if `cen.x = TRUE` and there are interacting terms
#'   in the model, all effects for lower order terms of the interaction are
#'   adjusted using an expression which ensures that each main effect or lower
#'   order term is estimated at the mean values of the terms they interact with
#'   (zero in a 'centred' model) – typically improving the interpretation of
#'   effects. The expression used comprises a weighted sum of all the effects
#'   that contain the lower order term, with the weight for the term itself
#'   being zero and those for 'containing' terms being the product of the means
#'   of the other variables involved in that term (i.e. those not in the lower
#'   order term itself). For example, for a three-way interaction (x1 * x2 *
#'   x3), the expression for main effect \eqn{\beta1} would be:
#'
#'   \deqn{\beta_{1} + \beta_{12}\bar{x}_{2} + \beta_{13}\bar{x}_{3} +
#'   \beta_{123}\bar{x}_{2}\bar{x}_{3}}{\beta1 + (\beta12 * x̄2) + (\beta13 *
#'   x̄3) + (\beta123 * x̄2 * x̄3)} (adapted from
#'   [here](https://stats.stackexchange.com/questions/65898/why-could-centering-independent-variables-change-the-main-effects-with-moderatio))
#'
#'   In addition, if `std.x = TRUE` or `unique.eff = TRUE` (see below), product
#'   terms for interactive effects will be recalculated using mean-centred
#'   variables, to ensure that standard deviations and variance inflation
#'   factors (VIF) for predictors are calculated correctly (the model must be
#'   refit for this latter purpose, to recalculate the variance-covariance
#'   matrix).
#'
#'   If `std.x = TRUE`, effects are scaled by multiplying by the standard
#'   deviations of predictor variables (or terms), while if `std.y = TRUE` they
#'   are divided by the standard deviation of the response variable (minus any
#'   offsets). If the model is a GLM, this latter is calculated using the
#'   link-transformed response (or an estimate of same) generated using the
#'   function [glt()]. If both arguments are true, the effects are regarded as
#'   'fully' standardised in the traditional sense, often referred to as
#'   'betas'.
#'
#'   If `unique.eff = TRUE` (default), effects are adjusted for
#'   multicollinearity among predictors by dividing by the square root of the
#'   VIFs (Dudgeon, 2016; Thompson et al., 2017; [RVIF()]). If they have also
#'   been scaled by the standard deviations of x and y, this converts them to
#'   semipartial correlations, i.e. the correlation between the unique
#'   components of predictors (residualised on other predictors) and the
#'   response variable. This measure of effect size is arguably much more
#'   interpretable and useful than the traditional standardised coefficient, as
#'   it always represents the unique effects of predictors and so can more
#'   readily be compared both within and across models. Values range from zero
#'   to +/- one rather than +/- infinity (as in the case of betas) – putting
#'   them on the same scale as the bivariate correlation between predictor and
#'   response. In the case of GLMs however, the measure is analogous but not
#'   exactly equal to the semipartial correlation, so its values may not always
#'   be bound between +/- one (such cases are likely rare). Importantly, for
#'   ordinary linear models, the square of the semipartial correlation equals
#'   the increase in R-squared when that variable is included last in the model
#'   – directly linking the measure to unique variance explained. See
#'   [here](https://www.daviddisabato.com/blog/2016/4/8/on-effect-sizes-in-multiple-regression)
#'   for additional arguments in favour of the use of semipartial correlations.
#'
#'   If `refit.x`, `cen.x`, and `unique.eff` are `TRUE` and there are
#'   interaction terms in the model, the model will be refit with any
#'   (newly-)centred continuous predictors, in order to calculate correct VIFs
#'   from the variance-covariance matrix. However, refitting may not be
#'   necessary in some circumstances, for example where predictors have already
#'   been mean-centred, and whose values will not subsequently be resampled
#'   (e.g. parametric bootstrap). Setting `refit.x = FALSE` in such cases will
#'   save time, especially with larger/more complex models and/or bootstrap
#'   runs.
#'
#'   If `incl.raw = TRUE`, raw (unstandardised) effects can also be appended,
#'   i.e. those with all centring and scaling options set to `FALSE` (though
#'   still adjusted for multicollinearity, where applicable). These may be of
#'   interest in some cases, for example to compare their bootstrapped
#'   distributions with those of standardised effects.
#'
#'   If `R.squared = TRUE`, model R-squared values are appended to effects via
#'   the [R2()] function, with any additional arguments passed via `R2.arg`.
#'
#'   Finally, if `weights` are specified, the function calculates a weighted
#'   average of standardised effects across a set (or sets) of different
#'   candidate models for a particular response variable(s) (Burnham & Anderson,
#'   2002), via the [avgEst()] function.
#' @return A numeric vector of the standardised effects, or a list or nested
#'   list of such vectors.
#' @references Burnham, K. P., & Anderson, D. R. (2002). *Model Selection and
#'   Multimodel Inference: A Practical Information-Theoretic Approach* (2nd
#'   ed.). Springer-Verlag. <https://link.springer.com/book/10.1007/b97636>
#'
#'   Dudgeon, P. (2016). A Comparative Investigation of Confidence Intervals for
#'   Independent Variables in Linear Regression. *Multivariate Behavioral
#'   Research*, *51*(2-3), 139-153. \doi{10/gfww3f}
#'
#'   Thompson, C. G., Kim, R. S., Aloe, A. M., & Becker, B. J. (2017).
#'   Extracting the Variance Inflation Factor and Other Multicollinearity
#'   Diagnostics from Typical Regression Results. *Basic and Applied Social
#'   Psychology*, *39*(2), 81-90. \doi{10/gfww2w}
#' @examples
#' library(lme4)
#'
#' # Standardised (direct) effects for SEM
#' m <- shipley.sem
#' stdEff(m)
#' stdEff(m, cen.y = FALSE, std.y = FALSE)  # x-only
#' stdEff(m, std.x = FALSE, std.y = FALSE)  # centred only
#' stdEff(m, cen.x = FALSE, cen.y = FALSE)  # scaled only
#' stdEff(m, unique.eff = FALSE)  # include multicollinearity
#' stdEff(m, R.squared = TRUE)  # add R-squared
#' stdEff(m, incl.raw = TRUE)  # add unstandardised
#'
#' # Demonstrate equality with effects from manually-standardised variables
#' # (gaussian models only)
#' m <- shipley.growth[[3]]
#' d <- data.frame(scale(na.omit(shipley)))
#' e1 <- stdEff(m, unique.eff = FALSE)
#' e2 <- coef(summary(update(m, data = d)))[, 1]
#' stopifnot(all.equal(e1, e2))
#'
#' # Demonstrate equality with square root of increment in R-squared
#' # (ordinary linear models only)
#' m <- lm(Growth ~ Date + DD + lat, data = shipley)
#' r2 <- summary(m)$r.squared
#' e1 <- stdEff(m)[-1]
#' en <- names(e1)
#' e2 <- sapply(en, function(i) {
#'   f <- reformulate(en[!en %in% i])
#'   r2i <- summary(update(m, f))$r.squared
#'   sqrt(r2 - r2i)
#' })
#' stopifnot(all.equal(e1, e2))
#'
#' # Model-averaged standardised effects
#' m <- shipley.growth  # candidate models
#' w <- runif(length(m), 0, 1)  # weights
#' stdEff(m, w)
#' @export
stdEff <- function(mod, weights = NULL, data = NULL, term.names = NULL,
                   unique.eff = TRUE, cen.x = TRUE, cen.y = TRUE, std.x = TRUE,
                   std.y = TRUE, refit.x = TRUE, incl.raw = FALSE,
                   R.squared = FALSE, R2.arg = NULL, env = NULL) {

  m <- mod; w <- weights; d <- data; en <- term.names

  # Function
  stdEff <- function(m) {

    # Refit model with any supplied data
    if (!is.null(d)) {
      m <- eval(update(m, data = d, evaluate = FALSE))
      env <- environment()
    }

    # Model coefficients
    b <- if (isMer(m)) lme4::fixef(m, add.dropped = TRUE) else coef(m)
    bn <- names(b)
    int <- any(isInt(bn))  # intercept?

    # Effects (subset only relevant parameters)
    e <- na.omit(b[!isPhi(bn)])
    en <- names(e)[!isInt(names(e))]
    k <- length(en)

    # Model weights
    w <- weights(m)
    if (is.null(w)) w <- rep(1, nobs(m))
    w <- w[w > 0 & !is.na(w)]

    # Response
    y <- getY(m, env = env)
    obs <- names(y)

    # Centre/scale (x)
    if (k > 0) {

      # Predictors
      x <- getX(m, add.data = TRUE, as.df = TRUE, env = env)
      xm <- colMeans(x)
      inx <- any(isInx(en))  # interactions?

      # Adjust intercept (and possibly effects/predictors)
      if (cen.x) {

        # Model offset
        d <- getData(m, subset = TRUE, env = env)
        mf <- suppressWarnings(model.frame(m, data = d))
        o <- model.offset(mf)
        if (is.null(o)) o <- 0

        # Adjust intercept (set to mean of predicted y)
        if (int) {
          f <- predict(m, re.form = NA)[obs] - o
          e[1] <- weighted.mean(f, w)
        }

        # For interactions, adjust effects and centre predictors
        if (inx) {

          # Adjust lower-order terms
          # (ti = terms containing term i; ni = non-i components of ti)
          sI <- function(x) {
            unlist(strsplit(x, "(?<!:):(?!:)", perl = TRUE))
          }
          EN <- sapply(en, sI)
          e[en] <- sapply(en, function(i) {
            ei <- e[[i]]
            ii <- EN[[i]]
            ti <- en[en != i & sapply(en, function(j) all(ii %in% EN[[j]]))]
            if (length(ti) > 0) {
              ei + sum(sapply(ti, function(j) {
                jj <- EN[[j]]
                ni <- jj[!jj %in% ii]
                prod(e[j], xm[ni])
              }))
            } else ei
          })

          # Centre predictors (for correct SDs for product terms)
          if (std.x) {
            x <- getX(m, centre = TRUE, as.df = TRUE, env = env)
          }

        }

      }

      # Scale effects (x)
      if (std.x) {
        xs <- sapply(x[en], sdW, w)
        e[en] <- e[en] * xs
      }

      # Calculate unique effects of predictors (adjust for multicollinearity)
      if (unique.eff && k > 1) {

        # Refit model with centred predictors
        # (to calculate correct VIFs for interacting terms)
        m2 <- if (cen.x && inx && refit.x) {

          # Update dataset
          # (add response, weights, offset; set sum contrasts for factors)
          d2 <- data.frame(y, w, o, d)
          d2 <- data.frame(sapply(d2, function(i) {
            if (!is.numeric(i)) {
              i <- factor(i)
              contrasts(i) <- contr.sum(levels(i))
            }
            return(i)
          }))

          # Update term names (scale numeric predictors)
          XN <- lapply(labels(terms(m)), sI)
          xn <- unlist(lapply(XN, function(i) {
            paste(sapply(i, function(j) {
              if (is.numeric(mf[, j])) paste0("scale(", j, ")") else j
            }), collapse = ":")
          }))

          # Add random effects
          if (isMer(m)) {
            re <- lme4::findbars(formula(m))
            re <- sapply(re, function(i) paste0("(", deparse(i), ")"))
            xn <- c(xn, re)
          }

          # Rename any existing terms named "y", "w", "o"
          if (any(c("y", "w", "o") %in% names(d))) {
            xn <- sapply(xn, function(i) {
              i <- gsub("([^\\w.])", "~\\1~", i, perl = TRUE)
              i <- unlist(strsplit(i, "~"))
              s <- i %in% c("y", "w", "o")
              i[s] <- paste0(i[s], ".1")
              paste(i, collapse = "")
            })
          }

          # Refit model
          eval(update(m, reformulate(xn, "y", int), data = d2, weights = w,
                      offset = o, contrasts = NULL, evaluate = FALSE))

        } else m

        # Divide effects by RVIFs
        rvif <- na.omit(RVIF(m2, env = environment()))
        e[en] <- e[en] / rvif
        if (incl.raw) b[en] <- b[en] / rvif

      }

    }

    # Centre intercept (y)
    if (cen.y && int) {
      ym <- weighted.mean(y, w)
      if (isGlm(m)) {
        f <- getFamily(m)
        ym <- f$linkfun(ym)
      }
      e[1] <- e[1] - ym
    }

    # Scale effects (y)
    if (std.y) {
      if (isGlm(m)) y <- getY(m, link = TRUE, env = env)
      ys <- sdW(y, w)
      e <- e / ys
    }

    # Return standardised effects (optionally append R-squared/raw effects)
    e <- sapply(bn, function(i) {
      if (i %in% names(e)) e[[i]] else b[[i]]
    })
    if (R.squared) {
      r2 <- do.call(R2, c(list(m, env = env), R2.arg))
      names(r2) <- paste0("(", names(r2), ")")
      e <- c(e, r2)
    }
    if (incl.raw) {
      b <- c(b, e[isR2(names(e))])
      names(b) <- paste0("(raw)_", names(b))
      c(e, b)
    } else e

  }

  # Apply recursively
  e <- rMapply(stdEff, m, SIMPLIFY = FALSE)

  # Output effects or weighted average
  if (isList(e) && !is.null(w)) avgEst(e, w, en) else {
    if (!is.null(en)) {
      rMapply(function(i) {
        i[en[en %in% names(i)]]
      }, e, SIMPLIFY = FALSE)
    } else e
  }

}


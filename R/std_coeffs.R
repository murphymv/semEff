

#' @title Weighted Variance
#' @description Calculate the weighted variance of \code{x}.
#' @param x A numeric vector.
#' @param w A numeric vector of weights of the same length as \code{x}.
#' @param na.rm Logical, should NA's in \code{x} be removed?
#' @param ... Not used.
#' @details This function calculates the weighted variance of \code{x} via the
#'   weighted covariance matrix (\code{cov.wt}). If no weights are supplied, the
#'   simple variance is returned instead. As in \code{weighted.mean},
#'   \code{NA}'s in \code{w} are not handled specially and will return \code{NA}
#'   as result.
#' @return A numeric value, the weighted variance of \code{x}.
#' @seealso \code{\link[stats]{var}}, \code{\link[stats]{cov.wt}},
#'   \code{\link[stats]{weighted.mean}}
#' @examples
#' varW(1:10, w = 1:10)
#' varW(1:10)  # simple var
#'
#' ## NA handling
#' varW(c(1:9, NA), w = 1:10, na.rm = TRUE)  # NA in x (removed)
#' varW(c(1:9, NA), w = 1:10, na.rm = FALSE)  # NA in x (NA returned)
#' varW(1:10, w = c(1:9, NA))  # NA in w (NA returned)
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
sdW <- function(...) {
  sqrt(varW(...))
}


#' @title Model Term Names
#' @description Extract term names from a fitted model object as a character
#'   vector or list.
#' @param m A fitted model object of class \code{lm}, \code{glm}, or
#'   \code{merMod}, or a list or nested list of such objects.
#' @param intercept Logical, whether the intercept term should be included.
#' @param aliased Logical, whether names of aliased terms should be included
#'   (see Details).
#' @param list Logical, whether names should be returned as a (named) list.
#' @param ... Not used.
#' @details This function can be used to extract term names from a fitted model.
#'   Names of terms for which coefficients cannot be estimated are also included
#'   if \code{aliased = TRUE} (default). These may be terms which are perfectly
#'   correlated with other terms in the model, so that the model design matrix
#'   is rank deficient.
#' @return A character vector or named list of term names.
#' @examples
xNam <- function(m, intercept = TRUE, aliased = TRUE, list = FALSE, ...) {

  ## Function
  xN <- function(m) {

    ## All names as list
    xn <- labels(terms(m))
    mf <- model.frame(m)
    if (ncol(mf) > 1) {
      x1 <- mf[, 2]
      if (is.matrix(x1)) xn <- paste0(xn, colnames(x1))
    }
    xn2 <- rownames(summary(m)$coef)
    xn <- c(xn2[isInt(xn2)], xn)
    XN <- sapply(xn, function(i) {
      if (i %in% names(mf)) {
        xi <- mf[, i]
        xic <- class(xi)[1]
        if (xic == "factor") {
          paste0(i, levels(xi)[-1])
        } else {
          if (xic == "poly") {
            paste0(i, colnames(xi))
          } else i
        }
      } else i
    }, simplify = F)

    ## Drop intercept?
    if (!intercept) {
      XN <- XN[!isInt(names(XN))]
      if (length(XN) < 1)
        stop("No valid (non-intercept) terms present.")
    }

    ## Drop aliased terms?
    if (!aliased) {
      XN <- XN[sapply(XN, function(i) all(i %in% xn2))]
      if (length(XN) < 1)
        stop("No valid (non-aliased) terms present.")
    }

    ## Return as list?
    if (!list) unlist(unname(XN)) else XN

  }

  ## Apply recursively
  rMapply(xN, m)

}


#' @title Model Data
#' @description Extract the data used to fit a model.
#' @param m A fitted model object of class \code{lm}, \code{glm}, or
#'   \code{merMod}, or a list or nested list of such objects.
#' @param subset Logical. If \code{TRUE}, rows with missing observations
#'   (\code{NA}) are removed (if those variable(s) were used to fit model(s)).
#' @param merge Logical. If \code{TRUE}, and \code{m} is a list or nested list,
#'   all unique data columns are merged into one data frame, in the order in
#'   which they are encountered.
#' @param ... Arguments to \code{eval}.
#' @details This is a simple function to extract the data used to fit a model
#'   by evaluating the data slot of the model call object. If the model was fit
#'   using separate variables (i.e. vectors) instead of a data frame, data is
#'   returned via \code{model.frame} instead (with a warning). If \code{m} is a
#'   list and \code{merge = TRUE}, a single dataset containing all variables
#'   used to fit models is returned.
#' @return A data frame of the variables used to fit the model(s).
#' @seealso \code{\link[base]{eval}}, \code{\link[stats]{getCall}},
#'   \code{\link[stats]{model.frame}}
#' @examples
#'
getData <- function(m, subset = FALSE, merge = FALSE, ...) {

  ## Function
  gD <- function(m) {
    d <- eval(getCall(m)$data, ...)
    mf <- model.frame(m)
    if (is.null(d)) {
      d <- mf
      wo <- as.character(getCall(m))[c("weights", "offset")]
      names(d)[names(d) %in% c("(weights)", "(offset)")] <- wo
      names(d) <- gsub("^offset\\((.*)\\)$", "\\1", names(d))
      warning("Model frame used for data. Column names may not match original variable names.")
    }
    if (subset) {
      w <- weights(m)
      if (!is.null(w)) mf <- mf[w > 0, ]
      d[rownames(mf), ]
    } else d
  }

  ## Apply recursively
  rgD <- function(m) {
    if (isList(m)) {
      d <- lapply(m, rgD)
      if (merge) {
        d <- do.call(cbind, unname(d))
        d[unique(names(d))]
      } else d
    } else gD(m)
  }
  rgD(m)

}


#' @title Model Response Variable
#' @description Extract the response variable from a fitted model in original or
#'   link scale (for GLM's).
#' @param m A fitted model object of class \code{lm}, \code{glm}, or
#'   \code{merMod}, or alternatively a numeric vector corresponding to a
#'   variable to be transformed. \code{m} can also be a list or nested list of
#'   such objects.
#' @param family A character string specifying the GLM error family (e.g.
#'   "binomial", "poisson"). Ignored if \code{m} is a fitted model(s).
#' @param data An optional dataset used to first re-fit the model(s) (if
#'   \code{m} is a model object or list).
#' @param link Logical. If \code{TRUE}, return the (GLM) response variable on
#'   the link scale (see Details).
#' @param ... Not used.
#' @details This function will reconstruct the response variable from a model by
#'   summing the fitted values and the response residuals. If \code{link =
#'   TRUE}, the response is then transformed using the GLM link function, and,
#'   if any values are undefined, the 'working' response variable of the model
#'   is returned instead (see below). The function can also be used to transform
#'   a variable (supplied to \code{m}) using the default link function of the
#'   specified \code{family} - in which case the \code{link} argument is
#'   ignored.
#'
#'   \strong{Use of the Working Response Variable}
#'
#'   A key challenge in generating fully standardised model coefficients for a
#'   generalised linear model (GLM) with a non-gaussian link function is the
#'   (in)ability to calculate appropriate standardised ranges (typically the
#'   standard deviation) for the response variable in the link scale. This is
#'   because the directly-transformed response will often have undefined values.
#'   Although methods for circumventing this issue by indirectly estimating the
#'   standard deviation of the response in the link scale have previously been
#'   proposed - including a general variance-based method using a pseudo
#'   R-squared (Menard 2011) and another using a latent-theoretic framework for
#'   binomial models (Grace \emph{et al.} 2018) - here another approach is used.
#'   Where any transformed values are undefined, the function will instead
#'   return the synthetic 'working' response from the final iteration of the
#'   iteratively reweighted least squares (IRLS) algorithm of the GLM (McCullagh
#'   & Nelder 1989). This is reconstructed as the sum of the linear predictor
#'   and the working residuals. The advantage of this approach is that a
#'   relatively straightforward 'transformation' of any non-gaussian response is
#'   readily attainable in all cases. The standard deviation (or other relevant
#'   range) can then be calculated using values of the transformed response and
#'   used to scale the coefficients. An additional benefit for piecewise SEM's
#'   is that the transformed rather than original response can be used as a
#'   predictor in other models, ensuring that standardised indirect and total
#'   effects are calculated correctly (i.e. using a common standard deviation).
#'
#'   As the working response from a well-fitting model will be more 'accurate' -
#'   in the sense that the inverse-transformed values will more closely resemble
#'   the original response - the function always calculates it from a new GLM
#'   fit with the original response as the sole predictor (and using the same
#'   error family). Using the response as the predictor is a convenient way to
#'   ensure a near-optimal model fit, one which is limited only by the extent of
#'   non-linearity in the relationship between the response and its
#'   link-transformed expected values. The inverse-transformation of the working
#'   response from this model consistently produces an almost perfect
#'   correlation with the original response (see examples). This approach also
#'   ensures that the transformed values, and hence the standard deviation, are
#'   the same for any GLM fitted to the same response variable (provided it uses
#'   the same link function), and so facilitates parameter comparisons across
#'   different models.
#'
#'   The solution here of using the working response variable from the GLM to
#'   scale coefficients is primarily a practical one, but would generally be
#'   classed as an 'observed-empirical' approach to standardisation by Grace
#'   \emph{et al.} (2018) - as it utilises model error variance (i.e. the
#'   working residuals) rather than theoretical distribution-specific variance.
#' @return A numeric vector corresponding to the response variable in the link
#'   scale, or an array, list or nested list of such vectors.
#' @references Grace, J.B., Johnson, D.J., Lefcheck, J.S. and Byrnes, J.E.K.
#'   (2018) Quantifying relative importance: computing standardized effects in
#'   models with binary outcomes. \emph{Ecosphere} \strong{9}, e02283.
#'   \url{https://doi.org/10.1002/ECS2.2283}
#'
#'   McCullagh P. and Nelder, J. A. (1989) \emph{Generalized Linear Models} (2nd
#'   Edition). London: Chapman and Hall.
#'
#'   Menard, S. (2011) Standards for Standardized Logistic Regression
#'   Coefficients. \emph{Social Forces} \strong{89}, 1409-1428.
#'   \url{https://doi.org/10.1093/SF/89.4.1409}
#' @examples
getY <- function(m, family = NULL, data = NULL, link = FALSE, ...) {

  f <- family; d <- data

  ## Function
  gY <- function(m) {

    ## Model response
    mod <- isMod(m)
    y <- if (mod) {
      if (!is.null(d)) m <- update(m, data = d)
      y <- fitted(m) + resid(m, "response")
      w <- weights(m)
      if (!is.null(w)) y[w > 0] else y
    } else {
      y <- as.numeric(m)
      setNames(y, names(m))
    }

    ## Return in original or link scale
    if (isGlm(m) && link || !mod && !is.null(f)) {

      ## Transform to link scale
      yl <- if (mod) {
        f <- family(m)$family
        family(m)$linkfun(y)
      }

      ## Return the transformed (or working) response
      # if (!mod || any(is.infinite(yl))) {
        y2 <- y
        m <- glm(y ~ y2, f)
        if (!mod) yl <- family(m)$linkfun(y)
        # if (any(is.infinite(yl))) {
          getZ <- function(m) predict(m) + resid(m, "working")
          z <- getZ(m)
          repeat {
            zi <- family(m)$linkinv(z)
            r <- cor(zi, y)
            if (isTRUE(all.equal(r, 1))) return(z) else {
              z <- getZ(glm(y ~ z, f))
            }
          }
      #   } else yl
      # } else yl

    } else y

  }

  ## Apply recursively
  rMapply(gY, m)

}

# ## Test dataset
# set.seed(1)
# dat <- data.frame(
#   y1 = rnorm(100, mean = rnorm(1), sd = runif(1, 0, 10)),
#   y2 = rbinom(100, size = 1, prob = 0.5),
#   y3 = sample(c(rpois(90, lambda = 5), rep(0, 10))),
#   x1 = rep(rnorm(50, mean = rnorm(1), sd = runif(1, 0, 10)), each = 2),
#   x2 = rep(rnorm(50, mean = rnorm(1), sd = runif(1, 0, 10)), each = 2),
#   x3 = rep(rnorm(50, mean = rnorm(1), sd = runif(1, 0, 10)), each = 2),
#   ran = rep(sprintf("r%02d", 1:50), each = 2)
# )
#
# ## Fit models
# m1 <- lm(y1 ~ x1, dat)
# m2 <- glm(y2 ~ x1, dat, family = binomial)
# m3 <- lme4::glmer(y3 ~ x1 + (1 | ran), dat, family = poisson)
#
# ## Get link-transformed responses
# getY(m1)  # gaussian (unchanged)
# getY(m2, link = T)  # binomial
# getY(m3, link = T)  # poisson
#
# ## Correlation between original and inverse-transformed working responses (mean
# ## of 1000 samples)
# f <- function(y, family, n = 1000) {
#   y <- replicate(n, {
#     y <- sample(y, replace = T)
#     y2 <- y; m <- glm(y ~ y2, family)
#     y2 <- family(m)$linkinv(predict(m) + resid(m, "working"))
#     r <- cor(y, y2)
#     rd <- var(y - y2) / var(y)
#     c(cor = r, rel.diff = rd)
#   })
#   apply(y, 1, mean)
# }
# f(dat$y2, "binomial")
# f(dat$y3, "poisson")



#' @title Generalised Variance Inflation Factors
#' @description Calculate (generalised) variance inflation factors for terms in
#'   a fitted model via the variance-covariance matrix of coefficients.
#' @param m A fitted model object of class \code{lm}, \code{glm}, or
#'   \code{merMod}, or a list or nested list of such objects.
#' @param data An optional dataset used to first re-fit the model(s).
#' @param ... Not used.
#' @details This function calculates generalised variance inflation factors
#'   (GVIF) as described in Fox & Monette (1992), and also implemented in the
#'   \code{vif} function in the \code{car} package. However, whereas \code{vif}
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
#' @seealso \code{\link[car]{vif}}
#' @examples
VIF <- function(m, data = NULL, ...) {

  d <- data

  ## Function
  VIF <- function(m) {

    ## Term names
    if (!is.null(d)) m <- update(m, data = d)
    xn <- xNam(m, intercept = F, aliased = F)

    ## VIF's
    if (length(xn) > 1) {

      ## All term names (incl. aliased) as list
      XN <- xNam(m, intercept = F, list = T)

      ## Calculate VIF's from var-cov matrix
      V <- as.matrix(vcov(m))[xn, xn]
      R <- cov2cor(V)
      det.R <- det(R)
      vif <- lapply(XN, function(i) {
        df <- length(i)
        vif <- if (all(i %in% xn)) {
          j <- !xn %in% i
          Ri <- R[i, i, drop = F]
          Rj <- R[j, j, drop = F]
          vif <- det(Ri) * det(Rj) / det.R
          (vif^(1 / (2 * df)))^2
        } else NA
        vif <- rep(vif, df)
        setNames(vif, i)
      })
      unlist(unname(vif))

    } else
      stop("Model must have 2 or more valid terms.")

  }

  ## Apply recursively
  rMapply(VIF, m, SIMPLIFY = F)

}
# VIF(models.sem$PA)
# car::vif(models.sem$PA)
# VIF(update(test.glmm, . ~ . + rep(1, 72) + blah + blah2 + poly(y1, 2)))
# car::vif(update(test.glmm, . ~ . + rep(1, 72) + blah + blah2 + poly(y1, 2)))


#' @title R-squared/Pseudo R-squared
#' @description Calculate R-squared or pseudo R-squared for a fitted model as
#'   the squared multiple correlation between the observed and fitted values for
#'   the response variable. 'Adjusted' and 'predicted' R-squared values can
#'   optionally also be calculated (see Details).
#' @param m A fitted model object of class \code{lm}, \code{glm}, or
#'   \code{merMod}, or a list or nested list of such objects.
#' @param data An optional dataset used to first re-fit the model(s).
#' @param adj Logical, should adjusted R squared also be calculated?
#' @param pred Logical, should predictive R squared also be calculated?
#' @param re.form For mixed models of class \code{merMod}, the formula for
#'   random effects to condition on when generating model fitted values used in
#'   the calculation of R-squared. Defaults to \code{NA}, meaning no random
#'   effects are included. See \code{?lme4::predict.merMod} for further
#'   information.
#' @param ... Not used.
#' @details Various approaches to the calculation of a goodness-of-fit measure
#'   for GLM's analogous to R-squared in the ordinary linear model have
#'   been proposed. Generally termed 'pseudo R-squared' measures, they include
#'   variance-based, likelihood-based, and distribution-specific approaches.
#'   Here however, a straightforward OLS definition is used, which can be
#'   applied to any type of model for which predicted values of the response
#'   variable are generated: R-squared is calculated as the squared (weighted)
#'   correlation between the observed and predicted values of the response (in
#'   the original scale). This is simply the squared version of the correlation
#'   measure advocated by Zheng & Agresti (2000), itself an intuitive measure of
#'   goodness-of-fit describing the predictive power of a model. As the measure
#'   does not depend on any specific error distibution or model estimating
#'   procedure, it is also generally comparable across many different types of
#'   model (Kvalseth 1985). In the case of the ordinary linear model, the
#'   measure equals the more 'traditional' R-squared based on sums of squares.
#'
#'   If argument \code{adj} is set to \code{TRUE}, the adjusted R-squared value
#'   is also returned, which provides an estimate of the population (as opposed
#'   to sample) R-squared via an analytical formula which adjusts R-squared for
#'   the degrees of freedom of the model. Here this is calculated here via the
#'   'Pratt' rather than standard ('Ezekiel/Wherry') formula, shown in a
#'   simulation study to be the most effective of a range of existing formulas
#'   at estimating the population R-squared, across a broad range of model
#'   specification scenarios (Yin & Fan 2001).
#'
#'   If \code{pred = TRUE}, then a 'predictive' R-squared is also returned,
#'   which is calculated via the same formula as for R-squared but using
#'   'cross-validated' rather than standard model predictions. These are
#'   obtained by dividing model residuals by the complement of the observation
#'   leverage values (diagonals of the hat matrix), then subtracting these
#'   inflated residuals from the response variable. This is essentially a short
#'   cut to obtaining 'out-of-sample' predictions, normally arising via a
#'   leave-one-out cross validation procedure (in a GLM however they will not be
#'   exactly equal to such predictions). The resulting R-squared is an estimate
#'   of the R-squared that would occur were the model to be fitted to new data,
#'   and will typically be lower than the original R-squared - highlighting the
#'   degree of overfitting in the original model.
#'
#'   For mixed models, the function will, by default, calculate all R-squared
#'   values using fitted values for the fixed effects only, meaning that random
#'   effects are averaged over. This is the measure likely of primary interest
#'   to most researchers, and is equivalent to the 'marginal' R-squared
#'   statistic of Nakagawa \emph{et al.} (2017). In order to incorporate some or
#'   all of the random effects however, simply set the appropriate formula using
#'   the argument \code{re.form} (see \code{\link[lme4]{predict.merMod}}).
#'   Setting this to \code{NULL} includes all random effects, making the measure
#'   equivalent to the 'conditional' R-squared of Nakagawa \emph{et al.} (2017).
#'
#'   R-squared values produced by this function will always be bounded between
#'   zero (no fit) and one (perfect fit), meaning that any negative values
#'   arising from calculations are rounded up to zero. Negative values for
#'   R-squared typically mean that the fit is 'worse' than the null expectation
#'   of no relationship between the variables, which is difficult to interpret
#'   in practice and in any case usually only occurs in rare situations, such as
#'   where the intercept is suppressed. Hence, for simplicity and ease of
#'   interpretation, values <= 0 here are presented as a complete lack of model
#'   fit.
#'
#'   PLEASE NOTE: caution must be exercised in interpreting the values of any
#'   (pseudo) R-squared measure calculated for a GLM or mixed model (including
#'   those produced by this function), as such measures do not hold all the
#'   properties of R-squared in the ordinary linear model and as such may not
#'   always behave as expected. They are, at best, approximations. However, the
#'   value of a standardised measure of fit for model assessment and comparison
#'   may outweigh such reservations.
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
#'   \url{https://doi.org/10.1098/RSIF.2017.0213}
#'
#'   Yin, P. and Fan, X. (2001) Estimating R2 Shrinkage in Multiple Regression:
#'   A Comparison of Different Analytical Methods. \emph{The Journal of
#'   Experimental Education} \strong{69}(2), 203-224.
#'   \url{https://doi.org/10.1080/00220970109600656}
#'
#'   Zheng, B. and Agresti, A. (2000) Summarizing the predictive power of a
#'   generalized linear model. \emph{Statistics in Medicine} \strong{19}(13),
#'   1771-1781.
#'   \url{https://doi.org/10.1002/1097-0258(20000715)19:13<1771::AID-SIM485>3.0
#'   .CO;2-P}
#' @examples
#'
R2 <- function(m, data = NULL, adj = FALSE, pred = FALSE, re.form = NA,
               ...) {

  d <- data

  ## Function
  R2 <- function(m) {

    ## R squared
    if (!is.null(d)) m <- update(m, data = d)
    i <- attr(terms(m), "intercept")
    k <- nrow(summary(m)$coef) - i
    R2 <- if (k > 0) {
      y <- getY(m)
      n <- length(y)
      w <- weights(m)
      if (is.null(w)) w <- rep(1, n)
      s <- w > 0; w <- w[s]
      f <- predict(m, type = "response", re.form = re.form)[s]
      R <- cov.wt(cbind(y, f), w, cor = T)$cor[1, 2]
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
        Rp <- cov.wt(cbind(y, f)[s, ], w[s], cor = T)$cor[1, 2]
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
#' @param ... Not used.
#' @details This function can be used to calculate a weighted average of model
#'   estimates such as coefficients, fitted values, or residuals, where models
#'   are typically competing candidate models fit to the same response variable
#'   and using the same data. Weights are typically a 'weight of evidence' type
#'   metric such as Akaike model weights (Burnham & Anderson 2002, Burnham
#'   \emph{et al.} 2011), which can be calculated in \code{R} using packages
#'   such as \pkg{MuMIn} or \pkg{AICcmodavg}. However, weights of any sort can
#'   be used, provided they are numeric. If none are supplied, a standard
#'   average is calculated instead.
#'
#'   Averaging is performed via the 'full/zero' rather than
#'   'subset/conditional/natural' method, meaning that zero is substituted for
#'   any missing estimats (coefficients) prior to calculations - providing a
#'   form of 'shrinkage' for estimates (Burnham & Anderson 2002, Grueber
#'   \emph{et al.} 2011).
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
#'   \strong{65}(1), 23-35. \url{https://doi.org/10.1007/S00265-010-1029-6}
#'
#'   Dormann, C. F., Calabrese, J. M., Guillera-Arroita, G., Matechou, E., Bahn,
#'   V., Barton, K., ... Hartig, F. (2018). Model averaging in ecology: a review
#'   of Bayesian, information-theoretic, and tactical approaches for predictive
#'   inference. \emph{Ecological Monographs}, \strong{88}(4), 485-504.
#'   \url{https://doi.org/10.1002/ECM.1309}
#'
#'   Grueber, C. E., Nakagawa, S., Laws, R. J., & Jamieson, I. G. (2011).
#'   Multimodel inference in ecology and evolution: challenges and solutions:
#'   Multimodel inference. \emph{Journal of Evolutionary Biology}, \code{24}(4),
#'   699-711. \url{https://doi.org/10.1111/j.1420-9101.2010.02210.x}
#' @seealso \code{\link[stats]{weighted.mean}}
#' @examples
avgEst <-  function(e, weights = "equal", est.names = NULL, ...) {

  w <- weights; en <- est.names

  ## Weights
  if (all(w == "equal")) {
    f <- function(i) rep(1, length(i))
    w <- if (any(sapply(e, isList))) lapply(e, f) else f(e)
  }

  ## Function
  aE <- function(e, w) {

    ## Sort names (for coefs)
    en2 <- unique(unlist(lapply(e, names)))
    en <- if (is.null(en)) {
      i <- isInt(en2)
      r2 <- isR2(en2)
      en <- sort(en2[!i & !r2])
      en <- names(sort(sapply(en, function(i) {
        lengths(regmatches(i, gregexpr(":", i)))
      })))
      c(en2[i], en, en2[r2])
    } else en[en %in% en2]

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
  rMapply(aE, e, w, SIMPLIFY = F)

}


#' @title Standardised Coefficients
#' @description Calculate fully standardised model coefficients in standard
#'   deviation units, adjusted for multicollinearity among predictors.
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
#' @param cen.x Logical, whether the intercept and coefficients should be
#'   calculated as if from mean-centred predictor variables.
#' @param cen.y Logical, whether the intercept should be calculated as if the
#'   response variable was mean-centred.
#' @param std.x Logical, whether coefficients should be scaled by the standard
#'   deviations of the predictors.
#' @param std.y Logical, whether coefficients should be scaled by the standard
#'   deviation of the response.
#' @param unique.x Logical, whether coefficients should be adjusted for
#'   multicollinearity among predictors.
#' @param r.squared Logical, whether R-squared values should also be returned.
#' @param ... Arguments to function \code{R2}.
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
#'   model where centring of x and y is specified, the intercept will equal zero
#'   (the (weighted) mean of y). In addition, if \code{cen.x = TRUE} and there
#'   are interacting terms in the model, all coefficients for lower order terms
#'   of the interation are adjusted using an expression which ensures that each
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
#'   Mx3) + (\beta123 * Mx2 * Mx3)} (adapted from here:
#'   \url{http://bit.ly/2FOJPk8})
#'
#'   In addition, if \code{std.x = TRUE} or \code{unique.x = TRUE} (see below),
#'   product terms for interactive effects will be re-calculated using
#'   mean-centred main effects, to ensure that standard deviations and variance
#'   inflation factors for predictors are calculated correctly (the model is
#'   re-fit for this latter purpose, to recalculate the variance-covariance
#'   matrix).
#'
#'   If \code{std.x = TRUE}, coefficients are standardised by multiplying by the
#'   standard deviations of predictor variables (or terms), while if \code{std.y
#'   = TRUE} they are divided by the standard deviation of the response. If the
#'   model is a GLM, this latter is calculated from the link-transformed
#'   response (or an estimate of same) generated using the function \code{getY}.
#'   If both arguments are true, the coefficients are regarded as 'fully'
#'   standardised in the traditional sense, often referred to as 'betas'.
#'
#'   If \code{unique.x = TRUE}, coefficients are adjusted for multicollinearity
#'   among predictors by dividing by the square root of the variance inflation
#'   factors (Dudgeon 2016, Thompson \emph{et al.} 2017). If they are also
#'   standardised by the standard deviations of x and y, this converts them to
#'   semipartial correlations, i.e. the correlation between the unique
#'   components of predictors (residualised on other predictors) and the
#'   response variable. This measure of effect size is arguably much more
#'   interpretable and useful than the traditional standardised coefficient, as
#'   it is always estimated independent of the effects of multicollinearity, and
#'   so values can more readily be compared both within and across models. The
#'   effect size ranges from zero (no effect) to +/-1 (perfect relationship),
#'   rather than from zero to +/- infinity (as in the case of betas) - putting
#'   it on the same scale as the bivariate correlation between predictor and
#'   response. In the case of GLM's however, the measure is analagous but not
#'   exactly equal to the semipartial correlation, so its value may not be
#'   always be bound between +/-1 (such cases are likely rare). Crucially, for
#'   ordinary linear models, the square of the semipartial correlation equals
#'   the increase in R-squared when that variable is added last in the model -
#'   directly linking the measure to model fit and 'variance explained'. See
#'   \url{http://bit.ly/2GmyrMA} for additional arguments for the use of
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
#'  and Multimodel Inference: A Practical Information-Theoretic Approach} (2nd
#'  ed.). New York: Springer-Verlag. Retrieved from
#'  \url{https://www.springer.com/gb/book/9780387953649}
#'
#'  Dudgeon, P. (2016). A Comparative Investigation of Confidence Intervals for
#'  Independent Variables in Linear Regression. \emph{Multivariate Behavioral
#'  Research}, \strong{51}(2-3), 139-153.
#'  \url{https://doi.org/10.1080/00273171.2015.1121372}
#'
#'  Thompson, C. G., Kim, R. S., Aloe, A. M., & Becker, B. J. (2017). Extracting
#'  the Variance Inflation Factor and Other Multicollinearity Diagnostics from
#'  Typical Regression Results. \emph{Basic and Applied Social Psychology},
#'  \strong{39}(2), 81-90. \url{https://doi.org/10.1080/01973533.2016.1277529}
#' @seealso \code{\link[stats]{coef}}, \code{\link[semEff]{VIF}},
#'   \code{\link[semEff]{getY}}, \code{\link[semEff]{R2}},
#'   \code{\link[semEff]{avgEst}}
#' @examples
stdCoeff <- function(m, weights = NULL, data = NULL, term.names = NULL,
                     cen.x = TRUE, cen.y = TRUE, std.x = TRUE, std.y = TRUE,
                     unique.x = TRUE, r.squared = FALSE, ...) {

  w <- weights; d <- data; bn <- term.names

  ## Function
  sC <- function(m) {

    ## Coefficients
    if (!is.null(d)) m <- update(m, data = d)
    bn <- xNam(m)
    b <- summary(m)$coef
    b <- setNames(b[, 1], rownames(b))
    xn <- names(b)

    ## Intercept?
    i <- isInt(xn)
    int <- if (any(i)) {xn <- xn[!i]; T} else F

    ## Model weights
    n <- nobs(m)
    w <- weights(m)
    if (is.null(w)) w <- rep(1, n)
    s <- w > 0; w <- w[s]

    ## Centre/standardise x
    k <- length(xn)
    if (k > 0) {

      ## Predictors
      x <- model.matrix(m)[s, , drop = F]
      x <- data.frame(x, check.names = F)[xn]
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
            x <- data.frame(x, row.names = obs, check.names = F)
          }

        }

        ## Adjust intercept (weighted mean of predicted y)
        if (int) b[1] <- weighted.mean(predict(m)[s], w)

      }

      ## Standardise by x
      if (std.x) b[xn] <- b[xn] * sapply(x, sdW, w)

      ## Calculate unique effects of predictors (adjust for multicollinearity)
      if (unique.x && k > 1) {

        ## Re-fit model with centred predictors
        ## (to calculate correct VIF's for interactions)
        if (cen.x && inx) {

          ## Add centred predictors to data (list)
          if (is.null(d)) d <- getData(m)
          d <- c(list(x = as.matrix(x)), as.list(d[obs, ]))

          ## Add any offset(s)
          o <- model.offset(model.frame(m)[s, ])
          if (is.null(o)) o <- rep(0, n)
          d <- c(list(o = o), d)

          ## New model formula
          ran <- if (isMerMod(m)) {
            sapply(lme4::findbars(formula(m)), function(i) {
              paste0("(", deparse(i), ")")
            })
          }  # rand. effects? (http://bit.ly/2V1yDeu)
          f <- reformulate(c("x", ran), response = ".")

          ## Update model
          m <- update(m, f, data = d, offset = o)

        }

        ## Divide coefs by square root of VIF's
        b[xn] <- b[xn] / sqrt(VIF(m))

      }

    }

    ## Centre/standardise y
    if (cen.y && int) {
      ym <- weighted.mean(getY(m), w)
      b[1] <- b[1] - family(m)$linkfun(ym)
    }
    if (std.y) b <- b / sdW(getY(m, link = T), w)

    ## Return standardised coefficients
    sapply(bn, function(i) unname(b[i]))

  }

  ## Add R-squared?
  sC2 <- if (r.squared) {
    function(m) c(sC(m), R2(m, d, ...))
  } else sC

  ## Apply recursively
  b <- rMapply(sC2, m, SIMPLIFY = F)

  ## Output coefs or weighted average
  if (!is.null(w) && isList(b)) avgEst(b, w, bn)
  else {
    if (!is.null(bn)) {
      f <- function(i) i[bn[bn %in% names(i)]]
      rMapply(f, b, SIMPLIFY = F)
    } else b
  }

}

# stdCoeff(models.sem.top$BD$`12`)
# stdCoeff(m, cen.x = F, cen.y = F, std.x = F, std.y = F)[-1]
# coef(m)[-1] / sqrt(VIF(m))
# stdCoeff(test.lmm)
# stdCoeff(test.lmm)[-1]; stdCoeff(test.lmm, unique.x = F)[-1] / sqrt(VIF(update(test.lmm, . ~ scale(x1) * scale(x2) * scale(x3) + (1 | ran1))))


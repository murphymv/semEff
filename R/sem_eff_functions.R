

#' @title Automatic Calculation of SEM Effects
#' @description Calculate direct, indirect, and total effects for endogenous
#'   variables in a piecewise Structural Equation Model (SEM).
#' @param sem A piecewise SEM, comprising a named list of boot objects
#'   containing bootstrapped estimates from fitted models, or, alternatively, a
#'   named list/nested list of fitted model objects (class \code{lm},
#'   \code{glm}, or \code{lmerMod}).
#' @param predictors A character vector of names of predictor variables for
#'   which to calculate effects. If \code{NULL} (default), all predictor
#'   variables will be used.
#' @param mediators A character vector of names of mediating variables through
#'   which indirect effects will be calculated. If \code{NULL} (default), all
#'   endogenous predictor variables will be used.
#' @param responses A character vector of names of response variables for which
#'   to calculate effects. If \code{NULL} (default), all endogenous variables
#'   will be used.
#' @param ci.conf A numeric value specifying the confidence level for confidence
#'   intervals on effects.
#' @param ci.type The type of confidence interval to return. Can be one of
#'   \code{c("norm","basic", "stud", "perc", "bca")} (see
#'   \code{\link[boot]{boot.ci}}). Defaults to \code{"bca"}.
#' @param digits The number of significant digits to return for numeric values.
#' @param bci.arg A named list of additional arguments to \code{boot.ci}, which
#'   should not include argument \code{index}.
#' @param ... Arguments to \code{bootSEM}.
#' @details The eponymous function of this package calculates all direct,
#'   indirect, and total effects for endogenous variables in a 'piecewise'
#'   Structural Equation Model (SEM), that is, one where parameter estimation is
#'   local rather than global (Shipley 2000, 2009; Lefcheck 2016). The SEM
#'   simply takes the form of a list of fitted models (or bootstrapped model
#'   estimates), describing hypothesised causal pathways from predictors to
#'   response ('endogenous') variables, which are either direct or which operate
#'   indirectly via other response variables ('mediators'). This list should
#'   represent a directed ('acyclic') causal model, which should be named
#'   (exactly) for each response variable and ordered from 'upstream' or
#'   'causal' variables through to 'downstream' (i.e. those at the end of the
#'   pathway).
#'
#'   Direct effects are calculated as fully standardised model coefficients for
#'   each response variable, with indirect effects being the product of these
#'   direct effects operating along causal pathways in the SEM. The total
#'   effects of any given predictor on a response are then the sum of its direct
#'   and (all) its indirect effects. 'Mediator' effects are also calculated, as
#'   the sum of all indirect paths which operate through each individual
#'   mediator (useful to assess the relative importance of mediators). All of
#'   these effect types are calculated automatically for all (default) or a
#'   subset of predictors, mediators, or response variables in the SEM.
#'
#'   Confidence intervals for effects are returned for each response, with
#'   BC\emph{a} intervals calculated by default using bootstrapped estimates for
#'   all effect types (MacKinnon \emph{et al.} 2004, Cheung 2009, Hayes &
#'   Scharkow 2013). As indirect, total and mediator effects are not directly
#'   bootstrapped using the fitted models for response variables (i.e. via
#'   \code{bootSEM}), their equivalent 'bootstrapped' estimates are generated
#'   instead using each bootstrapped direct effect.
#'
#'   Correlated errors (and confidence intervals) are also returned if their
#'   boostrapped estimates are present in \code{sem}, or, if \code{sem} is a
#'   list of fitted models, if specified to argument \code{cor.err} (see
#'   \code{\link[semEff]{bootSEM}}). These represent relationships among
#'   response variables unaccounted for by the SEM.
#'
#'   All effects and bootstrapped effects are returned in lists for each
#'   response variable, with all except mediator effects also including the
#'   model intercept(s) - required for prediction (these will be zero for
#'   ordinary linear models with fully standardised effects).
#' @return A list, comprising for each response variable: 1) all effects
#'   (numeric vectors), 2) all bootstrapped effects (matrices), and 3) summary
#'   tables of effects and confidence intervals.
#' @references Cheung, M. W. L. (2009). Comparison of methods for constructing
#'   confidence intervals of standardized indirect effects. \emph{Behavior
#'   Research Methods}, \strong{41}(2), 425–438.
#'   \url{https://doi.org/10.3758/brm.41.2.425}
#'
#'   Hayes, A. F., & Scharkow, M. (2013). The Relative Trustworthiness of
#'   Inferential Tests of the Indirect Effect in Statistical Mediation Analysis:
#'   Does Method Really Matter? \emph{Psychological Science}, \strong{24}(10),
#'   1918–1927. \url{https://doi.org/10.1177/0956797613480187}
#'
#'   Lefcheck, J. S. (2016). piecewiseSEM: Piecewise structural equation
#'   modelling in \code{R} for ecology, evolution, and systematics.
#'   \emph{Methods in Ecology and Evolution}, \strong{7}(5), 573–579.
#'   \url{https://doi.org/10.1111/2041-210x.12512}
#'
#'   MacKinnon, D. P., Lockwood, C. M., & Williams, J. (2004). Confidence Limits
#'   for the Indirect Effect: Distribution of the Product and Resampling
#'   Methods. \emph{Multivariate Behavioral Research}, \strong{39}(1), 99.
#'   \url{https://doi.org/10.1207/s15327906mbr3901_4}
#'
#'   Shipley, B. (2000). A New Inferential Test for Path Models Based on
#'   Directed Acyclic Graphs. \emph{Structural Equation Modeling: A
#'   Multidisciplinary Journal}, \strong{7}(2), 206–218.
#'   \url{https://doi.org/10.1207/s15328007sem0702_4}
#'
#'   Shipley, B. (2009). Confirmatory path analysis in a generalized multilevel
#'   context. \emph{Ecology}, 90(2), 363–368.
#'   \url{https://doi.org/10.1890/08-1034.1}
#' @seealso \code{\link[semEff]{bootSEM}}, \code{\link[semEff]{bootCI}},
#' @examples
semEff <- function(sem, predictors = NULL, mediators = NULL, responses = NULL,
                   # ci = FALSE
                   ci.conf = 0.95, ci.type = "bca", digits = 3, bci.arg = NULL,
                   ...) {


  ### Prep

  p <- predictors; m <- mediators; r <- responses

  ## Bootstrap SEM (if necessary)
  if (!all(sapply(sem, isBoot))) sem <- bootSEM(sem, ...)

  ## Function to (recursively) replace parts of names/colnames in an object/list
  subNam <- function(p, r, x) {

    ## Function
    sN <- function(x) {
      if (is.character(x)) gsub(p, r, x) else {
        if (!is.null(names(x))) {
          names(x) <- gsub(p, r, names(x)); x
        } else {
          if (is.matrix(x)) {
            colnames(x) <- gsub(p, r, colnames(x)); x
          } else x
        }
      }
    }

    ## Apply recursively
    rsN <- function(x, sN) {
      x <- sN(x)
      if (isList(x) || isBoot(x)) {
        for (i in 1:length(x)) {
          x[[i]] <- rsN(x[[i]], sN)
        }; x
      } else x
    }
    rsN(x, sN)

  }

  ## Replace any periods in SEM variable names
  ## (periods are used to separate var names for indirect effects)
  sem <- subNam("[.]", "_", sem)

  ## Response names (default to all endogenous)
  ce <- grepl("~~", names(sem))  # T/F for corr. errors
  if (is.null(r)) r <- names(sem)[!ce]

  ## Mediator names (default to endogenous predictors)
  all.p <- unique(unlist(lapply(sem[r], function(i) names(i$t0))))
  all.p <- all.p[!isInt(all.p) & !isR2(all.p)]
  en.p <- r[r %in% all.p]
  if (is.null(m)) m <- en.p

  ## Predictor names (default to all predictors)
  if (is.null(p)) {
    ex <- all.p[!all.p %in% r]  # exogenous
    ex <- names(sort(sapply(ex, function(i) {
      lengths(regmatches(i, gregexpr(":", i)))
    })))
    p <- c(ex, en.p)
  }


  ### Calculate all direct, total indirect, total, and mediator effects

  ## Function to extract direct effects for predictors
  dirEff <- function(D) {
    sapply(p, function(i) {
      D <- if (is.matrix(D)) D[, colnames(D) == i] else unname(D[names(D) == i])
      if (length(D) > 0) D else 0
    }, simplify = F)
  }

  ## Function to calculate all indirect effects for predictors
  indEff <- function(D) {

    ## Function to extract names from estimates
    eNam <- function(x) {
      eN <- function(i) if (is.matrix(i)) colnames(i) else names(i)
      unique(unlist(rMapply(eN, x)))
    }

    ## Function to multiply effects to calculate indirect effects
    ## (for each endogenous predictor on a response, multiply its effect by all
    ## direct effects on that predictor)
    mEffs <- function(x) {

      ## Function
      xE <- function(x) {
        if (is.matrix(x)) {
          x <- x[, colnames(x) %in% m, drop = F]
          Map(function(i, j) {
            eb <- sem[[j]]$t
            i * eb[, colnames(eb) %in% all.p, drop = F]
          }, data.frame(x), colnames(x))
        } else {
          x <- x[names(x) %in% m]
          Map(function(i, j) {
            e <- sem[[j]]$t0
            i * e[names(e) %in% all.p]
          }, x, names(x))
        }
      }

      ## Apply recursively
      rMapply(xE, x, SIMPLIFY = F)

    }

    ## Function to collapse a nested list of vectors/matrices into a single
    ## vector/list of vectors
    unl <- function(x) {
      if (all(rapply(x, is.matrix))) {
        x <- rapply(x, function(i) {
          as.list(data.frame(i, check.names = F))
        }, how = "list")
        lapply(rapply(x, enquote), eval)
      } else unlist(x)
    }

    ## Calculate indirect effects: start with last response and move backwards,
    ## repeat for n = no. of responses
    lapply(D, function(i) {
      if (any(m %in% eNam(i))) {
        I <- list()
        I[[1]] <- mEffs(i)
        for (j in 1:length(i)) {
          I[[j + 1]] <- if (any(m %in% eNam(I[j]))) mEffs(I[j])
        }
        unl(I)  # collapse results
      } else NA
    })

  }

  ## Function to sum all indirect effects for predictors
  totInd <- function(I) {
    sapply(p, function(i) {
      I <- I[grepl(paste0("[.]", i, "$"), names(I))]
      if (length(I) > 0) Reduce("+", I) else 0
    }, simplify = F)
  }

  ## Function to calculate total effects (direct + total indirect effects)
  totEff <- function(D, I) {
    rMapply(function(i, j) i + j, D, I, SIMPLIFY = F)
  }

  ## Function to sum indirect effects for predictors operating through each mediator
  totIndM <- function(I) {
    sapply(m, function(i) {
      M <- paste(paste0("^", i, "[.]"), paste0("[.]", i, "[.]"),
                 sep = "|")
      I <- I[grepl(M, names(I))]
      if (length(I) > 0) {
        P <- paste(sapply(p, function(j) paste0("[.]", j, "$")),
                   collapse = "|")
        I <- I[grepl(P, names(I))]
        if (length(I) > 0) Reduce("+", I) else 0
      } else 0
    }, simplify = F)
  }

  ## Calculate and compile all effects for original estimates
  D <- lapply(sem[r], "[[", 1)
  I <- indEff(D)
  ED <- lapply(D, function(i) list("Direct" = dirEff(i)))
  EI <- lapply(I, function(i) list("Indirect" = totInd(i)))
  ET <- lapply(totEff(ED, EI), function(i) setNames(i, "Total"))
  EM <- lapply(I, function(i) list("Mediators" = totIndM(i)))
  E <- Map(c, ED, EI, ET, EM)

  ## Calculate and compile all effects for bootstrapped estimates
  DB <- lapply(sem[r], "[[", 2)
  IB <- indEff(DB)
  EDB <- lapply(DB, function(i) list("Direct" = dirEff(i)))
  EIB <- lapply(IB, function(i) list("Indirect" = totInd(i)))
  ETB <- lapply(totEff(EDB, EIB), function(i) setNames(i, "Total"))
  EMB <- lapply(IB, function(i) list("Mediators" = totIndM(i)))
  EB <- Map(c, EDB, EIB, ETB, EMB)


  ### Calculate bootstrapped confidence intervals for all effects

  ## Function to calculate CI's
  bCI <- function(e, eb, r) {

    ## Boot object for response var
    B <- sem[[r]]

    ## Change default CI type for parametric bootstrapping
    if (B$sim == "parametric" && ci.type == "bca") {
      message("Percentile confidence intervals used for parametric bootstrap samples.")
      ci.type <- "perc"
    }

    ## Calculate CI's using boot object
    if (!is.na(e)) {
      if (e != 0) {
        B$t0 <- e; B$t <- t(t(eb))
        ci <- do.call(boot::boot.ci, c(list(B, ci.conf, ci.type), bci.arg))
        tail(as.vector(ci[[4]]), 2)
      } else c(0, 0)
    } else c(NA, NA)

  }

  ## Apply to lists of effects/bootstrapped effects
  CI <- rMapply(bCI, E, EB, r, SIMPLIFY = F)


  ### Compile and output results

  ## Extract all effects into lists of vectors/matrices
  extE <- function(E) {
    sapply(r, function(i) {
      sapply(names(E[[i]]), function(j) {
        e <- E[[i]][[j]]
        is.B <- any(sapply(e, length) > 1)
        e <- if (is.B) {
          e[sapply(e, function(k) length(k) > 1)]
        } else unlist(e)[unlist(e) != 0]
        if (length(e) > 0) {
          if (is.B) e <- as.matrix(data.frame(e, check.names = F))
          if (j != "Mediators") {
            B <- sem[[i]]
            if (is.B) {
              I <- B$t[, isInt(colnames(B$t)), drop = F]
              eb <- cbind(I, e)
              attr(eb, "sim") <- B$sim
              eb
            } else {
              I <- B$t0[isInt(names(B$t0))]
              c(I, e)
            }
          } else e
        } else NA
      }, simplify = F)
    }, simplify = F)
  }
  e <- extE(E); eb <- extE(EB)

  ## Generate summary tables of effects and CI's
  E <- rMapply(c, E, CI, SIMPLIFY = F)
  s <- lapply(E, function(i) {
    lapply(i, function(j) {
      e <- data.frame(
        t(do.call(rbind, j)),
        row.names = c("Estimate", "Lower CI", "Upper CI"),
        check.names = F
      )
      e <- e[, e[1, ] != 0, drop = F]
      if (ncol(e) > 0) {
        e <- round(e, digits)
        stars <- t(data.frame(sapply(e, function(k) {
          if (!any(is.na(k))) {
            if (all(k > 0) || all(k < 0)) "*" else ""
          } else ""
        })))
        e <- format(e, nsmall = digits)
        rbind(e, " " = stars)
      } else NA
    })
  })

  ## Add table for correlated errors
  if (any(ce)) {
    CE <- bootCI(sem[ce], ci.conf, ci.type, digits, bci.arg)
    s <- c(s, list("Correlated Errors" = do.call(cbind, CE)))
  }

  ## Reinstate periods to all variable names (if they were present)
  e <- subNam("_", ".", e)
  eb <- subNam("_", ".", eb)
  s <- subNam("_", ".", s)

  ## Output
  out <- list("Effects" = e, "Boot. Effects" = eb, "Summary" = s)
  class(out) <- c("semEff", class(out))
  out


}


#' @title Print SEM Effects
#' @description A print method for an object of class "semEff", which returns
#'   summary tables of effects and confidence intervals for all responses.
#' @param x An object of class "semEff".
## S3 method for class 'semEff'
print.semEff <- function(x) print(x$Summary)


#' @title Extract SEM Effects
#' @description Extract direct, indirect, or total effects for a piecewise SEM.
#' @param x An object of class "semEff".
#' @param type The type of effects to return. Must be one of \code{c("orig",
#'   "boot")} (defaults to \code{"orig"}).
#' @param responses A vector of names of response variables for which to return
#'   effects. If \code{NULL} (default), all responses are used.
#' @details These are simple extractor functions for effects calculated using
#'   \code{semEff}, intended for convenience.
#' @return A list containing the original or bootstrapped effects for each
#'   response variable, as numeric vectors or matrices (respectively).
#' @name Effects
NULL
#' @describeIn Effects Extract all effects.
effects <- function(x, type = "orig", responses = NULL) {
  e <- sem.eff; t <- type; r <- responses
  e <- if (t == "orig") e[[1]] else {
    if (t == "boot") e[[2]] else
      stop("Effect type must be either 'orig' (default) or 'boot'")
  }
  if (!is.null(r)) e[r] else e
}
#' @describeIn Effects Extract direct effects.
direct <- function(...) {
  lapply(effects(...), "[[", 1)
}
#' @describeIn Effects Extract indirect effects.
indirect <- function(...) {
  lapply(effects(...), "[[", 2)
}
#' @describeIn Effects Extract total effects.
total <- function(...) {
  lapply(effects(...), "[[", 3)
}


#' @title SEM Effect Predictions
#' @description Generate predicted values for SEM direct, indirect, or total
#'   effects on a response variable.
#' @param m A fitted model object of class \code{lm}, \code{glm}, or
#'   \code{merMod}, or a list or nested list of such objects.
#' @param newdata An optional data frame of new values to predict from, which
#'   should contain all the variables used to fit \code{m}, or all those named
#'   in \code{effects}.
#' @param effects A numeric vector of effects to predict, or a list or nested
#'   list of such vectors.
#' @param eff.boot A matrix of bootstrapped effects used to calculate confidence
#'   intervals for predictions, or a list or nested list of such matrices.
#' @param type The type of prediction to return (for GLM's). Can be either
#'   \code{"link"} (default) or \code{"response"}.
#' @param interaction An optional name of an interactive effect, for which to
#'   return standardised effects for the main (continuous) variable across
#'   different values or levels of interacting variables (see Details). The name
#'   should be of the form \code{"a:b:..."}, containing all the variables
#'   involved and matching the name of an interactive effect present in model(s)
#'   coefficients or in \code{effects}.
#' @param seed Integer, seed for the random number generator (for reproducible
#'   results). This should be set to the same value originally set for
#'   boostrapping effects.
#' @param ci.conf A numeric value specifying the confidence level for confidence
#'   intervals on predictions (and optional interactive effects).
#' @param ci.type The type of confidence interval to return. Can be one of
#'   \code{c("norm","basic", "stud", "perc", "bca")} (see
#'   \code{\link[boot]{boot.ci}}). Defaults to \code{"bca"} (see Details).
#' @param digits The number of significant digits to return for numeric values
#'   (for interactive effects only).
#' @param bci.arg A named list of additional arguments to \code{boot.ci}, which
#'   should not include argument \code{index}.
#' @param parallel The type of parallel processing to use for calculating
#'   confidence intervals on predictions. Can be one of \code{"snow"},
#'   \code{"multicore"}, or \code{"no"} (for none - the default).
#' @param ncpus Integer, number of system cores to use for parallel processing.
#'   If \code{NULL} (default), all available cores are used.
#' @param cl Optional cluster to use if \code{parallel = "snow"}. If \code{NULL}
#'   (default), a local cluster is created using the specified number of cores.
#' @param ... Arguments to \code{stdCoeff}.
#' @details This function can be used to generate predicted values for SEM
#'   direct, indirect, or total effects on a response variable (arising via
#'   function \code{semEff}), which should be supplied to argument
#'   \code{effects}. The standard prediction formula is used, except that
#'   effects are used in place of model coefficients, with values for predictors
#'   drawn from the original model data (or from \code{newdata}). Fully
#'   standardised effects will be assumed, but if this is not the case then the
#'   same standardisation options as were applied originally (via
#'   \code{stdCoeff}) should be specified to the function. If no effects are
#'   supplied, standardised coefficients will be calculated and used to generate
#'   predictions (these will equal the model(s) fitted values if \code{unique.x
#'   = F} and \code{newdata = NULL}). In addition, model averaged predictions
#'   can be generated if \code{m} is a list and \code{weights} are supplied
#'   (passed to \code{stdCoeff}) - but only if no effects are supplied.
#'
#'   If bootstrapped effects are also supplied (to \code{eff.boot}),
#'   bootstrapped predictions are calculated by predicting from each effect.
#'   Confidence intervals can then be returned, for which the \code{type} should
#'   be appropriate for the original form of bootstrap sampling (defaults to
#'   \code{"bca"}). If the number of observations to predict is very large,
#'   parallel processing may speed up the calculation of intervals.
#'
#'   Predictions are always returned in the original (typically unstandardised)
#'   units of the (link-)response variable. For GLM's, they can be converted to
#'   the response scale if \code{type = "response"}.
#'
#'   If the name of an interactive effect is supplied to \code{interaction},
#'   standardised effects (and confidence intervals) can be returned for
#'   predictions of a continuous 'main' variable across specified values or
#'   levels of interacting variable(s). The values of the variables to predict
#'   must be supplied to \code{newdata}, with the continuous variable being
#'   identified automatically as having the most unique values.
#' @return A numeric vector of the predictions; or, if bootstrapped effects are
#'   supplied, a list containing the predictions and the upper and lower
#'   confidence intervals. Optional interactive effects may also be appended to
#'   the output.
#' @seealso \code{\link[semEff]{semEff}}, \code{\link[semEff]{stdCoeff}},
#'   \code{\link[semEff]{bootCI}}, \code{\link[semEff]{pSapply}}
#'   \code{\link[stats]{predict}}
#' @examples
predEff <- function(m, newdata = NULL, effects = NULL, eff.boot = NULL,
                    type = "link", interaction = NULL, seed = 1, ci.conf = 0.95,
                    ci.type = "bca", digits = 3, bci.arg = NULL,
                    parallel = "no", ncpus = NULL, cl = NULL, ...) {

  nd <- newdata; e <- effects; eb <- eff.boot; ix <- interaction
  p <- parallel; nc <- ncpus

  ## Arguments to stdCoeff
  a <- list(...)

  ## Weights (for model averaging)
  w <- eval(a$weights); a$weights <- NULL
  if (is.null(w) && isList(m)) {
    w <- rMapply(function(i) NULL, m, SIMPLIFY = F)
    if (is.null(e)) e <- w
    if (is.null(eb)) eb <- w
  }

  ## Coef standardisation options (for back-transforming predictions)
  cen.x <- if (!is.null(a$cen.x)) a$cen.x else T
  cen.y <- if (!is.null(a$cen.y)) a$cen.y else T
  std.x <- if (!is.null(a$std.x)) a$std.x else T
  std.y <- if (!is.null(a$std.y)) a$std.y else T

  ## Function
  pE <- function(m, w, e, eb) {

    ## Effects
    if (is.null(e)) e <- do.call(stdCoeff, c(list(m, w), a))
    e <- e[!is.na(e) & !isR2(names(e))]
    en <- names(e)
    EN <- sapply(en, function(i) unlist(strsplit(i, ":")))

    ## Predictors (x)
    d <- getData(m, subset = T, merge = T)
    obs <- rownames(d)
    x <- sapply(en, function(i) {
      if (!isInt(i)) {
        pt <- function(j) parse(text = j)
        if (isInx(i)) {
          xi <- sapply(pt(EN[[i]]), eval, d)
          if (cen.x) xi <- sweep(xi, 2, colMeans(xi))
          apply(xi, 1, prod)
        } else eval(pt(i), d)
      } else 1
    })
    x <- data.frame(x, row.names = obs, check.names = F)
    n <- nrow(x)

    ## Weights
    if (isList(m)) m <- m[[1]]
    w <- weights(m)
    w <- if (is.null(w)) rep(1, n) else w[w > 0]

    ## Offset(s)
    o <- model.offset(model.frame(m)[obs, ])
    if (is.null(o)) o <- rep(0, n)
    om <- if (cen.x) weighted.mean(o, w) else 0
    if (!is.null(nd)) {
      tt <- terms(m)
      on <- attr(tt, "offset")
      on <- if (!is.null(on)) {
        tn <- attr(tt, "variables")
        sapply(on, function(i) tn[[i + 1]])
      }
      on <- c(on, getCall(m)$offset)
      o <- if (!is.null(on)) {
        rowSums(sapply(on, eval, nd))
      } else rep(0, nrow(nd))
    }
    o <- o - om

    ## Predictor means/sd's
    xm <- sapply(x, function(i) if (cen.x) mean(i) else 0)
    xmw <- sapply(x, function(i) if (cen.x) weighted.mean(i, w) else 0)
    xs <- sapply(x, function(i) if (std.x) sdW(i, w) else 1)

    ## Response mean/sd (link)
    lf <- family(m)$linkfun
    ym <- if (cen.y) lf(weighted.mean(getY(m), w)) else 0
    ys <- if (std.y) sdW(getY(m, link = T), w) else 1

    ## Data to predict (standardise using original means/SD's)
    d <- if (is.null(nd)) x else nd
    obs <- rownames(d)
    d <- sapply(en, function(i) {
      if (!isInt(i)) {
        pt <- function(j) parse(text = j)
        di <- if (isInx(i)) {
          ENi <- EN[[i]]
          di <- sapply(pt(ENi), eval, d)
          di <- sweep(di, 2, xm[ENi])
          apply(di, 1, prod)
        } else eval(pt(i), d)
        (di - xmw[i]) / xs[i]
      } else 1
    })
    d <- data.frame(d, row.names = obs, check.names = F)

    ## Predictions
    f <- colSums(e * t(d))
    f <- f * ys + ym + o
    f <- setNames(f, obs)
    if (type == "response") {
      li <- family(m)$linkinv
      f <- li(f)
    }

    ## Add CI's
    if (!is.null(eb)) {

      sim <- attr(eb, "sim")
      if (sim == "parametric" && ci.type == "bca") {
        message("Percentile confidence intervals used for parametric bootstrap samples.")
        ci.type <- "perc"
      }

      ## Bootstrapped predictions
      R <- nrow(eb)
      eb <- eb[, en, drop = F]
      fb <- t(sapply(1:R, function(i) {
        ei <- eb[i, ][!is.na(eb[i, ])]
        di <- d[names(ei)]
        fi <- colSums(ei * t(di))
        fi * ys + ym + o
      }))
      if (type == "response") fb <- li(fb)

      ## CI's (use dummy boot object)
      set.seed(seed)
      B <- list(
        t0 = f, t = fb, R = R, data = x, seed = .Random.seed, sim = sim,
        stype = "i", strata = rep(1, n)
      )
      class(B) <- "boot"
      ci <- pSapply(1:nrow(d), function(i) {
        ci <- do.call(boot::boot.ci, c(list(B, ci.conf, ci.type, i), bci.arg))
        tail(as.vector(ci[[4]]), 2)
      }, p, nc, cl)
      colnames(ci) <- obs
      f <- list(fit = f, ci.lower = ci[1, ], ci.upper = ci[2, ])

    }

    ## Add effects for interactions
    if (!is.null(ix) && !is.null(nd)) {

      ## Interaction variable names
      ## (a = main var, b = interacting vars, ab = all interactions)
      ixn <- EN[[ix]]  # all var names
      a <- ixn[which.max(sapply(nd[ixn], function(i) length(unique(i))))]
      b <- ixn[!ixn %in% a]
      ab <- en[sapply(en, function(i) {
        ENi <- EN[[i]]
        a %in% ENi && length(ENi) %in% 2:length(ixn)
      })]

      ## Values for interacting variable(s) ('b')
      db <- sweep(unique(nd[b]), 2, xm[b])  # centred
      db <- lapply(EN[ab], function(i) {
        apply(db[i[i %in% b]], 1, prod)
      })

      ## Effects
      e <- e * ys / xs
      e <- e[a] + rowSums(mapply("*", e[ab], db))
      e <- e * xs[a] / ys
      names(e) <- paste(ix, 1:length(e), sep = "_")

      ## Add CI's
      if (!is.null(eb)) {

        ## Bootstrapped effects
        eb <- t(sapply(1:R, function(i) {
          ei <- eb[i, ] * ys / xs
          ei <- ei[a] + rowSums(mapply("*", ei[ab], db))
          ei * xs[a] / ys
        }))

        ## CI's
        B$t0 <- e; B$t <- eb
        e <- bootCI(B, ci.conf, ci.type, digits, bci.arg)
        c(f, list(interactions = e))

      } else {
        list(fit = f, interactions = round(e, digits))
      }

    } else f

  }

  ## Apply recursively
  rMapply(pE, m, w, e, eb, SIMPLIFY = F)

}

# system.time(pred.D <- predEff(models.sem, newdata.D, tot, tot.boot, type = "response"))

# predEff(models.sem$BD, newdata.DA, stdCoeff(models.sem.top$BD$`12`))[1:10]
# logit(pred.BD.DA$fit)[1:10]

# predEff(models.sem$BD, newdata.DA)
# predEff(models.sem$BD, newdata.DA, interaction = c("A:D"))
# predEff(models.sem$BD, newdata.DAH, interaction = c("A:D:H"))
# predEff(models.sem$BD, newdata.DA, dir$BD, interaction = c("A:D"))
# predEff(models.sem$BD, newdata.DA, dir$BD, dir.boot$BD, interaction = c("A:D"))
# calcSlopes(models.sem$BD, newdata.DA, pred.BD.DA)

# predEff(models.sem$PA)
# predict(models.sem$PA)
# predEff(models.sem$PA, unique.x = F)
# predEff(models.sem.top$PA, weights = models.sem.wts$PA)
# predEff(models.sem.top$PA, weights = models.sem.wts$PA, effects = stdCoeff(models.sem$PA))
# predEff(models.sem.top$PA)
# predEff(models.sem$PA, effects = dir$PA, eff.boot = dir.boot$PA)
# predEff(models.sem[1:2], effects = tot[1:2], eff.boot = tot.boot[1:2])
# predict(lm(ND ~ A + D, dataset2, offset = NSr))
# predEff(lm(ND ~ A + D, dataset2, offset = NSr), unique.x = F)
# predEff(lm(ND ~ A + D, dataset2, offset = NSr), cen.x = F, unique.x = F)
# predict(lm(ND ~ A + D, dataset2, offset = NSr, weights = NC))
# predEff(lm(ND ~ A + D, dataset2, offset = NSr, weights = NC), unique.x = F)
# predEff(lm(ND ~ A + D, dataset2, offset = NSr, weights = NC), cen.x = F, unique.x = F)
# predict(update(models.sem$ND, offset = NSr))
# predEff(update(models.sem$ND, offset = NSr), unique.x = F)
# predict(update(models.sem$PA, offset = log(dataset$host.abun)))
# predEff(update(models.sem$PA, offset = log(dataset$host.abun)), unique.x = F)
# predict(update(models.sem$NS, offset = log(dataset$host.abun)))
# predEff(update(models.sem$NS, offset = log(dataset$host.abun)), unique.x = F)
# predict(update(models.sem$NS, . ~ . + offset(log(dataset$para.cell)), offset = log(dataset$host.abun)))
# predEff(update(models.sem$NS, . ~ . + offset(log(dataset$para.cell)), offset = log(dataset$host.abun)), unique.x = F)
#
# predEff(glm(cbind(stem.dead, stem.live) ~ D * A * H, x.quasibinomial, dataset))


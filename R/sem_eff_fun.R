

#' @title SEM Effects
#' @description Automatically calculate direct, indirect, total, and mediator
#'   effects for endogenous (response) variables in a 'piecewise' structural
#'   equation model (SEM).
#' @param sem A piecewise SEM, comprising a list of fitted model objects, or,
#'   alternatively, of boot objects (class \code{"boot"}), containing
#'   bootstrapped model effects.
#' @param predictors,mediators,responses Names of variables for/through which to
#'   calculate effects. If \code{NULL} (default), all predictors, endogenous
#'   predictors (mediators), and endogenous variables (responses) will be used.
#' @param ci.conf A numeric value specifying the confidence level for confidence
#'   intervals on effects.
#' @param ci.type The type of confidence interval to return (defaults to
#'   \code{"bca"} - see Details). See \code{\link[boot]{boot.ci}} for further
#'   specification details.
#' @param digits The number of significant digits to return for numeric values.
#' @param bci.arg A named list of any additional arguments to \code{boot.ci},
#'   excepting argument \code{index}.
#' @param ... Arguments to \code{bootEff}.
#' @details The eponymous function of this package calculates all direct,
#'   indirect, total, and mediator effects for endogenous variables in a
#'   'piecewise' structural equation model (SEM), that is, one where parameter
#'   estimation is local rather than global (Shipley 2000, 2009; Lefcheck 2016).
#'   The SEM simply takes the form of a list of fitted models, or bootstrapped
#'   estimates from such models, describing hypothesised causal pathways from
#'   predictors to response ('endogenous') variables. These are either direct,
#'   or operate indirectly via other response variables ('mediators'). This list
#'   should represent a directed ('acyclic') causal model, which should be named
#'   (exactly) for each response variable and ordered from 'upstream' or
#'   'causal' variables through to 'downstream' (i.e. those at the end of the
#'   pathway). If \code{sem} is a list of fitted models, effects will first be
#'   bootstrapped using \code{bootEff} (this may take a while!).
#'
#'   Direct effects are calculated as fully standardised model coefficients for
#'   each response variable, while indirect effects are the product of these
#'   direct effects operating along causal pathways in the SEM. The total
#'   effects of any given predictor on a response are then the sum of its direct
#'   and (all) its indirect effects. 'Mediator' effects are also calculated, as
#'   the sum of all indirect paths which operate through each individual
#'   mediator - useful to assess the relative importance of different mediators
#'   in affecting the response. All of these effect types are calculated
#'   automatically for all (default) or a subset of predictors, mediators, or
#'   response variables in the SEM.
#'
#'   Confidence intervals for effects are returned for each response, with
#'   BC\emph{a} intervals calculated by default using bootstrapped estimates for
#'   each effect type (MacKinnon \emph{et al.} 2004, Cheung 2009, Hayes &
#'   Scharkow 2013). As indirect, total, and mediator effects are not directly
#'   bootstrapped using the fitted models for response variables (i.e. via
#'   \code{bootEff}), their equivalent 'bootstrapped' estimates are calculated
#'   instead using each bootstrapped direct effect.
#'
#'   Correlated errors (and confidence intervals) are also returned if their
#'   bootstrapped values are present in \code{sem}, or, if \code{sem} is a list
#'   of fitted models, if specified to argument \code{cor.err} (see
#'   \code{\link[semEff]{bootEff}}). These represent residual relationships
#'   among response variables, unaccounted for by the SEM.
#'
#'   All effects and bootstrapped effects are returned in lists for each
#'   response variable, with all except mediator effects also including the
#'   model intercept(s) - required for prediction (this will be zero for
#'   ordinary linear models with fully standardised effects).
#' @return A list object of class \code{"semEff"}, comprising: \enumerate{\item
#'   all effects \item all bootstrapped effects \item summary tables of effects
#'   and confidence intervals}
#' @references Cheung, M. W. L. (2009). Comparison of methods for constructing
#'   confidence intervals of standardized indirect effects. \emph{Behavior
#'   Research Methods}, \strong{41}(2), 425-438. \url{https://doi.org/fnx7xk}
#'
#'   Hayes, A. F., & Scharkow, M. (2013). The Relative Trustworthiness of
#'   Inferential Tests of the Indirect Effect in Statistical Mediation Analysis:
#'   Does Method Really Matter? \emph{Psychological Science}, \strong{24}(10),
#'   1918-1927. \url{https://doi.org/bbhr}
#'
#'   Lefcheck, J. S. (2016). piecewiseSEM: Piecewise structural equation
#'   modelling in \code{R} for ecology, evolution, and systematics.
#'   \emph{Methods in Ecology and Evolution}, \strong{7}(5), 573-579.
#'   \url{https://doi.org/f8s8rb}
#'
#'   MacKinnon, D. P., Lockwood, C. M., & Williams, J. (2004). Confidence Limits
#'   for the Indirect Effect: Distribution of the Product and Resampling
#'   Methods. \emph{Multivariate Behavioral Research}, \strong{39}(1), 99.
#'   \url{https://doi.org/chqcnx}
#'
#'   Shipley, B. (2000). A New Inferential Test for Path Models Based on
#'   Directed Acyclic Graphs. \emph{Structural Equation Modeling: A
#'   Multidisciplinary Journal}, \strong{7}(2), 206-218.
#'   \url{https://doi.org/cqm32d}
#'
#'   Shipley, B. (2009). Confirmatory path analysis in a generalized multilevel
#'   context. \emph{Ecology}, \strong{90}(2), 363-368.
#'   \url{https://doi.org/bqd43d}
#' @seealso \code{\link[semEff]{bootEff}}, \code{\link[semEff]{bootCI}}
#' @examples
#' ## SEM effects
#' (Shipley.SEM.Eff <- semEff(Shipley.SEM.Boot))
#'
#' ## Effects for selected variables
#' # semEff(Shipley.SEM.Boot, predictors = "lat")
#' # semEff(Shipley.SEM.Boot, mediators = "DD")
#' # semEff(Shipley.SEM.Boot, responses = "Live")
#'
#' ## Effects calculated using original SEM (models)
#' ## (not typically recommended - better to use saved boot objects)
#' # system.time(
#' #  Shipley.SEM.Eff <- semEff(Shipley.SEM, ran.eff = "site", seed = 53908)
#' # )
#' @export
semEff <- function(sem, predictors = NULL, mediators = NULL, responses = NULL,
                   ci.conf = 0.95, ci.type = "bca", digits = 3, bci.arg = NULL,
                   ...) {


  ### Prep

  p <- predictors; m <- mediators; r <- responses

  ## Bootstrap SEM (if necessary)
  if (!all(sapply(sem, isBoot))) sem <- bootEff(sem, ...)
  if (!isList(sem)) sem <- list(sem)

  ## Function to (recursively) replace parts of names/colnames in an object/list
  subNam <- function(p, r, x) {

    ## Function
    subNam <- function(x) {
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
    rsubNam <- function(x, sN) {
      x <- sN(x)
      if (isList(x) || isBoot(x)) {
        for (i in 1:length(x)) {
          x[[i]] <- rsubNam(x[[i]], sN)
        }; x
      } else x
    }
    rsubNam(x, subNam)

  }

  ## Replace any periods in SEM variable names
  ## (periods are used to separate var. names for indirect effects)
  sem <- subNam("[.]", "_", sem)

  ## Response names (default to endogenous vars.)
  ce <- grepl("~~", names(sem))
  if (is.null(r)) r <- names(sem)[!ce]

  ## Mediator names (default to endogenous predictors)
  all.p <- unique(unlist(lapply(sem[r], function(i) names(i$t0))))
  all.p <- all.p[!isInt(all.p) & !isR2(all.p)]
  en.p <- r[r %in% all.p]
  m <- if (length(en.p) > 0) {
    if (is.null(m)) en.p else m
  }

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
    }, simplify = FALSE)
  }

  ## Function to calculate all indirect effects for predictors
  indEff <- function(D) {

    ## Function to extract effect names
    effNam <- function(x) {
      effNam <- function(x) if (is.matrix(x)) colnames(x) else names(x)
      unique(unlist(rMapply(effNam, x)))
    }

    ## Function to multiply effects to calculate indirect effects
    ## (for each endogenous predictor on a response, multiply its effect by all
    ## direct effects on that predictor)
    multEff <- function(x) {

      ## Function
      multEff <- function(x) {
        if (is.matrix(x)) {
          x <- x[, colnames(x) %in% m, drop = FALSE]
          Map(function(i, j) {
            eb <- sem[[j]]$t
            i * eb[, colnames(eb) %in% all.p, drop = FALSE]
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
      rMapply(multEff, x, SIMPLIFY = FALSE)

    }

    ## Function to collapse a nested list of vectors/matrices into a single
    ## vector/list of vectors
    unlist2 <- function(x) {
      if (all(rapply(x, is.matrix))) {
        x <- rapply(x, function(i) {
          as.list(data.frame(i, check.names = FALSE))
        }, how = "list")
        lapply(rapply(x, enquote), eval)
      } else unlist(x)
    }

    ## Calculate indirect effects: start with last response and move backwards,
    ## repeat for n = no. of responses
    lapply(D, function(i) {
      if (any(m %in% effNam(i))) {
        I <- list()
        I[[1]] <- multEff(i)
        for (j in 1:length(i)) {
          I[[j + 1]] <- if (any(m %in% effNam(I[j]))) multEff(I[j])
        }
        unlist2(I)  # collapse results
      } else NA
    })

  }

  ## Function to sum all indirect effects for predictors
  totInd <- function(I) {
    sapply(p, function(i) {
      I <- I[grepl(paste0("[.]", i, "$"), names(I))]
      if (length(I) > 0) Reduce("+", I) else 0
    }, simplify = FALSE)
  }

  ## Function to calculate total effects (direct + total indirect)
  totEff <- function(D, I) {
    rMapply(function(i, j) i + j, D, I, SIMPLIFY = FALSE)
  }

  ## Function to sum indirect effects operating through each mediator
  totIndM <- function(I) {
    if (length(m) > 0) {
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
      }, simplify = FALSE)
    } else 0
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
  bootCI2 <- function(e, eb, r) {

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
        B$t0 <- e; B$t <- matrix(eb)
        ci <- do.call(
          boot::boot.ci, c(list(B, ci.conf, ci.type), bci.arg)
        )
        tail(as.vector(ci[[4]]), 2)
      } else c(0, 0)
    } else c(NA, NA)

  }

  ## Apply to lists of effects/bootstrapped effects
  CI <- rMapply(bootCI2, E, EB, r, SIMPLIFY = FALSE)


  ### Compile and output effects

  ## Extract all effects into lists of vectors/matrices
  extE <- function(E) {
    sapply(r, function(i) {
      sapply(names(E[[i]]), function(j) {
        e <- E[[i]][[j]]
        is.B <- any(sapply(e, length) > 1)
        e <- if (is.B) {
          e[sapply(e, function(k) length(k) > 1)]
        } else {
          unlist(e)[unlist(e) != 0]
        }
        if (length(e) > 0) {
          if (is.B) e <- as.matrix(data.frame(e, check.names = FALSE))
          if (j != "Mediators") {
            B <- sem[[i]]
            if (is.B) {
              I <- B$t[, isInt(colnames(B$t)), drop = FALSE]
              eb <- cbind(I, e)
              a <- c("sim", "seed", "n")
              attributes(eb)[a] <- attributes(B$t)[a]
              eb
            } else {
              I <- B$t0[isInt(names(B$t0))]
              c(I, e)
            }
          } else e
        } else NA
      }, simplify = FALSE)
    }, simplify = FALSE)
  }
  e <- extE(E); eb <- extE(EB)

  ## Generate summary tables of effects and CI's
  E <- rMapply(c, E, CI, SIMPLIFY = FALSE)
  s <- lapply(E, function(i) {
    lapply(i, function(j) {
      e <- data.frame(
        if (isList(j)) t(do.call(rbind, j)) else j,
        row.names = c("Estimate", "Lower CI", "Upper CI"),
        check.names = FALSE
      )
      e <- e[, e[1, ] != 0, drop = FALSE]
      if (ncol(e) > 0) {
        e <- round(e, digits)
        stars <- t(data.frame(sapply(e, function(k) {
          if (!any(is.na(k))) {
            if (all(k > 0) || all(k < 0)) "*" else ""
          } else ""
        })))
        e <- format(e, nsmall = digits)
        e <- rbind(e, " " = stars)
        attr(e, "ci.conf") <- ci.conf
        attr(e, "ci.type") <- ci.type
        e
      } else NA
    })
  })

  ## Add correlated errors
  if (any(ce)) {
    CE <- bootCI(sem[ce], ci.conf, ci.type, digits, bci.arg)
    s <- c(s, list("Correlated Errors" = do.call(cbind, CE)))
  }

  ## Reinstate periods to all variable names (if they were present)
  e <- subNam("_", ".", e)
  eb <- subNam("_", ".", eb)
  s <- subNam("_", ".", s)

  ## Output effects
  e <- list("Effects" = e, "Boot. Effects" = eb, "Summary" = s)
  class(e) <- c("semEff", "list")
  e


}


#' @title Print SEM Effects
#' @description A print method for an object of class \code{"semEff"}, returning
#'   summary tables of effects and confidence intervals for all responses.
#' @param x An object of class \code{"semEff"}.
#' @param ... Further arguments passed to or from other methods.
## S3 method for class 'semEff'
#' @export
print.semEff <- function(x, ...) print(x$Summary)


#' @title Get SEM Effects
#' @description Extract SEM direct, indirect, and/or total effects from an
#'   object of class \code{"semEff"}.
#' @param eff An object of class \code{"semEff"}.
#' @param type The type of effects to return. Must be either \code{"orig"}
#'   (default) or \code{"boot"}.
#' @param ... Arguments (above) to be passed to \code{getEff} from other
#'   extractor functions.
#' @details These are simple extractor functions for effects calculated using
#'   \code{semEff}, intended for convenience (e.g. for use with \code{predEff}).
#' @return A list containing the original or bootstrapped effects for each
#'   response variable, as numeric vectors or matrices (respectively).
#' @name getEff
NULL
#' @describeIn getEff Extract all effects.
#' @export
getEff <- function(eff, type = "orig") {
  if (type == "orig") eff[[1]] else {
    if (type == "boot") eff[[2]] else
      stop("Effect type must be either 'orig' (default) or 'boot'")
  }
}
#' @describeIn getEff Extract direct effects.
#' @export
dirEff <- function(...) {
  lapply(getEff(...), "[[", 1)
}
#' @describeIn getEff Extract indirect effects.
#' @export
indEff <- function(...) {
  lapply(getEff(...), "[[", 2)
}
#' @describeIn getEff Extract total effects.
#' @export
totEff <- function(...) {
  lapply(getEff(...), "[[", 3)
}


#' @title Predict Effects
#' @description Generate predicted values for SEM direct, indirect, or total
#'   effects.
#' @param mod A fitted model object, or a list or nested list of such objects.
#' @param newdata An optional data frame of new values to predict, which should
#'   contain all the variables named in \code{effects} or all those used to fit
#'   \code{mod}.
#' @param effects A numeric vector of effects to predict, or a list or nested
#'   list of such vectors. These will typically have been calculated using
#'   \code{semEff}, \code{bootEff}, or \code{stdCoeff}. Alternatively, a boot
#'   object produced by \code{bootEff} can be supplied.
#' @param eff.boot A matrix of bootstrapped effects used to calculate confidence
#'   intervals for predictions, or a list or nested list of such matrices. These
#'   will have been calculated using \code{semEff} or \code{bootEff}.
#' @param re.form For mixed models of class \code{"merMod"}, the formula for
#'   random effects to condition on when predicting effects. Defaults to
#'   \code{NA}, meaning random effects are averaged over. See
#'   \code{\link[lme4]{predict.merMod}} for further specification details.
#' @param type The type of prediction to return (for GLMs). Can be either
#'   \code{"link"} (default) or \code{"response"}.
#' @param ci.conf A numeric value specifying the confidence level for confidence
#'   intervals on predictions (and any interactive effects).
#' @param ci.type The type of confidence interval to return (defaults to
#'   \code{"bca"} - see Details). See \code{\link[boot]{boot.ci}} for further
#'   specification details.
#' @param interaction An optional name of an interactive effect, for which to
#'   return standardised effects for a 'main' continuous variable across
#'   different values or levels of interacting variables (see Details).
#' @param digits The number of significant digits to return for interactive
#'   effects.
#' @param bci.arg A named list of any additional arguments to \code{boot.ci},
#'   excepting argument \code{index}.
#' @param parallel The type of parallel processing to use for calculating
#'   confidence intervals on predictions. Can be one of \code{"snow"},
#'   \code{"multicore"}, or \code{"no"} (for none - the default).
#' @param ncpus Number of system cores to use for parallel processing. If
#'   \code{NULL} (default), all available cores are used.
#' @param cl Optional cluster to use if \code{parallel = "snow"}. If \code{NULL}
#'   (default), a local cluster is created using the specified number of cores.
#' @param ... Arguments to \code{stdCoeff}.
#' @details Generate predicted values for SEM direct, indirect, or total effects
#'   on a response variable, which should be supplied to \code{effects}. These
#'   are used in place of model coefficients in the standard prediction formula,
#'   with values for predictors drawn either from the data used to fit the
#'   original model(s) (\code{mod}) or from \code{newdata}. It is assumed that
#'   effects are fully standardised; however, if this is not the case, then the
#'   same standardisation options originally specified to \code{stdCoeff} should
#'   be re-specified - which will then be used to standardise the data. If no
#'   effects are supplied, standardised model coefficients will be calculated
#'   and used to generate predictions. These will equal the model(s) fitted
#'   values if \code{newdata = NULL}, \code{unique.x = FALSE}, and \code{re.form
#'   = NULL} (where applicable).
#'
#'   Model-averaged predictions can be generated if averaged \code{effects} are
#'   supplied to the model in \code{mod}, or, alternatively, if \code{weights}
#'   are specified (passed to \code{stdCoeff}) and \code{mod} is a list of
#'   candidate models (\code{effects} can also be passed using this latter
#'   method). For mixed model predictions where random effects are included
#'   (e.g. \code{re.form = NULL}), the latter approach should be used, otherwise
#'   the contribution of random effects will be taken from the single model
#'   instead of (correctly) being averaged over a candidate set.
#'
#'   If bootstrapped effects are supplied to \code{eff.boot} (or to
#'   \code{effects}, as part of a boot object), bootstrapped predictions are
#'   calculated by predicting from each effect. Confidence intervals can then be
#'   returned, for which the \code{type} should be appropriate for the original
#'   form of bootstrap sampling (defaults to \code{"bca"}). If the number of
#'   observations to predict is very large, parallel processing may speed up the
#'   calculation of intervals.
#'
#'   Predictions are always returned in the original (typically unstandardised)
#'   units of the (link-)response variable. For GLMs, they can be returned in
#'   the response scale if \code{type = "response"}.
#'
#'   Additionally, if the name of an interactive effect is supplied to
#'   \code{interaction}, standardised effects (and confidence intervals) can be
#'   returned for effects of a continuous 'main' variable across different
#'   values or levels of interacting variable(s). The name should be of the form
#'   \code{"x1:x2..."}, containing all the variables involved and matching the
#'   name of an interactive effect in the model(s) terms or in \code{effects}.
#'   The values for all variables should be supplied in \code{newdata}, with the
#'   continuous variable being automatically identified as having the most
#'   unique values.
#' @return A numeric vector of the predictions, or, if bootstrapped effects are
#'   supplied, a list containing the predictions and the upper and lower
#'   confidence intervals. Optional interactive effects may also be appended.
#'   Predictions may also be returned in a list or nested list, depending on the
#'   structure of \code{mod} (and other arguments).
#' @seealso \code{\link[stats]{predict}}, \code{\link[semEff]{semEff}},
#'   \code{\link[semEff]{stdCoeff}}, \code{\link[semEff]{bootCI}},
#'   \code{\link[semEff]{pSapply}}
#' @examples
#' ## Predict effects (direct, total)
#' m <- Shipley.SEM
#' e <- Shipley.SEM.Eff
#' dir <- dirEff(e); tot <- totEff(e)
#' f.dir <- predEff(m, effects = dir, type = "response")
#' f.tot <- predEff(m, effects = tot, type = "response")
#'
#' ## Using new data for predictors
#' d <- na.omit(Shipley)
#' xn <- c("lat", "DD", "Date", "Growth")
#' seq100 <- function(x) seq(min(x), max(x), length = 100)
#' nd <- data.frame(sapply(d[xn], seq100))
#' f.dir <- predEff(m, nd, dir, type = "response")
#' f.tot <- predEff(m, nd, tot, type = "response")
#' ## Add CI's (can take a while...)
#' # dir.b <- dirEff(e, "boot"); tot.b <- totEff(e, "boot")
#' # f.dir <- predEff(m, nd, dir, dir.b, type = "response")
#' # f.tot <- predEff(m, nd, tot, tot.b, type = "response")
#'
#' ## Predict an interactive effect (e.g. Live ~ Growth * DD)
#' xn <- c("Growth", "DD")
#' d[xn] <- scale(d[xn])  # standardise predictors (improves fit)
#' m <- lme4::glmer(Live ~ Growth * DD + (1 | site) + (1 | tree),
#'                  family = binomial, data = d)
#' nd <- with(d, expand.grid(
#'   Growth = seq100(Growth),
#'   DD = mean(DD) + c(-sd(DD), sd(DD))  # two levels for DD
#' ))
#' f <- predEff(m, nd, type = "response", interaction = "Growth:DD")
#' ## Add CI's (need to bootstrap model...)
#' # system.time(B <- bootEff(m, ran.eff = "site", R = 1000))
#' # f <- predEff(m, nd, B, type = "response", interaction = "Growth:DD")
#'
#' ## Model-averaged predictions (several approaches)
#' m <- Shipley.Growth  # candidate models (list)
#' w <- runif(length(m), 0, 1)  # weights
#' e <- stdCoeff(m, w)  # averaged effects
#' f1 <- predEff(m[[1]], effects = e)  # pass avg. effects
#' f2 <- predEff(m, weights = w)  # pass weights argument
#' f3 <- avgEst(predEff(m), w)  # use avgEst function
#' stopifnot(all.equal(f1, f2))
#' stopifnot(all.equal(f2, f3))
#'
#' ## Compare model fitted values: predEff vs. predict
#' m <- Shipley.SEM$Live
#' f1 <- predEff(m, unique.x = FALSE, re.form = NULL)
#' f2 <- predict(m)
#' stopifnot(all.equal(f1, f2))
#' @export
predEff <- function(mod, newdata = NULL, effects = NULL, eff.boot = NULL,
                    re.form = NA, type = "link", ci.conf = 0.95, ci.type = "bca",
                    interaction = NULL, digits = 3, bci.arg = NULL,
                    parallel = "no", ncpus = NULL, cl = NULL, ...) {

  m <- mod; nd <- newdata; e <- effects; eb <- eff.boot; rf <- re.form;
  ix <- interaction; p <- parallel; nc <- ncpus

  ## Arguments to stdCoeff
  a <- list(...)

  ## Weights (for model averaging)
  w <- eval(a$weights); a$weights <- NULL
  if (is.null(w) && isList(m)) {
    w <- rMapply(function(i) NULL, m, SIMPLIFY = FALSE)
    if (is.null(e)) e <- w
    if (is.null(eb)) eb <- w
  }

  ## Coef standardisation options (for back-transforming predictions)
  cen.x <- if (!is.null(a$cen.x)) a$cen.x else TRUE
  cen.y <- if (!is.null(a$cen.y)) a$cen.y else TRUE
  std.x <- if (!is.null(a$std.x)) a$std.x else TRUE
  std.y <- if (!is.null(a$std.y)) a$std.y else TRUE

  ## Function
  predEff <- function(m, w, e, eb) {

    ## Effects
    if (is.null(e)) e <- do.call(stdCoeff, c(list(m, w), a))
    if (isBoot(e)) {eb <- e$t; e <- e$t0}
    e <- e[!is.na(e) & !isR2(names(e))]

    ## Effect names
    en <- names(e)
    EN <- sapply(en, function(i) {
      unlist(strsplit(i, "(?<!:):(?!:)", perl = TRUE))
    })

    ## Model data
    d <- getData(m, subset = TRUE, merge = TRUE)
    obs <- rownames(d)

    ## Extract a single model (if list)
    m1 <- if (isList(m)) m[[1]] else m

    ## Random effects
    is.re <- isMer(m1) && !identical(rf, NA) && !identical(rf, ~ 0)
    re <- if (is.re) {
      pRE <- function(x) predict(x, nd, re.form = rf, random.only = TRUE)
      re <- rMapply(pRE, m, SIMPLIFY = FALSE)
      if (isList(re)) re <- avgEst(re, w)
      if (is.null(nd)) re[obs] else re
    } else 0

    ## Model weights
    w <- weights(m1)
    w <- if (is.null(w)) rep(1, nrow(d)) else w[w > 0]

    ## Offset(s)
    mf <- model.frame(m1, data = d)
    o <- model.offset(mf[obs, ])
    om <- if (cen.x && !is.null(o)) weighted.mean(o, w) else 0
    if (!is.null(nd)) {
      nd <- data.frame(nd)
      tt <- terms(m1)
      on <- attr(tt, "offset")
      on <- if (!is.null(on)) {
        tn <- attr(tt, "variables")
        sapply(on, function(i) tn[[i + 1]])
      }
      on <- c(on, getCall(m1)$offset)
      if (!is.null(on)) {
        o <- rowSums(sapply(on, eval, nd))
      }
    }
    if (is.null(o)) o <- 0
    o <- o - om

    ## Predictors
    dF <- function(...) data.frame(..., check.names = FALSE)
    eT <- function(x, d) {
      eT <- function(x) eval(parse(text = x), d)
      if (grepl("poly\\(.*[0-9]$", x)) {
        n <- nchar(x)
        xd <- eT(substr(x, 1, n - 1))
        xd[, substr(x, n, n)]
      } else eT(x)
    }
    x <- dF(model.matrix(reformulate(names(d)), data = d))

    # head(sapply(en[-1], eT, x))

    x <- dF(sapply(en, function(i) {
      if (!isInt(i)) {
        if (isInx(i)) {
          xi <- sapply(EN[[i]], eT, x)
          if (cen.x) xi <- sweep(xi, 2, colMeans(xi))
          apply(xi, 1, prod)
        } else eT(i, x)
      } else 1
    }))

    head(x)

    # ## Predictor means/SDs
    # xm <- sapply(x, function(i) if (cen.x) mean(i) else 0)
    # xmw <- sapply(x, function(i) if (cen.x) weighted.mean(i, w) else 0)
    # xs <- sapply(x, function(i) if (std.x) sdW(i, w) else 1)
    #
    # ## Response mean/SD (link scale)
    # lF <- family(m1)$linkfun
    # ym <- if (cen.y) lF(weighted.mean(getY(m1), w)) else 0
    # ys <- if (std.y) sdW(getY(m1, link = TRUE), w) else 1
    #
    # ## Data to predict (standardise using original means/SDs)
    # if (!is.null(nd)) {
    #   obs <- rownames(nd)
    #   x <- dF(model.matrix(reformulate(names(nd)), data = nd))
    # }
    # d <- dF(sapply(en, function(i) {
    #   if (!isInt(i)) {
    #     xi <- if (isInx(i)) {
    #       xi <- sapply(EN[[i]], eT, x)
    #       xi <- sweep(xi, 2, xm[colnames(xi)])
    #       apply(xi, 1, prod)
    #     } else eT(i, x)
    #     (xi - xmw[i]) / xs[i]
    #   } else 1
    # }), row.names = obs)
    #
    # ## Predictions
    # f <- colSums(e * t(d))
    # f <- f * ys + ym + o + re
    # f <- setNames(f, obs)
    # if (type == "response") {
    #   lI <- family(m1)$linkinv
    #   f <- lI(f)
    # }
    #
    # ## Add CI's
    # if (!is.null(eb)) {
    #
    #   ## Bootstrap details
    #   sim <- attr(eb, "sim")
    #   seed <- attr(eb, "seed")
    #   n <- attr(eb, "n")
    #   R <- nrow(eb)
    #
    #   ## Change default CI type for parametric bootstrapping
    #   if (sim == "parametric" && ci.type == "bca") {
    #     message("Percentile confidence intervals used for parametric bootstrap samples.")
    #     ci.type <- "perc"
    #   }
    #
    #   ## Bootstrapped predictions
    #   eb <- eb[, en, drop = FALSE]
    #   fb <- t(sapply(1:R, function(i) {
    #     ei <- eb[i, ][!is.na(eb[i, ])]
    #     di <- d[names(ei)]
    #     fi <- colSums(ei * t(di))
    #     fi * ys + ym + o + re
    #   }))
    #   if (nrow(fb) != R) fb <- t(fb)
    #   if (type == "response") fb <- lI(fb)
    #
    #   ## Create dummy boot object (for CI's)
    #   set.seed(seed)
    #   dd <- dF(rep(1, n))  # dummy data
    #   B <- list(
    #     t0 = f, t = fb, R = R, data = dd, seed = .Random.seed, sim = sim,
    #     stype = "i", strata = dd[, 1]
    #   )
    #   class(B) <- "boot"
    #   attr(B, "boot_type") <- "boot"
    #
    #   ## Calculate CI's
    #   ci <- as.matrix(pSapply(1:nrow(d), function(i) {
    #     ci <- do.call(
    #       boot::boot.ci, c(list(B, ci.conf, ci.type, i), bci.arg)
    #     )
    #     tail(as.vector(ci[[4]]), 2)
    #   }, p, nc, cl))
    #   colnames(ci) <- obs
    #   f <- list(fit = f, ci.lower = ci[1, ], ci.upper = ci[2, ])
    #
    # }
    #
    # ## Add interactive effects
    # if (isTRUE(isInx(ix) && ix %in% en && !is.null(nd))) {
    #
    #   ## Names of variables involved in interaction
    #   ## (ab = all, a = main, b = interacting, a.b = interaction(s))
    #   ab <- EN[[ix]]; n <- length(ab)
    #   a <- ab[which.max(sapply(x[ab], function(i) length(unique(i))))]
    #   b <- ab[!ab %in% a]
    #   a.b <- if (n > 2) {
    #     a.b <- unlist(lapply(2:n, function(i) {
    #       combn(ab, i, paste, collapse = ":")
    #     }))
    #     a.b[sapply(a.b, function(i) a %in% EN[[i]])]
    #   } else ix
    #
    #   ## Values for interacting variable(s) (b)
    #   xb <- unique(dF(sapply(b, eT, x)))
    #   xb <- sweep(xb, 2, xm[b])
    #   xb <- lapply(EN[a.b], function(i) {
    #     apply(xb[i[i %in% b]], 1, prod)
    #   })
    #
    #   ## Effects
    #   e <- e * ys / xs
    #   e <- e[a] + rowSums(mapply("*", e[a.b], xb))
    #   e <- e * xs[a] / ys
    #   names(e) <- paste(ix, 1:length(e), sep = "_")
    #
    #   ## Add CI's
    #   f <- if (!is.null(eb)) {
    #
    #     ## Bootstrapped effects
    #     eb <- t(sapply(1:R, function(i) {
    #       ei <- eb[i, ] * ys / xs
    #       ei <- ei[a] + rowSums(mapply("*", ei[a.b], xb))
    #       ei * xs[a] / ys
    #     }))
    #     if (nrow(eb) != R) eb <- t(eb)
    #
    #     ## CI's
    #     B$t0 <- e; B$t <- eb
    #     e <- bootCI(B, ci.conf, ci.type, digits, bci.arg)
    #     c(f, list(interactions = e))
    #
    #   } else {
    #     list(fit = f, interactions = round(e, digits))
    #   }
    #
    # }
    #
    # ## Output
    # if (!is.null(eb)) set.seed(NULL); f

  }

  ## Apply recursively
  rMapply(predEff, m, w, e, eb, SIMPLIFY = FALSE)

}


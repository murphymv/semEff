

#' @title SEM Effects
#' @description Automatically calculate direct, indirect, total, and mediator
#'   effects for endogenous (response) variables in a 'piecewise' structural
#'   equation model (SEM).
#' @param sem A piecewise SEM, comprising a list of fitted model objects, or,
#'   alternatively, of boot objects (class `"boot"`), containing bootstrapped
#'   model effects.
#' @param predictors,mediators,responses Names of variables for/through which to
#'   calculate effects. If `NULL` (default), all predictors, endogenous
#'   predictors (mediators), and endogenous variables (responses) will be used.
#' @param use.raw Logical, whether to use 'raw' (unstandardised) effects for all
#'   calculations (if present in `sem` object).
#' @param ci.conf A numeric value specifying the confidence level for confidence
#'   intervals on effects.
#' @param ci.type The type of confidence interval to return (defaults to `"bca"`
#'   — see Details). See [boot.ci()] for further specification details.
#' @param digits The number of significant digits to return for numeric values
#'   in summary tables.
#' @param bci.arg A named list of any additional arguments to [boot.ci()],
#'   excepting argument `index`.
#' @param ... Arguments to [bootEff()].
#' @details The eponymous function of this package calculates all direct,
#'   indirect, total, and mediator effects for a 'piecewise' structural equation
#'   model (SEM), that is, one where parameter estimation is local rather than
#'   global (Shipley 2000, 2009; Lefcheck 2016). The SEM simply takes the form
#'   of a list of fitted models, or bootstrapped estimates from such models,
#'   describing hypothesised causal pathways from predictors to response
#'   ('endogenous') variables. These are either direct, or operate indirectly
#'   via other response variables ('mediators'). This list should represent a
#'   directed ('acyclic') causal model, which should be named (exactly) for each
#'   response variable and ordered from 'upstream' or 'causal' variables through
#'   to 'downstream' (i.e. those at the end of the pathway). If `sem` is a list
#'   of fitted models, effects will first be bootstrapped using [bootEff()]
#'   (this may take a while!).
#'
#'   Direct effects are calculated as fully standardised model coefficients for
#'   each response variable, while indirect effects are the product of these
#'   direct effects operating along causal pathways in the SEM. The total
#'   effects of any given predictor on a response are then the sum of its direct
#'   and (all) its indirect effects. 'Mediator' effects are also calculated, as
#'   the sum of all indirect paths which operate through each individual
#'   mediator — useful to assess the relative importance of different mediators
#'   in affecting the response. All of these effect types are calculated
#'   automatically for all (default) or a subset of predictors, mediators, or
#'   response variables in the SEM. As indirect, total, and mediator effects are
#'   not directly bootstrapped using the fitted models for response variables
#'   (i.e. via [bootEff()]), their equivalent 'bootstrapped' estimates are
#'   calculated instead using each bootstrapped direct effect.
#'
#'   Confidence intervals for all effects are returned in summary tables for
#'   each response (similarly to [bootCI()]), with BC*a* intervals calculated by
#'   default using the bootstrapped estimates for each effect type (MacKinnon
#'   *et al.* 2004, Cheung 2009, Hayes & Scharkow 2013). Effects for which the
#'   confidence intervals do not contain zero are highlighted with a star (i.e.
#'   'significant' at the `ci.conf` level). Bootstrap standard errors (standard
#'   deviations of the samples) and biases (sample means minus original
#'   estimates) are also included for reference.
#'
#'   Correlated errors (and confidence intervals) are also returned if their
#'   bootstrapped values are present in `sem`, or, if `sem` is a list of fitted
#'   models, if specified to argument `cor.err` (see [bootEff()]). These
#'   represent residual relationships among response variables, unaccounted for
#'   by the hypothesised SEM paths.
#'
#'   All calculated effects and bootstrapped effects are returned in lists for
#'   each response variable, with all except mediator effects also including the
#'   model intercept(s) — required for prediction (this will be zero for
#'   ordinary linear models with fully standardised effects).
#' @return A list object of class `"semEff"`, comprising: 1. All calculated
#'   effects 2. All calculated bootstrapped effects 3. Summary tables of effects
#'   and confidence intervals
#' @references Cheung, M. W. L. (2009). Comparison of methods for constructing
#'   confidence intervals of standardized indirect effects. *Behavior Research
#'   Methods*, **41**(2), 425-438. <https://doi.org/fnx7xk>
#'
#'   Hayes, A. F., & Scharkow, M. (2013). The Relative Trustworthiness of
#'   Inferential Tests of the Indirect Effect in Statistical Mediation Analysis:
#'   Does Method Really Matter? *Psychological Science*, **24**(10), 1918-1927.
#'   <https://doi.org/bbhr>
#'
#'   Lefcheck, J. S. (2016). piecewiseSEM: Piecewise structural equation
#'   modelling in `R` for ecology, evolution, and systematics. *Methods in
#'   Ecology and Evolution*, **7**(5), 573-579. <https://doi.org/f8s8rb>
#'
#'   MacKinnon, D. P., Lockwood, C. M., & Williams, J. (2004). Confidence Limits
#'   for the Indirect Effect: Distribution of the Product and Resampling
#'   Methods. *Multivariate Behavioral Research*, **39**(1), 99.
#'   <https://doi.org/chqcnx>
#'
#'   Shipley, B. (2000). A New Inferential Test for Path Models Based on
#'   Directed Acyclic Graphs. *Structural Equation Modeling: A Multidisciplinary
#'   Journal*, **7**(2), 206-218. <https://doi.org/cqm32d>
#'
#'   Shipley, B. (2009). Confirmatory path analysis in a generalized multilevel
#'   context. *Ecology*, **90**(2), 363-368. <https://doi.org/bqd43d>
#' @examples
#' # SEM effects
#' (Shipley.SEM.Eff <- semEff(Shipley.SEM.Boot))
#'
#' # Effects for selected variables
#' # semEff(Shipley.SEM.Boot, predictors = "lat")
#' # semEff(Shipley.SEM.Boot, mediators = "DD")
#' # semEff(Shipley.SEM.Boot, responses = "Live")
#'
#' # Effects calculated using original SEM (models)
#' # (not typically recommended — better to use saved boot objects)
#' # system.time(
#' #  Shipley.SEM.Eff <- semEff(Shipley.SEM, R = 10000, seed = 53908,
#' #                            ran.eff = "site")
#' # )
#' @export
semEff <- function(sem, predictors = NULL, mediators = NULL, responses = NULL,
                   use.raw = FALSE, ci.conf = 0.95, ci.type = "bca", digits = 3,
                   bci.arg = NULL, ...) {


  # Prep

  p <- predictors; m <- mediators; r <- responses

  # Bootstrap SEM (if necessary)
  is.B <- all(sapply(sem, isBoot))
  if (!is.B) sem <- bootEff(sem, ...)
  if (!isList(sem)) sem <- list(sem)

  # Use raw (unstandardised) effects (if present)?
  if (use.raw) {
    sem <- lapply(sem, function(i) {
      r <- isRaw(names(i$t0))
      if (any(r)) {
        i$t0[!r] <- i$t0[r]
        i$t[, !r] <- i$t[, r]
      }
      i
    })
  }

  # Function to (recursively) replace parts of names/colnames in object/list
  subNam <- function(p, r, x, ...) {

    # Function
    subNam <- function(x) {
      if (is.character(x)) {
        x <- gsub(p, r, x, ...)
      } else {
        if (!is.null(names(x))) {
          names(x) <- gsub(p, r, names(x), ...)
        } else if (is.matrix(x)) {
          colnames(x) <- gsub(p, r, colnames(x), ...)
        }
      }
      x
    }

    # Apply recursively
    rsubNam <- function(x, sN) {
      x <- sN(x)
      if (isList(x) || isBoot(x)) {
        for (i in 1:length(x)) {
          x[[i]] <- rsubNam(x[[i]], sN)
        }
      }
      x
    }
    rsubNam(x, subNam)

  }

  # Replace any periods in variable names
  # (periods are used to separate variable names for indirect effects)
  sem <- subNam("[.]", "_", sem)

  # Response names (default to endogenous variables)
  ce <- grepl("~~", names(sem))
  if (is.null(r)) r <- names(sem)[!ce]

  # Mediator names (default to endogenous predictors)
  all.p <- unique(unlist(lapply(sem[r], function(i) names(i$t0))))
  all.p <- all.p[!isInt(all.p) & !isPhi(all.p) & !isR2(all.p) & !isRaw(all.p)]
  en.p <- r[r %in% all.p]
  m <- if (length(en.p) > 0) {
    if (is.null(m)) en.p else m
  }

  # Predictor names (default to all predictors)
  if (is.null(p)) {
    ex <- all.p[!all.p %in% r]  # exogenous
    ex <- names(sort(sapply(ex, function(i) {
      lengths(regmatches(i, gregexpr(":", i)))
    })))
    p <- c(ex, en.p)
  }


  # Calculate all direct, total indirect, total, and mediator effects

  # Helper function to create data frames without name modification
  dF <- function(...) {
    data.frame(..., check.names = FALSE)
  }

  # Function to extract direct effects for predictors
  dirEff <- function(D) {
    sapply(p, function(i) {
      D <- if (is.matrix(D)) D[, colnames(D) == i] else unname(D[names(D) == i])
      if (length(D) > 0) D else 0
    }, simplify = FALSE)
  }

  # Function to calculate all indirect effects for predictors
  indEff <- function(D) {

    # Function to extract effect names
    effNam <- function(x) {
      effNam <- function(x) {
        if (is.matrix(x)) colnames(x) else names(x)
      }
      en <- rMapply(effNam, x)
      unique(unlist(en))
    }

    # Function to multiply effects to calculate indirect effects
    # (for each endogenous predictor on a response, multiply its effect by all
    # direct effects on that predictor)
    multEff <- function(x) {

      # Function
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

      # Apply recursively
      rMapply(multEff, x, SIMPLIFY = FALSE)

    }

    # Function to collapse a nested list of vectors/matrices into a single
    # vector/list of vectors
    unlist2 <- function(x) {
      if (all(rapply(x, is.matrix))) {
        x <- rapply(x, function(i) as.list(dF(i)), how = "list")
        lapply(rapply(x, enquote), eval)
      } else unlist(x)
    }

    # Calculate indirect effects: start with last response and move backwards,
    # repeat for n = no. of responses
    lapply(D, function(i) {
      en <- effNam(i)
      if (any(m %in% en)) {
        I <- list()
        I[[1]] <- multEff(i)
        for (j in 1:length(i)) {
          en <- effNam(I[j])
          I[[j + 1]] <- if (any(m %in% en)) multEff(I[j])
        }
        unlist2(I)
      } else NA
    })

  }

  # Function to sum all indirect effects for predictors
  totInd <- function(I) {
    sapply(p, function(i) {
      I <- I[grepl(paste0("[.]", i, "$"), names(I))]
      if (length(I) > 0) Reduce("+", I) else 0
    }, simplify = FALSE)
  }

  # Function to calculate total effects (direct + total indirect)
  totEff <- function(D, I) {
    rMapply(function(i, j) i + j, D, I, SIMPLIFY = FALSE)
  }

  # Function to sum indirect effects operating through each mediator
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

  # Calculate and compile all effects for original estimates
  D <- lapply(sem[r], "[[", 1)
  I <- indEff(D)
  ED <- lapply(D, function(i) list("Direct" = dirEff(i)))
  EI <- lapply(I, function(i) list("Indirect" = totInd(i)))
  ET <- lapply(totEff(ED, EI), function(i) setNames(i, "Total"))
  EM <- lapply(I, function(i) list("Mediators" = totIndM(i)))
  E <- Map(c, ED, EI, ET, EM)

  # Calculate and compile all effects for bootstrapped estimates
  DB <- lapply(sem[r], "[[", 2)
  IB <- indEff(DB)
  EDB <- lapply(DB, function(i) list("Direct" = dirEff(i)))
  EIB <- lapply(IB, function(i) list("Indirect" = totInd(i)))
  ETB <- lapply(totEff(EDB, EIB), function(i) setNames(i, "Total"))
  EMB <- lapply(IB, function(i) list("Mediators" = totIndM(i)))
  EB <- Map(c, EDB, EIB, ETB, EMB)


  # Calculate bootstrapped confidence intervals for all effects

  # Function to calculate CIs (& bias/standard errors)
  bootCI2 <- function(e, eb, r) {

    # Boot object for response var. (replace estimates)
    B <- sem[[r]]
    B$t0 <- e
    B$t <- matrix(eb)

    # Change default CI type for parametric bootstrapping
    if (B$sim == "parametric" && ci.type == "bca") {
      message("Percentile confidence intervals used for parametric bootstrap samples.")
      ci.type <- "perc"
    }

    # Calculate CIs (& bias/standard errors)
    if (!is.na(e)) {
      if (e != 0) {
        bi <- mean(eb, na.rm = TRUE) - e
        se <- sd(eb, na.rm = TRUE)
        ci <- do.call(
          boot::boot.ci,
          c(list(B, ci.conf, ci.type), bci.arg)
        )
        ci <- tail(as.vector(ci[[4]]), 2)
        c(bi, se, ci)
      } else rep(0, 4)
    } else rep(NA, 4)

  }

  # Apply to lists of effects/bootstrapped effects
  CI <- rMapply(bootCI2, E, EB, r, SIMPLIFY = FALSE)


  # Compile and output effects

  # Extract all effects into lists of vectors/matrices (add intercepts)
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
          if (is.B) e <- as.matrix(dF(e))
          if (j != "Mediators") {
            B <- sem[[i]]
            if (is.B) {
              I <- B$t[, isInt(colnames(B$t)), drop = FALSE]
              eb <- cbind(I, e)
              a <- attributes(B$t)[c("sim", "seed", "n")]
              attributes(eb)[names(a)] <- a
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
  e <- extE(E)
  eb <- extE(EB)

  # Generate summary tables of effects and CIs
  E <- rMapply(c, E, CI, SIMPLIFY = FALSE)
  s <- lapply(E, function(i) {

    # List of summary tables
    s <- lapply(names(i), function(j) {

      # Combine effects and CIs into table
      s <- dF(t(dF(i[[j]])))
      names(s) <- c("Effect", "Bias", "Std. Error", "Lower CI", "Upper CI")
      s <- s[s[, 1] != 0, ]
      s <- round(s, digits)
      if (nrow(s) < 1) {
        s[1, ] <- "-"
        rownames(s) <- "n/a"
      }

      # Add significance stars
      stars <- apply(s, 1, function(k) {
        if (is.numeric(k)) {
          k <- c(k[1], tail(k, 2))
          if (all(k > 0) || all(k < 0)) "*" else ""
        } else ""
      })
      s <- cbind(s, " " = stars)

      # Format table (add title columns & top space)
      s <- format(s, nsmall = digits)
      s <- cbind(" " = "", " " = rownames(s), s)
      s[1, 1] <- toupper(j)
      rbind("", s)

    })

    # Combine into one table
    s <- do.call(rbind, s)

    # Format table (text alignment, add top border)
    s[1:2] <- format(s[1:2], justify = "left")
    b <- mapply(function(j, k) {
      n1 <- nchar(k)
      n2 <- max(sapply(j, nchar), n1)
      b <- if (n1 > 1) rep("_", n2) else ""
      paste(b, collapse = "")
    }, s, names(s))
    s <- rbind(b, s)
    rownames(s) <- NULL

    # Set attributes and output
    attr(s, "ci.conf") <- ci.conf
    attr(s, "ci.type") <- ci.type
    s

  })

  # Add any correlated errors
  if (any(ce)) {
    CE <- bootCI(sem[ce], ci.conf, ci.type, digits, bci.arg)
    if (length(ce) > 1) {
      CE <- do.call(rbind, c(CE[1], lapply(CE[-1], "[", 3,)))
      CE[1] <- format(CE[1], justify = "left")
      rownames(CE) <- NULL
    }
    s <- c(s, list("Correlated Errors" = CE))
  }

  # Reinstate periods to all variable names (if they were present)
  e <- subNam("_", ".", e)
  eb <- subNam("_", ".", eb)
  s <- subNam("_", ".", s)

  # Output effects
  e <- list("Effects" = e, "Bootstrapped Effects" = eb, "Summary" = s)
  class(e) <- c("semEff", "list")
  e


}


#' @title Print SEM Effects
#' @description A print method for an object of class `"semEff"`, returning
#'   summary tables of effects and confidence intervals for all responses.
#' @param x An object of class `"semEff"`.
#' @param ... Further arguments passed to or from other methods.
# S3 method for class 'semEff'
#' @export
print.semEff <- function(x, ...) {
  print(x$Summary)
}


#' @title Summarise SEM Effects
#' @description A summary method for an object of class `"semEff"`, which acts
#'   the same as [print.semEff()] except that individual summary tables for
#'   different response variables can be extracted.
#' @param x An object of class `"semEff"`.
#' @param responses An optional character vector, the names of one or more SEM
#'   response variables for which to return summary tables. `"Correlated
#'   Errors"` can also be specified/included here (where applicable).
#' @param ... Further arguments passed to or from other methods.
# S3 method for class 'semEff'
#' @export
summary.semEff <- function(x, responses = NULL, ...) {
  s <- x$Summary
  if (!is.null(responses)) s <- s[responses]
  print(s)
}


#' @title Get SEM Effects
#' @description Extract SEM direct, indirect, and/or total effects from an
#'   object of class `"semEff"`.
#' @param eff An object of class `"semEff"`.
#' @param type The type of effects to return. Must be either `"orig"` (default)
#'   or `"boot"`.
#' @param ... Arguments (above) to be passed to [getEff()] from other extractor
#'   functions.
#' @details These are simple extractor functions for effects calculated using
#'   [semEff()], intended for convenience (e.g. for use with [predEff()]).
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
#'   contain all the variables named in `effects` or all those used to fit
#'   `mod`.
#' @param effects A numeric vector of effects to predict, or a list or nested
#'   list of such vectors. These will typically have been calculated using
#'   [semEff()], [bootEff()], or [stdEff()]. Alternatively, a boot object
#'   produced by [bootEff()] can be supplied.
#' @param eff.boot A matrix of bootstrapped effects used to calculate confidence
#'   intervals for predictions, or a list or nested list of such matrices. These
#'   will have been calculated using [semEff()] or [bootEff()].
#' @param re.form For mixed models of class `"merMod"`, the formula for random
#'   effects to condition on when predicting effects. Defaults to `NA`, meaning
#'   random effects are averaged over. See [predict.merMod()] for further
#'   specification details.
#' @param type The type of prediction to return (for GLMs). Can be either
#'   `"link"` (default) or `"response"`.
#' @param interaction An optional name of an interactive effect, for which to
#'   return standardised effects for a 'main' continuous variable across
#'   different values or levels of interacting variables (see Details).
#' @param use.raw Logical, whether to use raw (unstandardised) effects for all
#'   calculations (if present).
#' @param ci.conf A numeric value specifying the confidence level for confidence
#'   intervals on predictions (and any interactive effects).
#' @param ci.type The type of confidence interval to return (defaults to `"bca"`
#'   — see Details). See [boot.ci()] for further specification details.
#' @param digits The number of significant digits to return for interactive
#'   effects.
#' @param bci.arg A named list of any additional arguments to [boot.ci()],
#'   excepting argument `index`.
#' @param parallel The type of parallel processing to use for calculating
#'   confidence intervals on predictions. Can be one of `"snow"`, `"multicore"`,
#'   or `"no"` (for none — the default).
#' @param ncpus Number of system cores to use for parallel processing. If `NULL`
#'   (default), all available cores are used.
#' @param cl Optional cluster to use if `parallel = "snow"`. If `NULL`
#'   (default), a local cluster is created using the specified number of cores.
#' @param ... Arguments to [stdEff()].
#' @details Generate predicted values for SEM direct, indirect, or total effects
#'   on a response variable, which should be supplied to `effects`. These are
#'   used in place of model coefficients in the standard prediction formula,
#'   with values for predictors drawn either from the data used to fit the
#'   original model(s) (`mod`) or from `newdata`. It is assumed that effects are
#'   fully standardised; however, if this is not the case, then the same
#'   centring and scaling options originally specified to [stdEff()] should be
#'   re-specified — which will then be used to standardise the data. If no
#'   effects are supplied, standardised (direct) effects will be calculated from
#'   the model and used to generate predictions. These predictions will equal
#'   the model(s) fitted values if `newdata = NULL`, `unique.eff = FALSE`, and
#'   `re.form = NULL` (where applicable).
#'
#'   Model-averaged predictions can be generated if averaged `effects` are
#'   supplied to the model in `mod`, or, alternatively, if `weights` are
#'   specified (passed to [stdEff()]) and `mod` is a list of candidate models
#'   (`effects` can also be passed using this latter method). For mixed model
#'   predictions where random effects are included (e.g. `re.form = NULL`), the
#'   latter approach should be used, otherwise the contribution of random
#'   effects will be taken from the single model instead of (correctly) being
#'   averaged over a candidate set.
#'
#'   If bootstrapped effects are supplied to `eff.boot` (or to `effects`, as
#'   part of a boot object), bootstrapped predictions are calculated by
#'   predicting from each effect. Confidence intervals can then be returned via
#'   [bootCI()], for which the `type` should be appropriate for the original
#'   form of bootstrap sampling (defaults to `"bca"`). If the number of
#'   observations to predict is very large, parallel processing (via
#'   [pSapply()]) may speed up the calculation of intervals.
#'
#'   Predictions are always returned in the original (typically unstandardised)
#'   units of the (link-)response variable. For GLMs, they can be returned in
#'   the response scale if `type = "response"`.
#'
#'   Additionally, if the name of an interactive effect is supplied to
#'   `interaction`, standardised effects (and confidence intervals) can be
#'   returned for effects of a continuous 'main' variable across different
#'   values or levels of interacting variable(s). The name should be of the form
#'   `"x1:x2..."`, containing all the variables involved and matching the name
#'   of an interactive effect in the model(s) terms or in `effects`. The values
#'   for all variables should be supplied in `newdata`, with the continuous
#'   variable being automatically identified as having the most unique values.
#' @return A numeric vector of the predictions, or, if bootstrapped effects are
#'   supplied, a list containing the predictions and the upper and lower
#'   confidence intervals. Optional interactive effects may also be appended.
#'   Predictions may also be returned in a list or nested list, depending on the
#'   structure of `mod` (and other arguments).
#' @seealso [predict()]
#' @examples
#' # Predict effects (direct, total)
#' m <- Shipley.SEM
#' e <- Shipley.SEM.Eff
#' dir <- dirEff(e); tot <- totEff(e)
#' f.dir <- predEff(m, effects = dir, type = "response")
#' f.tot <- predEff(m, effects = tot, type = "response")
#'
#' # Using new data for predictors
#' d <- na.omit(Shipley)
#' xn <- c("lat", "DD", "Date", "Growth")
#' seq100 <- function(x) seq(min(x), max(x), length = 100)
#' nd <- data.frame(sapply(d[xn], seq100))
#' f.dir <- predEff(m, nd, dir, type = "response")
#' f.tot <- predEff(m, nd, tot, type = "response")
#' # Add CIs
#' # dir.b <- dirEff(e, "boot"); tot.b <- totEff(e, "boot")
#' # f.dir <- predEff(m, nd, dir, dir.b, type = "response")
#' # f.tot <- predEff(m, nd, tot, tot.b, type = "response")
#'
#' # Predict an interactive effect (e.g. Live ~ Growth * DD)
#' xn <- c("Growth", "DD")
#' d[xn] <- scale(d[xn])  # scale predictors (improves fit)
#' m <- lme4::glmer(Live ~ Growth * DD + (1 | site) + (1 | tree),
#'                  family = binomial, data = d)
#' nd <- with(d, expand.grid(
#'   Growth = seq100(Growth),
#'   DD = mean(DD) + c(-sd(DD), sd(DD))  # two levels for DD
#' ))
#' f <- predEff(m, nd, type = "response", interaction = "Growth:DD")
#' # Add CIs (need to bootstrap model...)
#' # system.time(B <- bootEff(m, ran.eff = "site", R = 1000))
#' # f <- predEff(m, nd, B, type = "response", interaction = "Growth:DD")
#'
#' # Model-averaged predictions (several approaches)
#' m <- Shipley.Growth  # candidate models (list)
#' w <- runif(length(m), 0, 1)  # weights
#' e <- stdEff(m, w)  # averaged effects
#' f1 <- predEff(m[[1]], effects = e)  # pass avg. effects
#' f2 <- predEff(m, weights = w)  # pass weights argument
#' f3 <- avgEst(predEff(m), w)  # use avgEst function
#' stopifnot(all.equal(f1, f2))
#' stopifnot(all.equal(f2, f3))
#'
#' # Compare model fitted values: predEff() vs. fitted()
#' m <- Shipley.SEM$Live
#' f1 <- predEff(m, unique.eff = FALSE, re.form = NULL, type = "response")
#' f2 <- fitted(m)
#' stopifnot(all.equal(f1, f2))
#'
#' # Compare predictions for standardised vs. raw effects
#' f1 <- predEff(m)
#' f2 <- predEff(m, use.raw = TRUE)
#' stopifnot(all.equal(f1, f2))
#' @export
predEff <- function(mod, newdata = NULL, effects = NULL, eff.boot = NULL,
                    re.form = NA, type = "link", interaction = NULL,
                    use.raw = FALSE, ci.conf = 0.95, ci.type = "bca",
                    digits = 3, bci.arg = NULL, parallel = "no", ncpus = NULL,
                    cl = NULL, ...) {

  m <- mod; nd <- newdata; e <- effects; eb <- eff.boot; rf <- re.form;
  ix <- interaction; p <- parallel; nc <- ncpus

  # Arguments to stdEff
  a <- list(...)

  # Weights (for model averaging)
  w <- a$weights; a$weights <- NULL
  if (isList(m) && (isList(w) || is.null(w))) {
    n <- function(x) {NULL}
    N <- rMapply(n, m, SIMPLIFY = FALSE)
    if (is.null(w)) w <- N else N <- lapply(m, n)
    if (is.null(e)) e <- N
    if (is.null(eb)) eb <- N
  }

  # Refit model(s) with any supplied data
  d <- a$data; a$data <- NULL
  if (!is.null(d)) {
    upd <- function(m) eval(update(m, data = d, evaluate = FALSE))
    m <- rMapply(upd, m, SIMPLIFY = FALSE)
  }

  # Environment to look for model data
  env <- a$env; a$env <- NULL
  if (is.null(env)) env <- parent.frame()
  if (!is.null(d)) env <- environment()

  # Effect standardisation options (for back-transforming predictions)
  a$cen.x <- !(isFALSE(a$cen.x) || use.raw)
  a$cen.y <- !(isFALSE(a$cen.y) || use.raw)
  a$std.x <- !(isFALSE(a$std.x) || use.raw)
  a$std.y <- !(isFALSE(a$std.y) || use.raw)
  a$incl.raw <- FALSE

  # Function
  predEff <- function(m, w, e, eb) {

    # Effects
    if (is.null(e)) e <- do.call(stdEff, c(list(m, w, env = env), a))
    if (isBoot(e)) {eb <- e$t; e <- e$t0}
    en <- names(e)

    # Use raw (unstandardised) effects (if present)?
    r <- isRaw(en)
    if (any(r) & use.raw) {
      e[!r] <- e[r]
      if (!is.null(eb)) eb[, !r] <- eb[, r]
    }

    # Subset only relevant parameters
    e <- na.omit(e[!isPhi(en) & !isR2(en) & !r])
    en <- names(e)
    EN <- sapply(en, function(i) {
      unlist(strsplit(i, "(?<!:):(?!:)", perl = TRUE))
    })

    # Data used to fit model(s)
    if (is.null(d)) d <- getData(m, subset = TRUE, merge = TRUE, env = env)
    obs <- rownames(d)

    # Extract the first model, if list
    # (model type and specification should be consistent)
    m1 <- if (isList(m)) m[[1]] else m

    # Model error family/link functions
    f <- if (isBet(m1)) m1$link$mean else family(m1)
    lF <- f$linkfun; lI <- f$linkinv

    # Random effects
    is.re <- isMer(m1) && !identical(rf, NA) && !identical(rf, ~ 0)
    re <- if (is.re) {
      pRE <- function(x) {
        predict(x, nd, re.form = rf, random.only = TRUE)
      }
      re <- rMapply(pRE, m, SIMPLIFY = FALSE)
      if (isList(re)) re <- avgEst(re, w)
      if (is.null(nd)) re[obs] else re
    } else 0

    # Model weights
    w <- weights(m1)
    if (is.null(w)) w <- rep(1, nobs(m1))
    w <- w[w > 0 & !is.na(w)]

    # Model offset(s)
    o <- if (!is.null(nd)) {
      nd <- data.frame(nd)
      tt <- terms(m1)
      on <- attr(tt, "offset")
      if (!is.null(on)) {
        tn <- attr(tt, "variables")
        on <- sapply(on, function(i) tn[[i + 1]])
      }
      on <- c(on, getCall(m1)$offset)
      if (!is.null(on)) rowSums(sapply(on, eval, nd))
    } else {
      mf <- model.frame(m1, data = d)
      model.offset(mf)
    }
    if (is.null(o)) o <- 0

    # Predictors
    dF <- function(...) {
      data.frame(..., check.names = FALSE)
    }
    dM <- function(d) {
      f <- reformulate(names(d))
      dF(model.matrix.lm(f, data = d, na.action = "na.pass"))
    }
    eT <- function(x, d) {
      eT <- function(x) {eval(parse(text = x), d)}
      if (grepl("poly\\(.*[0-9]$", x)) {
        n <- nchar(x)
        xd <- eT(substr(x, 1, n - 1))
        xd[, substr(x, n, n)]
      } else eT(x)
    }
    d2 <- dM(d)
    x <- dF(sapply(en, function(i) {
      if (!isInt(i)) {
        if (isInx(i)) {
          xi <- sapply(EN[[i]], eT, d2)
          if (a$cen.x) xi <- sweep(xi, 2, colMeans(xi))
          apply(xi, 1, prod)
        } else eT(i, d2)
      } else 1
    }))

    # Predictor means/SDs
    xm <- sapply(x, function(i) if (a$cen.x) mean(i) else 0)
    xmw <- sapply(x, function(i) if (a$cen.x) weighted.mean(i, w) else 0)
    xs <- sapply(x, function(i) if (a$std.x) sdW(i, w) else 1)

    # Response mean/SD (link scale)
    ym <- if (a$cen.y) lF(weighted.mean(getY(m1, env = env), w)) else 0
    ys <- if (a$std.y) sdW(getY(m1, link = TRUE, env = env), w) else 1

    # Data to predict (standardise using original means/SDs)
    if (!is.null(nd)) {obs <- rownames(nd); d2 <- dM(nd)}
    x <- dF(sapply(en, function(i) {
      if (!isInt(i)) {
        xi <- if (isInx(i)) {
          xi <- sapply(EN[[i]], eT, d2)
          xi <- sweep(xi, 2, xm[colnames(xi)])
          apply(xi, 1, prod)
        } else eT(i, d2)
        (xi - xmw[i]) / xs[i]
      } else 1
    }), row.names = obs)

    # Predictions
    f <- colSums(e * t(x))
    f <- f * ys + ym + re + o
    if (type == "response") f <- lI(f)
    names(f) <- obs

    # Add CIs
    if (!is.null(eb)) {

      # Bootstrap attributes
      sim <- attr(eb, "sim")
      seed <- attr(eb, "seed")
      n <- attr(eb, "n")
      R <- nrow(eb)

      # Change default CI type for parametric bootstrapping
      if (sim == "parametric" && ci.type == "bca") {
        message("Percentile confidence intervals used for parametric bootstrap samples.")
        ci.type <- "perc"
      }

      # Bootstrapped predictions
      eb <- eb[, en, drop = FALSE]
      fb <- t(sapply(1:R, function(i) {
        ei <- na.omit(eb[i, ])
        xi <- x[names(ei)]
        fi <- colSums(ei * t(xi))
        fi * ys + ym + re + o
      }))
      if (nrow(fb) != R) fb <- t(fb)
      if (type == "response") fb <- lI(fb)

      # Create dummy boot object (for CIs)
      set.seed(seed)
      dd <- dF(rep(1, n))  # dummy data
      B <- list(t0 = f, t = fb, R = R, data = dd, seed = .Random.seed,
                sim = sim, stype = "i", strata = dd[, 1])
      class(B) <- "boot"
      attr(B, "boot_type") <- "boot"

      # Calculate and add CIs
      ci <- pSapply(1:length(f), function(i) {
        ci <- do.call(boot::boot.ci, c(list(B, ci.conf, ci.type, i), bci.arg))
        tail(as.vector(ci[[4]]), 2)
      }, p, nc, cl)
      ci <- as.matrix(ci)
      colnames(ci) <- obs
      f <- list(fit = f, ci.lower = ci[1, ], ci.upper = ci[2, ])

    }

    # Add interactive effects
    if (isTRUE(isInx(ix) && ix %in% en && !is.null(nd))) {

      # Names of variables involved in interaction
      # (ab = all, a = main, b = interacting, a.b = interaction(s))
      ab <- EN[[ix]]; n <- length(ab)
      a <- ab[which.max(sapply(d2[ab], function(i) length(unique(i))))]
      b <- ab[!ab %in% a]
      a.b <- if (n > 2) {
        a.b <- unlist(lapply(2:n, function(i) {
          combn(ab, i, paste, collapse = ":")
        }))
        a.b[sapply(a.b, function(i) a %in% EN[[i]])]
      } else ix

      # Values for interacting variable(s)
      xb <- unique(dF(sapply(b, eT, d2)))
      xb <- sweep(xb, 2, xm[b])
      xb <- lapply(EN[a.b], function(i) {
        apply(xb[i[i %in% b]], 1, prod)
      })

      # Effects
      e <- e * ys / xs
      e <- e[a] + rowSums(mapply("*", e[a.b], xb))
      e <- e * xs[a] / ys
      names(e) <- paste(ix, 1:length(e), sep = "_")

      # Add CIs
      f <- if (!is.null(eb)) {

        # Bootstrapped effects
        eb <- t(sapply(1:R, function(i) {
          ei <- eb[i, ] * ys / xs
          ei <- ei[a] + rowSums(mapply("*", ei[a.b], xb))
          ei * xs[a] / ys
        }))
        if (nrow(eb) != R) eb <- t(eb)

        # CIs
        B$t0 <- e; B$t <- eb
        e <- bootCI(B, ci.conf, ci.type, digits, bci.arg)
        c(f, list(interactions = e))

      } else {
        list(fit = f, interactions = round(e, digits))
      }

    }

    # Output
    if (!is.null(eb)) set.seed(NULL)
    f

  }

  # Apply recursively
  rMapply(predEff, m, w, e, eb, SIMPLIFY = FALSE)

}


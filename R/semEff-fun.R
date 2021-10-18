

#' @title SEM Effects
#' @description Automatically calculate direct, indirect, total, and mediator
#'   effects for endogenous (response) variables in a 'piecewise' structural
#'   equation model (SEM).
#' @param sem A piecewise SEM, comprising a list of fitted model objects or of
#'   boot objects (containing bootstrapped model effects). Alternatively, a
#'   `"psem"` object from
#'   [`piecewiseSEM::psem()`](https://rdrr.io/cran/piecewiseSEM/man/psem.html).
#'   If list is unnamed, response variable names will be used.
#' @param predictors,mediators Names of variables for/through which to calculate
#'   effects. If `NULL` (default), all predictors/mediators in the SEM will be
#'   used.
#' @param use.raw Logical, whether to use 'raw' (unstandardised) effects for all
#'   calculations (if present in `sem`).
#' @param ci.conf A numeric value specifying the confidence level for confidence
#'   intervals on effects.
#' @param ci.type The type of confidence interval to return (defaults to `"bca"`
#'   – see Details). See [boot.ci()] for further specification details.
#' @param digits The number of significant digits to return for numeric values
#'   (for summary tables).
#' @param bci.arg A named list of any additional arguments to [boot.ci()],
#'   excepting argument `index`.
#' @param ... Arguments to [bootEff()].
#' @details The eponymous function of this package calculates all direct,
#'   indirect, total, and mediator effects for a 'piecewise' structural equation
#'   model (SEM), that is, one where parameter estimation is local rather than
#'   global (Lefcheck, 2016; Shipley, 2000, 2009). The SEM simply takes the form
#'   of a list of fitted models, or bootstrapped estimates from such models,
#'   describing hypothesised causal pathways from predictors to response
#'   ('endogenous') variables. These are either direct, or operate indirectly
#'   via other response variables ('mediators'). This list should represent a
#'   directed ('acyclic') causal model, which should be named exactly for each
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
#'   mediator – useful to assess the relative importance of different mediators
#'   in affecting the response. All of these effect types can be calculated
#'   automatically for all (default) or for a specified subset of predictors
#'   and/or mediators in the SEM. As indirect, total, and mediator effects are
#'   not directly bootstrapped using the fitted models for response variables
#'   (i.e. via [bootEff()]), their equivalent 'bootstrapped' estimates are
#'   calculated instead using each bootstrapped direct effect.
#'
#'   Confidence intervals for all effects are returned in summary tables for
#'   each response (see [bootCI()]), with BC*a* intervals calculated by default
#'   using the bootstrapped estimates for each effect type (Cheung, 2009; Hayes
#'   & Scharkow, 2013; MacKinnon et al., 2004). Effects for which the confidence
#'   intervals do not contain zero are highlighted with a star (i.e.
#'   'significant' at the `ci.conf` level). Bootstrap standard errors (standard
#'   deviations of the samples) and biases (sample means minus original
#'   estimates) are also included. Correlated errors (and confidence intervals)
#'   are also returned if their bootstrapped values are present in `sem`, or if
#'   they are specified to argument `cor.err` or as part of a `"psem"` object
#'   (see [bootEff()]). These represent residual relationships among response
#'   variables, unaccounted for by the hypothesised SEM paths. Use `summary()`
#'   for effect summary tables and `print()` to return a table of variable names
#'   and associated details.
#'
#'   All calculated effects and bootstrapped effects are also returned in lists
#'   for each response variable, with all except mediator effects also including
#'   the model intercept(s) – required for prediction (these will be zero for
#'   ordinary linear models with fully standardised effects). Effects can be
#'   conveniently extracted with [getEff()] and related functions.
#' @return A list object of class `"semEff"` for which several methods and
#'   extractor functions exist. Contains:
#'   1. Summary tables of effects and confidence intervals
#'   2. All effects
#'   3. All bootstrapped effects
#'   4. All indirect effects (individual, not summed)
#' @references Cheung, M. W. L. (2009). Comparison of methods for constructing
#'   confidence intervals of standardized indirect effects. *Behavior Research
#'   Methods*, *41*(2), 425-438. \doi{10/fnx7xk}
#'
#'   Hayes, A. F., & Scharkow, M. (2013). The Relative Trustworthiness of
#'   Inferential Tests of the Indirect Effect in Statistical Mediation Analysis:
#'   Does Method Really Matter? *Psychological Science*, *24*(10), 1918-1927.
#'   \doi{10/bbhr}
#'
#'   Lefcheck, J. S. (2016). piecewiseSEM: Piecewise structural equation
#'   modelling in `R` for ecology, evolution, and systematics. *Methods in
#'   Ecology and Evolution*, *7*(5), 573-579. \doi{10/f8s8rb}
#'
#'   MacKinnon, D. P., Lockwood, C. M., & Williams, J. (2004). Confidence Limits
#'   for the Indirect Effect: Distribution of the Product and Resampling
#'   Methods. *Multivariate Behavioral Research*, *39*(1), 99. \doi{10/chqcnx}
#'
#'   Shipley, B. (2000). A New Inferential Test for Path Models Based on
#'   Directed Acyclic Graphs. *Structural Equation Modeling: A Multidisciplinary
#'   Journal*, *7*(2), 206-218. \doi{10/cqm32d}
#'
#'   Shipley, B. (2009). Confirmatory path analysis in a generalized multilevel
#'   context. *Ecology*, *90*(2), 363-368. \doi{10/bqd43d}
#' @examples
#' # SEM effects
#' (shipley.sem.eff <- semEff(shipley.sem.boot))
#' summary(shipley.sem.eff)
#'
#' # Effects for selected variables
#' summary(shipley.sem.eff, response = "Live")
#' # summary(semEff(shipley.sem.boot, predictor = "lat"))
#' # summary(semEff(shipley.sem.boot, mediator = "DD"))
#'
#' # Effects calculated using original SEM (models)
#' # (not typically recommended – better to use saved boot objects)
#' # system.time(
#' #  shipley.sem.eff <- semEff(shipley.sem, R = 1000, seed = 13,
#' #                            ran.eff = "site")
#' # )
#' @export
semEff <- function(sem, predictors = NULL, mediators = NULL, use.raw = FALSE,
                   ci.conf = 0.95, ci.type = "bca", digits = 3, bci.arg = NULL,
                   ...) {

  p <- predictors; m <- mediators


  # Prep

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
      return(i)
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
      return(x)
    }

    # Apply recursively
    rsubNam <- function(x, sN) {
      x <- sN(x)
      if (isList(x) || isBoot(x)) {
        for (i in 1:length(x)) {
          x[[i]] <- rsubNam(x[[i]], sN)
        }
      }
      return(x)
    }
    rsubNam(x, subNam)

  }

  # Replace any periods in variable names with underscores
  # (periods are used to separate variable names for indirect effects)
  sem <- subNam("[.]", "_", sem)
  p <- subNam("[.]", "_", p)
  m <- subNam("[.]", "_", m)

  # Response names
  ce <- grepl("~~", names(sem))
  r <- names(sem)[!ce]

  # Mediator names (default to all endogenous predictors)
  ap <- unique(unlist(lapply(sem[r], function(i) names(i$t0))))
  ap <- ap[!isInt(ap) & !isPhi(ap) & !isR2(ap) & !isRaw(ap)]
  am <- r[r %in% ap]
  m <- if (length(am) > 0) {
    if (is.null(m)) m <- am
    if (!any(m %in% am))
      stop("Mediator(s) not in SEM.")
    am[am %in% m]
  }

  # Predictor names (default to all predictors)
  ex <- ap[!ap %in% r]
  ex <- names(sort(sapply(ex, function(i) {
    lengths(regmatches(i, gregexpr(":", i)))
  })))
  ap <- c(ex, am)
  if (is.null(p)) p <- ap
  if (!any(p %in% ap))
    stop("Predictor(s) not in SEM.")
  p <- ap[ap %in% p]


  # Calculate all direct, total indirect, total, and mediator effects

  # Helper function to create data frames without modifying names
  dF <- function(...) {
    data.frame(..., check.names = FALSE, fix.empty.names = FALSE)
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
      en <- rMapply(function(i) {
        if (is.matrix(i)) colnames(i) else names(i)
      }, x)
      unique(unlist(en))
    }

    # Function to multiply effects to calculate indirect effects
    # (for each endogenous predictor on a response, multiply its effect by all
    # direct effects on that predictor)
    multEff <- function(x) {

      # Function
      multEff <- function(x) {
        if (is.matrix(x)) {
          x <- x[, colnames(x) %in% am, drop = FALSE]
          Map(function(i, j) {
            eb <- sem[[j]]$t
            i * eb[, colnames(eb) %in% ap, drop = FALSE]
          }, data.frame(x), colnames(x))
        } else {
          x <- x[names(x) %in% am]
          Map(function(i, j) {
            e <- sem[[j]]$t0
            i * e[names(e) %in% ap]
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

    # Calculate indirect effects
    # (start with last response variable and move backwards through others)
    lapply(D, function(i) {
      en <- effNam(i)
      if (any(am %in% en)) {
        I <- list()
        I[[1]] <- multEff(i)
        for (j in 1:length(i)) {
          en <- effNam(I[j])
          I[[j + 1]] <- if (any(am %in% en)) multEff(I[j])
        }
        unlist2(I)
      } else NA
    })

  }

  # Function to sum all indirect effects for predictors
  totInd <- function(I) {
    sapply(p, function(i) {
      P <- paste0("[.]", i, "$")
      I <- I[grepl(P, names(I))]
      if (length(I) > 0) {
        M <- paste(sapply(m, function(j) {
          paste(paste0("^", j, "[.]"), paste0("[.]", j, "[.]"), sep = "|")
        }), collapse = "|")
        I <- I[grepl(M, names(I))]
        if (length(I) > 0) Reduce("+", I) else 0
      } else 0
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
        M <- paste(paste0("^", i, "[.]"), paste0("[.]", i, "[.]"), sep = "|")
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
    if (B$sim == "parametric" && ci.type[1] == "bca") {
      message("Percentile confidence intervals used for parametric bootstrap samples.")
      ci.type <- "perc"
    }

    # Calculate CIs (& bias/standard errors)
    if (!is.na(e)) {
      if (e != 0) {
        bi <- mean(eb, na.rm = TRUE) - e
        se <- sd(eb, na.rm = TRUE)
        ci <- suppressWarnings(
          do.call(
            boot::boot.ci,
            c(list(B, ci.conf, ci.type), bci.arg)
          )
        )
        ci <- tail(as.vector(ci[[4]]), 2)
        c(bi, se, ci)
      } else rep(0, 4)
    } else rep(NA, 4)

  }

  # Calculate CIs and append to effects
  CI <- rMapply(bootCI2, E, EB, r, SIMPLIFY = FALSE)
  ECI <- rMapply(c, E, CI, SIMPLIFY = FALSE)


  # Compile and output effects

  # Helper function to add a top border to a data frame
  tB <- function(d) {
    b <- mapply(function(i, j) {
      n1 <- nchar(j)
      n2 <- max(sapply(i, nchar), n1, 3)
      b <- if (n1 > 1) rep("-", n2) else ""
      paste(b, collapse = "")
    }, d, names(d))
    rbind(b, d)
  }

  # Extract all effects into lists of vectors/matrices
  # (remove zeros, add intercepts)
  extEff <- function(E) {
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
  e <- extEff(E)
  eb <- extEff(EB)

  # All indirect effects
  ai <- lapply(I, function(i) {
    P <- paste(sapply(p, function(j) paste0("[.]", j, "$")), collapse = "|")
    i <- i[grepl(P, names(i))]
    if (length(i) > 0) {
      M <- paste(sapply(m, function(j) {
        paste(paste0("^", j, "[.]"), paste0("[.]", j, "[.]"), sep = "|")
      }), collapse = "|")
      i <- i[grepl(M, names(i))]
      if (length(i) > 0) i else NA
    } else NA
  })

  # Summary tables of effects and CIs
  s <- lapply(ECI, function(i) {

    # List of summary tables
    s <- lapply(names(i), function(j) {

      # Combine effects and CIs into table
      s <- dF(t(dF(i[[j]])))
      names(s) <- c("Effect", "Bias", "Std. Err.", "Lower CI", "Upper CI")
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
      s <- dF(s, " " = stars)

      # Format table (add title columns, borders, top space)
      s <- format(s, nsmall = digits)
      s <- dF(" " = "", " " = rownames(s), "|",
              s[1], "|", s[2], "|", s[3], "|", s[4:5], "|", s[6])
      s[1, 1] <- toupper(j)
      rbind("", s)

    })

    # Combine into one and format (add borders, text alignment, etc.)
    s <- tB(do.call(rbind, s))[-2, ]
    s[, 2] <- subNam("_", ".", s[, 2])
    s[1:2] <- format(s[1:2], justify = "left")
    rownames(s) <- 1:nrow(s)

    # Set attributes and output
    class(s) <- c("semEff", class(s))
    attr(s, "ci.conf") <- ci.conf
    attr(s, "ci.type") <- ci.type
    s

  })

  # Add table of correlated errors
  if (any(ce)) {
    CE <- bootCI(sem[ce], ci.conf, ci.type, digits, bci.arg)
    if (length(ce) > 1) {
      CE <- c(CE[1], lapply(CE[-1], "[", 2,))
      CE <- do.call(rbind, CE)
      CE[1] <- format(CE[1], justify = "left")
      rownames(CE) <- 1:nrow(CE)
    }
    CE[, 1] <- subNam("_", ".", CE[, 1])
    class(CE) <- c("semEff", class(CE))
    s <- c(s, list("Correlated Errors" = CE))
  }

  # Add table of variables
  v <- c(ex[ex %in% p], r)
  # y <- "\u2713"; n <- "x"  # issues w/ unicode tick marks...
  y <- "Y"; n <- "N"
  v <- dF(
    " " = subNam("_", ".", v),
    "|",
    "Category" = ifelse(v %in% ex, "Exog.", "Endog."),
    "|",
    "Predictor" = ifelse(v %in% p, y, n),
    "Mediator" = ifelse(v %in% m, y, n),
    "Response" = sapply(v, function(i) {
      if (sum(E[[i]]$Total != 0)) y else n
    }),
    "|",
    "Dir. Eff." = sapply(v, function(i) {
      if (!i %in% ex) sum(E[[i]]$Direct != 0) else "-"
    }),
    "Ind. Eff." = sapply(v, function(i) {
      if (!i %in% ex) length(na.omit(ai[[i]])) else "-"
    }),
    "|"
  )
  if (any(ce)) {
    v <- dF(
      v,
      "Cor. Err." = sapply(v[, 1], function(i) {
        if (!i %in% ex) {
          cv <- unlist(lapply(CE[, 1], function(j) {
            gsub(" ", "", unlist(strsplit(j, "~~")))
          }))
          sum(cv == i)
        } else "-"
      }),
      "|"
    )
  }
  v <- tB(v)
  v[c(1:3)] <- format(v[c(1:3)], justify = "left")
  v[c(5:7)] <- format(v[c(5:7)], justify = "centre")
  rownames(v) <- 1:nrow(v)
  class(v) <- c("semEff", class(v))
  s <- c(list("Variables" = v), s)

  # Reinstate periods to variable names
  e <- subNam("_", ".", e)
  eb <- subNam("_", ".", eb)
  names(s) <- subNam("_", ".", names(s))
  names(ai) <- subNam("_", ".", names(ai))

  # Output effects
  e <- list("Summary" = s, "Effects" = e, "Bootstrapped Effects" = eb,
            "All Indirect Effects" = ai)
  class(e) <- c("semEff", class(e))
  e


}


#' @title Print `"semEff"` Objects
#' @description A [print()] method for an object of class `"semEff"`.
#' @param x An object of class `"semEff"`.
#' @param ... Further arguments passed to or from other methods. Not currently
#'   used.
#' @details This print method returns a summary table for the SEM variables,
#'   giving their status as exogenous or endogenous and as predictor, mediator
#'   and/or response. It also gives the number of direct vs. indirect paths
#'   leading to each variable, and the number of correlated errors (if
#'   applicable).
#'
#'   Printing of summary tables uses a custom version of `print.data.frame()`,
#'   facilitating correct rendering of unicode characters by bypassing
#'   [format.data.frame()] ([bug
#'   details](https://stat.ethz.ch/pipermail/r-devel/2015-May/071252.html),
#'   workaround adapted from
#'   [here](https://stat.ethz.ch/pipermail/r-devel/2015-May/071259.html)). Row
#'   names (numbers) are also suppressed by default.
#' @return A summary table for the SEM variables (data frame).
# S3 method for class 'semEff'
#' @export
print.semEff <- function(x, ...) {

  # Custom print.data.frame() for summary tables
  # (unicode support, rownames suppressed)
  print.semEff.table <- function(x, ..., digits = NULL, quote = FALSE,
                                 right = TRUE, row.names = FALSE, max = NULL) {

    n <- length(row.names(x))
    if (length(x) == 0L) {
      cat(sprintf(ngettext(n, "data frame with 0 columns and %d row",
                           "data frame with 0 columns and %d rows"), n),
          "\n", sep = "")
    }
    else if (n == 0L) {
      print.default(names(x), quote = FALSE)
      cat(gettext("<0 rows> (or 0-length row.names)\n"))
    }
    else {
      if (is.null(max))
        max <- getOption("max.print", 99999L)
      if (!is.finite(max))
        stop("invalid 'max' / getOption(\"max.print\"): ",
             max)
      omit <- (n0 <- max %/% length(x)) < n
      # m <- as.matrix(format.data.frame(if (omit)
      #   x[seq_len(n0), , drop = FALSE]
      #   else x, digits = digits, na.encode = FALSE))
      m <- as.matrix(if (omit)
        x[seq_len(n0), , drop = FALSE]
        else x)
      if (!isTRUE(row.names))
        dimnames(m)[[1L]] <- if (isFALSE(row.names))
          rep.int("", if (omit)
            n0
            else n)
      else row.names
      print(m, ..., quote = quote, right = right, max = max)
      if (omit)
        cat(" [ reached 'max' / getOption(\"max.print\") -- omitted",
            n - n0, "rows ]\n")
    }
    invisible(x)

  }

  # Print semEff object
  if ("list" %in% class(x)) {

    # SEM variable details
    v <- x$Summary$Variables
    ct <- v$Category
    di <- v[c("Dir. Eff.", "Ind. Eff.")]
    di <- suppressWarnings(
      na.omit(apply(di, 2, as.numeric))
    )
    n1 <- paste(sum(grepl("^Ex", ct)))
    n2 <- paste(sum(grepl("^En", ct)))
    n3 <- paste(sum(di[, 1]))
    n4 <- paste(sum(di[, 2]))

    # Correlated errors?
    ce <- "Cor. Err."
    ce <- if (ce %in% names(v)) {
      n5 <- suppressWarnings(
        sum(na.omit(as.numeric(v[, ce]))) / 2
      )
      paste0("  * ", n5, " correlated error(s)\n")
    }

    # Print variable table
    message("\nPiecewise SEM with:\n  * ", n1, " exogenous vs. ", n2,
            " endogenous variable(s)\n  * ", n3, " direct vs. ", n4,
            " indirect effect(s)\n", ce, "\nVariables:\n")
    print.semEff.table(v)
    message("\nUse summary() for effects and confidence intervals for endogenous variables.\n")

  }
  else print.semEff.table(x)

}


#' @title Summarise SEM Effects
#' @description A [summary()] method for an object of class `"semEff"`.
#' @param object An object of class `"semEff"`.
#' @param responses An optional character vector, the names of one or more SEM
#'   response variables for which to return summaries (and/or `"Correlated
#'   Errors"`, where applicable). Can also be a numeric vector of indices of
#'   `object`. If `NULL` (default), all summaries are returned.
#' @param ... Further arguments passed to or from other methods. Not currently
#'   used.
#' @details This summary method prints tables of effects and confidence
#'   intervals for SEM endogenous (response) variables.
#' @return A summary table or tables of effects for the endogenous variables
#'   (data frames).
# S3 method for class 'semEff'
#' @export
summary.semEff <- function(object, responses = NULL, ...) {

  # SEM response names
  s <- object$Summary[-1]
  r <- names(s)
  r2 <- responses
  if (is.null(r2)) r2 <- r
  if (is.numeric(r2)) r2 <- r[r2]
  if (!any(r2 %in% r))
    stop("Response(s) not in SEM.")

  # Print summary tables
  ce <- "Correlated Errors"
  if (length(r2[r2 != ce]) > 0) {
    message("\nSEM direct, summed indirect, total, and mediator effects:")
  }
  r <- r[r != ce]
  for (i in r2) {
    n <- which(r == i)
    ii <- if (length(n) > 0) {
      paste0(i, " (", n, "/", length(r), ")")
    } else i
    message("\n", ii, ":\n")
    print(s[[i]])
    cat("\n")
  }

}


#' @title Get SEM Effects
#' @description Extract SEM effects from an object of class `"semEff"`.
#' @param eff An object of class `"semEff"`.
#' @param responses An optional character vector, the names of one or more SEM
#'   response variables for which to return effects. Can also be a numeric
#'   vector of indices of `eff`. If `NULL` (default), all effects are returned.
#' @param type The type of effects to return. Can be `"orig"` (default) or
#'   `"boot"` (for bootstrapped).
#' @param ... Arguments (above) to be passed to `getEff()` from the other
#'   extractor functions (`type = "boot"` is not available for `getAllInd()`).
#' @details These are simple extractor functions for effects calculated using
#'   [semEff()], intended for convenience (e.g. for use with [predEff()]).
#' @return A list containing the original or bootstrapped effects for each
#'   response variable, as numeric vectors or matrices (respectively).
#' @name getEff
NULL
#' @describeIn getEff Extract effects.
#' @export
getEff <- function(eff, responses = NULL, type = c("orig", "boot")) {

  e <- eff; r <- responses; type <- match.arg(type)
  if (class(e)[1] != "semEff")
    stop("Object is not of class 'semEff'")

  # Extract effects
  e <- if (type == "boot") e[[3]] else e[[2]]
  en <- names(e)

  # Subset responses
  if (is.null(r)) r <- en
  if (is.numeric(r)) r <- en[r]
  if (!any(r %in% en))
    stop("Response(s) not in SEM.")
  e[en %in% r]

}
#' @describeIn getEff Extract direct effects.
#' @export
getDirEff <- function(...) {
  e <- lapply(getEff(...), "[[", 1)
  if (length(e) < 2) e[[1]] else e
}
#' @describeIn getEff Extract indirect effects.
#' @export
getIndEff <- function(...) {
  e <- lapply(getEff(...), "[[", 2)
  if (length(e) < 2) e[[1]] else e
}
#' @describeIn getEff Extract total effects.
#' @export
getTotEff <- function(...) {
  e <- lapply(getEff(...), "[[", 3)
  if (length(e) < 2) e[[1]] else e
}
#' @describeIn getEff Extract mediator effects.
#' @export
getMedEff <- function(...) {
  e <- lapply(getEff(...), "[[", 4)
  if (length(e) < 2) e[[1]] else e
}
#' @describeIn getEff Extract all indirect effects.
#' @export
getAllInd <- function(eff, ...) {
  e <- eff[c(1, 4)]
  class(e) <- class(eff)
  e <- getEff(e, type = "orig", ...)
  if (length(e) < 2) e[[1]] else e
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
#'   – see Details). See [boot.ci()] for further specification details.
#' @param digits The number of significant digits to return for interactive
#'   effects.
#' @param bci.arg A named list of any additional arguments to [boot.ci()],
#'   excepting argument `index`.
#' @param parallel The type of parallel processing to use for calculating
#'   confidence intervals on predictions. Can be one of `"snow"`, `"multicore"`,
#'   or `"no"` (for none – the default).
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
#'   re-specified – which will then be used to standardise the data. If no
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
#'   units of the (link-transformed) response variable. For GLMs, they can be
#'   returned in the response scale if `type = "response"`.
#'
#'   Additionally, if the name of an interactive effect is supplied to
#'   `interaction`, standardised effects (and confidence intervals) can be
#'   returned for effects of a continuous 'main' variable across different
#'   values or levels of interacting variable(s). The name should be of the form
#'   `"x1:x2..."`, containing all the variables involved and matching the name
#'   of an interactive effect in the model(s) terms or in `effects`. The values
#'   for all variables should be supplied in `newdata`, with the main continuous
#'   variable being automatically identified as that having the most unique
#'   values.
#' @return A numeric vector of the predictions, or, if bootstrapped effects are
#'   supplied, a list containing the predictions and the upper and lower
#'   confidence intervals. Optional interactive effects may also be appended.
#'   Predictions may also be returned in a list or nested list, depending on the
#'   structure of `mod` (and other arguments).
#' @seealso [predict()]
#' @examples
#' # Predict effects (direct, total)
#' m <- shipley.sem
#' e <- shipley.sem.eff
#' dir <- getDirEff(e)
#' tot <- getTotEff(e)
#' f.dir <- predEff(m, effects = dir, type = "response")
#' f.tot <- predEff(m, effects = tot, type = "response")
#' f.dir$Live[1:10]
#' f.tot$Live[1:10]
#'
#' # Using new data for predictors
#' d <- na.omit(shipley)
#' xn <- c("lat", "DD", "Date", "Growth")
#' seq100 <- function(x) seq(min(x), max(x), length = 100)
#' nd <- data.frame(sapply(d[xn], seq100))
#' f.dir <- predEff(m, nd, dir, type = "response")
#' f.tot <- predEff(m, nd, tot, type = "response")
#' f.dir$Live[1:10]
#' f.tot$Live[1:10]
#' # Add CIs
#' # dir.b <- getDirEff(e, "boot")
#' # tot.b <- getTotEff(e, "boot")
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
#' f$fit[1:10]
#' f$interaction
#' # Add CIs (need to bootstrap model...)
#' # system.time(B <- bootEff(m, R = 1000, ran.eff = "site"))
#' # f <- predEff(m, nd, B, type = "response", interaction = "Growth:DD")
#'
#' # Model-averaged predictions (several approaches)
#' m <- shipley.growth  # candidate models (list)
#' w <- runif(length(m), 0, 1)  # weights
#' e <- stdEff(m, w)  # averaged effects
#' f1 <- predEff(m[[1]], effects = e)  # pass avg. effects
#' f2 <- predEff(m, weights = w)  # pass weights argument
#' f3 <- avgEst(predEff(m), w)  # use avgEst function
#' stopifnot(all.equal(f1, f2))
#' stopifnot(all.equal(f2, f3))
#'
#' # Compare model fitted values: predEff() vs. fitted()
#' m <- shipley.sem$Live
#' f1 <- predEff(m, unique.eff = FALSE, re.form = NULL, type = "response")
#' f2 <- fitted(m)
#' stopifnot(all.equal(f1, f2))
#'
#' # Compare predictions using standardised vs. raw effects (same)
#' f1 <- predEff(m)
#' f2 <- predEff(m, use.raw = TRUE)
#' stopifnot(all.equal(f1, f2))
#' @export
predEff <- function(mod, newdata = NULL, effects = NULL, eff.boot = NULL,
                    re.form = NA, type = c("link", "response"),
                    interaction = NULL, use.raw = FALSE, ci.conf = 0.95,
                    ci.type = "bca", digits = 3, bci.arg = NULL,
                    parallel = "no", ncpus = NULL, cl = NULL, ...) {

  m <- mod; nd <- newdata; e <- effects; eb <- eff.boot; rf <- re.form;
  type <- match.arg(type); ix <- interaction; nc <- ncpus

  # Arguments to stdEff()
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
    d <- getData(m, subset = TRUE, merge = TRUE, env = env)
    obs <- rownames(d)

    # Extract the first model, if list
    # (model type and specification should be consistent)
    m1 <- if (isList(m)) m[[1]] else m

    # Model error family/link functions
    f <- if (isBet(m1)) m1$link$mean else family(m1)
    lF <- f$linkfun
    lI <- f$linkinv

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
    x <- getX(en, d, as.df = TRUE)
    xc <- getX(en, d, centre = TRUE, as.df = TRUE)
    x[isInx(names(x))] <- xc[isInx(names(xc))]

    # Predictor means/SDs
    xm <- sapply(x, function(i) if (a$cen.x) mean(i) else 0)
    xmw <- sapply(x, function(i) if (a$cen.x) weighted.mean(i, w) else 0)
    xs <- sapply(x, function(i) if (a$std.x) sdW(i, w) else 1)

    # Response mean/SD (link scale)
    ym <- if (a$cen.y) lF(weighted.mean(getY(m1, env = env), w)) else 0
    ys <- if (a$std.y) sdW(getY(m1, link = TRUE, env = env), w) else 1

    # Data to predict (standardise using original means/SDs)
    if (!is.null(nd)) {
      d <- nd
      obs <- rownames(d)
      x <- getX(en, d, as.df = TRUE)
      xc <- getX(en, d, centre = xm, as.df = TRUE)
      x[isInx(names(x))] <- xc[isInx(names(xc))]
    }
    x <- sweep(x, 2, xmw)
    x <- sweep(x, 2, xs, "/")
    x[isInt(en)] <- 1

    # Predictions
    f <- colSums(e * t(x))
    f <- f * ys + ym + re + o
    if (type == "response") f <- lI(f)
    names(f) <- obs

    # Add CIs (& bias/SE)
    if (!is.null(eb)) {

      # Bootstrap attributes
      sim <- attr(eb, "sim")
      seed <- attr(eb, "seed")
      n <- attr(eb, "n")
      R <- nrow(eb)

      # Change default CI type for parametric bootstrapping
      if (sim == "parametric" && ci.type[1] == "bca") {
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

      # Bootstrap bias/standard errors
      bi <- colMeans(fb, na.rm = TRUE) - f
      se <- apply(fb, 2, sd, na.rm = TRUE)
      se[is.na(f)] <- NA

      # Create dummy boot object (for CIs)
      set.seed(seed)
      dd <- data.frame(rep(1, n))  # dummy data
      B <- list(t0 = f, t = fb, R = R, data = dd, seed = .Random.seed,
                sim = sim, stype = "i", strata = dd[, 1])
      class(B) <- "boot"
      attr(B, "boot_type") <- "boot"

      # Calculate and add CIs
      ci <- pSapply(1:length(f), function(i) {
        ci <- do.call(
          boot::boot.ci,
          c(list(B, ci.conf, ci.type, i), bci.arg)
        )
        tail(as.vector(ci[[4]]), 2)
      }, parallel, nc, cl)
      ci <- as.matrix(ci)
      colnames(ci) <- obs
      f <- list(fit = f, bias = bi, se.fit = se, ci.lower = ci[1, ],
                ci.upper = ci[2, ])

    }

    # Add interactive effects
    if (isTRUE(isInx(ix) && ix %in% en && !is.null(nd))) {

      # Names of variables involved in interaction
      # (ab = all, a = main, b = interacting, a.b = interaction(s))
      ab <- EN[[ix]]
      xi <- getX(ab, d, as.df = TRUE)
      a <- ab[which.max(sapply(xi, function(i) length(unique(i))))]
      b <- ab[!ab %in% a]
      n <- length(ab)
      a.b <- if (n > 2) {
        a.b <- unlist(lapply(2:n, function(i) {
          combn(ab, i, paste, collapse = ":")
        }))
        a.b[sapply(a.b, function(i) a %in% EN[[i]])]
      } else ix

      # Values for interacting variable(s)
      xb <- sweep(unique(xi[b]), 2, xm[b])
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
        B$t0 <- e
        B$t <- eb
        e <- bootCI(B, ci.conf, ci.type, digits, bci.arg)
        c(f, list(interactions = e))

      } else {
        list(fit = f, interactions = round(e, digits))
      }

    }

    # Output
    if (!is.null(eb)) {
      set.seed(NULL)
      attr(f, "ci.conf") <- ci.conf
      attr(f, "ci.type") <- ci.type
    }
    f

  }

  # Apply recursively
  rMapply(predEff, m, w, e, eb, SIMPLIFY = FALSE)

}


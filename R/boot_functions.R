

#' @title Bootstrap a Piecewise SEM
#' @description Generate bootstrapped standardised coefficients and correlated
#'   errors for a list of fitted models, comprising a piecewise Structural
#'   Equation Model (SEM).
#' @param sem A piecewise SEM, comprising a named list/nested list of fitted
#'   model objects of class \code{lm}, \code{glm}, or \code{lmerMod}.
#' @param cor.err An optional character vector describing the correlated errors
#'   in the SEM. This should be of the form: \code{c("a ~~ b", "c ~~ d", ...)}
#'   (spaces optional), with the names corresponding to model names.
#' @param data An optional dataset used to first re-fit the model(s).
#' @param ran.eff For mixed models with nested random effects, the name of the
#'   variable containing the highest-level random effect. For non-nested random
#'   effects, specify \code{"crossed"}. If argument is not specified and
#'   \code{m} is a mixed model(s), an error will be thrown.
#' @param R Number of bootstrap replicates to generate.
#' @param seed Seed for the random number generator. If not provided, a random
#'   five-digit integer is used (see Details).
#' @param catch.err Logical, should errors generated during model fitting or
#'   estimation be caught and \code{NA} returned? If \code{FALSE}, any such
#'   errors will cause the function to exit.
#' @param parallel The type of parallel processing to use. Can be one of
#'   \code{"snow"}, \code{"multicore"}, or \code{"no"} (for none).
#' @param ncpus Number of system cores to use for parallel processing. If
#'   \code{NULL} (default), all available cores are used.
#' @param cl Optional cluster to use if \code{parallel = "snow"}. If \code{NULL}
#'   (default), a local cluster is created using the specified number of cores.
#' @param bM.arg A named list of additional arguments to \code{bootMer}.
#' @param ... Arguments to \code{stdCoeff}.
#' @details \code{bootSEM} uses the \code{boot} function (primarily) to
#'   bootstrap standardised coefficients from a fitted model or list of models
#'   (calculated using \code{stdCoeff}), where bootstrapping is typically
#'   nonparametric, i.e. coefficients are calculated using data where the rows
#'   have been randomly sampled with replacement. The number of replicates is
#'   set by default to 10,000, which should provide accurate coverage for
#'   confidence intervals in most situations.
#'
#'   To ensure that data is resampled in the same way across bootstrap
#'   operations, the same seed is set per individual operation, with the value
#'   then saved as an additional attribute to all boot objects (for
#'   reproducibility). The seed can either be user-supplied or a random
#'   five-digit number (default), and is always re-initialised on exit (i.e.
#'   \code{set.seed(NULL)}).
#'
#'   Where \code{weights} are specified, bootstrapped coefficients will be a
#'   weighted average across the set of candidate models for each response
#'   variable, calculated after each model is first refit to the resampled
#'   dataset (specifying \code{weights = "equal"} will use a simple average
#'   instead). If no weights are specified and \code{sem} is a nested list of
#'   models, the function will throw an error, as it will be expecting weights
#'   for a presumed model averaging scenario. If the user wishes to override
#'   this behaviour and instead return a boot object of estimates for each
#'   individual model, they should recursively apply the function using
#'   \code{rMapply}.
#'
#'   Where names of response variables with correlated errors are supplied to
#'   \code{cor.err}, the function will also return bootstrapped Pearson
#'   correlated errors (residuals) for the models fitting those variables
#'   (residuals used are the default type returned by \code{resid}). Where
#'   weights are supplied and \code{sem} is a nested list, errors are averaged
#'   across candidate models (as for coefficients) prior to correlating. If any
#'   two models (or candidate sets) with correlated errors were fit to different
#'   subsets of data observations, both models/sets are first refit to data
#'   containing only the common observations.
#'
#'   For mixed models with nested random effects, the highest-level random
#'   effect group (only) in the dataset is resampled, a procedure which best
#'   retains the hierarchical structure of the data (Davison & Hinkley 1997, Ren
#'   \emph{et al.} 2010). Lower-level groups or individual observations are not
#'   themselves resampled, as these are not independent. The name of this random
#'   effect must be supplied to \code{ran.eff}, matching the name in the data.
#'   This type of resampling will result in different sized datasets if
#'   observations are unbalanced across groups; however this should not be a
#'   problem as the number of independent units (groups), and hence the 'degrees
#'   of freedom', remains unchanged (see \url{http://bit.ly/2YFObSE}). For
#'   non-nested random effects (\code{ran.eff = "crossed"}), group resampling
#'   will not be appropriate, and (semi-)parametric bootstrapping is performed
#'   instead via \code{bootMer} in the \pkg{lme4} package. Users should think
#'   carefully about whether their random effects are nested or not (see
#'   \url{http://bit.ly/2K244jK}). (As \code{bootMer} takes only a fitted model
#'   as its first argument, any model averaging is calculated 'post-hoc' using
#'   the estimates in boot objects for each candidate model, rather than during
#'   the bootstrapping process itself (i.e. the default procedure via
#'   \code{boot}). Results are then returned in a new boot object for each
#'   response variable or correlated error estimate.)
#'
#'   Parallel processing is used by default via the \pkg{parallel} package and
#'   option \code{parallel = "snow"} (and is generally recommended), but users
#'   can specify the type of parallel processing to use, or none. If
#'   \code{"snow"}, a cluster of workers is created using \code{makeCluster},
#'   and the user can specify the number of system cores to incorporate in the
#'   cluster (defaults to all available). \code{bootSEM} then exports all
#'   required objects and functions to this cluster using \code{clusterExport},
#'   after performing a (rough) match of all objects and functions in the
#'   current global environment to those referenced in the model call(s). Users
#'   should load any required external packages prior to calling the function.
#'
#' @note Bootstrapping mixed (or any other) models may take a very long time
#'   when the number of replicates, observations, parameters, and/or models is
#'   high. To decrease processing time, it may be worth trying different
#'   optimizers and/or other options to generate faster estimates (always check
#'   results).
#' @return An object of class \code{boot} containing the bootstrapped
#'   coefficients, or a list/nested list of such objects.
#' @references Burnham, K. P., & Anderson, D. R. (2002). \emph{Model Selection
#'   and Multimodel Inference: A Practical Information-Theoretic Approach} (2nd
#'   ed.). New York: Springer-Verlag. Retrieved from \url{http://bit.ly/2MwlTHa}
#'
#'   Davison, A. C., & Hinkley, D. V. (1997). \emph{Bootstrap Methods and their
#'   Application}. Cambridge University Press.
#'
#'   Ren, S., Lai, H., Tong, W., Aminzadeh, M., Hou, X., & Lai, S. (2010).
#'   Nonparametric bootstrapping for hierarchical data. \emph{Journal of Applied
#'   Statistics}, \strong{37}(9), 1487–1498. \url{https://doi.org/dvfzcn}
#' @seealso \code{\link[boot]{boot}}, \code{\link[lme4]{bootMer}},
#'   \code{\link[semEff]{stdCoeff}}, \code{\link[semEff]{avgEst}},
#'   \code{\link[stats]{resid}}
#' @examples
#' ## Bootstrap Shipley SEM
#' ## (set 'site' as random effect group for resampling - highest-level)
#'
#' \dontrun{
#'
#' ## 45-60 mins with parallel processing and eight cores (YMMV)
#' system.time(
#'   Shipley.SEM.boot <- bootSEM(Shipley.SEM, ran.eff = "site", seed = 53908)
#' )
#' }
#'
#' ## Original estimates
#' lapply(Shipley.SEM.boot, "[[", 1)
#'
#' ## Bootstrapped estimates
#' lapply(Shipley.SEM.boot, function(i) head(i$t))
#' @export
bootSEM <- function(sem, cor.err = NULL, data = NULL, ran.eff = NULL, R = 10000,
                    seed = NULL, catch.err = TRUE, parallel = "snow",
                    ncpus = NULL, cl = NULL, bM.arg = NULL, ...) {

  m <- sem; ce <- cor.err; d <- data; re <- ran.eff; p <- parallel; nc <- ncpus

  ## Main function environment
  env <- environment()

  ## Arguments to stdCoeff
  a <- list(...)

  ## Weights (for model averaging)
  w <- eval(a$weights); a$weights <- NULL
  nw <- is.null(w)
  if (isList(m) && (nw || all(w == "equal"))) {
    if (nw && any(sapply(m, isList)))
      stop("'weights' must be supplied for model averaging (or specify 'equal').")
    w <- lapply(m, function(i) w)
  }

  ## Names of models with correlated errors (list)
  any.ce <- isList(m) && !is.null(ce)
  if (any.ce) {
    cv <- sapply(ce, function(i) {
      gsub(" ", "", unlist(strsplit(i, "~~")))
    }, simplify = FALSE)
    cv <- cv[sapply(cv, function(i) {
      all(i %in% names(m)) && all(i %in% names(w))
    })]
    if (length(cv) < 1)
      stop("Names of var(s) with correlated errors missing from list of models and/or weights.")
  }

  ## Mixed models?
  mer <- all(unlist(rMapply(isMerMod, m)))
  if (mer && is.null(re))
    stop("Name of highest-level random effect to sample must be specified to 'ran.eff' (or specify 'crossed').")
  mer2 <- mer && re == "crossed"
  if (mer2) {

    ## Modified bootMer function
    bootMer2 <- function(...) {

      ## Set up function call with specified arguments
      C <- match.call()
      n <- length(C)
      a <- c(list(FUN = s, nsim = R, parallel = p, ncpus = nc, cl = cl),
             bM.arg)
      for (i in 1:length(a)) {
        C[n + i] <- a[i]
        names(C)[n + i] <- names(a)[i]
      }
      C[[1]] <- as.name("bootMer")

      ## Run
      set.seed(seed)
      eval.parent(C)

    }

    ## Create boot statistic object for later assignment
    ## (avoids package check note: "no visible binding for global variable 's'")
    s <- NULL

  }

  ## Update models with any supplied data
  if (!is.null(d)) {
    m <- rMapply(function(i) update(i, data = d), m, SIMPLIFY = FALSE)
  }

  ## Set up parallel processing
  if (p != "no") {

    ## No. cores to use
    if (is.null(nc)) nc <- parallel::detectCores()

    if (p == "snow") {

      ## Create local cluster using system cores
      if (is.null(cl)) {
        cl <- parallel::makeCluster(getOption("cl.cores", nc))
      }

      ## Export required objects/functions to cluster
      ## (search global env for objects in model call(s))
      P <- function(...) paste(..., collapse = " ")
      mc <- P(unlist(rMapply(function(i) P(getCall(i)), m)))
      o <- unlist(lapply(search(), ls))
      o <- o[sapply(o, function(i) grepl(i, mc, fixed = TRUE))]
      o <- c(o, ls("package:semEff"))
      parallel::clusterExport(cl, o)

    }

  }

  ## Generate a seed for bootstrapping
  if (is.null(seed)) {
    set.seed(NULL)
    seed <- sample(10000:99999, 1)
  }

  ## Function to bootstrap coefficients
  bootCoeff <- function(m, w) {

    ## Data to resample (x)
    ## (mixed model: x = highest-level random effect)
    d <- getData(m, subset = TRUE, merge = TRUE, envir = env)
    if (!mer2) x <- if (mer) unique(d[re]) else d

    ## Bootstrap statistic (s)
    stat <- if (!mer2) {
      function(x, i) {
        xi <- if (mer) {
          do.call(rbind, lapply(x[i, ], function(j) {
            d[d[, re] == j, ]
          }))
        } else x[i, ]
        do.call(stdCoeff, c(list(m, w, xi), a))
      }
    } else {
      function(x) do.call(stdCoeff, c(list(x), a))
    }
    s <- if (catch.err) {
      function(...) tryCatch(stat(...), error = function(e) NA)
    } else stat
    if (mer2) assign("s", s, env)

    ## Perform bootstrap
    B <- if (!mer2) {
      set.seed(seed)
      boot::boot(x, s, R, parallel = p, ncpus = nc, cl = cl)
    } else {
      if (isList(m)) {

        ## Bootstrap
        B <- lapply(m, bootMer2)

        ## Weighted average of coefs
        bn <- a$term.names
        b <- lapply(B, "[[", 1)
        b <- avgEst(b, w, bn)

        ## Weighted average of bootstrapped coefs
        bb <- t(sapply(1:R, function(i) {
          bb <- lapply(B, function(j) j$t[i, ])
          avgEst(bb, w, bn)
        }))
        if (ncol(bb) == R) bb <- t(bb)

        ## Output new boot object
        B <- B[[1]]
        B$t0 <- b; B$t <- bb
        B$call <- NULL; B$statistic <- NULL; B$mle <- NULL
        B

      } else bootMer2(m)
    }

    ## Throw warning if any model fits produced errors
    n.err <- sum(apply(rbind(B$t0, B$t), 1, function(i) any(is.na(i))))
    if (n.err > 0)
      warning(paste(n.err, "or more model fit(s) or estimation(s) failed. NA's reported/generated."))

    ## Set attributes and output
    colnames(B$t) <- names(B$t0)
    attributes(B$t)[c("sim", "seed", "n")] <- c(B$sim, seed, nrow(B$data))
    B

  }

  ## Calculate bootstrapped coefficients
  BC <- rMapply(bootCoeff, m, w, SIMPLIFY = FALSE)

  # Add bootstrapped correlated errors
  if (any.ce) {

    ## Data used to fit models
    if (is.null(d)) d <- getData(m, merge = TRUE)
    obs <- rownames(d)

    ## Function to get resids/avg. resids from model/boot obj./list
    res <- function(x, w = NULL) {
      if (isList(x)) {
        if (all(sapply(x, isBoot))) {
          r <- lapply(x, "[[", 1)
          r <- avgEst(r, w)
          rb <- t(sapply(1:R, function(i) {
            rb <- lapply(x, function(j) j$t[i, ])
            avgEst(rb, w)
          }))
          list(r, rb)
        } else {
          r <- lapply(x, resid)
          avgEst(r, w)
        }
      } else {
        if (isBoot(x)) list(x$t0, x$t) else resid(x)
      }
    }

    ## Calculate bootstrapped correlated errors
    BCE <- Map(function(i, j) {

      ## Model/list 1
      m1 <- m[[i[1]]]; w1 <- w[[i[1]]]
      d1 <- getData(m1, subset = TRUE, merge = TRUE, envir = env)

      ## Model/list 2
      m2 <- m[[i[2]]]; w2 <- w[[i[2]]]
      d2 <- getData(m2, subset = TRUE, merge = TRUE, envir = env)

      ## Data to resample (x)
      ## (mixed model: x = highest-level random effect)
      o <- obs %in% rownames(d1) & obs %in% rownames(d2)
      if (!all(o)) d <- d[o, ]
      if (!mer2) {
        x <- if (mer) unique(d[re]) else d
      } else {
        if (!all(o)) {
          m1 <- rMapply(function(i) update(i, data = d), m1, SIMPLIFY = FALSE)
          m2 <- rMapply(function(i) update(i, data = d), m2, SIMPLIFY = FALSE)
        }
      }

      ## Bootstrap statistic (s)
      stat <- if (!mer2) {
        function(x, i) {
          xi <- if (mer) {
            do.call(rbind, lapply(x[i, ], function(j) {
              d[d[, re] == j, ]
            }))
          } else x[i, ]
          m1 <- rMapply(function(i) update(i, data = xi), m1, SIMPLIFY = FALSE)
          m2 <- rMapply(function(i) update(i, data = xi), m2, SIMPLIFY = FALSE)
          r1 <- res(m1, w1)
          r2 <- res(m2, w2)
          cor(r1, r2)
        }
      } else res
      s <- if (catch.err) {
        e <- function(e) if (mer2) rep(NA, nrow(d)) else NA
        function(...) tryCatch(stat(...), error = e)
      } else stat
      if (mer2) assign("s", s, env)

      ## Perform bootstrap
      B <- if (!mer2) {
        set.seed(seed)
        boot::boot(x, s, R, parallel = p, ncpus = nc, cl = cl)
      } else {

        ## Bootstrapped resids for model 1
        B1 <- rMapply(bootMer2, m1, SIMPLIFY = FALSE)
        R1 <- res(B1)
        r1 <- R1[[1]]; rb1 <- R1[[2]]

        ## Bootstrapped resids for model 2
        B2 <- rMapply(bootMer2, m2, SIMPLIFY = FALSE)
        R2 <- res(B2)
        r2 <- R2[[1]]; rb2 <- R2[[2]]

        ## Correlate (and add to new boot object)
        B <- if (isList(B1)) B1[[1]] else B1
        B$t0 <- cor(r1, r2)
        B$t <- matrix(sapply(1:R, function(i) {
          cor(rb1[i, ], rb2[i, ])
        }))
        B$call <- NULL; B$statistic <- NULL; B$mle <- NULL
        B

      }

      ## Throw a warning if any model fits produced errors
      n.err <- sum(apply(rbind(B$t0, B$t), 1, function(i) any(is.na(i))))
      if (n.err > 0)
        warning(paste(n.err, "or more model fit(s) or estimation(s) failed. NA's reported/generated."))

      ## Set attributes and output
      names(B$t0) <- j; colnames(B$t) <- j
      attributes(B$t)[c("sim", "seed", "n")] <- c(B$sim, seed, nrow(B$data))
      B

    }, cv, names(cv))

    ## Append
    if (!isList(BC)) BC <- list(BC)
    BC <- c(BC, BCE)

  }

  ## Output results
  if (p == "snow") parallel::stopCluster(cl)
  set.seed(NULL)
  BC

}


#' @title Bootstrap Confidence Intervals
#' @description Calculate confidence intervals using bootstrapped model
#'   estimates.
#' @param m A fitted model object of class \code{lm}, \code{glm}, or
#'   \code{merMod}, or a named list/nested list of such objects. \code{m} can
#'   also be a boot object(s)/list (class \code{"boot"}), corresponding to
#'   bootstrapped estimates from fitted models.
#' @param conf A numeric value specifying the confidence level for the
#'   intervals.
#' @param type The type of confidence interval to return (defaults to
#'   \code{"bca"} - see Details). See \code{\link[boot]{boot.ci}} for further
#'   specification details.
#' @param digits The number of significant digits to return for numeric values.
#' @param bci.arg A named list of additional arguments to \code{boot.ci}, which
#'   should not include argument \code{index}.
#' @param ... Arguments to \code{bootSEM}.
#' @details This function is essentially a wrapper for \code{boot.ci} from the
#'   \pkg{boot} package, which will return confidence intervals of the specified
#'   type and level calculated from bootstrapped model estimates. If a model or
#'   models are supplied to \code{m}, bootstrapping will first be performed via
#'   \code{bootSEM}. Estimates for which the confidence intervals do not contain
#'   zero are highlighted with a star.
#'
#'   Nonparametric bias-corrected and accelerated confidence intervals
#'   (BC\emph{a}, Efron 1987) are calculated by default, which provide the most
#'   accurate coverage across the widest variety of bootstrap sampling
#'   distributions (Puth \emph{et al.} 2015). They will, however, be
#'   inappropriate for parametric resampling - due to the requirement for
#'   empirical influence values calculated from the data
#'   (\url{http://bit.ly/2UzmA7Y}) - in which case the default will be set to
#'   the boostrap percentile method instead (\code{"perc"}).
#'
#'   NOTE: All bootstrapped confidence intervals will underestimate the true
#'   nominal coverage to some extent when sample size is small (Chernick &
#'   Labudde 2009), so the appropriate caution should be exercised in
#'   interpretation in these cases. Comparison of different interval types may
#'   be informative. For example, normal-theory based intervals may outperform
#'   bootstrap percentile methods when n < 34 (Hesterberg 2015). Ultimately
#'   though, the bootstrap is not a solution to small sample size
#'   (\url{http://bit.ly/2GpfLMn}).
#' @return A data frame of the original estimates and the bootstrapped
#'   confidence intervals, or a list of same.
#' @references Chernick, M. R., & Labudde, R. A. (2009). Revisiting Qualms about
#'   Bootstrap Confidence Intervals. \emph{American Journal of Mathematical and
#'   Management Sciences}, \strong{29}(3–4), 437–456. \url{https://doi.org/c8zv}
#'
#'   Efron, B. (1987). Better Bootstrap Confidence Intervals. \emph{Journal of
#'   the American Statistical Association}, \strong{82}(397), 171–185.
#'   \url{https://doi.org/gfww2z}
#'
#'   Hesterberg, T. C. (2015). What Teachers Should Know About the Bootstrap:
#'   Resampling in the Undergraduate Statistics Curriculum. \emph{The American
#'   Statistician}, \strong{69}(4), 371–386. \url{https://doi.org/gd85v5}
#'
#'   Puth, M.-T., Neuhäuser, M., & Ruxton, G. D. (2015). On the variety of
#'   methods for calculating confidence intervals by bootstrapping.
#'   \emph{Journal of Animal Ecology}, \strong{84}(4), 892–897.
#'   \url{https://doi.org/f8n9rq}
#' @seealso \code{\link[boot]{boot.ci}}, \code{\link[semEff]{bootSEM}}
#' @examples
#' ## CI's from bootstrapped SEM
#' bootCI(Shipley.SEM.boot)
#'
#' ## From original SEM
#' ## (not usually recommended - better to use saved boot objects)
#'
#' \dontrun{
#'
#' system.time(
#'   Shipley.SEM.CI <- bootCI(Shipley.SEM, ran.eff = "site", seed = 53908)
#' )
#' }
#' @export
bootCI <- function(m, conf = 0.95, type = "bca", digits = 3, bci.arg = NULL,
                   ...) {

  ## Create boot object(s) from model(s) if necessary
  is.B <- all(unlist(rMapply(isBoot, m)))
  B <- if (!is.B) bootSEM(m, ...) else m

  ## Function
  bootCI <- function(B) {

    ## Change default CI type for parametric bootstrapping
    if (B$sim == "parametric" && type == "bca") {
      message("Percentile confidence intervals used for parametric bootstrap samples.")
      type <- "perc"
    }

    ## Calculate confidence intervals
    e <- B$t0
    ci <- sapply(1:length(e), function(i) {
      if (!is.na(e[i])) {
        if (e[i] != 0) {
          ci <- do.call(
            boot::boot.ci, c(list(B, conf, type, i), bci.arg)
          )
          tail(as.vector(ci[[4]]), 2)
        } else c(0, 0)
      } else c(NA, NA)
    })

    ## Combine estimates and CI's into table (add significance stars)
    e <- data.frame(
      rbind(e, ci),
      row.names = c("Estimate", "Lower CI", "Upper CI"),
      check.names = FALSE
    )
    e <- round(e, digits)
    stars <- t(data.frame(sapply(e, function(i) {
      if (!any(is.na(i))) {
        if (all(i > 0) || all(i < 0)) "*" else ""
      } else ""
    })))
    e <- format(e, nsmall = digits)
    e <- rbind(e, " " = stars)
    attr(e, "seed") <- attr(B, "seed")
    e

  }

  ## Apply recursively
  rMapply(bootCI, B, SIMPLIFY = FALSE)

}


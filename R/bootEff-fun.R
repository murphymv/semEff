

#' @title Bootstrap Effects
#' @description Bootstrap model effects (standardised coefficients) and optional
#'   SEM correlated errors.
#' @param mod A fitted model object, or a list or nested list of such objects.
#' @param R Number of bootstrap samples to generate.
#' @param seed Seed for the random number generator. If not provided, a random
#'   five-digit integer is used (see Details).
#' @param type The type of bootstrapping to perform. Can be
#'   \code{"nonparametric"} (default), \code{"parametric"}, or
#'   \code{"semiparametric"} (the last two currently only for mixed models, via
#'   \code{bootMer}).
#' @param ran.eff For nonparametric bootstrapping of mixed models, the name of
#'   the (highest-level) random effect to resample (see Details).
#' @param cor.err Optional, names of SEM correlated errors to be bootstrapped.
#'   Should be of the form: \code{c("mod1 ~~ mod2", "mod3 ~~ mod4", ...)}
#'   (spaces optional), with names matching model names.
#' @param catch.err Logical, should errors generated during model fitting or
#'   estimation be caught and \code{NA} returned for estimates? If \code{FALSE},
#'   any such errors will cause the function to exit.
#' @param parallel The type of parallel processing to use. Can be one of
#'   \code{"snow"}, \code{"multicore"}, or \code{"no"} (for none).
#' @param ncpus Number of system cores to use for parallel processing. If
#'   \code{NULL} (default), all available cores are used.
#' @param cl Optional cluster to use if \code{parallel = "snow"}. If \code{NULL}
#'   (default), a local cluster is created using the specified number of cores.
#' @param bM.arg A named list of any additional arguments to \code{bootMer}.
#' @param ... Arguments to \code{\link[semEff]{stdEff}}.
#' @details \code{bootEff} uses the \code{\link[boot]{boot}} function
#'   (primarily) to bootstrap standardised effects from a fitted model or list
#'   of models (calculated using \code{stdEff}). Bootstrapping is typically
#'   nonparametric, i.e. model effects are calculated from data where the rows
#'   have been randomly sampled with replacement. For confidence intervals,
#'   10,000 such resamples should provide accurate coverage in most situations,
#'   with fewer sufficing in some cases. To ensure that data is resampled in the
#'   same way across individual bootstrap operations within the same run (e.g.
#'   models in a list), the same seed is set per operation, with the value saved
#'   as an attribute to the matrix of bootstrapped values (for reproducibility).
#'   The seed can either be user-supplied or a randomly-generated five-digit
#'   number (default), and is always re-initialised on exit (i.e.
#'   \code{set.seed(NULL)}).
#'
#'   Where \code{weights} are specified, bootstrapped effects will be a weighted
#'   average across the set of candidate models for each response variable,
#'   calculated after each model is first refit to the resampled dataset
#'   (specifying \code{weights = "equal"} will use a simple average instead -
#'   see \code{\link[semEff]{avgEst}}). If no weights are specified and
#'   \code{mod} is a nested list of models, the function will throw an error, as
#'   it will be expecting weights for a presumed model averaging scenario. If
#'   instead the user wishes to bootstrap each individual model, they should
#'   recursively apply the function using \code{rMapply} (remember to set a
#'   seed).
#'
#'   Where names of models with correlated errors are specified to
#'   \code{cor.err}, the function will also return bootstrapped Pearson
#'   correlated errors (weighted residuals) for those models. If \code{weights}
#'   are supplied and \code{mod} is a nested list, residuals will first be
#'   averaged across candidate models. If any two models (or candidate sets)
#'   with correlated errors were fit to different subsets of data observations,
#'   both models/sets are first refit to data containing only the observations
#'   in common.
#'
#'   For nonparametric bootstrapping of mixed models, resampling should occur at
#'   the group-level, as individual observations are not independent. The name
#'   of the random effect to resample must be supplied to \code{ran.eff}. For
#'   nested random effects, this should be the highest-level group (Davison &
#'   Hinkley 1997, Ren \emph{et al.} 2010). This form of resampling will result
#'   in differently-sized datasets if observations are unbalanced across groups;
#'   however this should not generally be an issue, as the number of independent
#'   units (groups), and hence the 'degrees of freedom', remains
#'   \href{https://stats.stackexchange.com/questions/46965/bootstrapping-unbalanced-clustered-data-non-parametric-bootstrap}{unchanged}.
#'
#'   For mixed models with
#'   \href{https://stats.stackexchange.com/questions/228800/crossed-vs-nested-random-effects-how-do-they-differ-and-how-are-they-specified}{non-nested
#'   random effects}, nonparametric resampling will not be appropriate. In these
#'   cases, (semi-)parametric bootstrapping can be performed instead via
#'   \code{\link[lme4]{bootMer}} in the \pkg{lme4} package (with additional
#'   arguments passed to that function as necessary). NOTE: As \code{bootMer}
#'   takes only a fitted model as its first argument, any model averaging is
#'   calculated 'post-hoc' using the estimates in boot objects for each
#'   candidate model, rather than during the bootstrapping process itself (i.e.
#'   the default procedure via \code{boot}). Results are then returned in a new
#'   boot object for each response variable or correlated error estimate.
#'
#'   If supplied a list containing both mixed and non-mixed models,
#'   \code{bootEff} with nonparametric bootstrapping will still work and will
#'   treat all models as mixed for resampling (with a warning). This is probably
#'   a relatively rare scenario, but may occur where the user decides that
#'   non-mixed models perform similarly and/or cause less issues than their
#'   mixed counterparts for at least some response variables (e.g. where random
#'   effect variance estimates are at or near zero). If nonparametric
#'   bootstrapping is not used however, an error will occur, as \code{bootMer}
#'   will only accept mixed models.
#'
#'   Parallel processing is used by default via the \pkg{parallel} package and
#'   option \code{parallel = "snow"} (and is generally recommended), but users
#'   can specify the type of parallel processing to use, or none. If
#'   \code{"snow"}, a cluster of workers is created using \code{makeCluster},
#'   and the user can specify the number of system cores to incorporate in the
#'   cluster (defaults to all available). \code{bootEff} then exports all
#'   required objects and functions to this cluster using \code{clusterExport},
#'   after performing a (rough) match of all objects and functions in the
#'   current global environment to those referenced in the model call(s). Users
#'   should load any required external packages prior to calling the function.
#'
#' @note Bootstrapping mixed (or indeed any other) models may take a very long
#'   time when the number of replicates, observations, parameters, and/or models
#'   is high. To decrease processing time, it may be worth trying different
#'   optimizers and/or other options to generate faster estimates (always check
#'   results).
#' @return An object of class \code{"boot"} containing the bootstrapped effects,
#'   or a list/nested list of such objects.
#' @references Burnham, K. P., & Anderson, D. R. (2002). \emph{Model Selection
#'   and Multimodel Inference: A Practical Information-Theoretic Approach} (2nd
#'   ed.). New York: Springer-Verlag. Retrieved from
#'   \url{https://www.springer.com/gb/book/9780387953649}
#'
#'   Davison, A. C., & Hinkley, D. V. (1997). \emph{Bootstrap Methods and their
#'   Application}. Cambridge University Press.
#'
#'   Ren, S., Lai, H., Tong, W., Aminzadeh, M., Hou, X., & Lai, S. (2010).
#'   Nonparametric bootstrapping for hierarchical data. \emph{Journal of Applied
#'   Statistics}, \strong{37}(9), 1487–1498. \url{https://doi.org/dvfzcn}
#' @examples
#' # Bootstrap Shipley SEM (test)
#' # (set 'site' as group for resampling - highest-level random effect)
#' bootEff(Shipley.SEM, ran.eff = "site", R = 1, parallel = "no")
#'
#' # Estimates (use saved boot object, 10000 resamples)
#' lapply(Shipley.SEM.Boot, "[[", 1)  # original
#' lapply(Shipley.SEM.Boot, function(i) head(i$t))  # bootstrapped
#' @export
bootEff <- function(mod, R, seed = NULL, type = "nonparametric", ran.eff = NULL,
                    cor.err = NULL, catch.err = TRUE, parallel = "snow",
                    ncpus = NULL, cl = NULL, bM.arg = NULL, ...) {

  m <- mod; re <- ran.eff; ce <- cor.err; p <- parallel; nc <- ncpus
  if (re == "crossed")
    warning("Use of 'crossed' for 'ran.eff' is deprecated; specify the bootstrap 'type' argument in future.")

  # Arguments to stdEff
  a <- list(...)

  # Weights (for model averaging)
  w <- a$weights; a$weights <- NULL
  if (isList(m) && (is.null(w) || all(w == "equal"))) {
    if (is.null(w) && any(sapply(m, isList)))
      stop("'weights' must be supplied for model averaging (or specify 'equal').")
    w <- lapply(m, function(i) w)
  }

  # Refit model(s) with any supplied data
  upd <- function(m, d) {
    upd <- function(m) eval(update(m, data = d, evaluate = FALSE))
    rMapply(upd, m, SIMPLIFY = FALSE)
  }
  d <- a$data; a$data <- NULL
  if (!is.null(d)) m <- upd(m, d)

  # Environment to look for model data
  main.env <- environment()
  env <- a$env; a$env <- NULL
  if (is.null(env)) env <- parent.frame()
  if (!is.null(d)) env <- main.env

  # Mixed models?
  mer <- unlist(rMapply(isMer, m))
  if (any(mer) && !all(mer))
    warning("Mixed and non-mixed models together in list. Resampling will treat all models as mixed.")
  mer <- any(mer)
  mer2 <- isTRUE(if (mer) {
    if (re == "crossed" && type == "nonparametric") type <- "parametric"
    pb <- type %in% c("parametric", "semiparametric")
    if (!pb && is.null(re))
      stop("Name of random effect to resample must be specified to 'ran.eff' (or use parametric bootstrapping).")
    mer && pb
  })
  if (mer2) {

    # Modified bootMer function
    bootMer2 <- function(...) {

      # Set up function call with specified arguments
      C <- match.call()
      n <- length(C)
      a <- c(list(FUN = s, nsim = R, seed = NULL, type = type, parallel = p,
                  ncpus = nc, cl = cl), bM.arg)
      for (i in 1:length(a)) {
        C[n + i] <- a[i]
        names(C)[n + i] <- names(a)[i]
      }
      C[[1]] <- as.name("bootMer")

      # Run
      set.seed(seed)
      eval.parent(C)

    }

    # Create boot statistic object for later assignment
    # (avoids package check note: "no visible binding for global variable 's'")
    s <- NULL

  }

  # Names of models with correlated errors (list)
  any.ce <- isList(m) && !is.null(ce)
  if (any.ce) {
    cv <- sapply(ce, function(i) {
      gsub(" ", "", unlist(strsplit(i, "~~")))
    }, simplify = FALSE)
    cv <- cv[sapply(cv, function(i) {
      all(i %in% names(m)) && all(i %in% names(w))
    })]
    if (length(cv) < 1)
      stop("Names of variable(s) with correlated errors missing from list of models and/or weights.")
  }

  # Set up parallel processing
  if (p != "no") {

    # No. cores to use
    if (is.null(nc)) nc <- parallel::detectCores()

    if (p == "snow") {

      # Create local cluster using system cores
      if (is.null(cl)) {
        cl <- parallel::makeCluster(getOption("cl.cores", nc))
      }

      # Export required objects/functions to cluster
      # (search global env for objects in model call(s))
      P <- function(...) paste(..., collapse = " ")
      mc <- P(unlist(rMapply(function(i) P(getCall(i)), m)))
      o <- unlist(lapply(search(), ls))
      o <- o[sapply(o, function(i) grepl(i, mc, fixed = TRUE))]
      o <- c(o, ls("package:semEff"))
      parallel::clusterExport(cl, o)

    }

  }

  # Generate a seed for bootstrapping
  if (is.null(seed)) {
    set.seed(NULL)
    seed <- sample(10000:99999, 1)
  }

  # Function to bootstrap effects
  bootEff <- function(m, w) {

    # Data to resample (x)
    # (mixed model: x = highest-level random effect)
    d <- getData(m, subset = TRUE, merge = TRUE, env = env)
    if (!mer2) x <- if (mer) unique(d[re]) else d

    # Bootstrap statistic (s)
    stat <- if (!mer2) {
      function(x, i) {
        xi <- if (mer) {
          do.call(rbind, lapply(x[i, ], function(j) {
            d[d[, re] == j, ]
          }))
        } else x[i, ]
        do.call(stdEff, c(list(m, w, xi), a))
      }
    } else {
      function(x) do.call(stdEff, c(list(x), a))
    }
    s <- if (catch.err) {
      function(...) tryCatch(stat(...), error = function(e) NA)
    } else stat
    if (mer2) assign("s", s, main.env)

    # Perform bootstrap
    B <- if (!mer2) {
      set.seed(seed)
      boot::boot(x, s, R, parallel = p, ncpus = nc, cl = cl)
    } else {
      if (isList(m)) {

        # Bootstrap
        B <- lapply(m, bootMer2)

        # Weighted average of effects
        bn <- a$term.names
        b <- lapply(B, "[[", 1)
        b <- avgEst(b, w, bn)

        # Weighted average of bootstrapped effects
        bb <- t(sapply(1:R, function(i) {
          bb <- lapply(B, function(j) j$t[i, ])
          avgEst(bb, w, bn)
        }))
        if (ncol(bb) == R) bb <- t(bb)

        # Output new boot object
        B <- B[[1]]
        B$t0 <- b; B$t <- bb
        B$call <- NULL; B$statistic <- NULL; B$mle <- NULL
        B

      } else bootMer2(m)
    }

    # Throw warning if any model fits produced errors
    n.err <- sum(apply(rbind(B$t0, B$t), 1, function(i) any(is.na(i))))
    if (n.err > 0)
      warning(paste(n.err, "model fit(s) or parameter estimation(s) failed. NAs reported/generated."))

    # Set attributes and output
    colnames(B$t) <- names(B$t0)
    attributes(B$t)[c("sim", "seed", "n")] <- c(B$sim, seed, nrow(B$data))
    B

  }

  # Bootstrapped effects
  BE <- rMapply(bootEff, m, w, SIMPLIFY = FALSE)

  # Add bootstrapped correlated errors
  if (any.ce) {

    # Data used to fit models
    if (is.null(d)) d <- getData(m, merge = TRUE, env = env)
    obs <- rownames(d)

    # Function to get (weighted) resids/avg. resids from model/boot obj./list
    res <- function(x, w = NULL) {
      f <- function(m) {
        if (!isGls(m) && !isGlm(m)) resid(m, "deviance") else resid(m)
      }
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
          r <- lapply(x, f)
          avgEst(r, w)
        }
      } else {
        if (isBoot(x)) list(x$t0, x$t) else f(x)
      }
    }

    # Bootstrapped correlated errors
    BCE <- Map(function(i, j) {

      # Model(s), weights, and data 1
      m1 <- m[[i[1]]]; w1 <- w[[i[1]]]
      d1 <- getData(m1, subset = TRUE, merge = TRUE, env = env)

      # Model(s), weights, and data 2
      m2 <- m[[i[2]]]; w2 <- w[[i[2]]]
      d2 <- getData(m2, subset = TRUE, merge = TRUE, env = env)

      # Data to resample (x)
      # (mixed model: x = highest-level random effect)
      o <- obs %in% rownames(d1) & obs %in% rownames(d2)
      if (!all(o)) d <- d[o, ]
      if (!mer2) {
        x <- if (mer) unique(d[re]) else d
      } else {
        if (!all(o)) {m1 <- upd(m1, d); m2 <- upd(m2, d)}
      }

      # Bootstrap statistic (s)
      stat <- if (!mer2) {
        function(x, i) {
          xi <- if (mer) {
            do.call(rbind, lapply(x[i, ], function(j) {
              d[d[, re] == j, ]
            }))
          } else x[i, ]
          r1 <- res(upd(m1, xi), w1)
          r2 <- res(upd(m2, xi), w2)
          cor(r1, r2)
        }
      } else res
      s <- if (catch.err) {
        e <- function(e) if (mer2) rep(NA, nrow(d)) else NA
        function(...) tryCatch(stat(...), error = e)
      } else stat
      if (mer2) assign("s", s, main.env)

      # Perform bootstrap
      B <- if (!mer2) {
        set.seed(seed)
        boot::boot(x, s, R, parallel = p, ncpus = nc, cl = cl)
      } else {

        # Bootstrapped resids for model 1
        B1 <- rMapply(bootMer2, m1, SIMPLIFY = FALSE)
        R1 <- res(B1)
        r1 <- R1[[1]]; rb1 <- R1[[2]]

        # Bootstrapped resids for model 2
        B2 <- rMapply(bootMer2, m2, SIMPLIFY = FALSE)
        R2 <- res(B2)
        r2 <- R2[[1]]; rb2 <- R2[[2]]

        # Correlate (and add to new boot object)
        B <- if (isList(B1)) B1[[1]] else B1
        B$t0 <- cor(r1, r2)
        B$t <- matrix(sapply(1:R, function(i) cor(rb1[i, ], rb2[i, ])))
        B$call <- NULL; B$statistic <- NULL; B$mle <- NULL
        B

      }

      # Throw a warning if any model fits produced errors
      n.err <- sum(apply(rbind(B$t0, B$t), 1, function(i) any(is.na(i))))
      if (n.err > 0)
        warning(paste(n.err, "model fit(s) or parameter estimation(s) failed. NAs reported/generated."))

      # Set attributes and output
      names(B$t0) <- j; colnames(B$t) <- j
      attributes(B$t)[c("sim", "seed", "n")] <- c(B$sim, seed, nrow(B$data))
      B

    }, cv, names(cv))

    # Append to bootstrapped effects
    if (!isList(BE)) BE <- list(BE)
    BE <- c(BE, BCE)

  }

  # Output results
  if (p == "snow") parallel::stopCluster(cl)
  set.seed(NULL)
  BE

}


#' @title Bootstrap Confidence Intervals
#' @description Calculate confidence intervals from bootstrapped model effects.
#' @param mod A fitted model object. Alternatively, a boot object (class
#'   \code{"boot"}), containing bootstrapped model effects. Can also be a list
#'   or nested list of such objects.
#' @param conf A numeric value specifying the confidence level for the
#'   intervals.
#' @param type The type of confidence interval to return (defaults to
#'   \code{"bca"} - see Details). See \code{\link[boot]{boot.ci}} for further
#'   specification details.
#' @param digits The number of significant digits to return for numeric values.
#' @param bci.arg A named list of any additional arguments to \code{boot.ci},
#'   excepting argument \code{index}.
#' @param ... Arguments to \code{\link[semEff]{bootEff}}.
#' @details This is essentially a wrapper for \code{boot.ci} from the \pkg{boot}
#'   package, returning confidence intervals of the specified type and level
#'   calculated from bootstrapped model effects. If a model or models is
#'   supplied, bootstrapping will first be performed via \code{bootEff}. Effects
#'   for which the confidence intervals do not contain zero are highlighted with
#'   an asterix.
#'
#'   Nonparametric bias-corrected and accelerated confidence intervals
#'   (BC\emph{a}, Efron 1987) are calculated by default, which should provide
#'   the most accurate coverage across a range of bootstrap sampling
#'   distributions (Puth \emph{et al.} 2015). They will, however, be
#'   \href{https://stackoverflow.com/questions/7588388/adjusted-bootstrap-confidence-intervals-bca-with-parametric-bootstrap-in-boot}{inappropriate}
#'   for parametric resampling - in which case the default will be set to the
#'   bootstrap percentile method instead (\code{"perc"}).
#'
#' @note All bootstrapped confidence intervals will tend to underestimate the
#'   true nominal coverage to some extent when sample size is small (Chernick &
#'   Labudde 2009), so the appropriate caution should be exercised in
#'   interpretation in such cases. Comparison of different interval types may be
#'   informative. For example, normal-theory based intervals may outperform
#'   bootstrap percentile methods when n < 34 (Hesterberg 2015). Ultimately
#'   however, the bootstrap is
#'   \href{https://stats.stackexchange.com/questions/112147/can-bootstrap-be-seen-as-a-cure-for-the-small-sample-size}{not
#'    a solution to small sample size}.
#' @return A data frame of the effects and bootstrapped confidence intervals, or
#'   a list or nested list of same.
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
#' @examples
#' # CIs from bootstrapped SEM
#' (Shipley.SEM.CI <- bootCI(Shipley.SEM.Boot))
#'
#' # From original SEM (models)
#' # (not typically recommended - better to use saved boot objects)
#' # system.time(
#' #   Shipley.SEM.CI <- bootCI(Shipley.SEM, ran.eff = "site", seed = 53908)
#' # )
#' @export
bootCI <- function(mod, conf = 0.95, type = "bca", digits = 3, bci.arg = NULL,
                   ...) {

  m <- mod

  # Create boot object(s) from model(s) if necessary
  is.B <- all(unlist(rMapply(isBoot, m)))
  B <- if (!is.B) bootEff(m, ...) else m

  # Function
  bootCI <- function(B) {

    # Change default CI type for parametric bootstrapping
    if (B$sim == "parametric" && type == "bca") {
      message("Percentile confidence intervals used for parametric bootstrap samples.")
      type <- "perc"
    }

    # Calculate confidence intervals
    e <- B$t0
    ci <- sapply(1:length(e), function(i) {
      if (!is.na(e[i])) {
        if (e[i] != 0) {
          ci <- do.call(boot::boot.ci, c(list(B, conf, type, i), bci.arg))
          tail(as.vector(ci[[4]]), 2)
        } else c(0, 0)
      } else c(NA, NA)
    })

    # Combine effects and CIs into table (add significance stars)
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
    attr(e, "ci.conf") <- conf
    attr(e, "ci.type") <- type
    e

  }

  # Apply recursively
  rMapply(bootCI, B, SIMPLIFY = FALSE)

}


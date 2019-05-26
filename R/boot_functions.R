

#' @title Bootstrap a Piecewise Structural Equation Model
#' @description Generate bootstrapped standardised coefficients and correlated
#'   errors for a list of models, comprising a piecewise Structural Equation
#'   Model (SEM).
#' @param sem A piecewise SEM, comprising a named list/nested list of fitted
#'   model objects of class \code{lm}, \code{glm}, or \code{lmerMod}.
#' @param cor.err An optional character vector of the correlated errors in the
#'   SEM. This should be of the form: \code{c("a ~~ b", "c ~~ d", ...)} (spaces
#'   optional), with names corresponding to model names.
#' @param data An optional dataset used to first re-fit the model(s).
#' @param ran.eff For mixed models with nested random effects, the name of the
#'   variable containing the highest-level random effect. If \code{NULL} and
#'   \code{m} is a mixed model(s), crossed random effects will be assumed
#'   instead.
#' @param R Integer, number of bootstrap replicates to generate.
#' @param seed Integer, seed for the random number generator (for reproducible
#'   results).
#' @param catch.err Logical, should errors generated during model fitting or
#'   estimation be caught and \code{NA} returned? If \code{FALSE}, any such
#'   errors will cause the function to exit.
#' @param parallel The type of parallel processing to use. Can be one of
#'   \code{"snow"}, \code{"multicore"}, or \code{"no"} (for none).
#' @param ncpus Integer, number of system cores to use for parallel processing.
#'   If \code{NULL} (default), all available cores are used.
#' @param cl Optional cluster to use if \code{parallel = "snow"}. If \code{NULL}
#'   (default), a local cluster is created using the specified number of cores.
#' @param bM.arg A named list of additional arguments to \code{bootMer}.
#' @param ... Arguments to \code{stdCoeff}.
#' @details \code{BootSEM} uses the \code{boot} function (primarily) to
#'   bootstrap standardised coefficients from a fitted model or list of models
#'   (calculated using function \code{stdCoeff}), where bootstrapping is
#'   nonparametric, i.e. coefficients are calculated from data where the rows
#'   have been randomly sampled with replacement. The number of replicates is
#'   set by default to 10,000, which should ensure accurate coverage for
#'   confidence intervals under most scenarios.
#'
#'   Where \code{weights} are specified, bootstrapped coefficients will be a
#'   weighted average across the set of candidate models for each response
#'   variable, calculated after each model is first refitted to the resampled
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
#'   weights are supplied and \code{sem} is a nested list, errors are first
#'   averaged across candidate models (as for coefficients). If any two models
#'   (or candidate sets) with correlated errors were fit to different subsets of
#'   data observations, both models/sets are first refit to data containing only
#'   the common observations.
#'
#'   For mixed models with nested random effects, the highest-level random
#'   effect (only) in the dataset is resampled, a procedure which best retains
#'   the hierarchical structure of the data (Davison & Hinkley 1997, Ren
#'   \emph{et al.} 2010). Lower-level groups or individual data rows are not
#'   themselves resampled, either within or across the higher groups. The name
#'   of the variable containing this random effect must be supplied to
#'   \code{ran.eff}. For non-nested ('crossed') random effects however (assumed
#'   when \code{ran.eff = NULL}), this form of resampling will not be
#'   appropriate, and (semi-)parametric bootstrapping is performed instead via
#'   \code{bootMer} in the \pkg{lme4} package. Users should think carefully
#'   about whether their random effects are truly nested or not (see
#'   \url{http://bit.ly/2K244jK}). NOTE: As \code{bootMer} takes only a fitted
#'   model as its first argument, any model averaging is calculated 'post-hoc'
#'   using the estimates in boot objects for each candidate model, rather than
#'   during the bootstrapping process itself (i.e. the default procedure using
#'   \code{boot} in \code{bootSEM}). Results are then returned in a new boot
#'   object for each response variable or correlated error estimate.
#'
#'   Parallel processing is used by default via the \pkg{parallel} package and
#'   option \code{parallel = "snow"} (and is generally recommended), but users
#'   can specify the type of parallel processing to use, or none. If
#'   \code{parallel = "snow"}, a cluster of workers is created using
#'   \code{makeCluster}, and the user can specify the number of system cores to
#'   incorporate in the cluster (defaults to all available). \code{BootSEM} then
#'   exports all required objects and functions to this cluster using
#'   \code{clusterExport}, after performing a (rough) match of all objects and
#'   functions in the current global environment to those referenced in the
#'   model call(s). Users should attach any required packages prior to calling
#'   the function.
#'
#'   NOTE: Bootstrapping mixed (or any other) models may take a very long time
#'   when the number of replicates, observations, parameters, and/or models is
#'   high. To decrease processing time, it may be worth trying different
#'   optimizers (e.g. \code{nloptwrap}) and/or other options to generate faster
#'   estimates (always check results).
#' @return An object of class \code{boot} containing the bootstrapped
#'   coefficients, or a list/nested list of such objects.
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
#'   Statistics}, \strong{37}(9), 1487–1498.
#'   \url{https://doi.org/10.1080/02664760903046102}
#' @seealso \code{\link[semEff]{stdCoeff}}, \code{\link[semEff]{avgEst}},
#'   \code{\link[stats]{resid}}, \code{\link[boot]{boot}},
#'   \code{\link[lme4]{bootMer}}
#' @examples
bootSEM <- function(sem, cor.err = NULL, data = NULL, ran.eff = NULL, R = 10000,
                    seed = 1, catch.err = TRUE, parallel = "snow", ncpus = NULL,
                    cl = NULL, bM.arg = NULL, ...) {

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
      stop("Argument 'weights' must be specified for model averaging.")
    w <- lapply(m, function(i) w)
  }

  ## Names of models with correlated errors (list)
  any.ce <- isList(m) && !is.null(ce)
  if (any.ce) {
    cv <- sapply(ce, function(i) {
      gsub(" ", "", unlist(strsplit(i, "~~")))
    }, simplify = F)
    cv <- cv[sapply(cv, function(i) {
      all(i %in% names(m)) && all(i %in% names(w))
    })]
    if (length(cv) < 1)
      stop("Names of var(s) with correlated errors missing from list of models and/or weights.")
  }

  ## Mixed model(s)?
  mer <- all(unlist(rMapply(isMerMod, m)))
  mer2 <- mer && is.null(re)
  if (mer2) {

    warning("Non-nested ('crossed') random effects assumed. Parametric bootstrapping used.")

    ## Modified bootMer
    bM <- function(...) {
      set.seed(seed)
      C <- match.call()
      n <- length(C)
      a <- c(list(FUN = s, nsim = R, parallel = p, ncpus = nc, cl = cl),
             bM.arg)
      for (i in 1:length(a)) {
        C[n + i] <- a[i]
        names(C)[n + i] <- names(a)[i]
      }
      C[[1]] <- as.name("bootMer")
      eval.parent(C)
    }

  }

  ## Update models with any supplied data
  if (!is.null(d)) {
    m <- rMapply(function(i) update(i, data = d), m, SIMPLIFY = F)
  }

  ## Function to bootstrap coefficients
  bC <- function(m, w) {

    ## Data to bootstrap (x)
    ## (mixed model: x = highest-level random effect)
    d <- getData(m, subset = T, merge = T, envir = env)
    if (!mer2) {
      x <- if (mer) {
        rn <- as.character(d[, re])
        unique(d[re])
      } else d
    }

    ## Bootstrap statistic (s)
    s <- if (!mer2) {
      if (mer) {
        function(x, i) {
          rni <- as.character(x[i, ])
          i <- unlist(lapply(rni, function(j) which(rn %in% j)))
          do.call(stdCoeff, c(list(m, w, d[i, ]), a))
        }
      } else {
        function(x, i) {
          do.call(stdCoeff, c(list(m, w, x[i, ]), a))
        }
      }
    } else {
      function(x) do.call(stdCoeff, c(list(x), a))
    }
    if (catch.err) s <- tryCatch(s, error = function(e) NA)
    if (mer2) assign("s", s, env)

    ## Perform bootstrap
    B <- if (!mer2) {
      set.seed(seed)
      boot::boot(x, s, R, parallel = p, ncpus = nc, cl = cl)
    } else {
      if (isList(m)) {

        ## Bootstrap
        B <- lapply(m, bM)

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
        B <- B[[1]]; B$t0 <- b; B$t <- bb
        B$call <- NULL; B$statistic <- NULL; B$mle <- NULL
        B

      } else bM(m)
    }

    ## Throw warning if any model fits produced errors
    est <- rbind(B$t0, B$t)
    n.err <- sum(apply(est, 1, function(i) any(is.na(i))))
    if (n.err > 0)
      warning(paste(n.err, "or more model fit(s) or estimation(s) failed. NA's reported/generated."))

    ## Output results
    colnames(B$t) <- names(B$t0)
    B

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
      o <- o[sapply(o, function(i) grepl(i, mc, fixed = T))]
      parallel::clusterExport(cl, c(o, ls("package:semEff")))

    }

  }

  ## Bootstrap coefficients
  BC <- rMapply(bC, m, w, SIMPLIFY = F)

  ## Bootstrap correlated errors
  BCE <- if (any.ce) {

    ## Data used to fit models
    if (is.null(d)) d <- getData(m, merge = T)
    obs <- rownames(d)

    ## Function to get resids/weighted avg from a model/boot obj/list
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
    Map(function(i, j) {

      ## Model/list 1
      m1 <- m[[i[1]]]; w1 <- w[[i[1]]]
      d1 <- getData(m1, subset = T, merge = T, envir = env)

      ## Model/list 2
      m2 <- m[[i[2]]]; w2 <- w[[i[2]]]
      d2 <- getData(m2, subset = T, merge = T, envir = env)

      ## Data to bootstrap (x)
      ## (mixed model: x = highest level random effect)
      z <- obs %in% rownames(d1) & obs %in% rownames(d2)
      if (!all(z)) d <- d[z, ]
      if (!mer2) {
        x <- if (mer) {
          rn <- as.character(d[, re])
          unique(d[re])
        } else d
      } else {
        if (!all(z)) {
          m1 <- rMapply(function(i) update(i, data = d), m1, SIMPLIFY = F)
          m2 <- rMapply(function(i) update(i, data = d), m2, SIMPLIFY = F)
        }
      }

      ## Bootstrap statistic (s)
      s <- if (!mer2) {
        function(x, i) {
          xi <- if (mer) {
            rni <- as.character(x[i, ])
            i <- unlist(lapply(rni, function(j) which(rn %in% j)))
            d[i, ]
          } else x[i, ]
          m1 <- rMapply(function(i) update(i, data = xi), m1, SIMPLIFY = F)
          m2 <- rMapply(function(i) update(i, data = xi), m2, SIMPLIFY = F)
          r1 <- res(m1, w1)
          r2 <- res(m2, w2)
          cor(r1, r2)
        }
      } else res
      if (catch.err) {
        e <- if (mer2) function(e) rep(NA, nrow(d)) else function(e) NA
        s <- tryCatch(s, error = e)
      }
      if (mer2) assign("s", s, env)

      ## Perform bootstrap
      B <- if (!mer2) {
        set.seed(seed)
        boot::boot(x, s, R, parallel = p, ncpus = nc, cl = cl)
      } else {

        ## Bootstrapped resids for model 1
        B1 <- rMapply(bM, m1, SIMPLIFY = F)
        R1 <- res(B1); r1 <- R1[[1]]; rb1 <- R1[[2]]

        ## Bootstrapped resids for model 2
        B2 <- rMapply(bM, m2, SIMPLIFY = F)
        R2 <- res(B2); r2 <- R2[[1]]; rb2 <- R2[[2]]

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
      est <- rbind(B$t0, B$t)
      n.err <- sum(apply(est, 1, function(i) any(is.na(i))))
      if (n.err > 0)
        warning(paste(n.err, "or more model fit(s) or estimation(s) failed. NA's reported/generated."))

      ## Set names and output results
      names(B$t0) <- j
      colnames(B$t) <- j
      B

    }, cv, names(cv))

  }

  ## Output results
  if (p == "snow") parallel::stopCluster(cl)
  if (any.ce) {
    if (!isList(BC)) BC <- list(BC)
    c(BC, BCE)
  } else BC

}

# bootCI(lm(ND ~ 1, dataset2))

# (blah <- bootSEM(models.sem, cor.err, dataset2, term.names = all.terms,
#                  r.squared = T, adj = T, pred = T, R = 10))
# (blah <- bootSEM(models.sem[c("PB", "SB")], cor.err, dataset2, term.names = all.terms,
#                  r.squared = T, adj = T, pred = T, R = 10))
# (blah <- bootSEM(models.sem[c("PB", "SB")], "PB ~~ SB", dataset2, term.names = all.terms,
#                  r.squared = T, adj = T, pred = T, R = 10))
# (blah <- bootSEM(models.sem.top, cor.err, dataset2, weights = models.sem.wts, term.names = all.terms,
#                  r.squared = T, adj = T, pred = T, R = 10))
# (blah <- bootSEM(models.sem.top, cor.err, dataset2, term.names = all.terms,
#                  r.squared = T, adj = T, pred = T, R = 10))
# (blah <- bootSEM(models.sem.top, data = dataset2, term.names = all.terms,
#                  r.squared = T, adj = T, pred = T, R = 10))
# (blah <- bootSEM(list("a" = test.lmm, "b" = test.glmm), "a ~~ b", ran.eff = "ran1",
#                  r.squared = T, adj = T, pred = T, R = 10))
# (blah <- bootSEM(list("a" = test.lmm, "b" = test.glmm), "a ~~ b",
#                  r.squared = T, adj = T, pred = T, R = 10))
# (blah <- bootSEM(list("a" = list(test.lmm, test.glmm), "b" = list(test.lmm, test.glmm)),
#                  "a ~~ b", r.squared = T, adj = T, pred = T, R = 10))
# (blah <- bootSEM(list("a" = test.lmm, "b" = test.glmm), "a ~~ b", weights = c("a" = 1, "b" = 2),
#                  r.squared = T, adj = T, pred = T, R = 10))
# (blah <- bootSEM(list("a" = list(test.lmm, test.glmm), "b" = list(test.lmm, test.glmm)),
#                  "a ~~ b", weights = list("a" = c(1, 2), "b" = c(3, 4)),
#                  r.squared = T, adj = T, pred = T, R = 10))#, term.names = "x1"))


#' @title Bootstrap Confidence Intervals
#' @description Calculate confidence intervals using bootstrapped model
#'   estimates.
#' @param x A fitted model object of class \code{lm}, \code{glm}, or
#'   \code{merMod}, or a named list/nested list of such objects. \code{x} can
#'   also be a boot object(s)/list (class \code{"boot"}), corresponding to
#'   bootstrapped estimates from fitted models.
#' @param conf A numeric value specifying the confidence level for the
#'   intervals.
#' @param type The type of confidence interval to return. Can be one of
#'   \code{c("norm","basic", "stud", "perc", "bca")} (see
#'   \code{\link[boot]{boot.ci}}). Defaults to \code{"bca"} (see Details).
#' @param digits The number of significant digits to return for numeric values.
#' @param bci.arg A named list of additional arguments to \code{boot.ci}, which
#'   should not include argument \code{index}.
#' @param ... Arguments to \code{bootSEM}.
#' @details This function is essentially a wrapper for \code{boot.ci} from the
#'   \pkg{boot} package, which will return confidence intervals of the specified
#'   type and confidence level calculated from bootstrapped model estimates. If
#'   a model or models are supplied to \code{x}, bootstrapping will first be
#'   performed via \code{bootSEM}. Estimates for which the confidence intervals
#'   do not contain zero are highlighted with a star.
#'
#'   Nonparametric bias-corrected and accelerated confidence intervals
#'   (BC\emph{a}, Efron 1987) are calculated by default, as these are likely the
#'   most accurate solution across the widest variety of bootstrap sampling
#'   distributions (Puth \emph{et al.} 2015). They will not be appropriate for
#'   parametric resampling however - due to the requirement for empirical
#'   influence values calculated from the data (see \url{http://bit.ly/2UzmA7Y})
#'   - in which case the default will be set to the boostrap percentile method
#'   instead (\code{"perc"}).
#'
#'   NOTE: All bootstrapped confidence intervals will tend to underestimate the
#'   true nominal coverage to some extent when sample size is small (Chernick &
#'   Labudde 2009), so the appropriate caution should be exercised in
#'   interpretation in these cases. Comparison of different interval types may
#'   be informative. For example, normal theory-based intervals may outperform
#'   bootstrap percentile methods when n < 34 (Hesterberg 2015). Ultimately
#'   however, the bootstrap is not a solution to small sample size
#'   (\url{http://bit.ly/2GpfLMn}).
#' @return A data frame of the original estimates and the bootstrapped
#'   confidence intervals, or a list of same.
#' @references Chernick, M. R., & Labudde, R. A. (2009). Revisiting Qualms about
#'   Bootstrap Confidence Intervals. \emph{American Journal of Mathematical and
#'   Management Sciences}, \strong{29}(3–4), 437–456.
#'   \url{https://doi.org/10.1080/01966324.2009.10737767}
#'
#'   Efron, B. (1987). Better Bootstrap Confidence Intervals. \emph{Journal of
#'   the American Statistical Association}, \strong{82}(397), 171–185.
#'   \url{https://doi.org/10.1080/01621459.1987.10478410}
#'
#'   Hesterberg, T. C. (2015). What Teachers Should Know About the Bootstrap:
#'   Resampling in the Undergraduate Statistics Curriculum. \code{The American
#'   Statistician}, \strong{69}(4), 371–386.
#'   \url{https://doi.org/10.1080/00031305.2015.1089789}
#'
#'   Puth, M.-T., Neuhäuser, M., & Ruxton, G. D. (2015). On the variety of
#'   methods for calculating confidence intervals by bootstrapping.
#'   \emph{Journal of Animal Ecology}, \strong{84}(4), 892–897.
#'   \url{https://doi.org/10.1111/1365-2656.12382}
#' @seealso \code{\link[boot]{boot.ci}}, \code{\link[semEff]{bootSEM}}
#' @examples
bootCI <- function(x, conf = 0.95, type = "bca", digits = 3, bci.arg = NULL,
                   ...) {

  ## Create boot object(s) from model(s) if necessary
  is.B <- all(unlist(rMapply(isBoot, x)))
  B <- if (!is.B) bootSEM(x, ...) else x

  ## Function
  bCI <- function(B) {

    ## Change default CI type for parametric bootstrapping
    if (B$sim == "parametric" && type == "bca") {
      message("Percentile confidence intervals used for parametric bootstrap samples.")
      type <- "perc"
    }

    ## Calculate confidence intervals
    e <- B$t0; eb <- B$t
    ci <- sapply(1:length(e), function(i) {
      if (!is.na(e[i])) {
        if (!all(eb[, i] == 0)) {
          ci <- do.call(boot::boot.ci, c(list(B, conf, type, i), bci.arg))
          tail(as.vector(ci[[4]]), 2)
        } else c(0, 0)
      } else c(NA, NA)
    })

    ## Combine into table (add significance stars)
    e <- data.frame(
      rbind(e, ci),
      row.names = c("Estimate", "Lower CI", "Upper CI"),
      check.names = F
    )
    e <- round(e, digits)
    stars <- t(data.frame(sapply(e, function(i) {
      if (!any(is.na(i))) {
        if (all(i > 0) || all(i < 0)) "*" else ""
      } else ""
    })))
    e <- format(e, nsmall = digits)
    rbind(e, " " = stars)

  }

  ## Apply recursively
  rMapply(bCI, B, SIMPLIFY = F)

}
# system.time(
# blah <- bootSEM(models.sem.top, weights = models.sem.wts, cor.err = cor.err,
#                 data = dataset2, term.names = all.terms, r.squared = T, adj = T,
#                 pred = T, R = 10000)
# )
# bootCI(blah[1])#, bci.arg = list(index = 1))


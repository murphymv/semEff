

#' @title Object Types
#' @keywords internal
#' @description Functions to determine the 'type' of an R object using classes.
#'   Intended largely for convenience and internal use.
#' @param x An R object.
#' @return A logical value.
#' @name Object.Type
NULL
#' @describeIn Object.Type Is object a list (class \code{"list"})?
isList <- function(x) class(x)[1] == "list"
#' @describeIn Object.Type Is object a boot object?
isBoot <- function(x) "boot" %in% class(x)
#' @describeIn Object.Type Is object a fitted model?
isMod <- function(x) {
  any(c("lm", "lmerMod", "glmerMod", "gls", "betareg") %in% class(x))
}
#' @describeIn Object.Type Is object a generalised least squares model?
isGls <- function(x) "gls" %in% class(x)
#' @describeIn Object.Type Is object a beta regression model?
isBet <- function(x) "betareg" %in% class(x)
#' @describeIn Object.Type Is object a generalised linear model (i.e. uses a
#'   link function)?
isGlm <- function(x) any(c("glm", "glmerMod", "betareg") %in% class(x))
#' @describeIn Object.Type Is object a mixed model (class \code{"merMod"})?
isMerMod <- function(x) any(c("lmerMod", "glmerMod") %in% class(x))


#' @title Parameter Types
#' @keywords internal
#' @description Functions to determine the presence/absence of certain model
#'   parameter types using their names. Intended largely for convenience and
#'   internal use.
#' @param x A character vector of parameter names (e.g. names of coefficients
#'   from \code{coef} or \code{stdCoeff}).
#' @return A logical vector of the same length as \code{x}.
#' @name Param.Type
NULL
#' @describeIn Param.Type Is parameter an intercept?
isInt <- function(x) x == "(Intercept)"
#' @describeIn Param.Type Is parameter a variable interaction (product term)?
isInx <- function(x) grepl("(?<!:):(?!:)", x, perl = TRUE)
#' @describeIn Param.Type Is parameter an R-squared value?
isR2 <- function(x) grepl("r.squared", x)


#' @title Recursive \code{mapply}
#' @description Recursively apply a function to a list or lists.
#' @param FUN Function to apply.
#' @param ... Object(s) to which \code{FUN} can be applied, or lists of such
#'   objects to iterate over (defined narrowly, as of class \code{"list"}).
#' @param MoreArgs A list of additional arguments to \code{FUN}.
#' @param SIMPLIFY Logical, whether to simplify the results to a vector or
#'   array.
#' @param USE.NAMES Logical, whether to use the names of the first list object
#'   in \code{...} for the output.
#' @details \code{rMapply} recursively applies \code{FUN} to the elements of the
#'   lists in \code{...} via \code{mapply}. If only a single list is supplied,
#'   the function acts like a recursive version of \code{sapply}. The particular
#'   condition that determines if the function should stop recursing is if
#'   either the first or second objects in \code{...} are not of class
#'   \code{"list"}. Thus, unlike \code{mapply}, it will not iterate over
#'   non-list elements in these objects, but instead returns the output of
#'   \code{f(...)}.
#'
#'   This is primarily a convenience function used internally to enable
#'   recursive application of functions to lists or nested lists. Its particular
#'   stop condition for recursing is also designed to either a) act as a wrapper
#'   for \code{FUN} if the first object in \code{...} is not a list, or b) apply
#'   a model averaging operation if the first object is a list and the second
#'   object is a numeric vector (of weights).
#' @return The output of \code{FUN} in a list or nested list, or simplified to a
#'   vector or array (or list of arrays).
#' @seealso \code{\link[base]{mapply}}
#' @export
rMapply <- function(FUN, ..., MoreArgs = NULL, SIMPLIFY = TRUE,
                    USE.NAMES = TRUE) {
  l <- list(...)
  n <- length(l)
  i <- if (n > 0) l[[1]] else l
  j <- if (n > 1) l[[2]] else i
  if (!isList(i) || !isList(j)) {
    do.call(FUN, c(l, MoreArgs))
  } else {
    mapply(
      rMapply, ...,
      MoreArgs = list(FUN = FUN, MoreArgs = MoreArgs, SIMPLIFY = SIMPLIFY,
                      USE.NAMES = USE.NAMES),
      SIMPLIFY = SIMPLIFY, USE.NAMES = USE.NAMES
    )
  }
}


#' @title Parallel \code{sapply}
#' @description Apply a function to a vector using parallel processing.
#' @param X A vector object (numeric, character, or list).
#' @param FUN Function to apply to the elements of \code{X}.
#' @param parallel The type of parallel processing to use. Can be one of
#'   \code{"snow"}, \code{"multicore"}, or \code{"no"} (for none). If none,
#'   \code{sapply} is used instead.
#' @param ncpus Number of system cores to use for parallel processing. If
#'   \code{NULL} (default), all available cores are used.
#' @param cl Optional cluster to use if \code{parallel = "snow"}. If \code{NULL}
#'   (default), a local cluster is created using the specified number of cores.
#' @param add.obj A character vector of any additional object names to be
#'   exported to the cluster for parallel processing. Use if a required object
#'   or function cannot be found.
#' @param ... Arguments to \code{parSapply} or \code{sapply}.
#' @details This is a wrapper for \code{parSapply} from the \pkg{parallel}
#'   package, enabling (potentially) faster processing of a function over a
#'   vector of objects. Parallel processing via option \code{"snow"} (default)
#'   is carried out using a cluster of workers, which is automatically set up
#'   via \code{makeCluster} using all available system cores or a user supplied
#'   number of cores. The function then exports the required objects and
#'   functions to this cluster using \code{clusterExport}, after performing a
#'   (rough) match of all objects and functions in the current global
#'   environment to those referenced in the call to \code{FUN} (and also any
#'   calls in \code{X}). Any additional required objects can be supplied using
#'   \code{add.obj}.
#' @return The output of \code{FUN} in a list, or simplified to a vector or
#'   array.
#' @seealso \code{\link[parallel]{parSapply}}, \code{\link[base]{sapply}}
#' @export
pSapply <- function(X, FUN, parallel = "snow", ncpus = NULL, cl = NULL,
                    add.obj = NULL, ...) {

  p <- parallel; nc <- ncpus; ao <- add.obj

  if (p != "no") {

    ## No. cores to use
    if (is.null(nc)) nc <- detectCores()

    ## Cluster
    if (p == "snow") {

      ## Create local cluster using system cores
      if (is.null(cl)) {
        cl <- parallel::makeCluster(getOption("cl.cores", nc))
      }

      ## Export required objects/functions to cluster
      ## (search global env. for objects in calls to X/FUN)
      P <- function(...) paste(..., collapse = " ")
      xc <- P(unlist(rMapply(function(i) {
        if (isMod(i) || isBoot(i)) P(getCall(i))
      }, X)))
      fa <- P(sapply(match.call(expand.dots = FALSE)$..., deparse))
      fc <- P(xc, enquote(FUN)[2], fa)
      o <- unlist(lapply(search(), ls))
      o <- o[sapply(o, function(i) grepl(i, fc, fixed = TRUE))]
      o <- c("X", o, ao)
      parallel::clusterExport(cl, o, environment())

    }

    ## Run parSapply using cluster and output results
    out <- parallel::parSapply(cl, X, FUN, ...)
    parallel::stopCluster(cl)
    out

  } else sapply(X, FUN, ...)

}


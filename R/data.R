

#' @title Simulated Data from Shipley (2009)
#' @format A data frame with 1900 observations and 9 variables:
#' \describe{
#' \item{site}{a numeric code giving the site from which the observation comes}
#' \item{tree}{a numeric code giving the tree from which the observation comes}
#' \item{lat}{the latitude of the site}
#' \item{year}{the year in which the observation was taken}
#' \item{Date}{the Julian date when the bud burst occurs}
#' \item{DD}{the number of degree days when bud burst occurs}
#' \item{Growth}{the increase in diameter growth of the tree}
#' \item{Survival}{the probability of survival until the next year (used only
#' for the simulation)}
#' \item{Live}{a binary value (1 = tree lived the following winter, 0 = tree
#' died the following winter)}
#' }
#' @source <https://doi.org/c886>
#' @references Shipley, B. (2009). Confirmatory path analysis in a generalized
#'   multilevel context. *Ecology*, **90**(2), 363-368. <https://doi.org/bqd43d>
"shipley"


#' @title Hypothesised SEM from Shipley (2009)
#' @format A list of fitted mixed models of class `"lmer"` and `"glmer"`,
#'   representing structured equations.
#' @references Shipley, B. (2009). Confirmatory path analysis in a generalized
#'   multilevel context. *Ecology*, **90**(2), 363-368. <https://doi.org/bqd43d>
#' @examples
#' # Specification
#' # shipley.sem <- list(
#' #   DD = lme4::lmer(DD ~ lat + (1 | site) + (1 | tree), data = shipley),
#' #   Date = lme4::lmer(Date ~ DD + (1 | site) + (1 | tree), data = shipley),
#' #   Growth = lme4::lmer(Growth ~ Date + (1 | site) + (1 | tree),
#' #                       data = shipley),
#' #   Live = lme4::glmer(Live ~ Growth + (1 | site) + (1 | tree), binomial,
#' #                      data = shipley)
#' # )
"shipley.sem"


#' @title Candidate Model Set from Shipley 'Growth' Model
#' @description A set of hypothetical competing models fit to the same response
#'   variable ('Growth') using the simulated data in Shipley (2009), for which
#'   model estimates can be compared and/or averaged.
#' @format A list of mixed models of class `"lmer"` and `"glmer"`, fit to the
#'   same response variable.
#' @references Shipley, B. (2009). Confirmatory path analysis in a generalized
#'   multilevel context. *Ecology*, **90**(2), 363-368. <https://doi.org/bqd43d>
#' @examples
#' # Specification
#' # shipley.growth <- list(
#' #   lme4::lmer(Growth ~ Date + (1 | site) + (1 | tree), data = shipley),
#' #   lme4::lmer(Growth ~ Date + DD + (1 | site) + (1 | tree), data = shipley),
#' #   lme4::lmer(Growth ~ Date + DD + lat + (1 | site) + (1 | tree),
#' #              data = shipley)
#' # )
"shipley.growth"


#' @title Bootstrapped Estimates for Shipley SEM
#' @description Bootstrapped estimates generated from the hypothesised SEM from
#'   Shipley (2009), using [bootEff()].
#' @format A list of objects of class `"boot"`, representing bootstrapped
#'   estimates from fitted mixed models.
#' @references Shipley, B. (2009). Confirmatory path analysis in a generalized
#'   multilevel context. *Ecology*, **90**(2), 363-368. <https://doi.org/bqd43d>
#' @examples
#' # Specification
#' # shipley.sem.boot <- bootEff(shipley.sem, R = 10000, seed = 53908,
#' #                             ran.eff = "site")
"shipley.sem.boot"


#' @title Effects for Shipley SEM
#' @description SEM effects calculated from bootstrapped estimates of the
#'   hypothesised SEM from Shipley (2009), using [semEff()].
#' @format A list object of class `"semEff"`, containing SEM effects and summary
#'   tables.
#' @references Shipley, B. (2009). Confirmatory path analysis in a generalized
#'   multilevel context. *Ecology*, **90**(2), 363-368. <https://doi.org/bqd43d>
#' @examples
#' # Specification
#' # shipley.sem.eff <- semEff(shipley.sem.boot)
"shipley.sem.eff"


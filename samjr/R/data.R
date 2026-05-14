##' North Sea cod example data
##'
##' Data, configuration and parameter list prepared from the North Sea cod
##' files shipped with the SAM \code{stockassessment} package's nscod test
##' case. The configuration matches the workflow used in
##' \code{SAM/testmore/nscod/script.R}: AR(1) F-correlation, an F-at-age key
##' of \code{0:5}, the \code{keyVarObs} and \code{keyLogFpar} matrices used
##' there, two F-variance groups via \code{keyVarF[1,] = c(0,1,1,1,1,1)},
##' and the catch-scaling block (\code{noScaledYears = 13},
##' \code{keyScaledYears = 1993:2005}). \code{fbarRange} is set to
##' \code{c(2, 4)}.
##'
##' @format
##' \describe{
##' \item{\code{nscodData}}{a list of class \code{sam_data} as returned by
##' \code{\link{setup.sam.data}}.}
##' \item{\code{nscodConf}}{a configuration list as returned by
##' \code{\link{defcon}} with the modifications above.}
##' \item{\code{nscodParameters}}{a list of initial parameter values from
##' \code{\link{defpar}}.}
##' }
##' @seealso \code{\link{sam.fit}}.
##' @name nscodData
##' @aliases nscodConf nscodParameters
##' @keywords datasets
NULL

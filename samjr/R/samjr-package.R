##' samjr: a minimal RTMB re-implementation of SAM
##'
##' Provides a small subset of the SAM \code{stockassessment} user interface
##' (\code{\link{read.ices}}, \code{\link{setup.sam.data}}, \code{\link{defcon}},
##' \code{\link{defpar}}, \code{\link{sam.fit}}) backed by a pure-R RTMB
##' implementation of the state-space age-structured assessment model. Plotting
##' helpers (\code{\link{ssbplot}}, \code{\link{fbarplot}}, \code{\link{recplot}},
##' \code{\link{catchplot}}) mirror the corresponding SAM functions.
##'
##' @name samjr-package
##' @aliases samjr
##' @keywords package
##' @import RTMB
##' @import Matrix
##' @importFrom methods as
##' @importFrom stats nlminb complete.cases median quantile rnorm
##' @importFrom utils read.table
##' @importFrom graphics lines points polygon grid par
##' @importFrom grDevices gray
"_PACKAGE"

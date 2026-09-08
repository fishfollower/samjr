##' samjr: a minimal RTMB re-implementation of SAM
##'
##' Provides a small subset of the SAM \code{stockassessment} user interface
##' (\code{\link{read.ices}}, \code{\link{setup.sam.data}}, \code{\link{defcon}},
##' \code{\link{defpar}}, \code{\link{sam.fit}}) backed by a pure-R RTMB
##' implementation of the state-space age-structured assessment model. Plotting
##' helpers (\code{\link{ssbplot}}, \code{\link{fbarplot}}, \code{\link{recplot}},
##' \code{\link{catchplot}}) mirror the corresponding SAM functions.
##'
##' Also provides deterministic per-recruit and equilibrium reference points
##' (\code{\link{referencepoints}}, \code{\link{perRecruitClosure}}) with
##' delta-method confidence intervals, plus stochastic harvest control rule
##' projections (\code{\link{hcr}}, \code{\link{icesAdviceRule}}), mirroring the
##' reference-point catalogue of SAM's \code{stockassessment} package.
##'
##' @name samjr-package
##' @aliases samjr
##' @keywords package
##' @import RTMB
##' @import Matrix
##' @importFrom methods as
##' @importFrom stats nlminb complete.cases median quantile rnorm na.omit
##' @importFrom utils read.table capture.output
##' @importFrom graphics lines points polygon grid par
##' @importFrom graphics plot.new
##' @importFrom grDevices gray rgb
"_PACKAGE"

## RTMB exposes the entries of the parameter / data list as bare names inside
## the likelihood closure - the static analyser in R CMD check cannot see
## that, so register them here.
utils::globalVariables(c(
  "PF", "PM", "Wc", "Wd", "Wp", "age", "aux", "bhpar",
  "catchWeightModel", "corList", "covType", "cwNobs", "fbarIdx", "fcormode",
  "fixVarToWeight", "fleetDim", "fleetTypes", "idxCor", "isTag", "itrans_rho",
  "keyCatchWeightMean", "keyCatchWeightObsVar", "keyF", "keyIGAR",
  "keyMatureMean", "keyMortalityMean", "keyMortalityObsVar", "keyQ", "keyQpow",
  "keySd", "keyStockWeightMean", "keyStockWeightObsVar", "keyVarFperState",
  "logCW", "logF", "logIGARdist", "logN", "logNM", "logPhiCW", "logPhiMO",
  "logPhiNM", "logPhiSW", "logQ", "logQpow", "logSW", "logScale",
  "logSdLogCW", "logSdLogN", "logSdLogNM", "logSdLogObs", "logSdLogSW",
  "logSdMO", "logSdProcLogCW", "logSdProcLogNM", "logSdProcLogSW",
  "logSdProcLogitMO", "logSRpar", "logXtraSd", "logitMO", "logitRecapturePhi",
  "logitReleaseSurvival", "logsdF", "matureModel", "meanLogCW", "meanLogNM",
  "meanLogSW", "meanLogitMO", "minAge", "minYear", "moNobs", "mortalityModel",
  "nmNobs", "parUS", "predVarObs", "predVarObsLink", "rickerpar",
  "sampleTimes", "scaleIdxByObs", "srmode", "stockWeightModel", "swNobs",
  "tagNscan", "tagR", "tagTypeIdx", "weight", "xtraSdIdxByObs", "year",
  "useCKMR", "usePOP", "useHSP", "ckmrPsi", "ckmrScale", "ckmrPrep",
  "ckmrPOPobs", "ckmrHSPobs", "ckmrNpop", "ckmrNhsp",
  "ckmrEstPsi", "logPsim1", "ckmrSWref"
))

## Internal utilities for the reference-point machinery.
##
## averageBio     - assemble averaged biological inputs over a year window,
##                  honouring the {stock,catch,mortality,mature}Model
##                  smoothing flags (tables.R:344-370).
## selectivityFromLogF - relative selectivity from fit$pl$logF averaged
##                       over selYears and normalised to fbarRange.

##' @keywords internal
##' @noRd
.aveYearsIndex <- function(fit, aveYears){
  yrs <- fit$data$years
  if(is.null(aveYears)) aveYears <- yrs[max(1L, length(yrs) - 9L):length(yrs)]
  idx <- match(aveYears, yrs)
  if(any(is.na(idx)))
    stop("aveYears outside fit$data$years: ",
         paste(aveYears[is.na(idx)], collapse = ", "))
  idx
}

##' Build averaged biology vectors for the equilibrium calculation.
##'
##' Each per-(year,age) input is picked either from the observation matrix
##' on \code{fit$data} (when the corresponding smoothing model is 0) or from
##' the fitted GMRF-smoothed state in \code{fit$pl} (when the smoothing
##' model is 1), then averaged over \code{aveYears}.
##'
##' Returns a list of length-nAge numeric vectors:
##'   stockMeanWeight, catchMeanWeight, natMor, propMat,
##'   landFrac, landMeanWeight, disMeanWeight, propF, propM.
##' @keywords internal
##' @noRd
averageBio <- function(fit, aveYears){
  data <- fit$data; conf <- fit$conf; pl <- fit$pl
  ageRange <- conf$minAge:conf$maxAge
  nY <- length(data$years)
  pickMat <- function(modelFlag, plMat, dataMat, link){
    if(isTRUE(modelFlag == 1L)){
      M <- link(plMat)
      M[seq_len(nY), , drop = FALSE]
    } else dataMat
  }
  sw <- pickMat(conf$stockWeightModel, pl$logSW,   data$stockMeanWeight, exp)
  cw <- pickMat(conf$catchWeightModel, pl$logCW,   data$catchMeanWeight, exp)
  nm <- pickMat(conf$mortalityModel,   pl$logNM,   data$natMor,          exp)
  mo <- pickMat(conf$matureModel,      pl$logitMO, data$propMat,         plogis)
  avg <- function(M) colMeans(M[aveYears, , drop = FALSE], na.rm = TRUE)
  defOrAvg <- function(M, def){
    if(is.null(M)) rep(def, length(ageRange))
    else            avg(M)
  }
  list(stockMeanWeight = avg(sw),
       catchMeanWeight = avg(cw),
       natMor          = avg(nm),
       propMat         = avg(mo),
       landFrac        = defOrAvg(data$landFrac,       1),
       landMeanWeight  = if(is.null(data$landMeanWeight)) avg(cw) else avg(data$landMeanWeight),
       disMeanWeight   = defOrAvg(data$disMeanWeight,  0),
       propF           = defOrAvg(data$propF,          0),
       propM           = defOrAvg(data$propM,          0))
}

##' Relative selectivity from a samjr fit
##'
##' Computes the unitless age-vector \code{aveFa / mean(aveFa, a in fbarRange)},
##' where \code{aveFa} is \code{faytable(fit)} averaged over \code{selYears}.
##' Scaling by \code{exp(logFbar)} reproduces absolute F-at-age. NA entries
##' (ages with no F state) are treated as zero.
##'
##' Matches SAM's \code{.logFtoSel}
##' (\code{deterministic_referencepoints.R:37-41}), using samjr's
##' user-facing F-at-age convention via \code{\link{faytable}}.
##'
##' @keywords internal
##' @noRd
selectivityFromLogF <- function(fit, selYears){
  conf <- fit$conf
  yrs <- fit$data$years
  idx <- match(selYears, yrs)
  if(any(is.na(idx)))
    stop("selYears outside fit$data$years: ",
         paste(selYears[is.na(idx)], collapse = ", "))
  Fmat <- faytable(fit)
  Fmat[is.na(Fmat)] <- 0
  ages <- conf$minAge:conf$maxAge
  aveFa <- colMeans(Fmat[idx, , drop = FALSE])
  fbarAges <- which(ages %in% conf$fbarRange[1]:conf$fbarRange[2])
  scale <- mean(aveFa[fbarAges])
  if(scale <= 0) stop("selectivityFromLogF: average F on fbarRange is zero")
  aveFa / scale
}

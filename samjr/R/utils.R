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

#' @title Make path matrices
#'
#' @description Constructs path matrices for dynamic structural equation model (DSEM)
#'        using a vector of parameters and specification of the DSEM
#'
#' @param beta_p vector parameters.
#' @param model matrix or data.frame with the following columns, and one row
#'   per one-headed or two-headed arrow in the dynamic structural model:
#' \describe{
#' \item{direction}{whether a path coefficient is one-headed (1) or two-headed (2)}
#' \item{lag}{whether the lag associated with a given coefficient}
#' \item{start}{starting value, used when \code{parameter=0}}
#' \item{parameter}{The parameter number from \code{beta_p} associated with a given path}
#' \item{first}{The variable at the tail of a given path}
#' \item{second}{The variable at the head of a given path}
#' }
#' @param times integer-vector of times to use when defining matrices
#' @param variables character-vector listing variables
#'
#' @importFrom Matrix solve Diagonal sparseMatrix drop0 kronecker bandSparse
#' @importFrom RTMB ADoverload AD
#'
#' @details
#' When \code{length(times)} is \eqn{T} and \code{length(variables)} is \eqn{J},
#' \code{make_matrices} returns matrices of dimension \eqn{TJ \times TJ} representing
#' paths among \eqn{vec(\mathbf{X})} where matrix \eqn{\mathbf{X}} has dimension
#' \eqn{T \times J} and \eqn{vec} stacks columns into a single long vector
#'
#' @return
#' A named list of matrices including:
#' \describe{
#' \item{P_kk}{The matrix of interactions, i.e., one-headed arrows}
#' \item{G_kk}{The matrix of exogenous covariance, i.e., two-headed arrows}
#' }
#' @export
dsem_make_matrices <-
function( beta_p,
          model,
          times,
          variables ){

  "c" <- ADoverload("c")
  "[<-" <- ADoverload("[<-")
  model <- as.data.frame(model)
  model[,'parameter'] = as.integer(model[,'parameter'])

  # Combine fixed, estimated, and mapped parameters into vector
  beta_i = rep(0, nrow(model))
  off = which(model[,'parameter'] == 0)
  if( length(off) > 0 ){
    beta_i[off] = as.numeric(model[off,'start'])
  }
  not_off = which(model[,'parameter'] > 0)
  if( length(not_off) > 0 ){
    beta_i[not_off] = beta_p[model[not_off,'parameter']]
  }

  # Loop through paths
  P_kk = drop0(sparseMatrix( i=1, j=1, x=0, dims=rep(length(variables)*length(times),2) ))   # Make with a zero
  #P_kk = AD(P_kk)
  G_kk = (P_kk)
  for( i in seq_len(nrow(model)) ){
    lag = as.numeric(model[i,'lag'])
    # Time-lag matrix
    L_tt = bandSparse( n = length(times),
                       m = length(times),
                       k = -1 * lag,  # positive lag means below diagonal
                       diagonals = list(rep(1,length(times))) )

    # Interaction matrix
    P_jj = sparseMatrix( i = match(model[i,'second'],variables),
                         j = match(model[i,'first'],variables),
                         x = 1,
                         dims = rep(length(variables),2) )

    # Assemble
    tmp_kk = kronecker(P_jj, L_tt)
    if(abs(as.numeric(model[i,'direction']))==1){
      P_kk = P_kk + beta_i[i] * tmp_kk # AD(tmp_kk)
    }else{
      G_kk = G_kk + beta_i[i] * tmp_kk # AD(tmp_kk)
    }
  }

  # Diagonal component
  I_kk = Diagonal(nrow(P_kk))

  # Assemble
  IminusP_kk = AD(I_kk - P_kk)
  invV_kk = AD(G_kk)
  invV_kk@x = AD(1 / G_kk@x^2)
  #Q_kk = t(IminusP_kk) %*% invV_kk %*% IminusP_kk

  out = list(
    "P_kk" = P_kk,
    "G_kk" = G_kk,
    "invV_kk" = invV_kk,
    #"Q_kk" = Q_kk,     # NOT USED
    "IminusP_kk" = IminusP_kk
  )
  return(out)
}

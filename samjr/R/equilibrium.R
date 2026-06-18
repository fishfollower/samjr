## Per-recruit / equilibrium core for samref.
##
## perRecruitClosure(fit, ...) returns a function of (logFbar, recPars) that
## evaluates the deterministic equilibrium for a fixed averaged biology and
## selectivity. The body uses only base-R arithmetic / exp / log / sum /
## cumsum so it is transparently AD-traceable when MakeADFun retapes it for
## the implicit-function delta method later.
##
## Mirrors SAM's perRecruit_D in equilibrium.hpp.

##' Build a per-recruit / equilibrium evaluator from a samjr fit
##'
##' Returns a function \code{perRec(logFbar, recPars)} that produces a list
##' of equilibrium per-recruit quantities matching SAM's \code{PERREC_t}:
##' \code{logYPR, logSPR, logSe, logRe, logYe, logLifeExpectancy,
##' logYearsLost, logDiscYPR, logDiscYe, dSR0}.
##'
##' Averaging of biology (\code{averageBio}) and selectivity
##' (\code{selectivityFromLogF}) happens once at closure-construction time.
##' The returned function depends on \code{logFbar} and \code{recPars} only.
##'
##' @param fit a samjr \code{sam} fit.
##' @param aveYears years over which to average biology (subset of
##'   \code{fit$data$years}). Default: last 10 years of the fit.
##' @param selYears years over which to average selectivity (subset of
##'   \code{fit$data$years}). Default: last year of the fit.
##' @param catchType one of \code{"catch"} (default), \code{"landing"},
##'   \code{"discard"}.
##' @param customSel optional numeric vector to override selectivity. Must
##'   sum to one over the fbar range.
##' @return a function \code{perRec(logFbar, recPars)}.
##' @examples
##' \donttest{
##' data(nscodData); data(nscodConf)
##' fit <- samjr::sam.fit(nscodData, nscodConf,
##'                       samjr::defpar(nscodData, nscodConf), silent = TRUE)
##' pr <- perRecruitClosure(fit, catchType = "landing")
##' pr(log(0.3), recPars(fit))   # per-recruit quantities at F-bar = 0.3
##' }
##' @export
perRecruitClosure <- function(fit,
                              aveYears  = NULL,
                              selYears  = NULL,
                              catchType = c("catch", "landing", "discard"),
                              customSel = NULL){
  catchType <- match.arg(catchType)
  yrs <- fit$data$years
  if(is.null(aveYears)) aveYears <- yrs[max(1L, length(yrs) - 9L):length(yrs)]
  if(is.null(selYears)) selYears <- yrs[length(yrs)]
  aveIdx <- .aveYearsIndex(fit, aveYears)
  selIdx <- .aveYearsIndex(fit, selYears)
  bio <- averageBio(fit, aveIdx)
  sel <- if(!is.null(customSel)) customSel
         else selectivityFromLogF(fit, selYears)
  conf <- fit$conf
  ages <- conf$minAge:conf$maxAge
  if(length(sel) != length(ages))
    stop("selectivity length (", length(sel), ") != nAge (", length(ages), ")")
  M  <- bio$natMor
  sw <- bio$stockMeanWeight
  cw <- bio$catchMeanWeight
  lw <- bio$landMeanWeight
  dw <- bio$disMeanWeight
  lf <- bio$landFrac
  pm <- bio$propMat
  PF <- bio$propF
  PM <- bio$propM
  nAge <- length(ages)
  maxAge <- conf$maxAge
  srCode <- conf$stockRecruitmentModelCode
  function(logFbar, recPars = NULL){
    Fa <- sel * exp(logFbar)
    Z  <- M + Fa
    cumZ <- c(0, cumsum(Z[-nAge]))
    N <- exp(-cumZ)
    N[nAge] <- N[nAge] / (1 - exp(-Z[nAge]))
    Ctot <- Fa / Z * (1 - exp(-Z)) * N
    yieldW <- switch(catchType,
                     catch   = cw,
                     landing = lf * lw,
                     discard = (1 - lf) * dw)
    YPR <- sum(Ctot * yieldW)
    SPRa <- N * exp(-Fa * PF - M * PM) * pm * sw
    SPR <- sum(SPRa)
    yearsLost <- sum(Fa / Z * (1 - exp(-Z)) * N * seq_len(nAge))
    lifeExp <- sum(N * (1 - exp(-Z)) / Z) + (conf$minAge - 1)
    DYPR <- sum(Ctot * ((1 - lf) * dw))
    eps <- 1e-300
    logYPR <- log(YPR + eps)
    logSPR <- log(SPR + eps)
    logYearsLost <- log(yearsLost + eps)
    logLifeExpectancy <- log(min(lifeExp, 10 * maxAge) + eps)
    logDiscYPR <- log(DYPR + eps)
    if(srCode %in% c(1L, 2L) && length(recPars) >= 2L){
      logSe <- srEquilibriumR(srCode, recPars, logSPR)
      logRe <- logSe - logSPR
      logYe <- logRe + logYPR
      logDiscYe <- logRe + logDiscYPR
      dSR0  <- srGradAt0(srCode, recPars)
    } else {
      logSe <- NA_real_
      logRe <- NA_real_
      logYe <- NA_real_
      logDiscYe <- NA_real_
      dSR0  <- NA_real_
    }
    list(logYPR = logYPR, logSPR = logSPR,
         logSe = logSe, logRe = logRe, logYe = logYe,
         logLifeExpectancy = logLifeExpectancy,
         logYearsLost = logYearsLost,
         logDiscYPR = logDiscYPR, logDiscYe = logDiscYe,
         dSR0 = dSR0)
  }
}

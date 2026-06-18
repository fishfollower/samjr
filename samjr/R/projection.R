## Self-contained per-year population projection for samref.
##
## Used by hcr.sam to apply a feedback rule. Independent of samjr's internal
## forecast machinery so samjr stays read-only.

##' @keywords internal
##' @noRd
.initState <- function(fit){
  N <- exp(fit$pl$logN)
  list(N = N[nrow(N), ],
       year0 = max(fit$data$years),
       sigmaR = exp(fit$pl$logSdLogN[1]))
}

##' Advance population by one year under a fixed Fbar
##'
##' Returns the new numbers-at-age (length nAge), the SSB at the start of
##' the year, the catch (in weight) and the predicted recruitment (before
##' adding noise).
##'
##' @keywords internal
##' @noRd
.projectOneYear <- function(N, fval, bio, sel, srCode, recPars, recNoise,
                            maxAgePlusGroup = TRUE){
  nAge <- length(N)
  Fa <- sel * fval
  Z  <- bio$natMor + Fa
  ssb <- sum(N * exp(-Fa * bio$propF - bio$natMor * bio$propM) *
             bio$propMat * bio$stockMeanWeight)
  catchN <- N * Fa / Z * (1 - exp(-Z))
  catchW <- sum(catchN * bio$landFrac * bio$landMeanWeight)
  survival <- N * exp(-Z)
  if(srCode %in% c(1L, 2L) && length(recPars) >= 2L){
    logRPred <- predictLogR(srCode, recPars, log(ssb))
  } else {
    logRPred <- log(N[1])
  }
  newN <- numeric(nAge)
  newN[1] <- exp(logRPred + recNoise)
  if(nAge > 1L) newN[2:nAge] <- survival[1:(nAge - 1L)]
  if(maxAgePlusGroup) newN[nAge] <- newN[nAge] + survival[nAge]
  list(N = newN, ssb = ssb, catch = catchW, logRPred = logRPred)
}

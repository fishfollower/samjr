## SR-model helpers for samref.
##
## samjr supports three recruitment codes (likelihood.R):
##   0  random walk        logR_y = logR_{y-1}
##   1  Ricker             logR   = a + log(S) - exp(b) * S
##   2  Beverton-Holt      logR   = a + log(S) - log(1 + exp(b) * S)
## with rec parameters stored in fit$pl$rickerpar / fit$pl$bhpar of length 2.

##' Properties of a recruitment model
##'
##' Returns the equilibrium / compensatory / finite-max-gradient flags used
##' by the reference-point machinery to decide which RP types are defined for
##' a given SR model code.
##'
##' @param srCode integer SR-model code (0 random walk, 1 Ricker,
##'   2 Beverton-Holt). Other codes are not supported.
##' @return list with components \code{name}, \code{hasEquilibrium},
##'   \code{isCompensatory}, \code{hasFiniteMaxGradient}.
##' @examples
##' srProperties(0L)   # random walk: no equilibrium
##' srProperties(1L)   # Ricker
##' srProperties(2L)   # Beverton-Holt
##' @export
srProperties <- function(srCode){
  if(srCode == 0L)
    return(list(name = "random walk",
                hasEquilibrium = FALSE,
                isCompensatory = FALSE,
                hasFiniteMaxGradient = FALSE))
  if(srCode == 1L)
    return(list(name = "Ricker",
                hasEquilibrium = TRUE,
                isCompensatory = FALSE,
                hasFiniteMaxGradient = TRUE))
  if(srCode == 2L)
    return(list(name = "Beverton-Holt",
                hasEquilibrium = TRUE,
                isCompensatory = TRUE,
                hasFiniteMaxGradient = TRUE))
  stop("srCode ", srCode, " not supported in samref (only 0, 1, 2)")
}

##' Equilibrium log-SSB from a per-recruit SPR
##'
##' Closed form. AD-safe.
##'
##' Ricker:        Se = (a + log(SPR)) * exp(-b)
##' Beverton-Holt: Se = (exp(a) * SPR - 1) * exp(-b)
##' The degenerate case (Se <= 0) returns \code{-Inf}.
##'
##' @param srCode 1 (Ricker) or 2 (Beverton-Holt).
##' @param recPars length-2 vector \code{c(a, b)} on samjr's
##'   parameterisation (b on log scale).
##' @param logSPR scalar log spawners-per-recruit at the current F.
##' @return scalar \code{logSe}.
##' @examples
##' ## Beverton-Holt log-equilibrium SSB with a=1.5, b=-3 at log(SPR)=0
##' srEquilibriumR(2L, c(1.5, -3), logSPR = 0)
##' ## Ricker
##' srEquilibriumR(1L, c(1.5, -3), logSPR = 0)
##' @export
srEquilibriumR <- function(srCode, recPars, logSPR){
  a <- recPars[1]
  b <- recPars[2]
  if(srCode == 1L){
    z <- a + logSPR
    if(z <= 0) return(-Inf)
    return(log(z) - b)
  }
  if(srCode == 2L){
    z <- exp(a + logSPR) - 1
    if(z <= 0) return(-Inf)
    return(log(z) - b)
  }
  if(srCode == 0L) stop("random walk SR has no equilibrium")
  stop("srCode ", srCode, " not supported in samref")
}

##' Slope of the SR curve at S = 0
##'
##' For both Ricker and Beverton-Holt in samjr's parameterisation, the
##' slope at the origin is \code{exp(a)}.
##'
##' @param srCode 1 or 2.
##' @param recPars length-2 vector.
##' @return scalar \code{dR/dS} at \code{S = 0}.
##' @examples
##' srGradAt0(1L, c(1.5, -3))   # exp(1.5) for Ricker
##' srGradAt0(2L, c(1.5, -3))   # exp(1.5) for Beverton-Holt
##' @export
srGradAt0 <- function(srCode, recPars){
  if(srCode %in% c(1L, 2L)) return(exp(recPars[1]))
  if(srCode == 0L) stop("random walk SR has no gradient at zero")
  stop("srCode ", srCode, " not supported in samref")
}

##' Predicted log-recruitment for an SR model
##'
##' Closed-form predictor matching the closures in
##' \code{samjr/R/likelihood.R}. Useful for SR-curve overlays.
##'
##' @param srCode 1 (Ricker) or 2 (Beverton-Holt).
##' @param recPars length-2 vector.
##' @param logSsb scalar or vector log-SSB.
##' @return scalar or vector log-recruitment.
##' @examples
##' ## Ricker R(S) at SSB = exp(10..12)
##' predictLogR(1L, c(1.5, -3), logSsb = 10:12)
##' @export
predictLogR <- function(srCode, recPars, logSsb){
  a <- recPars[1]
  b <- recPars[2]
  if(srCode == 1L) return(a + logSsb - exp(b + logSsb))
  if(srCode == 2L) return(a + logSsb - log1p(exp(b + logSsb)))
  stop("predictLogR: srCode ", srCode, " not supported")
}

##' Extract SR parameters from a samjr fit
##'
##' @param fit a samjr \code{sam} object.
##' @return length-2 numeric vector of SR parameters or
##'   \code{numeric(0)} for SR code 0.
##' @examples
##' \donttest{
##' data(nscodData); data(nscodConf)
##' fit <- samjr::sam.fit(nscodData, nscodConf,
##'                       samjr::defpar(nscodData, nscodConf), silent = TRUE)
##' recPars(fit)   # numeric(0) for random-walk recruitment
##' }
##' @export
recPars <- function(fit){
  cc <- fit$conf$stockRecruitmentModelCode
  if(cc == 1L) return(fit$pl$rickerpar)
  if(cc == 2L) return(fit$pl$bhpar)
  numeric(0)
}

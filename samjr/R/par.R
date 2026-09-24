##' Setup initial values for all samjr model parameters and random effects
##'
##' Builds a parameter list compatible with the samjr likelihood from a
##' data object and a configuration. Most parameters are initialised to zero
##' (in log-space); the catchabilities \code{logQ} are initialised to
##' \eqn{-5} so that they act roughly multiplicatively on small abundance
##' indices. Recruitment, F-process, and observation random effects are also
##' included.
##'
##' If \code{conf$noScaledYears > 0}, a \code{logScale} parameter vector of
##' length \code{max(keyParScaledYA) + 1} is allocated and used by
##' \code{\link{sam.fit}} to multiplicatively rescale catch predictions in
##' the configured years.
##'
##' @param data data object as returned by \code{\link{setup.sam.data}}.
##' @param conf configuration list as returned by \code{\link{defcon}}
##' (and possibly modified).
##' @param spinoutyear number of extra years to append to the GMRF biology
##' process matrices (\code{logSW}, \code{logCW}, \code{logitMO},
##' \code{logNM}) beyond the data span, to allow for projection / spin-out.
##' Default 10.
##' @return a named list of initial values for all model parameters and
##' random effects.
##' @export
defpar <- function(data, conf, spinoutyear = 10){
  dat <- toBabyDat(data, conf, spinoutyear = spinoutyear)
  par <- list()
  par$logSdLogN <- c(0, 0)
  par$logsdF <- numeric(max(dat$keyVarFperState))
  par$rickerpar <- if(dat$srmode == 1) c(1, 1) else numeric(0)
  par$itrans_rho <- if(dat$fcormode == 0) numeric(0) else 0.1
  par$bhpar <- if(dat$srmode == 2) c(1, 1) else numeric(0)
  par$logQ <- numeric(max(dat$keyQ, na.rm = TRUE)) - 5
  qpowMax <- if(is.null(dat$keyQpow)) -1 else suppressWarnings(max(dat$keyQpow, na.rm = TRUE))
  par$logQpow <- if(!is.finite(qpowMax) || qpowMax <= 0) numeric(0) else numeric(qpowMax)
  pvolMax <- if(is.null(dat$predVarObsLink)) -1 else suppressWarnings(max(dat$predVarObsLink, na.rm = TRUE))
  par$predVarObs <- if(!is.finite(pvolMax) || pvolMax <= 0) numeric(0) else numeric(pvolMax)
  xsdMax <- suppressWarnings(max(dat$xtraSdIdxByObs, na.rm = TRUE))
  par$logXtraSd <- if(!is.finite(xsdMax) || xsdMax <= 0) numeric(0) else numeric(xsdMax)
  par$logSdLogObs <- numeric(max(dat$keySd, na.rm = TRUE))
  if(!is.null(dat$tagTypeIdx) && any(!is.na(dat$tagTypeIdx))){
    nT <- max(dat$tagTypeIdx, na.rm = TRUE)
    par$logitReleaseSurvival <- numeric(nT)
    par$logitRecapturePhi    <- numeric(nT)
  }else{
    par$logitReleaseSurvival <- numeric(0)
    par$logitRecapturePhi    <- numeric(0)
  }
  if(dat$useCKMR == 1L && dat$ckmrEstPsi == 1L){
    if(dat$ckmrPsi <= 1)
      stop("conf$ckmrPsi must exceed 1 when conf$ckmrEstimatePsi is 1")
    par$logPsim1 <- log(dat$ckmrPsi - 1)
  }else{
    par$logPsim1 <- numeric(0)
  }
  if(dat$useCKMRL == 1L && dat$ckmrlEstPsi == 1L){
    if(dat$ckmrlPsi <= 0)
      stop("conf$ckmrlPsi must be positive when conf$ckmrlEstimatePsi is 1")
    par$logPsiL <- log(dat$ckmrlPsi)
  }else{
    par$logPsiL <- numeric(0)
  }
  if(dat$useCKMRL == 1L && dat$ckmrlEstOmega == 1L){
    if(dat$ckmrlOmega <= 0)
      stop("conf$ckmrlOmega must be positive when conf$ckmrlEstimateOmega is 1")
    par$logOmegaL <- log(dat$ckmrlOmega)
  }else{
    par$logOmegaL <- numeric(0)
  }
  par$logIGARdist <- if(sum(dat$covType == 1) == 0) numeric(0)
                     else numeric(max(dat$keyIGAR, na.rm = TRUE))
  par$parUS <- unlist(lapply(seq_along(dat$covType), function(f)
    if(dat$covType[f] == 2) unstructured(dat$fleetDim[f])$parms()))
  if(is.null(par$parUS)) par$parUS <- numeric(0)
  par$missing <- numeric(sum(is.na(dat$logobs)))
  par$logN <- matrix(0, nrow = length(dat$year), ncol = length(dat$age))
  par$logF <- matrix(0, nrow = length(dat$year), ncol = max(dat$keyF))
  if(dat$noScaledYears == 0){
    par$logScale <- numeric(0)
  }else{
    par$logScale <- numeric(max(dat$keyParScaledYA, na.rm = TRUE) + 1L)
  }
  if(dat$stockWeightModel == 0){
    par$meanLogSW      <- numeric(0)
    par$logPhiSW       <- numeric(0)
    par$logSdProcLogSW <- numeric(0)
    par$logSdLogSW     <- numeric(0)
    par$logSW          <- matrix(0, 0, 0)
  }else{
    par$meanLogSW      <- numeric(max(dat$keyStockWeightMean, na.rm = TRUE))
    par$logPhiSW       <- numeric(dat$stockWeightModel + 1L)
    par$logSdProcLogSW <- 0
    par$logSdLogSW     <- numeric(max(dat$keyStockWeightObsVar, na.rm = TRUE))
    par$logSW          <- matrix(0, nrow = length(dat$year) + spinoutyear,
                                  ncol = length(dat$age))
  }
  if(dat$catchWeightModel == 0){
    par$meanLogCW      <- numeric(0)
    par$logPhiCW       <- numeric(0)
    par$logSdProcLogCW <- numeric(0)
    par$logSdLogCW     <- numeric(0)
    par$logCW          <- matrix(0, 0, 0)
  }else{
    par$meanLogCW      <- numeric(max(dat$keyCatchWeightMean, na.rm = TRUE))
    par$logPhiCW       <- numeric(dat$catchWeightModel + 1L)
    par$logSdProcLogCW <- 0
    par$logSdLogCW     <- numeric(max(dat$keyCatchWeightObsVar, na.rm = TRUE))
    par$logCW          <- matrix(0, nrow = length(dat$year) + spinoutyear,
                                  ncol = length(dat$age))
  }
  if(dat$matureModel == 0){
    par$meanLogitMO        <- numeric(0)
    par$logPhiMO           <- numeric(0)
    par$logSdProcLogitMO   <- numeric(0)
    par$logSdMO            <- numeric(0)
    par$logitMO            <- matrix(0, 0, 0)
  }else{
    par$meanLogitMO        <- numeric(max(dat$keyMatureMean, na.rm = TRUE))
    par$logPhiMO           <- numeric(dat$matureModel + 1L)
    par$logSdProcLogitMO   <- 0
    par$logSdMO            <- 0
    par$logitMO            <- matrix(0, nrow = length(dat$year) + spinoutyear,
                                      ncol = length(dat$age))
  }
  if(dat$mortalityModel == 0){
    par$meanLogNM      <- numeric(0)
    par$logPhiNM       <- numeric(0)
    par$logSdProcLogNM <- numeric(0)
    par$logSdLogNM     <- numeric(0)
    par$logNM          <- matrix(0, 0, 0)
  }else{
    par$meanLogNM      <- numeric(max(dat$keyMortalityMean, na.rm = TRUE))
    par$logPhiNM       <- numeric(dat$mortalityModel + 1L)
    par$logSdProcLogNM <- 0
    par$logSdLogNM     <- numeric(max(dat$keyMortalityObsVar, na.rm = TRUE))
    par$logNM          <- matrix(0, nrow = length(dat$year) + spinoutyear,
                                  ncol = length(dat$age))
  }
  par
}

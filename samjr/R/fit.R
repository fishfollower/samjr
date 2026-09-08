##' Internal: build a tracer that records, optionally prints, and
##' optionally live-plots negative log-likelihood and max absolute
##' gradient at every gradient evaluation. Returns a list with
##' \code{wrap(obj)} (yielding side-effecting fn/gr wrappers for
##' \code{\link[stats]{nlminb}}) and \code{record(nll, g)} (used after
##' Newton steps which operate on the unwrapped fn/gr).
##' @keywords internal
##' @noRd
makeConvTracer <- function(silent = FALSE, convplot = FALSE){
  if(silent && !convplot){
    return(list(
      record = function(nll, g) invisible(),
      wrap   = function(obj) list(fn = obj$fn, gr = obj$gr)
    ))
  }
  state <- new.env(parent = emptyenv())
  state$nllHist   <- numeric(0)
  state$maxgrHist <- numeric(0)
  state$lastNll   <- NA_real_
  record <- function(nll, g){
    mg <- max(abs(as.numeric(g)))
    state$nllHist   <- c(state$nllHist,   nll)
    state$maxgrHist <- c(state$maxgrHist, mg)
    if(!silent) cat(sprintf("nll=%-12.6g  |maxgr|=%.6g\n", nll, mg))
    if(convplot) drawConvPlot(state$nllHist, state$maxgrHist)
  }
  wrap <- function(obj){
    fn <- function(par){
      v <- obj$fn(par)
      state$lastNll <- v
      v
    }
    gr <- function(par){
      g <- obj$gr(par)
      record(state$lastNll, g)
      g
    }
    list(fn = fn, gr = gr)
  }
  list(record = record, wrap = wrap)
}

##' Internal: draw a two-axis convergence plot. Left y-axis: nll
##' (linear). Right y-axis: \code{|maxgr|} on log10 scale. Called by
##' \code{makeConvTracer} when \code{convplot = TRUE} in
##' \code{\link{sam.fit}}; re-drawn from scratch on every update.
##' @keywords internal
##' @noRd
drawConvPlot <- function(nllVec, maxgrVec){
  n <- length(nllVec)
  if(n < 1) return(invisible())
  it <- seq_len(n)
  opar <- par(mar = c(4, 4, 1, 4))
  on.exit(par(opar))
  plot(it, nllVec, type = "o", pch = 16, lwd = 2,
       xlab = "Iteration", ylab = "nll", col = "black", las = 1)
  par(new = TRUE)
  plot(it, maxgrVec, type = "o", pch = 16, lwd = 2, col = "red",
       log = "y", axes = FALSE, xlab = "", ylab = "")
  axis(side = 4, col = "red", col.axis = "red", las = 1)
  mtext("|maxgr|", side = 4, line = 2.5, col = "red")
}

##' Internal: translate a (data, conf) pair into the flat \code{dat} list
##' consumed by the samjr RTMB likelihood.
##'
##' Converts the SAM-style configuration keys (\code{keyLogFsta},
##' \code{keyLogFpar}, \code{keyVarObs}, \code{keyVarF}, \code{obsCorStruct},
##' \code{keyCorObs}) to the integer indices used inside the likelihood,
##' precomputes the per-F-state variance index (\code{keyVarFperState}),
##' the per-observation logScale index (\code{scaleIdxByObs}), and the
##' \eqn{\bar F} age slice (\code{fbarIdx}).
##'
##' @param data data object from \code{\link{setup.sam.data}}.
##' @param conf configuration list from \code{\link{defcon}}.
##' @param spinoutyear extra years appended to the GMRF biology matrices.
##' @param warn if \code{TRUE}, warn about close-kin pairs that carry no
##' information and are excluded from the likelihood. Only \code{\link{sam.fit}}
##' sets this, so the message appears once per fit.
##' @return a list of arrays and integer keys.
##' @keywords internal
##' @noRd
toBabyDat <- function(data, conf, spinoutyear = 10, warn = FALSE){
  dat <- list()
  dat$logobs      <- data$logobs
  dat$aux         <- data$aux
  dat$minYear     <- min(data$years)
  dat$minAge      <- conf$minAge
  dat$fleetTypes  <- data$fleetTypes
  dat$sampleTimes <- data$sampleTimes
  dat$year        <- data$years
  dat$age         <- conf$minAge:conf$maxAge
  dat$M  <- data$natMor
  dat$SW <- data$stockMeanWeight
  dat$MO <- data$propMat
  dat$PF <- data$propF
  dat$PM <- data$propM
  dat$CW <- data$catchMeanWeight
  ages <- dat$age
  dat$fbarIdx <- which(ages >= conf$fbarRange[1] & ages <= conf$fbarRange[2])
  dat$srmode   <- conf$stockRecruitmentModelCode
  dat$fcormode <- conf$corFlag[1]
  dat$keyF  <- conf$keyLogFsta[1, ] + 1
  dat$keyQ  <- conf$keyLogFpar + 1
  dat$keyQpow <- if(is.null(conf$keyQpow)) NULL else conf$keyQpow + 1L
  dat$keySd <- conf$keyVarObs + 1
  dat$keySd[dat$keySd <= 0] <- NA
  dat$predVarObsLink <- if(is.null(conf$predVarObsLink)) NULL else conf$predVarObsLink + 1L
  nObsTot <- length(data$logobs)
  xtraSdIdxByObs <- rep(NA_integer_, nObsTot)
  if(!is.null(conf$keyXtraSd) && nrow(conf$keyXtraSd) > 0){
    kx <- conf$keyXtraSd
    for(i in seq_len(nObsTot)){
      m <- which(kx[, 1] == data$aux[i, 2] &
                 kx[, 2] == data$aux[i, 1] &
                 kx[, 3] == data$aux[i, 3])
      if(length(m) > 0) xtraSdIdxByObs[i] <- kx[m[1], 4] + 1L
    }
  }
  dat$xtraSdIdxByObs <- xtraSdIdxByObs

  hasTag <- any(data$fleetTypes == 5) && ncol(data$aux) >= 8
  if(hasTag){
    isTag <- data$fleetTypes[data$aux[, 2]] == 5
    typeRaw <- ifelse(isTag, data$aux[, 8], NA_integer_)
    typeUniq <- sort(unique(stats::na.omit(typeRaw)))
    dat$tagNscan   <- ifelse(isTag, data$aux[, 6], NA_real_)
    dat$tagR       <- ifelse(isTag, data$aux[, 7], NA_real_)
    dat$tagTypeIdx <- match(typeRaw, typeUniq)
    dat$isTag      <- isTag
  }else{
    dat$tagNscan   <- rep(NA_real_, length(data$logobs))
    dat$tagR       <- rep(NA_real_, length(data$logobs))
    dat$tagTypeIdx <- rep(NA_integer_, length(data$logobs))
    dat$isTag      <- rep(FALSE, length(data$logobs))
  }

  dat$weight  <- if(is.null(data$weight)) rep(NA_real_, length(data$logobs)) else data$weight
  dat$corList <- if(is.null(data$corList)) list() else data$corList
  dat$idxCor  <- if(is.null(data$idxCor))
                   matrix(NA_integer_, nrow = length(data$fleetTypes), ncol = length(data$years))
                 else data$idxCor
  dat$fixVarToWeight <- if(is.null(conf$fixVarToWeight)) rep(0L, length(data$fleetTypes))
                        else as.integer(rep(conf$fixVarToWeight, length.out = length(data$fleetTypes)))
  nFstate <- max(dat$keyF)
  dat$keyVarFperState <- sapply(seq_len(nFstate),
                                function(s) conf$keyVarF[1, which(dat$keyF == s)[1]] + 1L)
  dat$fleetDim <- apply(dat$keySd, 1, function(x) sum(!is.na(x)))
  dat$covType  <- as.integer(conf$obsCorStruct) - 1
  dat$keyIGAR  <- conf$keyCorObs + 1
  dat$keyIGAR[dat$keyIGAR == 0] <- NA
  dat$keyIGAR[is.na(conf$keyCorObs)] <- -1

  dat$noScaledYears   <- as.integer(conf$noScaledYears)
  dat$keyScaledYears  <- as.integer(conf$keyScaledYears)
  dat$keyParScaledYA  <- conf$keyParScaledYA
  nobs <- length(data$logobs)
  scaleIdxByObs <- rep(NA_integer_, nobs)
  if(dat$noScaledYears > 0){
    for(i in seq_len(nobs)){
      ff <- data$aux[i, 2]
      if(data$fleetTypes[ff] == 0){
        j <- match(data$aux[i, 1], dat$keyScaledYears)
        if(!is.na(j)){
          ai <- data$aux[i, 3] - dat$minAge + 1L
          sidx <- conf$keyParScaledYA[j, ai]
          if(!is.na(sidx) && sidx >= 0) scaleIdxByObs[i] <- sidx + 1L
        }
      }
    }
  }
  dat$scaleIdxByObs <- scaleIdxByObs

  swm <- if(is.null(conf$stockWeightModel)) 0L else as.integer(conf$stockWeightModel)
  cwm <- if(is.null(conf$catchWeightModel)) 0L else as.integer(conf$catchWeightModel)
  mom <- if(is.null(conf$matureModel)) 0L else as.integer(conf$matureModel)
  nmm <- if(is.null(conf$mortalityModel)) 0L else as.integer(conf$mortalityModel)
  dat$stockWeightModel      <- swm
  dat$keyStockWeightMean    <- if(swm == 0) integer(0) else (conf$keyStockWeightMean + 1L)
  dat$keyStockWeightObsVar  <- if(swm == 0) integer(0) else (conf$keyStockWeightObsVar + 1L)
  dat$catchWeightModel      <- cwm
  dat$keyCatchWeightMean    <- if(cwm == 0) integer(0) else (conf$keyCatchWeightMean[1, ] + 1L)
  dat$keyCatchWeightObsVar  <- if(cwm == 0) integer(0) else (conf$keyCatchWeightObsVar[1, ] + 1L)
  dat$matureModel       <- mom
  dat$keyMatureMean     <- if(mom == 0) integer(0) else (conf$keyMatureMean + 1L)
  dat$mortalityModel    <- nmm
  dat$keyMortalityMean    <- if(nmm == 0) integer(0) else (conf$keyMortalityMean + 1L)
  dat$keyMortalityObsVar  <- if(nmm == 0) integer(0) else (conf$keyMortalityObsVar + 1L)
  if(swm >= 1) dat$swNobs <- length(dat$year)
  if(cwm >= 1) dat$cwNobs <- length(dat$year)
  if(mom >= 1) dat$moNobs <- length(dat$year)
  if(nmm >= 1) dat$nmNobs <- length(dat$year)
  if(swm >= 1 || cwm >= 1 || mom >= 1 || nmm >= 1){
    nyr <- length(dat$year) + spinoutyear; nag <- length(dat$age); n <- nyr * nag
    Pshape <- matrix(0, nrow = nyr, ncol = nag)
    R <- as.vector(row(Pshape)); C <- as.vector(col(Pshape))
    Wc <- Matrix::spMatrix(n, n)
    Wc[(abs(outer(R, R, "-")) == 1) & (outer(C, C, "-") == 0)] <- 1
    diag(Wc) <- -Matrix::rowSums(Wc)
    Wd <- Matrix::spMatrix(n, n)
    Wd[((outer(R, R, "-") ==  1) & (outer(C, C, "-") ==  1)) |
       ((outer(R, R, "-") == -1) & (outer(C, C, "-") == -1))] <- 1
    diag(Wd) <- -Matrix::rowSums(Wd)
    Wp <- Matrix::spMatrix(n, n)
    Wp[(abs(outer(R, R, "-")) == 1) & (C == nag) & (outer(C, C, "-") == 0)] <- 1
    diag(Wp) <- -Matrix::rowSums(Wp)
    dat$Wc <- Wc; dat$Wd <- Wd; dat$Wp <- Wp
  }

  ## ---- close-kin mark-recapture ----
  ## Every index the close-kin probabilities need is a function of the pair
  ## table and the model dimensions only, so it is built once here and cached;
  ## the likelihood then does nothing but whole-vector gathers.
  set <- ckmrSettings(conf)
  dat$usePOP    <- set$usePOP
  dat$useHSP    <- set$useHSP
  dat$ckmrPsi   <- set$psi
  dat$ckmrScale <- set$scale
  dat$ckmrEstPsi <- set$estPsi
  hasCK <- !is.null(data$ckmr) && nrow(data$ckmr) > 0 &&
           (set$usePOP == 1L || set$useHSP == 1L)
  if(hasCK){
    ck   <- ckmrCheck(data$ckmr)
    prep <- ckmrPrep(ck, dat$year, dat$age,
                     propMat = if(mom == 0) data$propMat else NULL)
    nDropP <- nrow(ck) - length(prep$popRow)
    nDropH <- nrow(ck) - length(prep$hspRow)
    if(set$usePOP == 1L && length(prep$popRow) == 0L)
      stop("conf$usePOP is on but no pair contributes a parent-offspring probability")
    if(set$useHSP == 1L && length(prep$hspRow) == 0L)
      stop("conf$useHSP is on but no pair contributes a half-sibling probability")
    if(warn && set$usePOP == 1L && nDropP > 0L)
      warning(nDropP, " of ", nrow(ck), " CKMR pairs excluded from the POP ",
              "likelihood (birth year before the model starts, parent sampled ",
              "before the offspring was born, or no mature candidate parent age)")
    if(warn && set$useHSP == 1L && nDropH > 0L)
      warning(nDropH, " of ", nrow(ck), " CKMR pairs excluded from the HSP ",
              "likelihood (same cohort, or birth year before the model starts)")
    swPos <- data$stockMeanWeight[is.finite(data$stockMeanWeight) &
                                  data$stockMeanWeight > 0]
    dat$ckmrSWref  <- if(length(swPos) > 0) exp(mean(log(swPos))) else 1
    dat$useCKMR    <- 1L
    dat$ckmrPrep   <- prep
    dat$ckmrPOPobs <- ck$nPOP[prep$popRow]
    dat$ckmrHSPobs <- ck$nHSP[prep$hspRow]
    dat$ckmrNpop   <- ck$nComp[prep$popRow]
    dat$ckmrNhsp   <- ck$nComp[prep$hspRow]
  }else{
    dat$ckmrEstPsi <- 0L
    dat$ckmrSWref  <- 1
    dat$useCKMR    <- 0L
    dat$ckmrPrep   <- NULL
    dat$ckmrPOPobs <- numeric(0)
    dat$ckmrHSPobs <- numeric(0)
    dat$ckmrNpop   <- numeric(0)
    dat$ckmrNhsp   <- numeric(0)
  }
  dat
}

##' Fit a samjr state-space stock assessment model
##'
##' Fits the samjr state-space age-structured assessment model to data
##' prepared with \code{\link{setup.sam.data}} using the configuration from
##' \code{\link{defcon}} and the parameter list from \code{\link{defpar}}.
##' Internally the likelihood is built with RTMB; the optimiser is
##' \code{\link[stats]{nlminb}}.
##'
##' Supported features in v1: stock-recruitment modes 0 (random walk),
##' 1 (Ricker) and 2 (Beverton-Holt), or a user-supplied recruitment
##' function via \code{logSRfun} / \code{logSRinit}; F-process correlation modes 0
##' (independent), 1 (compound symmetry) and 2 (AR(1)-like decay);
##' observation correlation structures \code{ID}, \code{IGAR} (single AR-style
##' decay) and \code{US} (unstructured, Cholesky-parameterised); fleet types
##' 0 (catch) and 2 (survey); plus-group accumulation in the last age;
##' missing-observation imputation; and the \code{logScale}
##' catch-scaling facility (\code{noScaledYears}, \code{keyScaledYears},
##' \code{keyParScaledYA}).
##'
##' Out of scope (v1): plotting beyond \code{\link{ssbplot}} /
##' \code{\link{fbarplot}} / \code{\link{recplot}} / \code{\link{catchplot}},
##' retrospective analysis, forecast, reference points,
##' multi-residual-fleet support, and one-step-ahead residuals.
##'
##' @param data data object from \code{\link{setup.sam.data}}.
##' @param conf configuration list from \code{\link{defcon}} (and possibly
##' modified).
##' @param parameters parameter list from \code{\link{defpar}}.
##' @param map optional named list passed to
##' \code{\link[RTMB]{MakeADFun}} to fix or share parameters.
##' @param rm.unidentified if \code{TRUE}, parameters flagged as
##' unidentified by the configuration are mapped to \code{NA} and dropped
##' from estimation.
##' @param lower,upper optional named numeric vectors of lower / upper
##' bounds on fixed-effect parameters, passed through to
##' \code{\link[stats]{nlminb}}.
##' @param newtonsteps number of post-optimisation Newton steps applied to
##' polish the fit. Default 3.
##' @param run if \code{FALSE}, return the prepared RTMB object without
##' running the optimiser (useful for debugging).
##' @param silent if \code{FALSE} (default), print the negative log
##' likelihood and the maximum absolute gradient component on every
##' iteration; set to \code{TRUE} to suppress.
##' @param convplot if \code{TRUE}, live-update a two-axis convergence
##' plot (left: nll, right: \code{|maxgr|} on log10 scale) on every
##' iteration. \code{FALSE} by default.
##' @param logSRfun optional custom stock-recruitment function of two
##' arguments, \code{S} (the relevant lagged spawning stock biomass) and
##' \code{pv} (the recruitment parameter vector), returning predicted
##' log-recruitment. When non-\code{NULL} it overrides the model selected by
##' \code{conf$stockRecruitmentModelCode}, and any conf-based \code{rickerpar}
##' / \code{bhpar} are fixed (mapped out). Use only AD-compatible operations
##' (\code{log}, \code{exp}, arithmetic, indexing) so RTMB can tape it.
##' Defaults to \code{NULL}.
##' @param logSRinit initial values for the parameter vector \code{pv} passed
##' to \code{logSRfun}; its length sets the number of estimated recruitment
##' parameters. Required when \code{logSRfun} is set. Defaults to \code{NULL}.
##' @param ... currently ignored.
##' @return an object of class \code{sam} with elements \code{data},
##' \code{conf}, \code{parameters}, the flattened \code{dat}, the RTMB
##' objective \code{obj}, the optimiser result \code{opt}, and the
##' \code{sdrep} returned by \code{\link[RTMB]{sdreport}}.
##' @seealso \code{\link{ssbplot}}, \code{\link{fbarplot}},
##' \code{\link{recplot}}, \code{\link{catchplot}}.
##' @examples
##' \donttest{
##' data(nscodData)
##' data(nscodConf)
##' par <- defpar(nscodData, nscodConf)
##' fit <- sam.fit(nscodData, nscodConf, par)
##' fit$opt$objective
##'
##' opar <- par(mfrow = c(2, 2))
##' ssbplot(fit)
##' fbarplot(fit)
##' recplot(fit)
##' catchplot(fit)
##' par(opar)
##'
##' ## a custom recruitment function reproduces the built-in Ricker
##' ## (stockRecruitmentModelCode = 1): same negative log likelihood.
##' confR <- nscodConf; confR$stockRecruitmentModelCode <- 1
##' fitR <- sam.fit(nscodData, confR, defpar(nscodData, confR))
##' fitC <- sam.fit(nscodData, confR, defpar(nscodData, confR),
##'                 logSRfun  = function(S, pv) pv[1] + log(S) - exp(pv[2]) * S,
##'                 logSRinit = c(1, 1))
##' c(ricker = fitR$opt$objective, custom = fitC$opt$objective)
##' }
##' @export
sam.fit <- function(data, conf, parameters, map = list(),
                    rm.unidentified = FALSE, lower = NULL, upper = NULL,
                    newtonsteps = 3, run = TRUE, silent = FALSE,
                    convplot = FALSE, logSRfun = NULL, logSRinit = NULL, ...){
  dat <- toBabyDat(data, conf, warn = TRUE)
  if(!is.null(logSRfun)){
    if(!is.function(logSRfun)) stop("logSRfun must be a function of (S, pv)")
    if(is.null(logSRinit)) stop("logSRinit must be supplied when logSRfun is set")
    parameters$logSRpar <- as.numeric(logSRinit)
    for(nm in c("rickerpar", "bhpar")){
      if(length(parameters[[nm]]) > 0 && is.null(map[[nm]]))
        map[[nm]] <- factor(rep(NA, length(parameters[[nm]])))
    }
  }
  f <- makeBabyLikelihood(dat, logSRfun)
  nmissing <- sum(is.na(dat$logobs))
  parameters$missing <- numeric(nmissing)
  randomVars <- c("logN", "logF", "missing")
  if(dat$stockWeightModel >= 1) randomVars <- c(randomVars, "logSW")
  if(dat$catchWeightModel >= 1) randomVars <- c(randomVars, "logCW")
  if(dat$matureModel       >= 1) randomVars <- c(randomVars, "logitMO")
  if(dat$mortalityModel    >= 1) randomVars <- c(randomVars, "logNM")
  obj <- MakeADFun(f, parameters,
                   random = randomVars,
                   map = map,
                   silent = TRUE)
  if(!run){
    fit <- list(data = data, conf = conf, parameters = parameters,
                dat = dat, obj = obj, logSRfun = logSRfun, logSRinit = logSRinit)
    class(fit) <- "sam"
    return(fit)
  }
  expandBounds <- function(L, default){
    out <- rep(default, length(obj$par))
    if(is.null(L) || length(L) == 0) return(out)
    nms <- names(obj$par)
    for(nm in names(L)){
      idx <- which(nms == nm)
      if(length(idx) > 0)
        out[idx] <- rep(L[[nm]], length.out = length(idx))
    }
    out
  }
  lowVec <- expandBounds(lower, -Inf)
  hiVec  <- expandBounds(upper,  Inf)
  tracer  <- makeConvTracer(silent = silent, convplot = convplot)
  wrapped <- tracer$wrap(obj)
  opt <- nlminb(obj$par, wrapped$fn, wrapped$gr,
                lower = lowVec, upper = hiVec,
                control = list(eval.max = 2000, iter.max = 2000,
                               rel.tol = 1e-10))
  for(i in seq_len(newtonsteps)){
    step <- try({
      g <- as.numeric(obj$gr(opt$par))
      h <- stats::optimHess(opt$par, obj$fn, obj$gr)
      solve(h, g)
    }, silent = TRUE)
    if(inherits(step, "try-error")) break
    opt$par <- opt$par - step
    tracer$record(obj$fn(opt$par), obj$gr(opt$par))
  }
  opt$objective <- obj$fn(opt$par)
  sdr <- try(sdreport(obj, getJointPrecision = TRUE), silent = TRUE)
  rep <- try(obj$report(obj$env$last.par.best), silent = TRUE)
  if(!inherits(rep, "try-error") && !is.null(rep$Slist))
    rep$obsCov <- rep$Slist
  fit <- list(data = data, conf = conf, parameters = parameters,
              dat = dat, obj = obj, opt = opt, sdrep = sdr, rep = rep,
              map = map, low = lower, hig = upper,
              logSRfun = logSRfun, logSRinit = logSRinit)
  class(fit) <- "sam"
  fit$pl <- fitParList(fit)
  fit
}

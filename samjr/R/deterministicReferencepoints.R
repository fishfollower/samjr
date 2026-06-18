## Deterministic reference points for a samjr fit.
##
## Strategy: locate each RP's F numerically via stats::optimize / uniroot
## on the per-recruit closure (samref/R/equilibrium.R); propagate
## uncertainty by the delta method over the SR-parameter vector. For
## per-recruit-only RPs (no SR dependency) we hold biology+sel fixed and
## report zero CIs in v1.

##' Locate deterministic reference points
##'
##' @param fit a samjr \code{sam} fit.
##' @param refpoints character vector of RP specifications, e.g.
##'   \code{c("MSY", "0.35SPR", "Fmax")}. See \code{\link{parseRefpoint}}.
##' @param catchType one of \code{"catch"}, \code{"landing"}, \code{"discard"}.
##' @param Fsequence F-bar grid used for diagnostic graphs and for
##'   bracketing the solver. Default \code{seq(0, 2, length = 50)}.
##' @param aveYears years to average biology over (default: last 10 years).
##' @param selYears years to average selectivity over (default: last year).
##' @param customSel optional override of selectivity (numeric vector of
##'   length nAge).
##' @param mdyRate discount rate used by the \code{MDY} reference point
##'   (default 0.05).
##' @param ... unused.
##' @return list of class \code{sam_referencepoints} with components
##'   \code{tables}, \code{graphs}, \code{fbarlabel}, \code{stochastic}.
##' @examples
##' \donttest{
##' data(nscodData); data(nscodConf)
##' fit <- samjr::sam.fit(nscodData, nscodConf,
##'                       samjr::defpar(nscodData, nscodConf), silent = TRUE)
##' ## per-recruit reference points - no SR model needed
##' rp <- referencepoints(fit, c("Max", "0.1dYPR", "0.35SPR"),
##'                       catchType = "landing")
##' rp$tables$F
##' plot(rp, type = "ypr")
##'
##' ## SR-equilibrium points need Ricker or Beverton-Holt
##' fit$conf$stockRecruitmentModelCode <- 2L
##' fitBH <- samjr::runwithout(fit)
##' rpSR <- referencepoints(fitBH, c("MSY", "0.2B0", "Ext"),
##'                         catchType = "landing")
##' summary(rpSR)
##' }
##' @export
referencepoints <- function(fit, refpoints, ...) UseMethod("referencepoints")

##' @rdname referencepoints
##' @method referencepoints sam
##' @export
referencepoints.sam <- function(fit, refpoints,
                                catchType = c("catch", "landing", "discard"),
                                Fsequence = seq(0, 2, length = 50),
                                aveYears  = NULL,
                                selYears  = NULL,
                                customSel = NULL,
                                mdyRate   = 0.05,
                                ...){
  catchType <- match.arg(catchType)
  yrs <- fit$data$years
  if(is.null(aveYears)) aveYears <- yrs[max(1L, length(yrs) - 9L):length(yrs)]
  if(is.null(selYears)) selYears <- yrs[length(yrs)]

  perRec <- perRecruitClosure(fit, aveYears = aveYears, selYears = selYears,
                              catchType = catchType, customSel = customSel)
  srCode <- fit$conf$stockRecruitmentModelCode
  rpHat  <- recPars(fit)
  base0  <- perRec(log(1e-12), rpHat)

  parsed <- lapply(refpoints, parseRefpoint)
  parsed <- Filter(function(r) r$rpType != -99L, parsed)
  merged <- Reduce(mergeRefpoint, parsed, list())

  ## Pre-compute joint covariance once (slow on first call) and cache on fit.
  if(inherits(fit$sdrep$jointPrecision, "Matrix") &&
     is.null(attr(fit, ".samref_jointCov"))){
    attr(fit, ".samref_jointCov") <-
      as.matrix(Matrix::solve(fit$sdrep$jointPrecision))
  }
  ctx <- list(aveYears = aveYears, selYears = selYears,
              catchType = catchType, customSel = customSel,
              Fsequence = Fsequence, mdyRate = mdyRate)

  rows <- list()
  fmsyCached <- NA_real_
  for(rp in merged){
    .validateRP(rp$rpType, srCode)
    if(rp$rpType == 2L){
      if(is.na(fmsyCached)){
        loc <- .locateRP(1L, NA_real_, perRec, rpHat, base0, Fsequence, fit,
                         srCode, mdyRate)
        fmsyCached <- loc$Fbar
        ymaxCached <- exp(loc$out$logYe)
      }
      for(xv in rp$xVal){
        locLow <- .locateMSYRange(xv, perRec, rpHat, Fsequence, fmsyCached,
                                  ymaxCached, side = "low")
        locUp  <- .locateMSYRange(xv, perRec, rpHat, Fsequence, fmsyCached,
                                  ymaxCached, side = "high")
        labs <- refpointName(rp)
        idxLow <- 2L * match(xv, rp$xVal) - 1L
        ## For MSYRange we don't re-locate the bracket under each FD step.
        rows[[labs[idxLow]]]     <- .rowFromLoc(locLow, perRec, rpHat, fit, srCode,
                                                 rpType = NA_integer_, ctx = ctx)
        rows[[labs[idxLow + 1L]]] <- .rowFromLoc(locUp, perRec, rpHat, fit, srCode,
                                                 rpType = NA_integer_, ctx = ctx)
      }
      next
    }
    if(length(rp$xVal) == 0L){
      loc <- .locateRP(rp$rpType, NA_real_, perRec, rpHat, base0,
                       Fsequence, fit, srCode, mdyRate)
      rows[[refpointName(rp)]] <- .rowFromLoc(loc, perRec, rpHat, fit, srCode,
                                              rpType = rp$rpType,
                                              xVal = NA_real_, ctx = ctx)
    } else {
      labs <- refpointName(rp)
      for(i in seq_along(rp$xVal)){
        loc <- .locateRP(rp$rpType, rp$xVal[i], perRec, rpHat, base0,
                         Fsequence, fit, srCode, mdyRate)
        rows[[labs[i]]] <- .rowFromLoc(loc, perRec, rpHat, fit, srCode,
                                       rpType = rp$rpType, xVal = rp$xVal[i],
                                       ctx = ctx)
      }
    }
  }

  tables <- .buildTables(rows)
  graphs <- .buildGraphs(perRec, rpHat, Fsequence, srCode)
  fbarRange <- fit$conf$fbarRange
  ret <- list(tables = tables, graphs = graphs,
              fbarlabel = substitute(bar(F)[X - Y],
                                     list(X = fbarRange[1], Y = fbarRange[2])),
              stochastic = FALSE,
              aveYears = aveYears, selYears = selYears,
              catchType = catchType)
  attr(ret, "fit") <- fit
  class(ret) <- "sam_referencepoints"
  ret
}

##' @rdname referencepoints
##' @export
deterministicReferencepoints <- function(fit, refpoints, ...){
  referencepoints(fit, refpoints, ...)
}

## -------- internal helpers --------

.validateRP <- function(rpType, srCode){
  prop <- srProperties(srCode)
  needsEq    <- rpType %in% c(1L, 2L, 6L, 7L, 8L, 9L, 11L)
  needsGrad  <- rpType %in% c(10L, 11L)
  needsCmp   <- rpType == 10L
  if(needsEq && !prop$hasEquilibrium)
    stop("Reference point requires SR equilibrium; not available for ",
         prop$name)
  if(needsGrad && !prop$hasFiniteMaxGradient)
    stop("Reference point requires finite max SR gradient; not available for ",
         prop$name)
  if(needsCmp && !prop$isCompensatory)
    stop("Reference point (Crash) requires compensatory SR; not available for ",
         prop$name)
  if(rpType == 12L) stop("Lim reference point not implemented")
  invisible(NULL)
}

## Numerical first derivative of YPR w.r.t. F at logF (via central FD).
.dYPRdF <- function(perRec, rpHat, logF, h = 1e-5){
  yphi <- perRec(logF + h, rpHat)$logYPR
  ylo  <- perRec(logF - h, rpHat)$logYPR
  (exp(yphi) - exp(ylo)) / (2 * h * exp(logF))
}

.locateRP <- function(rpType, xVal, perRec, rpHat, base0, Fsequence,
                      fit, srCode, mdyRate){
  ## Returns list(Fbar = locatedF, logF = locatedLogF, out = perRec(out))
  logFseq <- log(pmax(Fsequence, 1e-12))
  if(rpType == -1L){
    Fbar <- xVal
    logF <- log(max(Fbar, 1e-12))
    return(list(Fbar = Fbar, logF = logF, out = perRec(logF, rpHat)))
  }
  if(rpType == 0L){
    lag <- if(length(xVal) == 0 || is.na(xVal)) 0L else as.integer(xVal)
    nY <- length(fit$data$years)
    fbarVec <- fbartable(fit)[, "Estimate"]
    Fbar <- fbarVec[nY - lag]
    logF <- log(max(Fbar, 1e-12))
    return(list(Fbar = Fbar, logF = logF, out = perRec(logF, rpHat)))
  }
  ## Grid scan to bracket
  vals <- vapply(logFseq, function(lf){
    o <- perRec(lf, rpHat)
    switch(as.character(rpType),
           "1" = if(is.na(o$logYe)) -Inf else o$logYe,
           "3" = o$logYPR,
           "4" = .dYPRdF(perRec, rpHat, lf) /
                 .dYPRdF(perRec, rpHat, log(1e-12)) - xVal,
           "5" = o$logSPR - log(xVal) - base0$logSPR,
           "6" = exp(o$logSe) - xVal * exp(base0$logSe),
           "7" = if(is.na(o$logYe)) -Inf
                 else o$logYe - log(1 + exp(o$logYearsLost) / fit$conf$maxAge),
           "8" = if(is.na(o$logYe)) -Inf
                 else o$logYe + log(1 - exp(xVal * (o$logYearsLost - log(fit$conf$maxAge)))),
           "9" = if(is.na(o$logYe)) -Inf
                 else o$logYe - log(1 + mdyRate) * (1 + exp(o$logYearsLost)) / 2,
           "10" = srGradAt0(srCode, rpHat) * exp(o$logSPR) - 1,
           "11" = exp(o$logSe) - 1,
           NA_real_)
  }, numeric(1))
  isOptim <- rpType %in% c(1L, 3L, 7L, 8L, 9L)
  isRoot  <- rpType %in% c(4L, 5L, 6L, 10L, 11L)
  ## Locate
  if(isOptim){
    idx <- which.max(vals)
    lo <- max(1L, idx - 1L); hi <- min(length(logFseq), idx + 1L)
    obj <- function(lf){
      o <- perRec(lf, rpHat)
      v <- switch(as.character(rpType),
                  "1" = o$logYe,
                  "3" = o$logYPR,
                  "7" = o$logYe - log(1 + exp(o$logYearsLost) / fit$conf$maxAge),
                  "8" = o$logYe + log(1 - exp(xVal * (o$logYearsLost - log(fit$conf$maxAge)))),
                  "9" = o$logYe - log(1 + mdyRate) * (1 + exp(o$logYearsLost)) / 2)
      -v
    }
    sol <- stats::optimize(obj, interval = c(logFseq[lo], logFseq[hi]),
                            tol = 1e-10)
    logF <- sol$minimum
  } else if(isRoot){
    signs <- sign(vals)
    cross <- which(diff(signs) != 0)
    if(length(cross) == 0L){
      ## Fall back to brackets at endpoints
      cross <- if(vals[1] * vals[length(vals)] < 0) 1L else integer(0)
    }
    if(length(cross) == 0L)
      stop("Could not bracket root for RP type ", rpType)
    j <- if(rpType == 11L) cross[1L] else cross[length(cross)]
    objR <- function(lf){
      o <- perRec(lf, rpHat)
      switch(as.character(rpType),
             "4" = .dYPRdF(perRec, rpHat, lf) /
                   .dYPRdF(perRec, rpHat, log(1e-12)) - xVal,
             "5" = o$logSPR - log(xVal) - base0$logSPR,
             "6" = exp(o$logSe) - xVal * exp(base0$logSe),
             "10" = srGradAt0(srCode, rpHat) * exp(o$logSPR) - 1,
             "11" = exp(o$logSe) - 1)
    }
    sol <- stats::uniroot(objR, interval = c(logFseq[j], logFseq[j + 1L]),
                          tol = 1e-12)
    logF <- sol$root
  } else {
    stop("Unhandled RP type ", rpType)
  }
  list(Fbar = exp(logF), logF = logF, out = perRec(logF, rpHat))
}

.locateMSYRange <- function(xVal, perRec, rpHat, Fsequence, fmsy, msyY, side){
  target <- xVal * msyY
  obj <- function(lf) exp(perRec(lf, rpHat)$logYe) - target
  logFseq <- log(pmax(Fsequence, 1e-12))
  if(side == "low"){
    lo <- logFseq[1]; hi <- log(fmsy)
  } else {
    lo <- log(fmsy); hi <- logFseq[length(logFseq)]
  }
  ## Walk inward to find a sign change
  vlo <- obj(lo); vhi <- obj(hi)
  if(vlo * vhi > 0)
    return(list(Fbar = NA_real_, logF = NA_real_, out = perRec(log(fmsy), rpHat)))
  sol <- stats::uniroot(obj, interval = c(lo, hi), tol = 1e-12)
  list(Fbar = exp(sol$root), logF = sol$root, out = perRec(sol$root, rpHat))
}

## ---- Helpers for the full delta method over (rec_pars + selectivity logF) ----

.paramSpec <- function(fit, selYears){
  H <- fit$sdrep$jointPrecision
  empty <- list(recIdx = integer(0), selLogFidx = integer(0),
                recPars = numeric(0), selLogF = numeric(0))
  if(!inherits(H, "Matrix") || is.null(rownames(H))) return(empty)
  nm <- rownames(H)
  srCode <- fit$conf$stockRecruitmentModelCode
  recNm <- if(srCode == 1L) "rickerpar"
           else if(srCode == 2L) "bhpar" else NA_character_
  recIdx <- if(!is.na(recNm)) which(nm == recNm) else integer(0)
  logFblock <- which(nm == "logF")
  if(length(logFblock) == 0L) return(empty)
  nLogFstates <- ncol(fit$pl$logF)
  nYears <- nrow(fit$pl$logF)
  selYearIdx <- match(selYears, fit$data$years)
  posInBlock <- as.integer(outer(selYearIdx,
                                  (0:(nLogFstates - 1L)) * nYears, "+"))
  selLogFidx <- logFblock[posInBlock]
  lpb <- fit$obj$env$last.par.best
  list(recIdx = recIdx, selLogFidx = selLogFidx,
       recPars = if(length(recIdx) > 0) as.numeric(lpb[recIdx]) else numeric(0),
       selLogF = as.numeric(lpb[selLogFidx]))
}

## Marginal covariance of the chosen subset from the joint precision.
## Cached on the fit object via attribute to avoid resolving 672x672 per RP.
.paramCov <- function(fit, spec){
  cached <- attr(fit, ".samref_jointCov")
  if(is.null(cached))
    cached <- as.matrix(Matrix::solve(fit$sdrep$jointPrecision))
  idx <- c(spec$recIdx, spec$selLogFidx)
  cached[idx, idx, drop = FALSE]
}

## Normalised selectivity vector from a length-nState logF row.
.selFromLogFRow <- function(logFvec, conf){
  keyF <- conf$keyLogFsta[1, ]
  ages <- conf$minAge:conf$maxAge
  fa <- numeric(length(ages))
  for(a in seq_along(ages)){
    k <- keyF[a]
    if(is.na(k) || k < 0) next
    fa[a] <- exp(logFvec[k + 1L])
  }
  fbarAges <- which(ages %in% conf$fbarRange[1]:conf$fbarRange[2])
  scale <- mean(fa[fbarAges])
  if(scale <= 0) return(rep(0, length(fa)))
  fa / scale
}

## Output vector on the log scale. We do the FD Jacobian and the
## delta-method variance on this scale (matching SAM), then build the
## CIs as exp(logEst +/- 1.96 * SE_log). This gives the SAM-style
## asymmetric CIs for positive quantities (F, Yield, SPR, ...) and is
## also numerically better behaved.
.outVec <- function(logF, perRec, recPars){
  o <- perRec(logF, recPars)
  c(F   = logF,
    Y   = if(is.na(o$logYe))  NA_real_ else o$logYe,
    YPR = o$logYPR,
    SPR = o$logSPR,
    SSB = if(is.na(o$logSe))  NA_real_ else o$logSe,
    R   = if(is.na(o$logRe))  NA_real_ else o$logRe,
    LifeExp   = o$logLifeExpectancy,
    YearsLost = o$logYearsLost)
}

## Convert a located RP into an Estimate/Low/High triple per output quantity.
## Delta method over (rec_pars + logF at selYears) when fit$sdrep$jointPrecision
## is present; falls back to rec_pars-only when it isn't.
.rowFromLoc <- function(loc, perRec, rpHat, fit, srCode,
                        rpType = NA_integer_, xVal = numeric(0),
                        ctx = NULL){
  out <- loc$out
  base <- .outVec(loc$logF, perRec, rpHat)
  vars <- rep(0, length(base)); names(vars) <- names(base)
  spec <- .paramSpec(fit, if(is.null(ctx)) NULL else ctx$selYears)
  nP <- length(spec$recIdx) + length(spec$selLogFidx)

  if(nP > 0 && !is.null(ctx)){
    Sigma <- .paramCov(fit, spec)
    params0 <- c(spec$recPars, spec$selLogF)
    nRec <- length(spec$recIdx)
    h <- 1e-5
    decompose <- function(p){
      list(rec = if(nRec > 0) p[seq_len(nRec)] else numeric(0),
           selLogF = p[(nRec + 1L):length(p)])
    }
    evalAt <- function(rec, selLogF){
      newSel <- .selFromLogFRow(selLogF, fit$conf)
      perRec2 <- perRecruitClosure(fit, aveYears = ctx$aveYears,
                                   selYears = ctx$selYears,
                                   catchType = ctx$catchType,
                                   customSel = newSel)
      newLogF <- loc$logF
      if(!is.na(rpType) && !(rpType %in% c(-1L, 0L))){
        base0New <- perRec2(log(1e-12), rec)
        locNew <- tryCatch(.locateRP(rpType, xVal, perRec2, rec, base0New,
                                      ctx$Fsequence, fit, srCode, ctx$mdyRate),
                            error = function(e) NULL)
        if(!is.null(locNew)) newLogF <- locNew$logF
      }
      .outVec(newLogF, perRec2, rec)
    }
    J <- matrix(0, length(base), nP)
    for(k in seq_len(nP)){
      dp <- rep(0, nP); dp[k] <- h
      pl <- decompose(params0 + dp); mi <- decompose(params0 - dp)
      fhi <- evalAt(pl$rec, pl$selLogF)
      flo <- evalAt(mi$rec, mi$selLogF)
      J[, k] <- (fhi - flo) / (2 * h)
    }
    V <- J %*% Sigma %*% t(J)
    vars <- pmax(diag(V), 0); names(vars) <- names(base)
  } else if(srCode %in% c(1L, 2L) && !is.null(fit$sdrep$cov.fixed)){
    ## Fallback path: rec-pars-only delta method.
    cv <- fit$sdrep$cov.fixed
    nm <- if(srCode == 1L) "rickerpar" else "bhpar"
    keep <- names(diag(cv)) %in% nm
    if(any(keep)){
      Sigma <- cv[keep, keep, drop = FALSE]
      h <- 1e-5
      J <- matrix(0, length(base), length(rpHat))
      for(k in seq_along(rpHat)){
        dp <- rep(0, length(rpHat)); dp[k] <- h
        J[, k] <- (.outVec(loc$logF, perRec, rpHat + dp) -
                   .outVec(loc$logF, perRec, rpHat - dp)) / (2 * h)
      }
      V <- J %*% Sigma %*% t(J)
      vars <- pmax(diag(V), 0); names(vars) <- names(base)
    }
  }
  ## Convert log-scale (Estimate, Var) into linear CIs the SAM way:
  ## Estimate = exp(logEst); Low/High = exp(logEst +/- 1.96 * SE_log).
  ci <- function(logEst, vr){
    if(is.na(logEst)) return(c(Estimate = NA_real_, Low = NA_real_, High = NA_real_))
    est <- exp(logEst)
    if(is.na(vr) || vr <= 0) return(c(Estimate = est, Low = est, High = est))
    se <- sqrt(vr)
    c(Estimate = est,
      Low      = exp(logEst - 2 * se),
      High     = exp(logEst + 2 * se))
  }
  mapply(ci, base, vars, SIMPLIFY = TRUE)
}

.buildTables <- function(rows){
  if(length(rows) == 0L)
    return(list())
  qnames <- rownames(rows[[1]])
  fields <- c(F = "F", Yield = "Y", YieldPerRecruit = "YPR",
              SpawnersPerRecruit = "SPR", Biomass = "SSB", Recruitment = "R",
              LifeExpectancy = "LifeExp", LifeYearsLost = "YearsLost")
  out <- lapply(fields, function(field){
    mat <- t(vapply(rows, function(m) m[, field], numeric(3)))
    colnames(mat) <- qnames
    mat
  })
  out
}

.buildGraphs <- function(perRec, rpHat, Fsequence, srCode){
  rows <- lapply(log(pmax(Fsequence, 1e-12)), function(lf) perRec(lf, rpHat))
  ext <- function(nm){
    vapply(rows, function(r){
      v <- r[[nm]]
      if(is.null(v) || is.na(v)) NA_real_ else exp(v)
    }, numeric(1))
  }
  list(F                  = Fsequence,
       Yield              = ext("logYe"),
       YieldPerRecruit    = ext("logYPR"),
       SpawnersPerRecruit = ext("logSPR"),
       Biomass            = ext("logSe"),
       Recruitment        = ext("logRe"),
       YearsLost          = ext("logYearsLost"),
       LifeExpectancy     = ext("logLifeExpectancy"))
}

##' Internal: GMRF process negative log-density on a year-by-age latent matrix.
##' @keywords internal
##' @noRd
gmrfProcNll <- function(latent, mean, keyMean, logPhi, logSdProc,
                        Wc, Wd, Wp, mode){
  phi   <- exp(logPhi)
  nA    <- ncol(latent); nFull <- nrow(latent)
  Q <- Matrix::Diagonal(nFull * nA) - phi[1] * Wc - phi[2] * Wd
  if(mode == 2) Q <- Q - phi[3] * Wp
  mP <- outer(rep(1, nFull), mean[keyMean])
  -dgmrf(as.vector(latent - mP), mu = 0, Q = Q, log = TRUE, scale = exp(logSdProc))
}

##' Internal: per-cell log-normal observation negative log-likelihood
##' for a positive observation matrix \code{Wobs} against the first
##' \code{nObs} rows of a log-scale latent matrix \code{logW}.
##' @keywords internal
##' @noRd
logNormalObsNll <- function(logW, Wobs, nObs, logSdObs, keyObsVar){
  logWobs <- logW[1:nObs, , drop = FALSE]
  sdMat <- matrix(exp(logSdObs[keyObsVar]), nObs, ncol(Wobs), byrow = TRUE)
  ok <- !is.na(Wobs)
  -sum(dnorm(log(Wobs[ok]), logWobs[ok], sdMat[ok], log = TRUE))
}

##' Internal: per-cell beta observation negative log-likelihood for a
##' proportion observation matrix \code{Pobs} (values in [0,1]) against
##' the first \code{nObs} rows of a logit-scale latent matrix
##' \code{logitP}. Uses the standard mean / precision parameterisation
##' with \code{prec = plogis(logSdP) * 1000} and SAM's "squash" trick to
##' keep observations strictly inside (0, 1).
##' @keywords internal
##' @noRd
betaObsNll <- function(logitP, Pobs, nObs, logSdP){
  squash <- function(u) (1 - 1e-6) * (u - 0.5) + 0.5
  logitObs <- logitP[1:nObs, , drop = FALSE]
  m    <- plogis(logitObs)
  prec <- plogis(logSdP[1]) * 1000
  a    <- m * prec
  b    <- (1 - m) * prec
  P0   <- ifelse(is.na(Pobs), 0.5, Pobs)
  hi   <- P0 > 0.5
  contrib <- ifelse(hi, dbeta(squash(1 - P0), b, a, log = TRUE), dbeta(squash(P0), a, b, log = TRUE))
  maskNum <- as.numeric(!is.na(Pobs))
  -sum(contrib * maskNum)
}

##' Internal: GMRF + log-normal weight-process contribution.
##' Used for stockWeightModel, catchWeightModel, mortalityModel.
##' @keywords internal
##' @noRd
gmrfWeightContrib <- function(logW, Wobs, nObs, mean, keyMean, logPhi, logSdProc, logSdObs, keyObsVar, Wc, Wd, Wp, mode){
  nll <- gmrfProcNll(logW, mean, keyMean, logPhi, logSdProc, Wc, Wd, Wp, mode) + logNormalObsNll(logW, Wobs, nObs, logSdObs, keyObsVar)
  list(nll = nll, smoothed = exp(logW[1:nObs, , drop = FALSE]))
}

##' Internal: GMRF + beta proportion-process contribution.
##' Used for matureModel.
##' @keywords internal
##' @noRd
gmrfMatureContrib <- function(logitP, Pobs, nObs, mean, keyMean, logPhi, logSdProc, logSdP, Wc, Wd, Wp, mode){
  nll <- gmrfProcNll(logitP, mean, keyMean, logPhi, logSdProc, Wc, Wd, Wp, mode) + betaObsNll(logitP, Pobs, nObs, logSdP)
  list(nll = nll, smoothed = plogis(logitP[1:nObs, , drop = FALSE]))
}

##' Internal: spawning stock biomass per year
##'
##' Vectorised SSB computation. Given N, F, M, SW, MO, PF, PM as year-by-age
##' matrices, returns the per-year SSB on the natural scale.
##' @keywords internal
##' @noRd
ssbFun <- function(N, FF, M, SW, MO, PF, PM){
  X <- SW * MO * N * exp(-PF * FF - PM * M)
  as.vector(X %*% rep(1, ncol(X)))
}

##' Internal: build the samjr RTMB likelihood closure
##'
##' Returns a function of the parameter list that computes the negative
##' joint log likelihood (jnll), reports the per-fleet observation covariance
##' matrices and predicted observations, and ADREPORTs \code{logssb},
##' \code{logfbar}, \code{logR}, \code{logCatch}.
##' @param dat the flattened \code{dat} list from \code{toBabyDat}.
##' @keywords internal
##' @noRd
makeBabyLikelihood <- function(dat){
  function(par){
    getAll(par, dat)
    resFlag <- if(is.null(dat$resFlag)) 0L else as.integer(dat$resFlag)
    logobs <- OBS(logobs)
    osaMode <- inherits(logobs, "osa")
    if(length(missing) > 0){
      logobs <- AD(logobs)
      logobs[is.na(logobs)] <- missing
    }
    nobs <- length(logobs)
    nrow <- nrow(M)
    ncol <- ncol(M)
    sdR <- exp(logSdLogN[1])
    sdS <- exp(logSdLogN[2])
    sdF <- exp(logsdF[keyVarFperState])
    sdO <- exp(logSdLogObs)
    hasF  <- keyF > 0
    logFF <- AD(matrix(0, nrow = nrow(logF), ncol = length(keyF)))
    logFF[, hasF] <- logF[, keyF[hasF]]
    F <- exp(logFF)
    F[, !hasF] <- 0

    jnll <- 0
    if(osaMode){
      jnll <- jnll - sum(dnorm(logN[1, ], 0, 10, log = TRUE))
      jnll <- jnll - sum(dnorm(logF[1, ], 0, 10, log = TRUE))
    }

    if(stockWeightModel >= 1){
      sw <- gmrfWeightContrib(logSW, SW, swNobs, meanLogSW, keyStockWeightMean, logPhiSW, logSdProcLogSW, logSdLogSW, keyStockWeightObsVar, Wc, Wd, Wp, stockWeightModel)
      jnll <- jnll + sw$nll
      SW <- sw$smoothed
      ADREPORT(logSW)
    }
    if(catchWeightModel >= 1){
      cw <- gmrfWeightContrib(logCW, CW, cwNobs, meanLogCW, keyCatchWeightMean, logPhiCW, logSdProcLogCW, logSdLogCW, keyCatchWeightObsVar, Wc, Wd, Wp, catchWeightModel)
      jnll <- jnll + cw$nll
      CW <- cw$smoothed
      ADREPORT(logCW)
    }
    if(matureModel >= 1){
      mo <- gmrfMatureContrib(logitMO, MO, moNobs, meanLogitMO, keyMatureMean, logPhiMO, logSdProcLogitMO, logSdMO, Wc, Wd, Wp, matureModel)
      jnll <- jnll + mo$nll
      MO <- mo$smoothed
      ADREPORT(logitMO)
    }
    if(mortalityModel >= 1){
      nmRes <- gmrfWeightContrib(logNM, M, nmNobs, meanLogNM, keyMortalityMean, logPhiNM, logSdProcLogNM, logSdLogNM, keyMortalityObsVar, Wc, Wd, Wp, mortalityModel)
      jnll <- jnll + nmRes$nll
      M <- nmRes$smoothed
      ADREPORT(logNM)
    }

    ssb <- ssbFun(exp(logN), exp(logFF), M, SW, MO, PF, PM)

    yIdx <- 2:nrow
    ssbLag <- ssb[pmax(yIdx - minAge, 1)]
    if(srmode == 0){
      pred <- logN[yIdx - 1, 1]
    }
    if(srmode == 1){
      pred <- rickerpar[1] + log(ssbLag) - exp(rickerpar[2]) * ssbLag
    }
    if(srmode == 2){
      pred <- bhpar[1] + log(ssbLag) - log(1.0 + exp(bhpar[2]) * ssbLag)
    }
    if(!(srmode %in% c(0, 1, 2))){
      stop(paste("srmode", srmode, "not implemented"))
    }
    jnll <- jnll - sum(dnorm(logN[yIdx, 1], pred, sdR, log = TRUE))
    if(resFlag == 1L){
      resN <- AD(matrix(0, nrow = ncol, ncol = nrow - 1L))
      resN[1, yIdx - 1L] <- (logN[yIdx, 1] - pred) / sdR
    }

    for(y in 2:nrow){
      for(a in 2:ncol){
        pred <- logN[y - 1, a - 1] - F[y - 1, a - 1] - M[y - 1, a - 1]
        if(a == ncol){
          pred <- log(exp(pred) + exp(logN[y - 1, a] - F[y - 1, a] - M[y - 1, a]))
        }
        jnll <- jnll - dnorm(logN[y, a], pred, sdS, log = TRUE)
        if(resFlag == 1L){
          resN[a, y - 1L] <- (logN[y, a] - pred) / sdS
        }
      }
    }

    nF <- ncol(logF)
    eyeF <- diag(nF)
    if(fcormode == 0){
      corMatF <- eyeF
    }
    if(fcormode == 1){
      rhoF <- 2 * plogis(itrans_rho[1]) - 1
      onesF <- matrix(1, nF, nF)
      corMatF <- onesF * rhoF + eyeF * (1 - rhoF)
    }
    if(fcormode == 2){
      rhoF <- 2 * plogis(itrans_rho[1]) - 1
      distF <- abs(outer(seq_len(nF), seq_len(nF), "-"))
      corMatF <- rhoF^distF
    }
    SigmaF <- outer(sdF, sdF) * corMatF
    if(resFlag == 1L){
      resF <- AD(matrix(0, nrow = nF, ncol = nrow - 1L))
      Lf <- t(chol(SigmaF))
    }
    for(y in 2:nrow){
      jnll <- jnll - dmvnorm(logF[y, ], logF[y - 1, ], SigmaF, log = TRUE)
      if(resFlag == 1L){
        resF[, y - 1L] <- solve(Lf, logF[y, ] - logF[y - 1, ])
      }
    }

    logPred <- AD(numeric(nobs))
    for(i in 1:nobs){
      ff <- aux[i, 2]
      if(fleetTypes[ff] == 5) next
      y <- aux[i, 1] - minYear + 1
      if(fleetTypes[ff] == 3){
        logPred[i] <- log(ssb[y]) + logQ[keyQ[ff, 1]]
        next
      }
      a <- aux[i, 3] - minAge + 1
      Z <- F[y, a] + M[y, a]
      if(fleetTypes[ff] == 0){
        logPred[i] <- logN[y, a] - log(Z) + log(-expm1(-Z)) + logFF[y, a]
        if(!is.na(scaleIdxByObs[i])){
          logPred[i] <- logPred[i] - logScale[scaleIdxByObs[i]]
        }
      }else if(fleetTypes[ff] == 2){
        survPred <- logN[y, a] - Z * sampleTimes[ff]
        if(!is.null(keyQpow) && !is.na(keyQpow[ff, a]) && keyQpow[ff, a] >= 1){
          survPred <- survPred * exp(logQpow[keyQpow[ff, a]])
        }
        logPred[i] <- logQ[keyQ[ff, a]] + survPred
      }else if(fleetTypes[ff] == 4){
        survPred <- logN[y, a] - Z * sampleTimes[ff]
        if(!is.null(keyQpow) && !is.na(keyQpow[ff, a]) && keyQpow[ff, a] >= 1){
          survPred <- survPred * exp(logQpow[keyQpow[ff, a]])
        }
        logPred[i] <- logQ[keyQ[ff, a]] + log(MO[y, a]) + survPred
      }else{
        stop("Fleet type not implemented")
      }
    }

    Slist <- list()
    for(ff in unique(aux[, 2])){
      if(fleetTypes[ff] == 5) next
      sdv <- sdO[na.omit(keySd[ff, ])]
      n <- length(sdv)
      eye <- diag(n)
      if(covType[ff] == 0){
        corMat <- eye
      }
      if(covType[ff] == 1){
        distExt <- exp(logIGARdist[na.omit(keyIGAR[ff, ])])
        dist <- cumsum(c(AD(0), distExt))
        corMat <- AD(diag(n))
        for(i in 2:n){
          for(j in 1:(i - 1)){
            v <- 0.5^(dist[i] - dist[j])
            corMat[i, j] <- v
            corMat[j, i] <- v
          }
        }
      }
      if(covType[ff] == 2){
        idx <- which(ff == unlist(sapply(seq_along(covType), function(g) if(covType[g] == 2) unstructured(fleetDim[g])$parms() + g)))
        corMat <- unstructured(n)$corr(parUS[idx])
      }
      if(!(covType[ff] %in% c(0, 1, 2))){
        stop("Covariance type not implemented")
      }
      S <- outer(sdv, sdv) * corMat
      Slist[[length(Slist) + 1]] <- S
      hasIdxCor <- any(!is.na(idxCor[ff, ]))
      for(y in unique(aux[, 1])){
        idx <- which((aux[, 2] == ff) & (aux[, 1] == y))
        if(length(idx) == 0) next
        thisCor <- NULL
        if(hasIdxCor){
          yIdx <- which(year == y)
          if(length(yIdx) > 0 && !is.na(idxCor[ff, yIdx]))
            thisCor <- corList[[idxCor[ff, yIdx]]]
        }
        R <- if(!is.null(thisCor)) thisCor
             else if(covType[ff] == 0) diag(length(idx)) else corMat
        effSd <- AD(numeric(length(idx)))
        for(j in seq_along(idx)){
          i <- idx[j]
          a <- if(fleetTypes[ff] == 3) 1L else aux[i, 3] - minAge + 1L
          base <- logSdLogObs[keySd[ff, a]]
          baseSd <- exp(base)
          xtraIdx <- xtraSdIdxByObs[i]
          link <- if(!is.null(predVarObsLink)) predVarObsLink[ff, a] else NA_integer_
          w <- weight[i]
          if(!is.na(w) && fixVarToWeight[ff] == 1){
            effSd[j] <- sqrt(w)
          }else if(!is.na(w)){
            effSd[j] <- baseSd / sqrt(w)
          }else if(!is.na(xtraIdx)){
            effSd[j] <- baseSd * exp(logXtraSd[xtraIdx])
          }else if(!is.na(link) && link >= 1){
            k <- base + (exp(predVarObs[link]) - 1) * logPred[i]
            effSd[j] <- sqrt(log(1 + exp(k)))
          }else{
            effSd[j] <- baseSd
          }
        }
        if(is.null(thisCor) && covType[ff] == 0){
          jnll <- jnll - sum(dnorm(logobs[idx], logPred[idx], effSd, log = TRUE))
        }else{
          Sigma <- outer(effSd, effSd) * R
          jnll <- jnll - dmvnorm(logobs[idx], logPred[idx], Sigma, log = TRUE)
        }
      }
    }
    if(any(isTag) && length(logitReleaseSurvival) > 0){
      maxAgeIdx <- length(age)
      for(i in which(isTag)){
        y     <- aux[i, 1] - minYear + 1
        a_raw <- aux[i, 3] - minAge + 1
        a     <- if(a_raw > maxAgeIdx) maxAgeIdx else a_raw
        ti    <- tagTypeIdx[i]
        log_pred <- log(tagR[i]) + log(tagNscan[i]) - logN[y, a] - log(1000) + log(plogis(logitReleaseSurvival[ti]))
        log_var_minus_mu <- log_pred - logitRecapturePhi[ti]
        jnll <- jnll - dnbinom_robust(logobs[i], log_pred, log_var_minus_mu, log = TRUE)
      }
    }

    N <- exp(logN)
    Z <- F + M
    cAA <- N * (F / Z) * (1 - exp(-Z))
    CW0 <- ifelse(is.na(CW), 0, CW)
    catchInWeight <- as.vector((cAA * CW0) %*% rep(1, ncol(CW0)))
    logCatch <- log(catchInWeight)
    logssb <- log(ssbFun(N, F, M, SW, MO, PF, PM))
    logR <- logN[, 1]
    fbarMat <- F[, fbarIdx, drop = FALSE]
    fbar <- as.vector(fbarMat %*% rep(1 / length(fbarIdx), length(fbarIdx)))
    logfbar <- log(fbar)
    SW0 <- ifelse(is.na(SW), 0, SW)
    tsbVec <- as.vector((N * SW0) %*% rep(1, ncol(SW0)))
    logtsb <- log(tsbVec)

    REPORT(Slist)
    REPORT(logPred)
    ADREPORT(logssb)
    ADREPORT(logfbar)
    ADREPORT(logR)
    ADREPORT(logCatch)
    ADREPORT(logtsb)
    if(resFlag == 1L){
      ADREPORT(resN)
      ADREPORT(resF)
    }
    jnll
  }
}

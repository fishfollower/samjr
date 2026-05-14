##' Short-term stochastic forecast for a samjr fit
##'
##' Runs forward projections from the last assessment year using simple
##' age-structured population dynamics. F is determined per year by exactly
##' one of \code{fscale}, \code{fval}, \code{catchval} or \code{landval};
##' biology (M, SW, MO, PF, PM, LF, LW) is averaged over \code{ave.years}
##' (or taken from observations when present); recruitment is sampled
##' (or averaged in \code{deterministic} mode) from the assessment-estimated
##' recruits in \code{rec.years}.
##'
##' Implemented options: \code{fscale}, \code{fval}, \code{catchval},
##' \code{landval}; \code{deterministic}; \code{processNoiseF};
##' \code{useSWmodel} (uses the GMRF-smoothed stock weights from the fit
##' when \code{conf$stockWeightModel >= 1}); \code{useCWmodel} (likewise for
##' catch weights). Out of scope: \code{useMOmodel} / \code{useNMmodel}
##' (no maturity / mortality process model in samjr yet),
##' \code{nextssb}, \code{cwF}, \code{customWeights}, \code{customSel},
##' \code{overwriteSelYears}, multi-fleet, exact-match scenarios.
##'
##' @param fit fitted \code{sam} object from \code{\link{sam.fit}}.
##' @param fscale,fval,catchval,landval scenario specifications; for each
##' forecast year exactly one of them must be non-\code{NA}. The scenario
##' length determines the forecast horizon.
##' @param nosim number of stochastic simulations (default 1000). Ignored when
##' \code{deterministic = TRUE}.
##' @param year.base last observed year used as the starting state.
##' @param ave.years years over which biology is averaged.
##' @param rec.years years over which recruitment is sampled.
##' @param label optional label.
##' @param deterministic if \code{TRUE} all stochastic noise is turned off
##' (recruitment becomes the mean of \code{rec.years}, no F process noise).
##' @param processNoiseF if \code{TRUE} add F-process noise from the fitted
##' \code{logsdF} between forecast years.
##' @param useSWmodel,useCWmodel use the smoothed stock/catch weights from
##' the assessment model instead of the year-by-year averages.
##' @param useMOmodel,useNMmodel currently must be \code{FALSE}; a
##' helpful error is raised otherwise.
##' @param splitLD include landing/discard split columns in the table.
##' @param estimate summary function applied per-year (default
##' \code{\link[stats]{median}}).
##' @param ... unused.
##' @return an object of class \code{samforecast} - a list of per-year
##' simulation outputs. Two key attributes:
##' \describe{
##'   \item{\code{tab}}{matrix with Estimate / low / high triples per year
##'     for fbar, rec, ssb and catch (and Land/Discard if \code{splitLD}).}
##'   \item{\code{shorttab}}{the Estimate-only block of \code{tab},
##'     transposed; same format as SAM's.}
##' }
##' @export
forecast <- function(fit, ...) UseMethod("forecast")

##' Internal: re-run the inner Laplace approximation with the GMRF
##' biology process matrices extended by \code{ns} extra rows beyond the
##' fit's existing spin-out window. Pushes the GMRF boundary \code{ns} rows
##' past the forecast horizon so smoothed values for the forecast years are
##' not pulled by edge effects (mirrors SAM's
##' \code{forecast.R::if(useSWmodel)\{...MakeADFun...sdreport...\}} block).
##'
##' Returns a parList-shaped list with the (extended) smoothed latent
##' matrices. The fixed effects are held at \code{fit$opt$par}; only the
##' random effects (logSW, logCW, logitMO, logNM, plus the original
##' \code{logN}, \code{logF}, \code{missing}) are re-optimised by the
##' inner Laplace step.
##' @keywords internal
##' @noRd
##' Internal: extract a parList-shaped list from a samjr fit using the
##' parameter NAMES in \code{last.par.best} rather than positional unpacking.
##' RTMB's \code{parList(vec)} indexes its input by lfixed positions, which
##' for the FULL fixed+random vector \code{last.par.best} reads the wrong
##' elements when many random effects exist.
##' @keywords internal
##' @noRd
fitParList <- function(fit){
  template <- fit$parameters
  optpar <- fit$opt$par
  lp     <- fit$obj$env$last.par.best
  for(nm in names(template)){
    if(length(template[[nm]]) == 0) next
    v <- if(nm %in% names(optpar)) optpar[names(optpar) == nm]
         else                       lp[names(lp) == nm]
    if(length(v) == length(template[[nm]])) template[[nm]][] <- v
  }
  template
}

samRmvnorm <- function(n, mu, Sigma){
  p <- length(mu)
  idx <- diag(Sigma) > .Machine$double.xmin
  L <- matrix(0, p, p)
  if(any(idx)) L[idx, idx] <- chol(Sigma[idx, idx, drop = FALSE])
  X <- matrix(stats::rnorm(p * n), n, p)
  X %*% L + matrix(mu, n, p, byrow = TRUE)
}

extendForBoundary <- function(fit, ns){
  oldPl <- fitParList(fit)
  active <- list(SW = fit$conf$stockWeightModel >= 1,
                 CW = fit$conf$catchWeightModel >= 1,
                 MO = fit$conf$matureModel       >= 1,
                 NM = fit$conf$mortalityModel    >= 1)
  if(!any(unlist(active))) return(list(par = oldPl, obj = NULL))

  origSpin <- if(active$SW) nrow(oldPl$logSW)
              else if(active$CW) nrow(oldPl$logCW)
              else if(active$MO) nrow(oldPl$logitMO)
              else nrow(oldPl$logNM)
  origSpin <- origSpin - length(fit$data$years)
  newSpin  <- origSpin + ns

  newDat <- toBabyDat(fit$data, fit$conf, spinoutyear = newSpin)

  newPar <- oldPl
  pad <- function(M, k) rbind(M, matrix(0, nrow = k, ncol = ncol(M)))
  if(active$SW) newPar$logSW   <- pad(oldPl$logSW,   ns)
  if(active$CW) newPar$logCW   <- pad(oldPl$logCW,   ns)
  if(active$MO) newPar$logitMO <- pad(oldPl$logitMO, ns)
  if(active$NM) newPar$logNM   <- pad(oldPl$logNM,   ns)

  newF <- makeBabyLikelihood(newDat)
  randomVars <- c("logN", "logF", "missing")
  if(active$SW) randomVars <- c(randomVars, "logSW")
  if(active$CW) randomVars <- c(randomVars, "logCW")
  if(active$MO) randomVars <- c(randomVars, "logitMO")
  if(active$NM) randomVars <- c(randomVars, "logNM")
  newObj <- MakeADFun(newF, newPar, random = randomVars, silent = TRUE)
  newObj$fn(fit$opt$par)
  lp <- newObj$env$last.par
  out <- newPar
  for(nm in names(out)){
    if(length(out[[nm]]) == 0) next
    v <- if(nm %in% names(fit$opt$par)) fit$opt$par[names(fit$opt$par) == nm]
         else                            lp[names(lp) == nm]
    if(length(v) == length(out[[nm]])) out[[nm]][] <- v
  }
  list(par = out, obj = newObj)
}

##' @rdname forecast
##' @method forecast sam
##' @export
forecast.sam <- function(fit,
                         fscale = NULL, fval = NULL,
                         catchval = NULL, landval = NULL,
                         nosim = 1000,
                         year.base = max(fit$data$years),
                         ave.years = year.base + (-4:0),
                         rec.years = year.base + (-9:0),
                         label = NULL,
                         deterministic = FALSE,
                         processNoiseF = TRUE,
                         useSWmodel = (fit$conf$stockWeightModel >= 1),
                         useCWmodel = (fit$conf$catchWeightModel >= 1),
                         useMOmodel = (fit$conf$matureModel    >= 1),
                         useNMmodel = (fit$conf$mortalityModel >= 1),
                         splitLD = FALSE,
                         estimate = stats::median,
                         ...){

  ns <- max(c(length(fscale), length(fval), length(catchval), length(landval)))
  if(ns == 0) stop("No scenario specified")
  if(is.null(fscale))   fscale   <- rep(NA_real_, ns)
  if(is.null(fval))     fval     <- rep(NA_real_, ns)
  if(is.null(catchval)) catchval <- rep(NA_real_, ns)
  if(is.null(landval))  landval  <- rep(NA_real_, ns)
  spec <- cbind(fscale, fval, catchval, landval)
  if(!all(rowSums(!is.na(spec)) == 1))
    stop("For each forecast year exactly one of fscale, fval, catchval, landval must be specified")

  conf <- fit$conf
  data <- fit$data
  dat  <- fit$dat
  keyF <- dat$keyF
  fbarIdx <- dat$fbarIdx
  ages <- dat$age
  nages <- length(ages)

  yIdx <- which(data$years == year.base)
  if(length(yIdx) == 0) stop("year.base must be in fit$data$years")
  pl <- fitParList(fit)

  baseN <- exp(pl$logN[yIdx, ])
  baseF <- exp(pl$logF[yIdx, ])[keyF]

  doAve <- function(x){
    if(is.null(rownames(x))) return(rep(NA_real_, ncol(x)))
    rs <- as.integer(rownames(x)) %in% ave.years
    if(!any(rs)) return(rep(NA_real_, ncol(x)))
    colMeans(x[rs, , drop = FALSE], na.rm = TRUE)
  }
  ave_sw <- doAve(data$stockMeanWeight)
  ave_cw <- doAve(data$catchMeanWeight)
  ave_mo <- doAve(data$propMat)
  ave_nm <- doAve(data$natMor)
  ave_pf <- doAve(data$propF)
  ave_pm <- doAve(data$propM)
  ave_lf <- doAve(data$landFrac)
  ave_lw <- doAve(data$landMeanWeight)
  getY <- function(x, y, ave){
    if(is.null(rownames(x))) return(ave)
    if(!(y %in% as.integer(rownames(x)))) return(ave)
    ret <- x[as.integer(rownames(x)) == y, ]
    ret[is.na(ret)] <- ave[is.na(ret)]
    ret
  }

  needExtend <- useSWmodel || useCWmodel || useMOmodel || useNMmodel
  ext   <- if(needExtend) extendForBoundary(fit, ns) else list(par = pl, obj = NULL)
  extPl <- ext$par
  yearStart <- as.integer(min(data$years))
  yearTag <- function(M) yearStart:(yearStart + nrow(M) - 1)
  if(useSWmodel){
    if(conf$stockWeightModel == 0)
      stop("useSWmodel cannot be used: stockWeightModel was not part of the fit")
    smoothSW <- exp(extPl$logSW); rownames(smoothSW) <- yearTag(smoothSW)
  }
  if(useCWmodel){
    if(conf$catchWeightModel == 0)
      stop("useCWmodel cannot be used: catchWeightModel was not part of the fit")
    smoothCW <- exp(extPl$logCW); rownames(smoothCW) <- yearTag(smoothCW)
  }
  if(useMOmodel){
    if(conf$matureModel == 0)
      stop("useMOmodel cannot be used: matureModel was not part of the fit")
    smoothMO <- plogis(extPl$logitMO); rownames(smoothMO) <- yearTag(smoothMO)
  }
  if(useNMmodel){
    if(conf$mortalityModel == 0)
      stop("useNMmodel cannot be used: mortalityModel was not part of the fit")
    smoothNM <- exp(extPl$logNM); rownames(smoothNM) <- yearTag(smoothNM)
  }

  simLogSW <- simLogCW <- simLogitMO <- simLogNM <- NULL
  sampleBio <- needExtend && !deterministic && (as.integer(nosim) > 1) &&
               !is.null(ext$obj)
  if(sampleBio){
    sdr_bio <- try(RTMB::sdreport(ext$obj, par.fixed = fit$opt$par,
                                   ignore.parm.uncertainty = TRUE,
                                   getReportCovariance = TRUE), silent = TRUE)
    ok <- !inherits(sdr_bio, "try-error") && !is.null(sdr_bio$cov)
    if(ok){
      val_full <- sdr_bio$value
      cov_full <- sdr_bio$cov
      nm_full  <- names(val_full)
      Nsim_pre <- as.integer(nosim)
      sampleOne <- function(name){
        idx <- which(nm_full == name)
        if(length(idx) == 0) return(NULL)
        Sigma <- as.matrix(cov_full[idx, idx, drop = FALSE])
        Sigma <- (Sigma + t(Sigma)) / 2
        samRmvnorm(Nsim_pre, val_full[idx], Sigma)
      }
      if(useSWmodel) simLogSW   <- sampleOne("logSW")
      if(useCWmodel) simLogCW   <- sampleOne("logCW")
      if(useNMmodel) simLogNM   <- sampleOne("logNM")
      if(useMOmodel) simLogitMO <- sampleOne("logitMO")
    }
  }

  recAll <- exp(pl$logN[, 1])
  recpool <- recAll[as.integer(data$years) %in% rec.years]
  if(length(recpool) == 0) stop("No recruitment available for the requested rec.years")

  nFstate <- length(baseF[!duplicated(keyF)])
  if(processNoiseF && !deterministic){
    sdF_state <- exp(pl$logsdF)[dat$keyVarFperState]
    nFstate <- length(sdF_state)
    rhoF <- 2 * stats::plogis(if(length(pl$itrans_rho) == 0) 0 else pl$itrans_rho[1]) - 1
    if(dat$fcormode == 0){
      corMatF <- diag(nFstate)
    }else if(dat$fcormode == 1){
      corMatF <- matrix(rhoF, nFstate, nFstate); diag(corMatF) <- 1
    }else if(dat$fcormode == 2){
      d <- abs(outer(seq_len(nFstate), seq_len(nFstate), "-"))
      corMatF <- rhoF^d
    }else corMatF <- diag(nFstate)
    SigmaF <- outer(sdF_state, sdF_state) * corMatF
    cholF <- tryCatch(chol(SigmaF), error = function(e) NULL)
  }else cholF <- NULL

  Nsim <- if(deterministic) 1L else as.integer(nosim)
  sampleStart <- !deterministic && Nsim > 1 &&
                 !inherits(fit$sdrep, "try-error") &&
                 !is.null(fit$sdrep$jointPrecision)
  Fmat_state <- NULL
  if(sampleStart){
    Q  <- fit$sdrep$jointPrecision
    nm <- rownames(Q)
    posLogN <- which(nm == "logN")
    posLogF <- which(nm == "logF")
    nyear <- length(data$years)
    idxN <- posLogN[yIdx + (seq_len(nages)   - 1L) * nyear]
    idxF <- posLogF[yIdx + (seq_len(nFstate) - 1L) * nyear]
    idx  <- c(idxN, idxF)
    mu   <- c(pl$logN[yIdx, ], pl$logF[yIdx, ])
    e <- Matrix::sparseMatrix(i = idx, j = seq_along(idx), x = 1,
                              dims = c(nrow(Q), length(idx)))
    Sigma <- as.matrix(Matrix::solve(Q, e)[idx, , drop = FALSE])
    Sigma <- (Sigma + t(Sigma)) / 2
    draws <- samRmvnorm(Nsim, mu, Sigma)
    Nmat       <- exp(draws[, seq_len(nages), drop = FALSE])
    Fmat_state <- exp(draws[, nages + seq_len(nFstate), drop = FALSE])
    Fmat       <- Fmat_state[, keyF, drop = FALSE]
  }else{
    Nmat <- matrix(rep(baseN, each = Nsim), nrow = Nsim, ncol = nages)
    Fmat <- matrix(rep(baseF, each = Nsim), nrow = Nsim, ncol = nages)
  }

  resampleN <- function(n){
    if(deterministic) rep(mean(recpool), n) else sample(recpool, n, replace = TRUE)
  }
  applyScenario <- function(F_mat, scenario, sw, cw, nm, lf, lw, Nmat){
    if(!is.na(scenario["fscale"])) return(F_mat * scenario["fscale"])
    if(!is.na(scenario["fval"])){
      curfbar <- mean(F_mat[1, fbarIdx])
      return(F_mat * (scenario["fval"] / curfbar))
    }
    Nmean <- colMeans(Nmat)
    Fbase <- F_mat[1, ]
    if(!is.na(scenario["catchval"])){
      f_root <- function(s){
        Fnew <- Fbase * s
        Z <- Fnew + nm
        sum(Nmean * Fnew / Z * (1 - exp(-Z)) * cw) - scenario["catchval"]
      }
      s <- stats::uniroot(f_root, c(1e-3, 100))$root
      return(F_mat * s)
    }
    if(!is.na(scenario["landval"])){
      f_root <- function(s){
        Fnew <- Fbase * s
        Z <- Fnew + nm
        sum(Nmean * Fnew / Z * (1 - exp(-Z)) * cw * lf) - scenario["landval"]
      }
      s <- stats::uniroot(f_root, c(1e-3, 100))$root
      return(F_mat * s)
    }
    F_mat
  }

  simlist <- list()
  prevZ <- NULL
  for(i in seq_len(ns)){
    y <- year.base + (i - 1)
    sw <- if(useSWmodel) getY(smoothSW, y, ave_sw) else getY(data$stockMeanWeight, y, ave_sw)
    cw <- if(useCWmodel) getY(smoothCW, y, ave_cw) else getY(data$catchMeanWeight, y, ave_cw)
    mo <- if(useMOmodel) getY(smoothMO, y, ave_mo) else getY(data$propMat, y, ave_mo)
    nm <- if(useNMmodel) getY(smoothNM, y, ave_nm) else getY(data$natMor, y, ave_nm)
    pf <- getY(data$propF, y, ave_pf)
    pm <- getY(data$propM, y, ave_pm)
    lf <- getY(data$landFrac, y, ave_lf)
    lw <- getY(data$landMeanWeight, y, ave_lw)

    if(i > 1){
      Nnew <- matrix(0, Nsim, nages)
      if(deterministic){
        Nnew[, 1] <- mean(recpool)
      }else{
        recIdx <- vapply(seq_len(Nsim),
                         function(.) sample.int(length(recpool), 1L),
                         integer(1))
        Nnew[, 1] <- recpool[recIdx]
      }
      survN <- Nmat * exp(-prevZ)
      if(nages >= 3){
        for(a in 2:(nages - 1)) Nnew[, a] <- survN[, a - 1]
      }else if(nages == 2){
        Nnew[, 2] <- survN[, 1]
      }
      Nnew[, nages] <- survN[, nages - 1] + survN[, nages]
      Nmat <- Nnew

      if(!deterministic){
        nvar <- diag(c(0, rep(exp(2 * pl$logSdLogN[2]), nages - 1)),
                     nrow = nages, ncol = nages)
        if(processNoiseF && !is.null(cholF)){
          fvar <- crossprod(cholF)
        }else{
          fvar <- matrix(0, nFstate, nFstate)
        }
        procVar <- as.matrix(Matrix::bdiag(nvar, fvar))
        eps <- samRmvnorm(Nsim, rep(0, nages + nFstate), procVar)
        Nmat       <- exp(log(Nmat) + eps[, seq_len(nages), drop = FALSE])
        if(!is.null(Fmat_state)){
          Fmat_state <- exp(log(Fmat_state) +
                            eps[, nages + seq_len(nFstate), drop = FALSE])
          Fmat <- Fmat_state[, keyF, drop = FALSE]
        }else if(processNoiseF && !is.null(cholF)){
          Fmat <- exp(log(Fmat) +
                      eps[, nages + dat$keyVarFperState[keyF], drop = FALSE])
        }
      }
    }

    swMat <- matrix(rep(sw, each = Nsim), nrow = Nsim)
    cwMat <- matrix(rep(cw, each = Nsim), nrow = Nsim)
    moMat <- matrix(rep(mo, each = Nsim), nrow = Nsim)
    nmMat <- matrix(rep(nm, each = Nsim), nrow = Nsim)
    if(sampleBio){
      pickRow <- function(simMat){
        if(is.null(simMat)) return(NULL)
        nyTraj <- ncol(simMat) %/% nages
        i_y <- y - yearStart + 1L
        if(i_y < 1L || i_y > nyTraj) return(NULL)
        simMat[, i_y + (0:(nages - 1L)) * nyTraj, drop = FALSE]
      }
      sw_traj <- pickRow(simLogSW);   if(!is.null(sw_traj)) swMat <- exp(sw_traj)
      cw_traj <- pickRow(simLogCW);   if(!is.null(cw_traj)) cwMat <- exp(cw_traj)
      mo_traj <- pickRow(simLogitMO); if(!is.null(mo_traj)) moMat <- stats::plogis(mo_traj)
      nm_traj <- pickRow(simLogNM);   if(!is.null(nm_traj)) nmMat <- exp(nm_traj)
    }
    pfMat <- matrix(rep(pf, each = Nsim), nrow = Nsim)
    pmMat <- matrix(rep(pm, each = Nsim), nrow = Nsim)
    lfMat <- matrix(rep(lf, each = Nsim), nrow = Nsim)

    scen <- c(fscale = fscale[i], fval = fval[i],
              catchval = catchval[i], landval = landval[i])
    Fmat <- applyScenario(Fmat, scen, sw, cw, nm, lf, lw, Nmat)
    if(!is.null(Fmat_state)){
      firstIdxByState <- match(seq_len(nFstate), keyF)
      Fmat_state <- Fmat[, firstIdxByState, drop = FALSE]
    }

    Z <- Fmat + nmMat
    prevZ <- Z

    fbarsim  <- rowMeans(Fmat[, fbarIdx, drop = FALSE])
    catchsim <- rowSums(Nmat * Fmat / Z * (1 - exp(-Z)) * cwMat)
    landsim  <- rowSums(Nmat * Fmat / Z * (1 - exp(-Z)) * cwMat * lfMat)
    Zssb     <- pmMat * nmMat + Fmat * pfMat
    ssbsim   <- rowSums(Nmat * exp(-Zssb) * moMat * swMat)
    recsim   <- Nmat[, 1]

    simlist[[i]] <- list(year = y, fbar = fbarsim, rec = recsim,
                         ssb = ssbsim, catch = catchsim, land = landsim,
                         sim = cbind(Nmat, Fmat))
  }

  collect <- function(x){
    if(length(x) == 1) return(c(Estimate = x, low = x, high = x))
    quan <- stats::quantile(x, c(0.025, 0.975), names = FALSE)
    c(Estimate = unname(estimate(x)), low = quan[1], high = quan[2])
  }
  fbar  <- round(do.call(rbind, lapply(simlist, function(xx) collect(xx$fbar))), 3)
  rec   <- round(do.call(rbind, lapply(simlist, function(xx) collect(xx$rec))))
  ssb   <- round(do.call(rbind, lapply(simlist, function(xx) collect(xx$ssb))))
  catch <- round(do.call(rbind, lapply(simlist, function(xx) collect(xx$catch))))
  land  <- round(do.call(rbind, lapply(simlist, function(xx) collect(xx$land))))
  tab <- cbind(fbar, rec, ssb, catch)
  basenames <- c("fbar:", "rec:", "ssb:", "catch:")
  if(splitLD){
    tab <- cbind(tab, land, catch - land)
    basenames <- c(basenames, "Land:", "Discard:")
  }
  colnames(tab) <- paste0(rep(basenames, each = 3), c("Estimate", "low", "high"))
  rownames(tab) <- vapply(simlist, function(xx) as.character(xx$year), character(1))
  shorttab <- t(tab[, grepl("Estimate$", colnames(tab)), drop = FALSE])
  rownames(shorttab) <- sub(":Estimate$", "",
                            paste0(if(!is.null(label)) paste0(label, ":") else "",
                                   rownames(shorttab)))

  attr(simlist, "fit")      <- fit
  attr(simlist, "tab")      <- tab
  attr(simlist, "shorttab") <- shorttab
  attr(simlist, "label")    <- label
  class(simlist) <- "samforecast"
  simlist
}

##' Print a samjr forecast
##' @param x a \code{samforecast} object.
##' @param ... unused.
##' @method print samforecast
##' @export
print.samforecast <- function(x, ...){
  cat("samjr forecast over", length(x), "year(s)\n")
  print(attr(x, "shorttab"))
  invisible(x)
}

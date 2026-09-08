##' Reduce a samjr data list by removing observations
##'
##' Removes observations matching any of the supplied
##' \code{(year, fleet, age)} triples from a \code{sam_data} list and
##' rebuilds the per-(fleet, year) index matrices, biology arrays and
##' fleet metadata. The optional \code{conf} argument is reindexed in
##' parallel and attached as \code{attr(., "conf")}.
##'
##' Mirrors SAM's \code{stockassessment::reduce}, simplified for
##' samjr's single-residual-fleet 2D layout.
##'
##' @param data a \code{sam_data} list as returned by \code{\link{setup.sam.data}}.
##' @param year vector of years to drop.
##' @param fleet vector of fleets to drop. When both \code{year} and
##' \code{fleet} are given they are treated as paired - only the
##' \code{(year, fleet)} cells listed get dropped.
##' @param age vector of ages to drop (paired in the same way).
##' @param conf optional configuration list to be reindexed alongside.
##' @return the reduced data list (with optional \code{attr(., "conf")}).
##' @export
reduce <- function(data, year = NULL, fleet = NULL, age = NULL, conf = NULL){
  nam <- c("year", "fleet", "age")[c(length(year) > 0,
                                       length(fleet) > 0,
                                       length(age) > 0)]
  if(length(nam) == 0){
    idx <- rep(TRUE, nrow(data$aux))
  }else{
    keyDat <- do.call(paste, as.data.frame(data$aux[, nam, drop = FALSE]))
    drop   <- do.call(paste, as.data.frame(cbind(year = year, fleet = fleet, age = age)))
    idx    <- !(keyDat %in% drop)
  }
  data$aux    <- data$aux[idx, , drop = FALSE]
  data$logobs <- data$logobs[idx]
  if(!is.null(data$weight)) data$weight <- data$weight[idx]
  suf <- sort(unique(data$aux[, "fleet"]))
  data$noFleets    <- length(suf)
  data$fleetTypes  <- data$fleetTypes[suf]
  data$sampleTimes <- data$sampleTimes[suf]
  oldYears <- data$years
  data$years <- min(as.numeric(data$aux[, "year"])):max(as.numeric(data$aux[, "year"]))
  ages <- min(as.numeric(data$aux[, "age"])):max(as.numeric(data$aux[, "age"]))
  data$noYears <- length(data$years)
  mmfun <- function(f, y, ff){
    idx <- which(data$aux[, "year"] == y & data$aux[, "fleet"] == f)
    if(length(idx) == 0) NA_integer_ else ff(idx) - 1L
  }
  data$idx1 <- outer(suf, data$years, Vectorize(mmfun, c("f", "y")), ff = min)
  data$idx2 <- outer(suf, data$years, Vectorize(mmfun, c("f", "y")), ff = max)
  if(!is.null(data$idxCor))
    data$idxCor <- data$idxCor[suf, match(data$years, oldYears), drop = FALSE]
  data$nobs <- length(data$logobs)

  cutMat <- function(M){
    if(is.null(M)) return(NULL)
    M[rownames(M) %in% data$years, colnames(M) %in% ages, drop = FALSE]
  }
  data$propMat         <- cutMat(data$propMat)
  data$stockMeanWeight <- cutMat(data$stockMeanWeight)
  data$natMor          <- cutMat(data$natMor)
  data$propF           <- cutMat(data$propF)
  data$propM           <- cutMat(data$propM)
  data$catchMeanWeight <- cutMat(data$catchMeanWeight)
  data$disMeanWeight   <- cutMat(data$disMeanWeight)
  data$landMeanWeight  <- cutMat(data$landMeanWeight)
  data$landFrac        <- cutMat(data$landFrac)

  ## close-kin pairs referring to sampling years that no longer exist must go;
  ## pairs whose birth year falls before the (possibly later) model start are
  ## dropped downstream by ckmrPrep.
  if(!is.null(data$ckmr)){
    ck <- data$ckmr
    keep <- (ck$year1 %in% data$years) & (ck$year2 %in% data$years)
    data$ckmr <- if(any(keep)) ck[keep, , drop = FALSE] else NULL
  }

  data$aux[, "fleet"] <- match(data$aux[, "fleet"], suf)
  data$minAgePerFleet <- tapply(as.integer(data$aux[, "age"]),
                                INDEX = data$aux[, "fleet"], FUN = min)
  data$maxAgePerFleet <- tapply(as.integer(data$aux[, "age"]),
                                INDEX = data$aux[, "fleet"], FUN = max)
  attr(data, "fleetNames") <- attr(data, "fleetNames")[suf]

  if(!is.null(conf)){
    .reidx <- function(x){
      if(any(x >= -0.5 & !is.na(x))){
        xx <- x[x >= -0.5 & !is.na(x)]
        x[x >= -0.5 & !is.na(x)] <- match(xx, sort(unique(xx))) - 1L
      }
      x
    }
    if(!is.null(conf$keyLogFsta))    conf$keyLogFsta    <- .reidx(conf$keyLogFsta[suf, , drop = FALSE])
    if(!is.null(conf$keyLogFpar))    conf$keyLogFpar    <- .reidx(conf$keyLogFpar[suf, , drop = FALSE])
    if(!is.null(conf$keyQpow))       conf$keyQpow       <- .reidx(conf$keyQpow[suf, , drop = FALSE])
    if(!is.null(conf$keyVarF))       conf$keyVarF       <- .reidx(conf$keyVarF[suf, , drop = FALSE])
    if(!is.null(conf$keyVarObs))     conf$keyVarObs     <- .reidx(conf$keyVarObs[suf, , drop = FALSE])
    if(!is.null(conf$predVarObsLink))conf$predVarObsLink<- .reidx(conf$predVarObsLink[suf, , drop = FALSE])
    if(!is.null(conf$obsCorStruct))  conf$obsCorStruct  <- conf$obsCorStruct[suf]
    if(!is.null(conf$maxAgePlusGroup))conf$maxAgePlusGroup <- conf$maxAgePlusGroup[suf]
    if(!is.null(conf$keyCorObs))     conf$keyCorObs     <- .reidx(conf$keyCorObs[suf, , drop = FALSE])
    if(!is.null(conf$keyScaledYears) && length(conf$keyScaledYears) > 0){
      yidx <- conf$keyScaledYears %in% data$aux[data$aux[, "fleet"] == 1, "year"]
      conf$noScaledYears  <- sum(yidx)
      conf$keyScaledYears <- as.vector(conf$keyScaledYears)[yidx]
      if(!is.null(conf$keyParScaledYA))
        conf$keyParScaledYA <- .reidx(conf$keyParScaledYA[yidx, , drop = FALSE])
    }
    if(!is.null(conf$fixVarToWeight))   conf$fixVarToWeight  <- conf$fixVarToWeight[suf]
    attr(data, "conf") <- conf
  }
  data
}

##' Re-fit a samjr model with some observations excluded
##'
##' Drops the supplied \code{(year, fleet)} pairs and re-fits, carrying
##' over the previous fit's fixed-effect estimates as starting values.
##' Used by \code{\link{retro}} and \code{\link{leaveout}}.
##'
##' @param fit a fitted \code{sam} object.
##' @param year,fleet vectors of years / fleets to drop (paired).
##' @param map optional \code{map} list forwarded to \code{\link{sam.fit}};
##' defaults to \code{fit$map}.
##' @param lower,upper optional named numeric vectors of parameter bounds
##' forwarded to \code{\link{sam.fit}}; default to \code{fit$low} and
##' \code{fit$hig}.
##' @param silent forwarded to \code{\link{sam.fit}}; defaults to
##' \code{TRUE} so re-fits driven by \code{\link{retro}} or
##' \code{\link{leaveout}} stay quiet.
##' @param ... extra arguments forwarded to \code{\link{sam.fit}}.
##' @return a new \code{sam} fit.
##' @export
runwithout <- function(fit, year, fleet, ...) UseMethod("runwithout")

##' @rdname runwithout
##' @method runwithout sam
##' @export
runwithout.sam <- function(fit, year = NULL, fleet = NULL,
                            map = fit$map, lower = fit$low, upper = fit$hig,
                            silent = TRUE, ...){
  data <- reduce(fit$data, year = year, fleet = fleet, conf = fit$conf)
  conf <- attr(data, "conf")
  fakefile <- file()
  sink(fakefile); saveConf(conf, file = ""); sink()
  conf <- loadConf(data, fakefile, patch = TRUE)
  close(fakefile)
  par <- defpar(data, conf)
  carry <- intersect(names(par), names(fit$pl))
  carry <- setdiff(carry, c("logN", "logF", "logSW", "logCW",
                            "logitMO", "logNM", "missing"))
  for(nm in carry){
    if(length(par[[nm]]) == length(fit$pl[[nm]])) par[[nm]] <- fit$pl[[nm]]
  }
  sam.fit(data, conf, par, map = map, lower = lower, upper = upper,
          silent = silent, ...)
}

##' Retrospective analysis of a samjr fit
##'
##' Refits the assessment after peeling 1, 2, ..., \code{year} years
##' from the right edge of the time series. Returns a \code{samset}.
##'
##' @param fit a fitted \code{sam} object.
##' @param year either a single integer \code{n} (peel 1..n years) or a
##' vector / matrix of years to drop (see SAM's \code{retro} for details).
##' @param ncores ignored in samjr (always serial).
##' @param ... extra arguments forwarded to \code{\link{sam.fit}}.
##' @return a list of fits with class \code{samset}.
##' @export
retro <- function(fit, year = NULL, ncores = 1, ...) UseMethod("retro")

##' @rdname retro
##' @method retro sam
##' @export
retro.sam <- function(fit, year = NULL, ncores = 1, ...){
  data <- fit$data
  y <- data$aux[, "year"]; f <- data$aux[, "fleet"]
  suf <- sort(unique(f))
  maxy <- sapply(suf, function(ff) max(y[f == ff]))
  if(length(year) == 1){
    mat <- sapply(suf, function(ff){ my <- maxy[ff]; my:(my - year + 1) })
    if(year == 1) mat <- matrix(mat, nrow = 1)
  }else if(is.vector(year) && length(year) > 1){
    mat <- sapply(suf, function(ff) year)
  }else if(is.matrix(year)){
    mat <- year
  }else stop("retro: 'year' must be scalar, vector or matrix")
  if(nrow(mat) > length(unique(y))) stop("retro: too many runs")
  if(ncol(mat) != length(suf)) stop("retro: wrong number of fleet columns")
  setup <- lapply(seq_len(nrow(mat)), function(i)
    do.call(rbind, lapply(suf, function(ff)
      if(mat[i, ff] <= maxy[ff]) cbind(mat[i, ff]:maxy[ff], ff))))
  runs <- lapply(setup, function(s) runwithout(fit, year = s[, 1], fleet = s[, 2], ...))
  conv <- vapply(runs, function(x) x$opt$convergence, integer(1))
  if(any(conv != 0))
    warning("retro run(s) ", paste(which(conv != 0), collapse = ","), " did not converge.")
  attr(runs, "fit") <- fit
  class(runs) <- "samset"
  runs
}

##' Leave-one-fleet-out analysis of a samjr fit
##' @param fit a fitted \code{sam} object.
##' @param fleet a list of fleet vectors; one run per element.
##' @param ncores ignored in samjr.
##' @param ... extra arguments forwarded to \code{\link{sam.fit}}.
##' @return a \code{samset}.
##' @export
leaveout <- function(fit, fleet = as.list(2:fit$data$noFleets),
                     ncores = 1, ...){
  runs <- lapply(fleet, function(f) runwithout(fit, fleet = f, ...))
  names(runs) <- paste0("w.o. ",
                        vapply(fleet,
                               function(x) paste(attr(fit$data, "fleetNames")[x], collapse = " and "),
                               character(1)))
  conv <- vapply(runs, function(x) x$opt$convergence, integer(1))
  if(any(conv != 0))
    warning("leaveout run(s) ", paste(which(conv != 0), collapse = ","), " did not converge.")
  attr(runs, "fit") <- fit
  class(runs) <- "samset"
  runs
}

##' Mohn's rho on a samjr retrospective set
##' @param fits a \code{samset} as returned by \code{\link{retro}}.
##' @param what optional function returning a (year x quantity) matrix
##' for a single fit; default uses recruitment, SSB and Fbar.
##' @param lag lag applied to the comparison year.
##' @param ... unused.
##' @return numeric vector of Mohn's rho values.
##' @export
mohn <- function(fits, what = NULL, lag = 0, ...) UseMethod("mohn")

##' @rdname mohn
##' @method mohn samset
##' @export
mohn.samset <- function(fits, what = NULL, lag = 0, ...){
  if(is.null(what)){
    what <- function(fit){
      ret <- cbind(rectable(fit)[, 1], ssbtable(fit)[, 1], fbartable(fit)[, 1])
      colnames(ret) <- c(paste0("R(age ", fit$conf$minAge, ")"), "SSB",
                          paste0("Fbar(", fit$conf$fbarRange[1], "-",
                                 fit$conf$fbarRange[2], ")"))
      ret
    }
  }
  ref <- what(attr(fits, "fit"))
  ret <- lapply(fits, what)
  bias <- lapply(ret, function(x){
    y <- rownames(x)[nrow(x) - lag]
    (x[rownames(x) == y, ] - ref[rownames(ref) == y, ]) /
      ref[rownames(ref) == y, ]
  })
  colMeans(do.call(rbind, bias))
}

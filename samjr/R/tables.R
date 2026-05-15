##' Generic helper that returns an Estimate / Low / High table for a named
##' quantity from a samjr fit (or samjr forecast). Mirrors
##' \code{stockassessment::tableit}.
##' @param fit a fitted \code{sam} object or a \code{samforecast}.
##' @param what name of the quantity (e.g. \code{"logssb"}, \code{"logfbar"},
##' \code{"logR"}, \code{"logCatch"}).
##' @param x row names for the returned table (defaults to \code{fit$data$years}).
##' @param trans transformation applied before tabulation (e.g. \code{exp}).
##' @param ... unused.
##' @return matrix with columns Estimate / Low / High.
##' @export
tableit <- function(fit, what, x = fit$data$years, trans = function(x) x, ...) UseMethod("tableit")

##' @rdname tableit
##' @method tableit sam
##' @export
tableit.sam <- function(fit, what, x = fit$data$years, trans = function(x) x, ...){
  idx <- names(fit$sdrep$value) == what
  if(!any(idx)) stop(paste("Quantity", what, "not found in sdrep"))
  y  <- fit$sdrep$value[idx]
  ci <- y + fit$sdrep$sd[idx] %o% c(-2, 2)
  ret <- trans(cbind(y, ci))
  rownames(ret) <- x
  colnames(ret) <- c("Estimate", "Low", "High")
  ret
}

##' @rdname tableit
##' @method tableit samforecast
##' @export
tableit.samforecast <- function(fit, what, x = fit$data$years, trans = function(x) x, ...){
  tab <- attr(fit, "tab")
  cn  <- switch(what,
                "logssb"   = c("ssb:Estimate",   "ssb:low",   "ssb:high"),
                "logfbar"  = c("fbar:Estimate",  "fbar:low",  "fbar:high"),
                "logCatch" = c("catch:Estimate", "catch:low", "catch:high"),
                "logR"     = c("rec:Estimate",   "rec:low",   "rec:high"),
                stop(paste0("tableit.samforecast: 'what' = '", what,
                            "' not supported")))
  ret <- trans(log(tab[, cn, drop = FALSE]))
  colnames(ret) <- c("Estimate", "Low", "High")
  ret
}

##' Internal: extract a samjr ADREPORT time series and form an
##' Estimate / Low / High table on the natural scale.
##' @keywords internal
##' @noRd
adrepTable <- function(fit, what){
  v <- fit$sdrep$value; s <- fit$sdrep$sd
  idx <- names(v) == what
  if(!any(idx)) stop(paste("Quantity", what, "not found in sdrep"))
  est <- exp(v[idx])
  lo  <- exp(v[idx] - 2 * s[idx])
  hi  <- exp(v[idx] + 2 * s[idx])
  ret <- cbind(Estimate = est, Low = lo, High = hi)
  rownames(ret) <- fit$data$years
  ret
}

##' Time series of estimated spawning stock biomass
##' @param fit a fitted \code{sam} object.
##' @param ... unused.
##' @return matrix with Estimate, Low and High columns and years as row names.
##' @export
ssbtable <- function(fit, ...) UseMethod("ssbtable")

##' @rdname ssbtable
##' @method ssbtable sam
##' @export
ssbtable.sam <- function(fit, ...) adrepTable(fit, "logssb")

##' Time series of estimated mean fishing mortality
##' @param fit a fitted \code{sam} object.
##' @param ... unused.
##' @return matrix with Estimate, Low and High columns and years as row names.
##' @export
fbartable <- function(fit, ...) UseMethod("fbartable")

##' @rdname fbartable
##' @method fbartable sam
##' @export
fbartable.sam <- function(fit, ...) adrepTable(fit, "logfbar")

##' Time series of estimated recruitment
##' @param fit a fitted \code{sam} object.
##' @param ... unused.
##' @return matrix with Estimate, Low and High columns and years as row names.
##' @export
rectable <- function(fit, ...) UseMethod("rectable")

##' @rdname rectable
##' @method rectable sam
##' @export
rectable.sam <- function(fit, ...) adrepTable(fit, "logR")

##' Time series of estimated total catch in weight
##' @param fit a fitted \code{sam} object.
##' @param ... unused.
##' @return matrix with Estimate, Low and High columns and years as row names.
##' @export
catchtable <- function(fit, ...) UseMethod("catchtable")

##' @rdname catchtable
##' @method catchtable sam
##' @export
catchtable.sam <- function(fit, ...){
  ret <- adrepTable(fit, "logCatch")
  CW  <- fit$data$catchMeanWeight
  if(!is.null(CW)){
    keep <- !apply(is.na(CW), 1, any)
    ret  <- ret[keep, , drop = FALSE]
  }
  ret
}

##' Time series of estimated total stock biomass
##' @param fit a fitted \code{sam} object.
##' @param ... unused.
##' @return matrix with Estimate, Low and High columns.
##' @export
tsbtable <- function(fit, ...) UseMethod("tsbtable")

##' @rdname tsbtable
##' @method tsbtable sam
##' @export
tsbtable.sam <- function(fit, ...) adrepTable(fit, "logtsb")

##' Numbers-at-age table from a samjr fit
##' @param fit a fitted \code{sam} object.
##' @param ... unused.
##' @return matrix with year row names and age column names; values are N (linear scale).
##' @export
ntable <- function(fit, ...) UseMethod("ntable")

##' @rdname ntable
##' @method ntable sam
##' @export
ntable.sam <- function(fit, ...){
  ret <- exp(fit$pl$logN)
  rownames(ret) <- fit$data$years
  colnames(ret) <- fit$conf$minAge:fit$conf$maxAge
  ret
}

##' F-at-age table from a samjr fit
##' @param fit a fitted \code{sam} object.
##' @param fleet integer; which catch fleet(s) to include (default = all
##' residual-catch fleets).
##' @param ... unused.
##' @return matrix with year row names and age column names.
##' @export
faytable <- function(fit, ...) UseMethod("faytable")

##' @rdname faytable
##' @method faytable sam
##' @export
faytable.sam <- function(fit, fleet = which(fit$data$fleetTypes == 0), ...){
  getfleet <- function(f){
    idx <- fit$conf$keyLogFsta[f, ] + 2L
    ret <- cbind(NA, exp(fit$pl$logF))[, idx]
    ret[is.na(ret)] <- 0
    ret
  }
  ret <- Reduce("+", lapply(fleet, getfleet))
  rownames(ret) <- fit$data$years
  colnames(ret) <- fit$conf$minAge:fit$conf$maxAge
  ret
}

##' Catch-at-age in numbers from a samjr fit
##' @param fit a fitted \code{sam} object.
##' @param fleet integer; which catch fleet(s) to include (default = all
##' residual-catch fleets).
##' @return matrix with year row names and age column names.
##' @export
caytable <- function(fit, fleet = which(fit$data$fleetTypes == 0)){
  getfleet <- function(f){
    idx <- fit$conf$keyLogFsta[f, ] + 2L
    F <- cbind(NA, exp(fit$pl$logF))[, idx]
    F[is.na(F)] <- 0
    M <- fit$data$natMor
    N <- exp(fit$pl$logN)
    F / (F + M) * N * (1 - exp(-F - M))
  }
  ret <- Reduce("+", lapply(fleet, getfleet))
  rownames(ret) <- fit$data$years
  colnames(ret) <- fit$conf$minAge:fit$conf$maxAge
  ret
}

##' Parameter table from a samjr fit
##' @param fit a fitted \code{sam} object.
##' @param ... unused.
##' @return matrix with par, sd(par), exp(par), Low, High columns.
##' @export
partable <- function(fit, ...) UseMethod("partable")

##' @rdname partable
##' @method partable sam
##' @export
partable.sam <- function(fit, ...){
  param <- coef(fit)
  nam <- names(param)
  dup <- duplicated(nam)
  namadd <- rep(0L, length(nam))
  for(i in 2:length(dup)) if(dup[i]) namadd[i] <- namadd[i - 1] + 1L
  nam <- paste(nam, namadd, sep = "_")
  ret <- cbind(param, attr(param, "sd"))
  ex <- exp(ret[, 1])
  lo <- exp(ret[, 1] - 2 * ret[, 2])
  hi <- exp(ret[, 1] + 2 * ret[, 2])
  ret <- cbind(ret, ex, lo, hi)
  colnames(ret) <- c("par", "sd(par)", "exp(par)", "Low", "High")
  rownames(ret) <- nam
  ret
}

##' Survey catchability (logQ) table from a samjr fit
##' @param fit a fitted \code{sam} object.
##' @param ... unused.
##' @return a matrix with one row per non-residual fleet and one column per
##' age, holding the estimated log-catchability for each (fleet, age) cell
##' that uses one. Cells without a catchability are \code{NA}. The \code{sd}
##' attribute carries the per-cell standard error.
##' @export
qtable <- function(fit, ...) UseMethod("qtable")

##' @rdname qtable
##' @method qtable sam
##' @export
qtable.sam <- function(fit, ...){
  nonres <- fit$data$fleetTypes != 0
  key <- fit$conf$keyLogFpar[nonres, , drop = FALSE] + 1L
  key[key == 0] <- NA
  cf <- coef(fit)
  cfsd <- attr(cf, "sd")
  qvec <- cf[names(cf) == "logQ"]
  qsd  <- cfsd[names(cf) == "logQ"]
  qt <- matrix(qvec[key], nrow = nrow(key), ncol = ncol(key))
  rownames(qt) <- attr(fit$data, "fleetNames")[nonres]
  colnames(qt) <- fit$conf$minAge:fit$conf$maxAge
  sds <- qt; sds[] <- qsd[key]
  attr(qt, "sd") <- sds
  class(qt) <- "samqtable"
  qt
}

##' Print method for \code{samqtable}
##' @param x a \code{samqtable}.
##' @param ... unused.
##' @method print samqtable
##' @export
print.samqtable <- function(x, ...){
  y <- as.matrix(x); attr(y, "sd") <- NULL; class(y) <- "matrix"; print(y)
}

##' Compare a set of fitted samjr models
##'
##' Returns a matrix with one row per fit, columns log-likelihood, degrees
##' of freedom, AIC and (when comparing exactly two nested models) the
##' likelihood-ratio test p-value. Mirrors SAM's \code{modeltable}.
##' @param fits a single \code{sam} fit, a list of fits, or a \code{samset}.
##' @param ... unused.
##' @return a matrix.
##' @importFrom stats AIC
##' @export
modeltable <- function(fits, ...) UseMethod("modeltable")

##' @rdname modeltable
##' @method modeltable sam
##' @export
modeltable.sam <- function(fits, ...) modeltable(c(fits))

##' @rdname modeltable
##' @method modeltable default
##' @export
modeltable.default <- function(fits, ...){
  fits <- fits[!vapply(fits, is.null, logical(1))]
  nam <- if(is.null(names(fits))) paste0("M", seq_along(fits))
         else ifelse(names(fits) == "", paste0("M", seq_along(fits)), names(fits))
  logL <- vapply(fits, function(f) as.numeric(logLik(f)), numeric(1))
  npar <- vapply(fits, function(f) attr(logLik(f), "df"),  integer(1))
  aic  <- vapply(fits, AIC, numeric(1))
  res  <- cbind("log(L)" = logL, "#par" = npar, "AIC" = aic)
  rownames(res) <- nam
  o <- seq_along(fits)
  if(length(fits) == 2){
    o <- order(npar, decreasing = TRUE)
    if(npar[o[1]] > npar[o[2]]){
      df <- npar[o[1]] - npar[o[2]]
      D  <- 2 * (logL[o[1]] - logL[o[2]])
      P  <- 1 - stats::pchisq(D, df)
      cnam <- paste0("Pval( ", nam[o[1]], " -> ", nam[o[2]], " )")
      res <- cbind(res, c(NA, P)[o])
      colnames(res)[ncol(res)] <- cnam
    }
  }
  res[o, , drop = FALSE]
}

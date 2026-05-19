##' Plot a samjr ADREPORT time series
##'
##' Generic for plotting a named ADREPORT quantity from a samjr fit, with
##' an optional two-standard-deviation confidence band. Used internally by
##' \code{\link{ssbplot}}, \code{\link{fbarplot}}, \code{\link{recplot}} and
##' \code{\link{catchplot}}, but exported so it can be called directly on
##' any reported quantity.
##'
##' @param fit an object of class \code{sam} from \code{\link{sam.fit}}.
##' @param what the name of an ADREPORT entry (e.g. \code{"logssb"},
##' \code{"logfbar"}, \code{"logR"}, \code{"logCatch"}).
##' @param x x-axis values (defaults to \code{fit$data$years}).
##' @param ylab,xlab axis labels.
##' @param trans transformation applied before plotting (e.g. \code{exp}).
##' @param add if \code{TRUE}, overlay on an existing plot.
##' @param ci if \code{TRUE}, draw a 2-sd confidence polygon.
##' @param cicol colour for the confidence polygon.
##' @param drop number of years at the right edge to omit.
##' @param xlim,ylim optional axis ranges.
##' @param unnamed.basename legend label used by the \code{samset} method for
##' the base fit when the set has no \code{names()}.
##' @param ... further arguments passed to \code{plot}/\code{lines}.
##' @return invisibly a list with the plotted x, y and CI bands.
##' @export
plotit <- function(fit, what, ...) UseMethod("plotit")

##' @rdname plotit
##' @method plotit sam
##' @export
plotit.sam <- function(fit, what, x = fit$data$years, ylab = what, xlab = "Years",
                       trans = function(x) x, add = FALSE, ci = TRUE,
                       cicol = gray(.5, alpha = .5), drop = 0,
                       xlim = NULL, ylim = NULL, ...){
  idx <- names(fit$sdrep$value) == what
  if(!any(idx)) stop(paste("Quantity", what, "not found in sdrep"))
  y <- fit$sdrep$value[idx]
  s <- fit$sdrep$sd[idx]
  lowhig <- y + s %o% c(-2, 2)
  didx <- 1:(length(x) - drop)
  xr <- if(is.null(xlim)) range(x) else xlim
  x <- x[didx]; y <- y[didx]; lowhig <- lowhig[didx, , drop = FALSE]
  if(add){
    lines(x, trans(y), lwd = 3, ...)
  }else{
    yr <- if(is.null(ylim)) range(c(trans(lowhig), trans(y), 0), na.rm = TRUE) else ylim
    plot(x, trans(y), xlab = xlab, ylab = ylab, type = "n",
         lwd = 3, xlim = xr, ylim = yr, las = 1, ...)
    grid(col = "black")
    lines(x, trans(y), lwd = 3, ...)
  }
  if(ci){
    ok <- !is.na(lowhig[, 1]) & !is.na(lowhig[, 2])
    if(any(ok)){
      polygon(c(x[ok], rev(x[ok])),
              c(trans(lowhig[ok, 1]), rev(trans(lowhig[ok, 2]))),
              border = gray(.5, alpha = .5), col = cicol)
    }
    lines(x, trans(y), lwd = 3, col = cicol)
    lines(x, trans(y), lwd = 2, col = "black", lty = "dotted")
  }
  invisible(list(x = x, y = trans(y), lo = trans(lowhig[, 1]), hi = trans(lowhig[, 2])))
}

##' samjr SSB plot
##'
##' Plot of estimated spawning stock biomass with two-standard-deviation
##' confidence band.
##' @param fit object returned from \code{\link{sam.fit}}.
##' @param ... extra arguments passed to \code{\link{plotit}}.
##' @return invisibly a list returned by \code{plotit}.
##' @export
ssbplot <- function(fit, ...) UseMethod("ssbplot")

##' @rdname ssbplot
##' @method ssbplot sam
##' @export
ssbplot.sam <- function(fit, ...) plotit(fit, "logssb", ylab = "SSB", trans = exp, ...)

##' samjr total stock biomass plot
##'
##' Plot of estimated total stock biomass (\eqn{\sum_a N_{y,a} \cdot SW_{y,a}})
##' with two-standard-deviation confidence band.
##' @param fit object returned from \code{\link{sam.fit}}.
##' @param ... extra arguments passed to \code{\link{plotit}}.
##' @return invisibly a list returned by \code{plotit}.
##' @export
tsbplot <- function(fit, ...) UseMethod("tsbplot")

##' @rdname tsbplot
##' @method tsbplot sam
##' @export
tsbplot.sam <- function(fit, ...) plotit(fit, "logtsb", ylab = "TSB", trans = exp, ...)

##' @rdname tsbplot
##' @method tsbplot samforecast
##' @export
tsbplot.samforecast <- function(fit, ...){
  plotit(fit, "logtsb", ylab = "TSB", trans = exp, ...)
  addforecast(fit, "tsb")
}

##' @rdname tsbplot
##' @method tsbplot samset
##' @export
tsbplot.samset <- function(fit, ...) plotit(fit, "logtsb", ylab = "TSB", trans = exp, ...)

##' samjr stock-recruitment plot
##'
##' Plots recruitment in year \eqn{t} against spawning stock biomass in
##' year \eqn{t - \mathrm{minAge}}, with year labels and (optionally)
##' the fitted Ricker / Beverton-Holt curve overlaid. Mirrors
##' \code{stockassessment::srplot} but omits the confidence ellipses
##' (which require a joint-covariance ADREPORT samjr does not currently
##' expose).
##'
##' @param fit object returned from \code{\link{sam.fit}}.
##' @param textcol colour for the year labels.
##' @param years if \code{TRUE} (default), label points with the year.
##' @param linetype type for the connecting line (default \code{"l"}).
##' @param linecol colour of the connecting line.
##' @param xlim,ylim optional axis ranges.
##' @param add if \code{TRUE}, overlay on an existing plot.
##' @param addCurve if \code{TRUE} (default), overlay the fitted
##' stock-recruitment curve when \code{conf$stockRecruitmentModelCode}
##' is 1 (Ricker) or 2 (Beverton-Holt).
##' @param curveCol colour of the fitted SR curve.
##' @param ... further arguments passed to \code{plot}/\code{lines}.
##' @return invisibly a list with the plotted \code{S} (SSB) and
##' \code{R} (recruitment) vectors.
##' @export
srplot <- function(fit, ...) UseMethod("srplot")

##' @rdname srplot
##' @method srplot sam
##' @export
srplot.sam <- function(fit, textcol = "red", years = TRUE,
                       linetype = "l", linecol = "black",
                       xlim = NULL, ylim = NULL, add = FALSE,
                       addCurve = TRUE, curveCol = "black", ...){
  X   <- summary(fit)
  n   <- nrow(X)
  lag <- fit$conf$minAge
  idxR <- (lag + 1):n
  idxS <- 1:(n - lag)
  R <- X[idxR, 1]
  S <- X[idxS, 4]
  Snam <- colnames(X)[4]
  Rnam <- colnames(X)[1]
  yLab <- rownames(X)[idxR]
  xr <- if(is.null(xlim)) range(0, S) else xlim
  yr <- if(is.null(ylim)) range(0, R) else ylim
  if(!add){
    plot(S, R, xlab = Snam, ylab = Rnam, type = "n",
         xlim = xr, ylim = yr, las = 1, ...)
  }
  srm <- fit$conf$stockRecruitmentModelCode
  if(addCurve && srm %in% c(1L, 2L)){
    sgrid <- seq(max(xr[1], 1e-8), xr[2], length.out = 200)
    if(srm == 1L){
      a <- fit$pl$rickerpar[1]; b <- fit$pl$rickerpar[2]
      logRcurve <- a + log(sgrid) - exp(b) * sgrid
    }else{
      a <- fit$pl$bhpar[1];     b <- fit$pl$bhpar[2]
      logRcurve <- a + log(sgrid) - log(1 + exp(b) * sgrid)
    }
    lines(sgrid, exp(logRcurve), col = curveCol, lwd = 2)
  }
  lines(S, R, col = linecol, type = linetype, ...)
  if(years) text(S, R, labels = yLab, cex = 0.7, col = textcol)
  invisible(list(S = S, R = R, year = yLab))
}

##' samjr Fbar plot
##'
##' Plot of estimated average fishing mortality (\eqn{\bar F}) over the
##' age range given by \code{fit$conf$fbarRange}, with two-standard-deviation
##' confidence band.
##' @param fit object returned from \code{\link{sam.fit}}.
##' @param ... extra arguments passed to \code{\link{plotit}}.
##' @return invisibly a list returned by \code{plotit}.
##' @export
fbarplot <- function(fit, ...) UseMethod("fbarplot")

##' @rdname fbarplot
##' @method fbarplot sam
##' @export
fbarplot.sam <- function(fit, ...){
  fr <- fit$conf$fbarRange
  lab <- substitute(bar(F)[X - Y], list(X = fr[1], Y = fr[2]))
  plotit(fit, "logfbar", ylab = lab, trans = exp, ...)
}

##' samjr recruitment plot
##'
##' Plot of estimated recruitment (numbers in the lowest age class) with
##' two-standard-deviation confidence band.
##' @param fit object returned from \code{\link{sam.fit}}.
##' @param ... extra arguments passed to \code{\link{plotit}}.
##' @return invisibly a list returned by \code{plotit}.
##' @export
recplot <- function(fit, ...) UseMethod("recplot")

##' @rdname recplot
##' @method recplot sam
##' @export
recplot.sam <- function(fit, ...){
  lab <- paste("Recruits (age ", fit$conf$minAge, ")", sep = "")
  plotit(fit, "logR", ylab = lab, trans = exp, ...)
}

##' samjr total catch plot
##'
##' Plot of estimated total catch in weight with two-standard-deviation
##' confidence band, optionally overlaid with the sum-of-products of the
##' observed catch in numbers and the catch mean weight.
##' @param fit object returned from \code{\link{sam.fit}}.
##' @param obs.show if \code{TRUE} (default) overlay the observed
##' sum-of-products catch in weight as crosses.
##' @param ... extra arguments passed to \code{\link{plotit}}.
##' @return invisibly a list returned by \code{plotit}.
##' @export
catchplot <- function(fit, obs.show = TRUE, ...) UseMethod("catchplot")

##' Internal: sum-of-products of observed catch numbers and catch mean weight.
##' @keywords internal
##' @noRd
catchSOP <- function(fit){
  CW <- fit$data$catchMeanWeight
  if(is.null(CW)) return(NULL)
  aux    <- fit$data$aux
  logobs <- fit$data$logobs
  catchFleet <- which(fit$data$fleetTypes == 0)[1]
  isCatch <- aux[, "fleet"] == catchFleet
  yrs <- as.integer(rownames(CW))
  sopY <- sapply(yrs, function(y){
    ii <- which(isCatch & aux[, "year"] == y)
    if(length(ii) == 0) return(NA_real_)
    ages <- aux[ii, "age"]
    acol <- match(as.character(ages), colnames(CW))
    sum(exp(logobs[ii]) * CW[as.character(y), acol], na.rm = TRUE)
  })
  list(x = yrs, y = sopY)
}

##' @rdname catchplot
##' @method catchplot sam
##' @export
catchplot.sam <- function(fit, obs.show = TRUE, ...){
  CW <- fit$data$catchMeanWeight
  if(!is.null(CW)){
    bad <- apply(is.na(CW), 1, any)
    if(any(bad)){
      idx <- which(names(fit$sdrep$value) == "logCatch")
      fit$sdrep$value[idx[bad]] <- NA_real_
      fit$sdrep$sd[idx[bad]]    <- NA_real_
    }
  }
  ret <- plotit(fit, "logCatch", ylab = "Catch", trans = exp, ...)
  if(obs.show){
    sop <- catchSOP(fit)
    if(!is.null(sop)) points(sop$x, sop$y, pch = 4, lwd = 2, cex = 1.2)
  }
  invisible(ret)
}

##' Default 2x2 panel plot for a samjr fit
##'
##' Produces a 2-by-2 layout showing SSB, Fbar, recruitment and catch.
##' @param x a \code{sam} object returned by \code{\link{sam.fit}}.
##' @param ... extra arguments passed to each panel.
##' @return invisible \code{NULL}.
##' @method plot sam
##' @export
plot.sam <- function(x, ...){
  oldpar <- par(mfrow = c(2, 2))
  on.exit(par(oldpar))
  ssbplot(x, ...)
  fbarplot(x, ...)
  recplot(x, ...)
  catchplot(x, ...)
  invisible(NULL)
}

##' Internal: add forecast point estimates and 95% intervals to an existing plot.
##' Mirrors \code{stockassessment::addforecast}.
##' @keywords internal
##' @noRd
##' @importFrom graphics arrows legend
addforecast <- function(fit, what, ...) UseMethod("addforecast")

##' @rdname addforecast
##' @method addforecast samforecast
##' @keywords internal
##' @noRd
addforecast.samforecast <- function(fit, what, dotcol = "black", dotpch = 19, dotcex = 1.5,
                                    intervalcol = gray(.5, alpha = .5), ...){
  tab <- attr(fit, "tab")
  yr  <- as.numeric(rownames(tab))
  for(i in seq_along(yr)){
    xx <- c(tab[i, paste(what, "low",  sep = ":")],
            tab[i, paste(what, "high", sep = ":")])
    units <- par(c("usr", "pin"))
    xx_to_inches <- units$pin[2L] / diff(units$usr[3:4])
    if(abs(xx_to_inches * diff(xx)) > 0.01){
      arrows(yr[i], xx[1], yr[i], xx[2], lwd = 3,
             col = intervalcol, angle = 90, code = 3, length = .1)
    }
  }
  points(yr, tab[, paste(what, "Estimate", sep = ":")],
         pch = dotpch, cex = dotcex, col = dotcol)
}

##' @rdname plotit
##' @method plotit samset
##' @param addCI logical (length 1 or \code{length(fit)}) indicating whether
##' confidence intervals should be drawn for the added fits in a
##' \code{samset}. Defaults to \code{FALSE} for all added fits.
##' @export
plotit.samset <- function(fit, what, x = NULL, ylab = what, xlab = "Years",
                          trans = function(x) x, add = FALSE, ci = TRUE,
                          cicol = gray(.5, alpha = .5), drop = 0,
                          unnamed.basename = "current", xlim = NULL,
                          addCI = rep(FALSE, length(fit)), ...){
  if(is.logical(addCI) && length(addCI) == 1) addCI <- rep(addCI, length(fit))
  colSet <- c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933",
              "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499")
  idxfrom <- 1
  leg <- names(fit)
  if(is.null(attr(fit, "fit"))){
    attr(fit, "fit") <- fit[[1]]
    idxfrom <- 2
  }else{
    leg <- c(unnamed.basename, leg)
  }
  if(is.null(x)) x <- attr(fit, "fit")$data$years
  xr <- if(is.null(xlim)) range(x) else xlim
  plotit(attr(fit, "fit"), what = what, x = x, ylab = ylab, xlab = xlab,
         trans = trans, add = add, ci = ci, cicol = cicol, drop = drop, xlim = xr, ...)
  invisible(lapply(idxfrom:length(fit), function(i)
    plotit(fit[[i]], what = what, trans = trans, add = TRUE, ci = addCI[i],
           col = colSet[(i - 1) %% length(colSet) + 1],
           cicol = paste0(colSet[(i - 1) %% length(colSet) + 1], "80"),
           drop = drop, ...)))
  if(!is.null(names(fit))){
    legend("bottom", legend = leg, lwd = 3,
           col = c(par("col"),
                   colSet[((idxfrom:length(fit)) - 1) %% length(colSet) + 1]),
           ncol = 3, bty = "n")
    legend("bottom", legend = leg, lwd = 2,
           col = rep("black", length(leg)), lty = "dotted",
           ncol = 3, bty = "n")
  }
}

##' @rdname plotit
##' @method plotit samforecast
##' @export
plotit.samforecast <- function(fit, what, x = NULL, ylab = what, xlab = "Years",
                               trans = function(x) x, add = FALSE, ci = TRUE,
                               cicol = gray(.5, alpha = .5), drop = 0,
                               xlim = NULL, ylim = NULL, ...){
  thisfit <- attr(fit, "fit")
  if(is.null(x)) x <- thisfit$data$years
  xy <- vapply(fit, function(xx) xx$year, numeric(1))
  xr <- if(is.null(xlim)) range(thisfit$data$years, xy) else xlim
  if(is.null(ylim)){
    keep <- seq_len(length(x) - drop)
    v1 <- tableit(fit, what = what, trans = trans, x = x[keep])
    v2 <- tableit(thisfit, what = what, trans = trans, x = x[keep])
    ylim <- range(v1, v2, na.rm = TRUE)
  }
  plotit(thisfit, what = what, ylab = ylab, xlab = xlab, trans = trans,
         add = add, ci = ci, cicol = cicol, drop = drop, xlim = xr, ylim = ylim, ...)
}

##' @rdname ssbplot
##' @method ssbplot samforecast
##' @export
ssbplot.samforecast <- function(fit, ...){
  plotit(fit, "logssb", ylab = "SSB", trans = exp, ...)
  addforecast(fit, "ssb")
}

##' @rdname ssbplot
##' @method ssbplot samset
##' @export
ssbplot.samset <- function(fit, ...) plotit(fit, "logssb", ylab = "SSB", trans = exp, ...)

##' @rdname fbarplot
##' @method fbarplot samforecast
##' @export
fbarplot.samforecast <- function(fit, ...){
  fr  <- attr(fit, "fit")$conf$fbarRange
  lab <- substitute(bar(F)[X - Y], list(X = fr[1], Y = fr[2]))
  plotit(fit, "logfbar", ylab = lab, trans = exp, ...)
  addforecast(fit, "fbar")
}

##' @rdname fbarplot
##' @method fbarplot samset
##' @export
fbarplot.samset <- function(fit, ...){
  fitlocal <- attr(fit, "fit"); if(is.null(fitlocal)) fitlocal <- fit[[1]]
  fr  <- fitlocal$conf$fbarRange
  lab <- substitute(bar(F)[X - Y], list(X = fr[1], Y = fr[2]))
  plotit(fit, "logfbar", ylab = lab, trans = exp, ...)
}

##' @rdname recplot
##' @method recplot samforecast
##' @export
recplot.samforecast <- function(fit, ...){
  fitlocal <- attr(fit, "fit")
  lab <- paste0("Recruits (age ", fitlocal$conf$minAge, ")")
  plotit(fit, "logR", ylab = lab, trans = exp, ...)
  addforecast(fit, "rec")
}

##' @rdname recplot
##' @method recplot samset
##' @export
recplot.samset <- function(fit, ...){
  fitlocal <- attr(fit, "fit"); if(is.null(fitlocal)) fitlocal <- fit[[1]]
  lab <- paste0("Recruits (age ", fitlocal$conf$minAge, ")")
  plotit(fit, "logR", ylab = lab, trans = exp, ...)
}

##' @rdname catchplot
##' @method catchplot samforecast
##' @export
catchplot.samforecast <- function(fit, obs.show = TRUE, ...){
  fitlocal <- attr(fit, "fit")
  ret <- plotit(fit, "logCatch", ylab = "Catch", trans = exp, ...)
  if(obs.show){
    sop <- catchSOP(fitlocal)
    if(!is.null(sop)) points(sop$x, sop$y, pch = 4, lwd = 2, cex = 1.2)
  }
  addforecast(fit, "catch")
  invisible(ret)
}

##' @rdname catchplot
##' @method catchplot samset
##' @export
catchplot.samset <- function(fit, obs.show = TRUE, ...){
  scrub <- function(f){
    CW <- f$data$catchMeanWeight
    if(!is.null(CW)){
      bad <- apply(is.na(CW), 1, any)
      if(any(bad)){
        idx <- which(names(f$sdrep$value) == "logCatch")
        f$sdrep$value[idx[bad]] <- NA_real_
        f$sdrep$sd[idx[bad]]    <- NA_real_
      }
    }
    f
  }
  ref <- attr(fit, "fit")
  if(!is.null(ref)) attr(fit, "fit") <- scrub(ref)
  for(i in seq_along(fit)) fit[[i]] <- scrub(fit[[i]])
  fitlocal <- attr(fit, "fit"); if(is.null(fitlocal)) fitlocal <- fit[[1]]
  ret <- plotit(fit, "logCatch", ylab = "Catch", trans = exp, ...)
  if(obs.show){
    sop <- catchSOP(fitlocal)
    if(!is.null(sop)) points(sop$x, sop$y, pch = 4, lwd = 2, cex = 1.2)
  }
  invisible(ret)
}

##' Default 2x2 panel plot for a samjr forecast
##' @param x a \code{samforecast} object returned by \code{\link{forecast}}.
##' @param ... extra arguments passed to each panel.
##' @return invisible \code{NULL}.
##' @method plot samforecast
##' @export
plot.samforecast <- function(x, ...){
  oldpar <- par(mfrow = c(2, 2))
  on.exit(par(oldpar))
  ssbplot(x, ...)
  fbarplot(x, ...)
  recplot(x, ...)
  catchplot(x, ...)
  invisible(NULL)
}

##' Default 2x2 panel plot for a samjr samset
##' @param x a \code{samset} object (as returned by \code{\link{retro}},
##' \code{\link{leaveout}}, \code{\link{jit}} or \code{\link{c.sam}}).
##' @param ... extra arguments passed to each panel.
##' @return invisible \code{NULL}.
##' @method plot samset
##' @export
plot.samset <- function(x, ...){
  oldpar <- par(mfrow = c(2, 2))
  on.exit(par(oldpar))
  ssbplot(x, ...)
  fbarplot(x, ...)
  recplot(x, ...)
  catchplot(x, ...)
  invisible(NULL)
}

##' samjr data overview plot
##'
##' Plots the data available for a samjr fit: one row per fleet, one tick
##' per (year, age) observation, with a side panel showing the fleet type.
##' Mirrors \code{stockassessment::dataplot}.
##'
##' @param fit object of class \code{sam} from \code{\link{sam.fit}}.
##' @param col colour per fleet; defaults to two alternating colours.
##' @param fleet_type character vector of fleet-type labels. Default is
##' derived from \code{fit$data$fleetTypes}: 0 \dQuote{Catch at age},
##' 1 \dQuote{Catch at age with effort}, 2 or 6 \dQuote{Index at age},
##' 3 \dQuote{Biomass or catch index}, 5 \dQuote{Tagging data},
##' 7 \dQuote{Sum of fleets}.
##' @param fleet_names character vector of fleet names; defaults to
##' \code{attr(fit$data, "fleetNames")}.
##' @return invisible \code{NULL}.
##' @importFrom graphics layout mtext axis
##' @export
dataplot <- function(fit, col = NULL, fleet_type = NULL, fleet_names = NULL)
  UseMethod("dataplot")

##' @rdname dataplot
##' @method dataplot sam
##' @export
dataplot.sam <- function(fit, col = NULL, fleet_type = NULL, fleet_names = NULL){
  years <- fit$data$years
  nf <- fit$data$noFleets
  for(k in 1:nf){
    if(fit$data$minAgePerFleet[k] == -1) fit$data$minAgePerFleet[k] <- NA
    if(fit$data$maxAgePerFleet[k] == -1) fit$data$maxAgePerFleet[k] <- NA
  }
  yspace <- 2
  noage <- length(min(fit$data$minAgePerFleet, na.rm = TRUE):
                  max(fit$data$maxAgePerFleet, na.rm = TRUE))

  dat <- cbind(fit$data$aux, fit$data$logobs)

  if(is.null(col)) col <- c("#67a9cf", "#ef8a62")
  col <- rep(col, length.out = nf)

  ynum <- noage * nf + (nf - 1) * 2
  ylab <- seq((noage / 2) - 1, ynum, (noage + yspace))

  if(is.null(fleet_names)) fleet_names <- attr(fit$data, "fleetNames")
  for(k in 1:nf){
    fleet_names[k] <- abbreviate(fleet_names[k], minlength = 19,
                                 dot = TRUE, use.classes = FALSE)
  }

  if(is.null(fleet_type)){
    fleet_type <- fit$data$fleetTypes
    for(i in 1:nf){
      if(fleet_type[i] == 0) fleet_type[i] <- "Catch at age"
      if(fleet_type[i] == 1) fleet_type[i] <- "Catch at age with effort"
      if(fleet_type[i] == 2 || fleet_type[i] == 6) fleet_type[i] <- "Index at age"
      if(fleet_type[i] == 3) fleet_type[i] <- "Biomass or catch index"
      if(fleet_type[i] == 5) fleet_type[i] <- "Tagging data"
      if(fleet_type[i] == 7) fleet_type[i] <- "Sum of fleets"
    }
  }

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  layout(matrix(c(rep(1, 3), 2), nrow = 1))
  par(oma = c(2, 10, 2, 0), mar = c(0, 0, 0, 0), xpd = NA)
  plot(x = rep(years, (ynum + 1)), y = rep(0:ynum, length(years)),
       type = "n", xlab = "", ylab = "", yaxt = "n")
  axis(2, labels = fleet_names, at = ylab, las = 1, cex = 0.8)
  x <- 0 - min(fit$data$minAgePerFleet, na.rm = TRUE)
  for(i in 1:nf){
    if(!is.na(fit$data$minAgePerFleet[i])){
      for(a in min(fit$data$minAgePerFleet, na.rm = TRUE):
              max(fit$data$maxAgePerFleet, na.rm = TRUE)){
        lines(x = years, y = rep(x + a, length(years)), col = "grey87")
      }
      for(a in fit$data$minAgePerFleet[i]:fit$data$maxAgePerFleet[i]){
        yval <- dat[which(dat[, 2] == i & dat[, 3] == a), 4]
        for(k in seq_along(yval)){
          if(!is.na(yval[k]) && yval[k] != 0) yval[k] <- x + a else yval[k] <- NA
        }
        lines(x = dat[which(dat[, 2] == i & dat[, 3] == a), 1],
              y = yval, lwd = 1, col = col[i])
        points(x = dat[which(dat[, 2] == i & dat[, 3] == a), 1],
               y = yval, lwd = 1, col = col[i], pch = 16)
        text(x = years[length(years)] + 0.5, y = x + a, labels = a, cex = 0.7)
      }
    }else{
      a <- -1
      yval <- dat[which(dat[, 2] == i & dat[, 3] == a), 4]
      for(k in seq_along(yval)){
        if(!is.na(yval[k]) && yval[k] != 0) yval[k] <- x + noage / 2 else yval[k] <- NA
      }
      lines(x = dat[which(dat[, 2] == i & dat[, 3] == a), 1],
            y = yval, lwd = 1, col = col[i])
      points(x = dat[which(dat[, 2] == i & dat[, 3] == a), 1],
             y = yval, lwd = 1, col = col[i], pch = 16)
    }
    x <- x + noage + yspace
  }
  plot(x = rep(1, (ynum + 1)), y = rep(0:ynum), type = "n",
       xlab = "", ylab = "", yaxt = "n", bty = "n", xaxt = "n")
  for(i in 1:nf) text(x = 1, y = ylab[i], labels = fleet_type[i])
  mtext(text = "Available data", side = 3, line = 0.5, at = 3 / 8, outer = TRUE)
  mtext(text = "Data type",      side = 3, line = 0.5, at = 7 / 8, outer = TRUE)
  invisible(NULL)
}

##' samjr F-selectivity plot
##'
##' Stacked-bar plot of normalised \eqn{F_{y,a}} (each year sums to 1)
##' so the bar height per age shows the relative F-selectivity. Mirrors
##' \code{stockassessment::fselectivityplot}; renamed for brevity.
##' @param fit a fitted \code{sam} object.
##' @param cexAge size multiplier for the age labels overlaid on the bars.
##' @param ... further arguments passed to \code{barplot}.
##' @return invisible \code{NULL}.
##' @importFrom graphics barplot text
##' @export
selplot <- function(fit, cexAge = 1, ...) UseMethod("selplot")

##' @rdname selplot
##' @method selplot sam
##' @export
selplot.sam <- function(fit, cexAge = 1, ...){
  fmat <- faytable(fit)
  fmat[is.na(fmat)] <- 0
  P <- t(fmat / rowSums(fmat))
  barplot(P, border = NA, space = 0, xlab = "Year",
          main = "Selectivity in F", ...)
  text(1, cumsum(P[, 1]) - 0.5 * P[, 1],
       labels = as.character(seq_len(ncol(fmat))),
       adj = c(0, 0.2), cex = cexAge)
  invisible(NULL)
}

##' samjr observation-SD bar plot
##'
##' Bar plot of estimated per-fleet, per-age observation standard
##' deviations in increasing order, with each fleet a different colour.
##' Mirrors \code{stockassessment::sdplot}.
##' @param fit a fitted \code{sam} object.
##' @param barcol vector of bar colours (one per fleet); a default is
##' derived from \code{colors()}.
##' @param marg \code{mar} setting for the plot.
##' @param ylim y-axis range.
##' @param show.rel.w if \code{TRUE}, plot the relative weight
##' \eqn{(1/\sigma^2)/\max(1/\sigma^2)} instead of the SD.
##' @param ... unused.
##' @return invisible \code{NULL}.
##' @importFrom graphics barplot box
##' @importFrom grDevices colors
##' @export
sdplot <- function(fit, barcol = NULL, marg = NULL, ylim = NULL, ...)
  UseMethod("sdplot")

##' @rdname sdplot
##' @method sdplot sam
##' @export
sdplot.sam <- function(fit, barcol = NULL, marg = NULL, ylim = NULL,
                       show.rel.w = FALSE, ...){
  cf <- fit$conf$keyVarObs
  fn <- attr(fit$data, "fleetNames")
  ages <- fit$conf$minAge:fit$conf$maxAge
  pt <- partable(fit)
  sd <- unname(exp(pt[grep("logSdLogObs", rownames(pt)), 1]))
  v <- cf
  v[] <- c(NA, sd)[cf + 2]
  res <- data.frame(
    fleet = fn[as.vector(row(v))],
    name  = paste0(fn[as.vector(row(v))], " age ", ages[as.vector(col(v))]),
    sd    = as.vector(v),
    stringsAsFactors = FALSE)
  res <- res[stats::complete.cases(res), ]
  res <- res[order(res$sd), ]
  res$rel.w <- (1 / res$sd^2) / max(1 / res$sd^2)
  if(is.null(barcol))
    barcol <- grDevices::colors()[as.integer(as.factor(res$fleet)) * 10]
  if(is.null(marg))
    marg <- c(max(nchar(fn)) * 0.8, 4, 2, 1)
  oldpar <- par(mar = marg)
  on.exit(par(oldpar))
  if(show.rel.w){
    barplot(res$rel.w, names.arg = res$name, las = 2, col = barcol,
            ylab = "Relative weight", ylim = ylim)
  }else{
    barplot(res$sd,   names.arg = res$name, las = 2, col = barcol,
            ylab = "SD", ylim = ylim)
  }
  box()
  invisible(NULL)
}

##' samjr fixed-effect parameter plot
##'
##' One sub-panel per fixed-effect parameter showing the estimate (line)
##' with its 2-SD confidence interval (grey polygon). Pairs of
##' parameters with \code{|correlation| > cor.report.limit} are listed
##' in the top-right (positive, blue) and bottom-right (negative, red)
##' legends. Mirrors \code{stockassessment::parplot}.
##' @param fit a fitted \code{sam} object.
##' @param cor.report.limit threshold above which a parameter
##' correlation is annotated on the panel.
##' @param ... further arguments passed to \code{plot}.
##' @return invisible \code{NULL}.
##' @importFrom graphics box layout legend lines polygon
##' @export
parplot <- function(fit, cor.report.limit = 0.95, ...)
  UseMethod("parplot")

##' @rdname parplot
##' @method parplot sam
##' @export
parplot.sam <- function(fit, cor.report.limit = 0.95, ...){
  param <- coef(fit)
  nam   <- names(param)
  dup   <- duplicated(nam)
  namadd <- rep(0L, length(nam))
  for(i in seq_along(dup)[-1]) if(dup[i]) namadd[i] <- namadd[i - 1] + 1L
  nam <- paste(nam, namadd, sep = "_")
  Sigma <- attr(param, "cov")
  corrs <- stats::cov2cor(Sigma) - diag(length(param))
  rownames(corrs) <- nam; colnames(corrs) <- nam
  higcor <- lapply(seq_len(nrow(corrs)),
                   function(i) corrs[i, ][corrs[i, ] >  cor.report.limit] * 100)
  lowcor <- lapply(seq_len(nrow(corrs)),
                   function(i) corrs[i, ][corrs[i, ] < -cor.report.limit] * 100)
  names(higcor) <- nam; names(lowcor) <- nam
  sds <- attr(param, "sd")
  est <- as.numeric(param)
  bands <- cbind(est - 2 * sds, est, est + 2 * sds)
  plotOne <- function(i){
    nm <- nam[i]
    y  <- bands[i, ]
    x  <- c(-0.5, 0, 0.5)
    plot(x, y, ylim = range(y), xlab = "", ylab = nm,
         axes = FALSE, type = "n", ...)
    box(); axis(2, las = 1)
    polygon(c(x, rev(x)), c(rep(y[1], 3), rep(y[3], 3)),
            border = gray(0.5, alpha = 0.5), col = gray(0.5, alpha = 0.5))
    lines(x, rep(y[2], 3), lwd = 3)
    if(length(higcor[[nm]]))
      legend("topright",
             legend = paste0(names(higcor[[nm]]), ": ",
                              round(higcor[[nm]]), "%"),
             bty = "n", text.col = "blue", cex = 0.75)
    if(length(lowcor[[nm]]))
      legend("bottomright",
             legend = paste0(names(lowcor[[nm]]), ": ",
                              round(lowcor[[nm]]), "%"),
             bty = "n", text.col = "red", cex = 0.75)
  }
  np  <- length(nam)
  div <- rep(ceiling(sqrt(np)), 2)
  if(div[1] * (div[2] - 1) >= np) div[2] <- div[2] - 1
  oldpar <- par(mar = c(0.2, par("mar")[2], 0.2, par("mar")[4]))
  on.exit(par(oldpar))
  layout(matrix(seq_len(div[1] * div[2]), nrow = div[1], ncol = div[2]))
  for(i in seq_len(np)) plotOne(i)
  invisible(NULL)
}

##' samjr between-age observation-correlation plot
##'
##' One panel per fleet showing the estimated between-age correlation
##' matrix from \code{fit$rep$obsCov} (\code{corplot.sam}) or the
##' empirical correlation of one-step residuals
##' (\code{corplot.samres}), drawn as correlation ellipses via
##' \code{ellipse::plotcorr}. Mirrors \code{stockassessment::corplot}
##' (requires the suggested \pkg{ellipse} package).
##' @param x a fitted \code{sam} object or a \code{samres} object from
##' \code{\link{residuals.sam}}.
##' @param ... further arguments passed to \code{ellipse::plotcorr}.
##' @return invisible \code{NULL}.
##' @export
corplot <- function(x, ...) UseMethod("corplot")

##' Internal: per-fleet correlation-ellipse layout used by
##' \code{corplot.sam} and \code{corplot.samres}.
##' @keywords internal
##' @noRd
corplotCommon <- function(cmats, fleetNames, ...){
  if(!requireNamespace("ellipse", quietly = TRUE))
    stop("Package 'ellipse' is required for corplot; install it.")
  ccolors <- c("#A50F15", "#DE2D26", "#FB6A4A", "#FCAE91", "#FEE5D9",
               "white",
               "#EFF3FF", "#BDD7E7", "#6BAED6", "#3182BD", "#08519C")
  nf <- length(cmats)
  div <- if(nf == 3) c(3, 1) else {
    d <- rep(ceiling(sqrt(nf)), 2)
    if(d[1] * (d[2] - 1) >= nf) d[2] <- d[2] - 1
    d
  }
  oldpar <- par(no.readonly = TRUE)
  on.exit({ try(par(oldpar), silent = TRUE); layout(1) })
  layout(matrix(seq_len(div[1] * div[2]), nrow = div[1], ncol = div[2]))
  for(i in seq_along(cmats)){
    xx <- cmats[[i]]
    if(any(is.na(xx))){ graphics::plot.new(); title(main = substr(fleetNames[i], 1, 20)); next }
    ellipse::plotcorr(xx, col = ccolors[5 * xx + 6],
                      mar = 0.1 + c(2, 2, 2, 2),
                      main = substr(fleetNames[i], 1, 20), ...)
  }
  invisible(NULL)
}

##' @rdname corplot
##' @method corplot sam
##' @export
corplot.sam <- function(x, ...){
  if(is.null(x$rep$obsCov))
    stop("fit$rep$obsCov is missing; nothing to plot.")
  fn <- attr(x$data, "fleetNames")
  cmats <- lapply(seq_along(x$rep$obsCov), function(i){
    S <- x$rep$obsCov[[i]]
    if(any(is.na(S))) return(matrix(NA, nrow(S), ncol(S)))
    M <- stats::cov2cor(S)
    a0 <- x$data$minAgePerFleet[i]
    ages <- seq(a0, a0 + nrow(M) - 1L)
    rownames(M) <- ages; colnames(M) <- ages
    M
  })
  corplotCommon(cmats, fn, ...)
}

##' @rdname corplot
##' @method corplot samres
##' @export
corplot.samres <- function(x, ...){
  dat <- data.frame(resid = x$residual, age = x$age,
                    year = x$year, fleet = x$fleet)
  fleets <- sort(unique(dat$fleet))
  fn <- attr(x, "fleetNames")
  cmats <- lapply(fleets, function(ff){
    tab <- stats::xtabs(resid ~ age + year, data = dat[dat$fleet == ff, ])
    M <- stats::cor(t(tab))
    rownames(M) <- rownames(tab); colnames(M) <- rownames(tab)
    M
  })
  fnsub <- if(!is.null(fn)) fn[fleets] else paste("Fleet", fleets)
  corplotCommon(cmats, fnsub, ...)
}

##' samjr observed-vs-predicted fit plot
##'
##' One small panel per (fleet, age) showing the observations (points)
##' and the predicted values from \code{fit$rep$logPred} (line) over
##' time. Mirrors \code{stockassessment::fitplot} but uses a simple
##' base-R grid layout instead of SAM's \code{plotby} helper.
##' @param fit a fitted \code{sam} object.
##' @param log if \code{TRUE} (default), plot on the log scale.
##' @param fleets integer fleet ids to include (default: all).
##' @param ... further arguments passed to \code{plot}.
##' @return invisible \code{NULL}.
##' @export
fitplot <- function(fit, log = TRUE, ...) UseMethod("fitplot")

##' @rdname fitplot
##' @method fitplot sam
##' @export
fitplot.sam <- function(fit, log = TRUE,
                         fleets = unique(fit$data$aux[, "fleet"]), ...){
  trans <- if(log) identity else exp
  aux  <- fit$data$aux
  obs  <- trans(fit$data$logobs)
  pre  <- trans(fit$rep$logPred)
  fn   <- attr(fit$data, "fleetNames")
  keep <- aux[, "fleet"] %in% fleets
  aux  <- aux[keep, , drop = FALSE]
  obs  <- obs[keep]
  pre  <- pre[keep]
  fleetsSorted <- sort(unique(as.integer(fleets)))
  ages <- sort(unique(as.integer(aux[, "age"])))
  ages <- ages[ages >= 0]
  if(length(ages) == 0L) ages <- sort(unique(as.integer(aux[, "age"])))
  nF <- length(fleetsSorted); nA <- length(ages)
  ylab <- if(log) "log obs" else "obs"
  oldpar <- par(mfrow = c(nA, nF), mar = c(2, 3, 1.5, 0.5),
                mgp = c(1.8, 0.6, 0), oma = c(2, 0, 0, 0))
  on.exit(par(oldpar))
  for(a in ages){
    for(ff in fleetsSorted){
      ii <- which(aux[, "fleet"] == ff & aux[, "age"] == a)
      if(length(ii) == 0L){
        graphics::plot.new(); next
      }
      yr <- aux[ii, "year"]; o <- order(yr)
      main <- if(a == ages[1]) strtrim(fn[ff], 30) else ""
      ylabHere <- if(ff == fleetsSorted[1]) paste0(ylab, " a=", a) else ""
      plot(yr[o], obs[ii][o], pch = 19, cex = 0.7,
           xlab = "", ylab = ylabHere, main = main, ...)
      lines(yr[o], pre[ii][o], col = "blue", lwd = 2)
    }
  }
  mtext("Year", side = 1, line = 0.5, outer = TRUE)
  invisible(NULL)
}

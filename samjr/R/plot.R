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
    polygon(c(x, rev(x)),
            c(trans(lowhig[, 1]), rev(trans(lowhig[, 2]))),
            border = gray(.5, alpha = .5), col = cicol)
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
##' @export
plotit.samset <- function(fit, what, x = NULL, ylab = what, xlab = "Years",
                          trans = function(x) x, add = FALSE, ci = TRUE,
                          cicol = gray(.5, alpha = .5), drop = 0,
                          unnamed.basename = "current", xlim = NULL, ...){
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
    plotit(fit[[i]], what = what, trans = trans, add = TRUE, ci = FALSE,
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
